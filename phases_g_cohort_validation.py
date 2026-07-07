"""
Thin cohort validation orchestrator: merged phenotype + PRS TSV → metrics + plots.

Does NOT reimplement GWAS. Discrimination, calibration, and DCA use scikit-learn
and matplotlib (analogous to common R pROC / rmda workflows, not byte-identical).
LD clumping is delegated to PLINK 2 via the ``clump-plink2`` subcommand. For
PRS-CS, lassosum, ancestry PCs, and similar, use external tools and put the
resulting PRS column (plus outcome and PHASES covariates) into the merged TSV this
script reads.

Typical input: one TSV from merging phases_g_prs_from_genotypes.py output with
outcome and PHASES columns (same sample_id).

  python phases_g_cohort_validation.py \\
      --cohort results/phases_g/merged_cohort.tsv \\
      --outcome-col y \\
      --prs-col PRS_z \\
      --baseline-cols phases_sum,size_mm \\
      --out-dir results/phases_g/cohort_validation_run1

Optional LD clumping (requires PLINK 2 on PATH and genotype + association files):

  python phases_g_cohort_validation.py clump-plink2 \\
      --plink2-exec plink2 \\
      --bfile /path/prefix \\
      --clump-assoc /path/gwas_or_cohort.assoc.logistic \\
      --out-dir results/phases_g/clump_out

Clumping uses PLINK's implementation; supply an association file with variant ID
and P column as required by your PLINK build (--clump-* args vary; we pass
standard flags).

See data/genotypes/README_GENOTYPE_INPUT.txt (cohort validation section).
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent


def _net_benefit(y: np.ndarray, pred_prob: np.ndarray, pt: float) -> float:
    """Vickers–Elkin style net benefit at threshold probability pt."""
    y = np.asarray(y).astype(int)
    p = np.asarray(pred_prob, dtype=float)
    pred = (p >= pt).astype(int)
    n = len(y)
    if n == 0:
        return 0.0
    tp = np.sum((pred == 1) & (y == 1))
    fp = np.sum((pred == 1) & (y == 0))
    denom = 1.0 - pt
    if denom < 1e-12:
        return 0.0
    return tp / n - (fp / n) * (pt / denom)


def _decision_curve(y: np.ndarray, pred_prob: np.ndarray, thresholds: np.ndarray) -> dict:
    y = np.asarray(y).astype(int)
    prev = float(np.mean(y))
    nb_model = [_net_benefit(y, pred_prob, t) for t in thresholds]
    nb_all = [prev - (1.0 - prev) * (t / (1.0 - t + 1e-12)) for t in thresholds]
    nb_none = [0.0 for _ in thresholds]
    return {"threshold": thresholds.tolist(), "model": nb_model, "treat_all": nb_all, "treat_none": nb_none}


def run_validation(args) -> dict:
    from sklearn.calibration import calibration_curve
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import brier_score_loss, roc_auc_score, roc_curve
    from sklearn.model_selection import train_test_split
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cohort_path: Path = args.cohort
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(cohort_path, sep="\t", low_memory=False)
    if args.sample_id_col not in df.columns:
        raise SystemExit(f"Missing column {args.sample_id_col!r}")
    if args.outcome_col not in df.columns:
        raise SystemExit(f"Missing outcome column {args.outcome_col!r}")
    if args.prs_col not in df.columns:
        raise SystemExit(f"Missing PRS column {args.prs_col!r}")

    d = df.copy()
    d[args.outcome_col] = pd.to_numeric(d[args.outcome_col], errors="coerce")
    d = d.dropna(subset=[args.outcome_col, args.prs_col])
    y = d[args.outcome_col].astype(int).values
    if len(np.unique(y)) < 2:
        raise SystemExit("Outcome must have at least two classes after cleaning.")

    prs = pd.to_numeric(d[args.prs_col], errors="coerce").values.astype(float)
    prs_mean, prs_std = float(np.nanmean(prs)), float(np.nanstd(prs))
    prs_z = (prs - prs_mean) / prs_std if prs_std > 1e-12 else np.zeros_like(prs)
    d["_prs_scaled"] = prs_z

    baseline_cols = [c.strip() for c in args.baseline_cols.split(",") if c.strip()]
    for c in baseline_cols:
        if c not in d.columns:
            raise SystemExit(f"Baseline column {c!r} not in cohort TSV")

    Xb = d[baseline_cols].apply(pd.to_numeric, errors="coerce") if baseline_cols else pd.DataFrame(index=d.index)
    Xb = Xb.fillna(Xb.median(numeric_only=True))
    Xf = pd.concat([Xb.reset_index(drop=True), pd.Series(prs_z, name="prs_z")], axis=1)

    summary: dict = {
        "n": int(len(d)),
        "n_events": int(np.sum(y == 1)),
        "n_nonevents": int(np.sum(y == 0)),
        "outcome_col": args.outcome_col,
        "prs_col": args.prs_col,
        "baseline_cols": baseline_cols,
        "prs_rescaled_in_script": True,
    }

    def build_lr_pipe(*, intercept_only: bool) -> Pipeline:
        lr = LogisticRegression(
            max_iter=4000,
            class_weight="balanced",
            solver="lbfgs",
            fit_intercept=not intercept_only,
        )
        if intercept_only:
            return Pipeline([("clf", lr)])
        return Pipeline([("sc", StandardScaler()), ("clf", lr)])

    def fit_auc(X: pd.DataFrame, yv: np.ndarray, name: str) -> tuple[Pipeline, float]:
        intercept_only = X.shape[1] == 0
        Xmat = np.ones((len(yv), 1), dtype=float) if intercept_only else X.values.astype(float)
        pipe = build_lr_pipe(intercept_only=intercept_only)
        pipe.fit(Xmat, yv)
        prob = pipe.predict_proba(Xmat)[:, 1]
        auc = float(roc_auc_score(yv, prob))
        summary[f"auc_train_{name}"] = auc
        summary[f"brier_train_{name}"] = float(brier_score_loss(yv, prob))
        return pipe, auc

    pipe_base, auc_base = fit_auc(Xb, y, "baseline")
    pipe_full, auc_full = fit_auc(Xf, y, "prs")

    summary["delta_auc_train"] = auc_full - auc_base

    # Held-out split if enough samples
    use_split = len(y) >= max(40, args.min_split_n) and min(np.sum(y == 1), np.sum(y == 0)) >= args.min_class_test
    if use_split:
        Xb_tr, Xb_te, Xf_tr, Xf_te, y_tr, y_te = train_test_split(
            Xb.values,
            Xf.values,
            y,
            test_size=args.test_fraction,
            stratify=y,
            random_state=args.seed,
        )
        io = Xb.shape[1] == 0
        pipe_b = build_lr_pipe(intercept_only=io)
        pipe_f = build_lr_pipe(intercept_only=False)
        if io:
            Xb_tr = np.ones((len(y_tr), 1), dtype=float)
            Xb_te = np.ones((len(y_te), 1), dtype=float)
        pipe_b.fit(Xb_tr, y_tr)
        pipe_f.fit(Xf_tr, y_tr)
        p_base = pipe_b.predict_proba(Xb_te)[:, 1]
        p_full = pipe_f.predict_proba(Xf_te)[:, 1]
        summary["auc_test_baseline"] = float(roc_auc_score(y_te, p_base))
        summary["auc_test_prs"] = float(roc_auc_score(y_te, p_full))
        summary["delta_auc_test"] = summary["auc_test_prs"] - summary["auc_test_baseline"]
        summary["n_test"] = int(len(y_te))
        y_plot, prob_plot = y_te, p_full
        roc_y, roc_p_base, roc_p_full = y_te, p_base, p_full
    else:
        summary["note_split"] = (
            f"Skipped train/test split (need n>={max(40, args.min_split_n)} and "
            f"min class count>={args.min_class_test}); metrics are in-sample."
        )
        if Xb.shape[1] == 0:
            p_base = pipe_base.predict_proba(np.ones((len(y), 1), dtype=float))[:, 1]
        else:
            p_base = pipe_base.predict_proba(Xb.values.astype(float))[:, 1]
        p_full = pipe_full.predict_proba(Xf.values.astype(float))[:, 1]
        summary["auc_test_baseline"] = None
        summary["auc_test_prs"] = None
        summary["delta_auc_test"] = None
        y_plot, prob_plot = y, p_full
        roc_y, roc_p_base, roc_p_full = y, p_base, p_full

    # ROC
    fpr_b, tpr_b, _ = roc_curve(roc_y, roc_p_base)
    fpr_f, tpr_f, _ = roc_curve(roc_y, roc_p_full)
    plt.figure(figsize=(6, 5))
    plt.plot(fpr_b, tpr_b, label=f"Baseline AUC={summary.get('auc_test_baseline') or summary['auc_train_baseline']:.3f}")
    plt.plot(fpr_f, tpr_f, label=f"+ PRS AUC={summary.get('auc_test_prs') or summary['auc_train_prs']:.3f}")
    plt.plot([0, 1], [0, 1], "k--", alpha=0.4)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title("ROC (test set if split; else train)")
    plt.legend(loc="lower right")
    plt.tight_layout()
    roc_path = out_dir / f"{args.prefix}roc.png"
    plt.savefig(roc_path, dpi=150)
    plt.close()
    summary["plot_roc"] = str(roc_path.resolve())

    # Calibration (full model)
    prob_plot = np.clip(prob_plot, 1e-6, 1 - 1e-6)
    frac_pos, mean_pred = calibration_curve(y_plot, prob_plot, n_bins=args.calibration_bins, strategy="quantile")
    plt.figure(figsize=(6, 5))
    plt.plot(mean_pred, frac_pos, "s-", label="Model")
    plt.plot([0, 1], [0, 1], "k--", alpha=0.4, label="Ideal")
    plt.xlabel("Mean predicted probability")
    plt.ylabel("Fraction positives")
    plt.title("Calibration (quantile bins)")
    plt.legend()
    plt.tight_layout()
    cal_path = out_dir / f"{args.prefix}calibration.png"
    plt.savefig(cal_path, dpi=150)
    plt.close()
    summary["plot_calibration"] = str(cal_path.resolve())

    # Decision curve
    t_min, t_max = 0.01, 0.99
    thresholds = np.linspace(t_min, t_max, 80)
    dc = _decision_curve(y_plot, prob_plot, thresholds)
    plt.figure(figsize=(6, 5))
    plt.plot(dc["threshold"], dc["model"], label="Model (PRS + baseline)")
    plt.plot(dc["threshold"], dc["treat_all"], "--", label="Treat all", alpha=0.7)
    plt.plot(dc["threshold"], dc["treat_none"], "k-", label="Treat none", alpha=0.5)
    plt.xlabel("Threshold probability")
    plt.ylabel("Net benefit")
    plt.title("Decision curve (simple Vickers–Elkin)")
    plt.xlim(t_min, t_max)
    plt.legend(loc="upper right")
    plt.tight_layout()
    dca_path = out_dir / f"{args.prefix}decision_curve.png"
    plt.savefig(dca_path, dpi=150)
    plt.close()
    summary["plot_decision_curve"] = str(dca_path.resolve())

    (out_dir / f"{args.prefix}decision_curve_data.json").write_text(
        json.dumps(dc, indent=2), encoding="utf-8"
    )

    (out_dir / f"{args.prefix}validation_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    print(f"Wrote plots and {args.prefix}validation_summary.json under {out_dir}")
    return summary


def run_clump_plink2(args) -> int:
    """Invoke PLINK 2 clumping; user supplies association output compatible with --clump."""
    exe = args.plink2_exec
    if not shutil.which(exe) and not Path(exe).is_file():
        print(f"plink2 not found: {exe}", file=sys.stderr)
        return 1
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / "clump_run"
    cmd = [
        exe,
        "--bfile",
        str(args.bfile),
        "--clump",
        str(args.clump_assoc),
        "--clump-p1",
        str(args.clump_p1),
        "--clump-r2",
        str(args.clump_r2),
        "--clump-kb",
        str(args.clump_kb),
        "--out",
        str(out_prefix),
    ]
    if args.clump_snp_field:
        cmd.extend(["--clump-snp-field", args.clump_snp_field])
    if args.clump_p_field:
        cmd.extend(["--clump-p-field", args.clump_p_field])
    log_path = out_dir / "clump_plink2.log.txt"
    print("RUN:", " ".join(cmd), flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    log_path.write_text(r.stdout + "\n--- stderr ---\n" + r.stderr, encoding="utf-8")
    if r.returncode != 0:
        print(f"PLINK2 exited {r.returncode}; see {log_path}", file=sys.stderr)
        return r.returncode
    print(f"OK; log {log_path}; outputs prefixed {out_prefix}")
    return 0


def main() -> None:
    import argparse

    raw_argv = sys.argv[1:]
    if not raw_argv or raw_argv in (["-h"], ["--help"]):
        raw_argv = ["validate", "--help"]
    elif raw_argv[0] not in ("validate", "clump-plink2"):
        raw_argv = ["validate"] + raw_argv

    ap = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = ap.add_subparsers(dest="command", help="Subcommand (default: validate if omitted)")

    p_val = sub.add_parser("validate", help="Merged cohort TSV → JSON + ROC/calibration/DCA plots")
    p_val.add_argument("--cohort", type=Path, required=True)
    p_val.add_argument("--out-dir", type=Path, default=BASE / "results/phases_g/cohort_validation")
    p_val.add_argument("--sample-id-col", default="sample_id")
    p_val.add_argument("--outcome-col", default="y", help="Binary 0/1 outcome")
    p_val.add_argument("--prs-col", default="PRS_z")
    p_val.add_argument(
        "--baseline-cols",
        default="",
        help="Comma-separated PHASES/clinical columns (empty = intercept-only baseline)",
    )
    p_val.add_argument("--test-fraction", type=float, default=0.25)
    p_val.add_argument("--min-split-n", type=int, default=40)
    p_val.add_argument("--min-class-test", type=int, default=5)
    p_val.add_argument("--seed", type=int, default=42)
    p_val.add_argument("--calibration-bins", type=int, default=10)
    p_val.add_argument("--prefix", default="", help="Prefix for output filenames")

    p_cl = sub.add_parser("clump-plink2", help="Call PLINK2 --clump (needs bfile + assoc file)")
    p_cl.add_argument("--plink2-exec", default="plink2")
    p_cl.add_argument("--bfile", type=Path, required=True, help="PLINK 1 binary prefix (no .bed)")
    p_cl.add_argument("--clump-assoc", type=Path, required=True)
    p_cl.add_argument("--out-dir", type=Path, default=BASE / "results/phases_g/clump_plink2")
    p_cl.add_argument("--clump-p1", default="1")
    p_cl.add_argument("--clump-r2", default="0.1")
    p_cl.add_argument("--clump-kb", default="250")
    p_cl.add_argument("--clump-snp-field", default="", help="Optional PLINK2 field name for variant ID")
    p_cl.add_argument("--clump-p-field", default="", help="Optional PLINK2 field name for p-value")

    args = ap.parse_args(raw_argv)
    if args.command == "validate":
        run_validation(args)
        return
    if args.command == "clump-plink2":
        sys.exit(run_clump_plink2(args))


if __name__ == "__main__":
    main()
