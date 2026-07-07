"""
Build a combined external-validation summary table and figure for the manuscript.

The output intentionally preserves the evidentiary hierarchy:
1. direct rupture support (GSE54083, GSE13353, GSE15629 rupture subset)
2. orthogonal tissue context (GSE26969, GSE15629 aneurysm-vs-control)
3. blood bridge cohorts for reduced-panel derivation (GSE159610, GSE106520)
4. post-rupture / SAH blood context only (GSE36791, GSE73378)
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from statsmodels.stats.multitest import multipletests

try:
    from project_paths import get_project_root
    BASE = get_project_root()
except ImportError:
    BASE = Path(__file__).resolve().parent
OUT_TSV = BASE / "results" / "combined_validation_summary.tsv"
OUT_FIG = BASE / "results" / "figures" / "ultra_hd" / "combined_validation_summary_ultra_hd.png"


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _bh_sig_count(pvals: pd.Series) -> int:
    vals = pd.to_numeric(pvals, errors="coerce").dropna()
    if vals.empty:
        return 0
    return int((multipletests(vals.values, method="fdr_bh")[1] < 0.05).sum())


def _shortlist_sig_counts(path: Path) -> tuple[int, int]:
    df = pd.read_csv(path, sep="\t")
    nominal = int((pd.to_numeric(df["P_value"], errors="coerce") < 0.05).sum())
    fdr = int((pd.to_numeric(df["FDR_BH"], errors="coerce") < 0.05).sum())
    return nominal, fdr


def build_combined_validation_dataframe(base: Path = BASE) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    g54083 = _load_json(base / "results" / "external_validation_gse54083" / "gse54083_validation_summary.json")
    rows.append(
        {
            "Order": 1,
            "Dataset": "GSE54083",
            "Role": "Direct rupture replication",
            "Context_group": "Direct rupture support",
            "Comparison": "ruptured aneurysm vs unruptured aneurysm",
            "Samples": f"{g54083['n_ruptured']} R / {g54083['n_unruptured']} U",
            "Measured_or_tested_genes": int(g54083["n_shortlist_genes_present"]),
            "Direction_matched_genes": int(g54083["n_shortlist_direction_match"]),
            "Direction_match_pct": 100.0 * float(g54083["n_shortlist_direction_match"]) / max(1, float(g54083["n_shortlist_genes_present"])),
            "Nominal_P_lt_0_05": int(g54083["n_shortlist_p_lt_0_05"]),
            "FDR_lt_0_05": int(g54083["n_shortlist_fdr_lt_0_05"]),
            "Interpretation": "Best matched public rupture replication, but very small and no FDR survivors.",
        }
    )

    g13353_summary = pd.read_csv(
        base / "results" / "external_validation_gse13353" / "gse13353_validation_summary.tsv",
        sep="\t",
        header=None,
        names=["key", "value"],
    )
    g13353_map = dict(zip(g13353_summary["key"], g13353_summary["value"]))
    g13353_tbl = pd.read_csv(
        base / "results" / "external_validation_gse13353" / "gse13353_validation_table_with_gse122897.tsv",
        sep="\t",
    )
    g13353_nominal = int((pd.to_numeric(g13353_tbl["P_value"], errors="coerce") < 0.05).sum())
    g13353_fdr = _bh_sig_count(g13353_tbl["P_value"])
    g13353_tested = int(float(g13353_map["genes_tested"]))
    g13353_match = int(float(g13353_map["direction_matches"]))
    rows.append(
        {
            "Order": 2,
            "Dataset": "GSE13353",
            "Role": "Targeted rupture support",
            "Context_group": "Direct rupture support",
            "Comparison": "ruptured aneurysm vs unruptured aneurysm",
            "Samples": f"{int(float(g13353_map['n_ruptured']))} R / {int(float(g13353_map['n_unruptured']))} U",
            "Measured_or_tested_genes": g13353_tested,
            "Direction_matched_genes": g13353_match,
            "Direction_match_pct": 100.0 * g13353_match / max(1, g13353_tested),
            "Nominal_P_lt_0_05": g13353_nominal,
            "FDR_lt_0_05": g13353_fdr,
            "Interpretation": "Targeted seven-gene support only; useful but narrower than genome-wide replication.",
        }
    )

    g15629 = _load_json(base / "results" / "external_validation_gse15629" / "gse15629_context_summary.json")
    g15629_rupture_nominal, g15629_rupture_fdr = _shortlist_sig_counts(
        base / "results" / "external_validation_gse15629" / "gse15629_shortlist_overlap_rupture.tsv"
    )
    rows.append(
        {
            "Order": 3,
            "Dataset": "GSE15629 (rupture subset)",
            "Role": "Secondary direct rupture support",
            "Context_group": "Direct rupture support",
            "Comparison": "ruptured aneurysm vs unruptured aneurysm (subset of mixed-design study)",
            "Samples": f"{g15629['n_ruptured']} R / {g15629['n_unruptured']} U",
            "Measured_or_tested_genes": int(g15629["n_shortlist_genes_present_rupture_subset"]),
            "Direction_matched_genes": int(g15629["n_shortlist_direction_match_rupture_subset"]),
            "Direction_match_pct": 100.0
            * float(g15629["n_shortlist_direction_match_rupture_subset"])
            / max(1, float(g15629["n_shortlist_genes_present_rupture_subset"])),
            "Nominal_P_lt_0_05": g15629_rupture_nominal,
            "FDR_lt_0_05": g15629_rupture_fdr,
            "Interpretation": "Useful secondary rupture subset, but embedded within a mixed ruptured, unruptured, and control tissue design.",
        }
    )

    g26969 = _load_json(base / "results" / "external_validation_gse26969" / "gse26969_context_summary.json")
    rows.append(
        {
            "Order": 4,
            "Dataset": "GSE26969",
            "Role": "Orthogonal tissue context",
            "Context_group": "Tissue context",
            "Comparison": "unruptured aneurysm wall vs normal superficial temporal artery",
            "Samples": f"{g26969['n_aneurysm']} A / {g26969['n_control']} C",
            "Measured_or_tested_genes": int(g26969["n_shortlist_genes_present"]),
            "Direction_matched_genes": int(g26969["n_shortlist_direction_match"]),
            "Direction_match_pct": 100.0 * float(g26969["n_shortlist_direction_match"]) / max(1, float(g26969["n_shortlist_genes_present"])),
            "Nominal_P_lt_0_05": int(g26969["n_shortlist_p_lt_0_05"]),
            "FDR_lt_0_05": int(g26969["n_shortlist_fdr_lt_0_05"]),
            "Interpretation": "Shows aneurysm-versus-normal wall context, not rupture replication.",
        }
    )

    g15629_avc_nominal, g15629_avc_fdr = _shortlist_sig_counts(
        base / "results" / "external_validation_gse15629" / "gse15629_shortlist_overlap_aneurysm_vs_control.tsv"
    )
    rows.append(
        {
            "Order": 5,
            "Dataset": "GSE15629 (aneurysm vs control)",
            "Role": "Mixed-design tissue context",
            "Context_group": "Tissue context",
            "Comparison": "aneurysm wall (ruptured or unruptured) vs control vessel",
            "Samples": f"{int(g15629['n_ruptured'] + g15629['n_unruptured'])} A / {g15629['n_control']} C",
            "Measured_or_tested_genes": int(g15629["n_shortlist_genes_present_aneurysm_vs_control"]),
            "Direction_matched_genes": int(g15629["n_shortlist_direction_match_aneurysm_vs_control"]),
            "Direction_match_pct": 100.0
            * float(g15629["n_shortlist_direction_match_aneurysm_vs_control"])
            / max(1, float(g15629["n_shortlist_genes_present_aneurysm_vs_control"])),
            "Nominal_P_lt_0_05": g15629_avc_nominal,
            "FDR_lt_0_05": g15629_avc_fdr,
            "Interpretation": "Provides aneurysm-versus-vessel context and is not interpreted as direct rupture replication.",
        }
    )

    g36791 = _load_json(base / "results" / "external_validation_gse36791" / "gse36791_context_summary.json")
    rows.append(
        {
            "Order": 8,
            "Dataset": "GSE36791",
            "Role": "Acute post-rupture blood context",
            "Context_group": "Blood context only",
            "Comparison": "ruptured aneurysm with SAH vs headache controls",
            "Samples": f"{g36791['n_ruptured']} R / {g36791['n_control']} C",
            "Measured_or_tested_genes": int(g36791["n_shortlist_genes_present"]),
            "Direction_matched_genes": int(g36791["n_shortlist_direction_match"]),
            "Direction_match_pct": 100.0 * float(g36791["n_shortlist_direction_match"]) / max(1, float(g36791["n_shortlist_genes_present"])),
            "Nominal_P_lt_0_05": int(g36791["n_shortlist_p_lt_0_05"]),
            "FDR_lt_0_05": int(g36791["n_shortlist_fdr_lt_0_05"]),
            "Interpretation": "Acute post-rupture / SAH blood response only; not used to derive the reduced blood panel.",
        }
    )

    g159610_meta = pd.read_csv(
        base / "data" / "blood_rnaseq" / "GSE159610" / "sample_metadata_with_labels.tsv",
        sep="\t",
    )
    g159610_shortlist = pd.read_csv(base / "results" / "phases_g" / "shortlist_evidence.tsv", sep="\t")
    g159610_has = g159610_shortlist["GSE159610_Log2FC_aneurysm_vs_control"].notna()
    g159610_measured = int(g159610_has.sum())
    g159610_match = int((g159610_has & (g159610_shortlist["direction_match_rupture_and_aneurysm"] == True)).sum())  # noqa: E712
    g159610_nominal = int((g159610_has & (pd.to_numeric(g159610_shortlist["GSE159610_P"], errors="coerce") < 0.05)).sum())
    g159610_fdr = _bh_sig_count(g159610_shortlist.loc[g159610_has, "GSE159610_P"])
    n_aneurysm = int((g159610_meta["Group"] == "Aneurysm").sum())
    n_control = int((g159610_meta["Group"] == "Control").sum())
    rows.append(
        {
            "Order": 6,
            "Dataset": "GSE159610",
            "Role": "Primary blood bridge",
            "Context_group": "Blood bridge",
            "Comparison": "aneurysm vs control blood RNA-seq",
            "Samples": f"{n_aneurysm} A / {n_control} C",
            "Measured_or_tested_genes": g159610_measured,
            "Direction_matched_genes": g159610_match,
            "Direction_match_pct": 100.0 * g159610_match / max(1, g159610_measured),
            "Nominal_P_lt_0_05": g159610_nominal,
            "FDR_lt_0_05": g159610_fdr,
            "Interpretation": "Primary blood bridge cohort used to derive the reduced blood-testable panel; still an aneurysm-presence rather than rupture comparison.",
        }
    )

    g106520 = _load_json(base / "results" / "external_validation_gse106520" / "gse106520_context_summary.json")
    rows.append(
        {
            "Order": 7,
            "Dataset": "GSE106520",
            "Role": "Neutrophil blood bridge",
            "Context_group": "Blood bridge",
            "Comparison": "unruptured aneurysm vs control neutrophil RNA-seq",
            "Samples": f"{g106520['n_aneurysm']} A / {g106520['n_control']} C",
            "Measured_or_tested_genes": int(g106520["n_shortlist_genes_present"]),
            "Direction_matched_genes": int(g106520["n_shortlist_direction_match"]),
            "Direction_match_pct": 100.0 * float(g106520["n_shortlist_direction_match"]) / max(1, float(g106520["n_shortlist_genes_present"])),
            "Nominal_P_lt_0_05": int(g106520["n_shortlist_p_lt_0_05"]),
            "FDR_lt_0_05": int(g106520["n_shortlist_fdr_lt_0_05"]),
            "Interpretation": "Independent neutrophil bridge cohort used with GSE159610 to derive the reduced blood-testable panel.",
        }
    )

    g73378 = _load_json(base / "results" / "external_validation_gse73378" / "gse73378_context_summary.json")
    rows.append(
        {
            "Order": 9,
            "Dataset": "GSE73378",
            "Role": "Long-horizon post-SAH blood context",
            "Context_group": "Blood context only",
            "Comparison": "past aneurysmal SAH history vs control whole blood",
            "Samples": f"{g73378['n_case']} past SAH / {g73378['n_control']} C",
            "Measured_or_tested_genes": int(g73378["n_shortlist_genes_present"]),
            "Direction_matched_genes": int(g73378["n_shortlist_direction_match"]),
            "Direction_match_pct": 100.0 * float(g73378["n_shortlist_direction_match"]) / max(1, float(g73378["n_shortlist_genes_present"])),
            "Nominal_P_lt_0_05": int(g73378["n_shortlist_p_lt_0_05"]),
            "FDR_lt_0_05": int(g73378["n_shortlist_fdr_lt_0_05"]),
            "Interpretation": "Long-horizon blood context after aneurysmal SAH history; retained as context only rather than as a bridge or rupture replicate.",
        }
    )

    df = pd.DataFrame(rows).sort_values("Order").reset_index(drop=True)
    df["Direction_match_pct"] = df["Direction_match_pct"].round(1)
    return df


def save_combined_validation_figure(df: pd.DataFrame, out_path: Path = OUT_FIG) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    color_map = {
        "Direct rupture support": "#c44e52",
        "Tissue context": "#4c72b0",
        "Blood bridge": "#8172b2",
        "Blood context only": "#55a868",
    }
    bar_colors = [color_map.get(v, "#808080") for v in df["Context_group"]]
    y = list(range(len(df)))[::-1]
    plot_df = df.iloc[::-1].reset_index(drop=True)
    plot_colors = [color_map.get(v, "#808080") for v in plot_df["Context_group"]]
    labels = [f"{d}\n{s}" for d, s in zip(plot_df["Dataset"], plot_df["Samples"])]

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), gridspec_kw={"width_ratios": [1.1, 1.2]})

    axes[0].barh(y, plot_df["Direction_match_pct"], color=plot_colors, alpha=0.9)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(labels, fontsize=10)
    axes[0].set_xlim(0, 100)
    axes[0].set_xlabel("Direction-matched genes (%)")
    axes[0].set_title("A. Concordance with GSE122897 shortlist effects", loc="left", fontsize=12, fontweight="bold")
    axes[0].grid(axis="x", alpha=0.25)
    for yi, pct, matched, measured in zip(
        y,
        plot_df["Direction_match_pct"],
        plot_df["Direction_matched_genes"],
        plot_df["Measured_or_tested_genes"],
    ):
        axes[0].text(min(pct + 1.5, 98), yi, f"{int(matched)}/{int(measured)}", va="center", fontsize=9)

    axes[1].barh(y, plot_df["Measured_or_tested_genes"], color="none", edgecolor="#c7c7c7", linewidth=1.2)
    nominal_only = plot_df["Nominal_P_lt_0_05"] - plot_df["FDR_lt_0_05"]
    nominal_only = nominal_only.clip(lower=0)
    axes[1].barh(y, nominal_only, color="#bdbdbd", alpha=0.9, label="Nominal only")
    axes[1].barh(y, plot_df["FDR_lt_0_05"], color=plot_colors, alpha=0.95, label="FDR < 0.05")
    axes[1].set_yticks(y)
    axes[1].set_yticklabels([])
    axes[1].set_xlabel("Genes among measured/tested shortlist set")
    axes[1].set_title("B. Significant genes within each cohort-specific context", loc="left", fontsize=12, fontweight="bold")
    axes[1].grid(axis="x", alpha=0.25)
    for yi, nominal, fdr, measured in zip(
        y,
        plot_df["Nominal_P_lt_0_05"],
        plot_df["FDR_lt_0_05"],
        plot_df["Measured_or_tested_genes"],
    ):
        axes[1].text(min(measured + 3, max(plot_df["Measured_or_tested_genes"]) + 10), yi, f"{int(nominal)} nominal / {int(fdr)} FDR", va="center", fontsize=9)

    from matplotlib.patches import Patch

    legend_handles = [
        Patch(facecolor=color_map["Direct rupture support"], label="Direct rupture support"),
        Patch(facecolor=color_map["Tissue context"], label="Orthogonal tissue context"),
        Patch(facecolor=color_map["Blood bridge"], label="Blood bridge cohorts"),
        Patch(facecolor=color_map["Blood context only"], label="Post-rupture / SAH blood context"),
        Patch(facecolor="#bdbdbd", label="Nominal-only hits"),
    ]
    axes[1].legend(handles=legend_handles, frameon=False, loc="lower right", fontsize=9)

    fig.suptitle("Combined external-validation hierarchy across public datasets", fontsize=14, fontweight="bold")
    fig.text(
        0.5,
        0.01,
        "R = ruptured, U = unruptured, A = aneurysm, C = control. GSE159610 and GSE106520 act as blood bridges, whereas GSE36791 and GSE73378 remain post-rupture / SAH context only.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def build_combined_validation_assets(base: Path = BASE) -> pd.DataFrame:
    df = build_combined_validation_dataframe(base)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    save_combined_validation_figure(df, OUT_FIG)
    return df


if __name__ == "__main__":
    out = build_combined_validation_assets()
    print(out.to_string(index=False))
