"""
Initialize git, commit analysis pipeline code, and push to GitHub.

Usage:
  python publish_github_repo.py
  python publish_github_repo.py --repo-name aneurysm-rupture-genomics --dry-run
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO_NAME_DEFAULT = "aneurysm-rupture-genomics"
GITHUB_USER = "vsanker21"
DESCRIPTION = (
    "End-to-end rupture genomics analysis: GSE122897 tissue DE, 300-gene shortlist, "
    "blood-bridge 8-gene panel, external validation (analysis code only; no manuscript docx)."
)


def find_git() -> str:
    git = shutil.which("git")
    if git:
        return git
    candidates = [
        Path(r"C:\Program Files\Git\bin\git.exe"),
        Path(r"C:\Program Files (x86)\Git\bin\git.exe"),
    ]
    for c in candidates:
        if c.is_file():
            return str(c)
    raise RuntimeError("git not found; install Git for Windows")


def run(cmd: list[str], cwd: Path, *, check: bool = True) -> subprocess.CompletedProcess:
    print(">>>", " ".join(cmd), flush=True)
    return subprocess.run(cmd, cwd=str(cwd), check=check, text=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-name", default=REPO_NAME_DEFAULT)
    ap.add_argument("--private", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    root = Path(__file__).resolve().parent
    git = find_git()
    os.environ["PATH"] = str(Path(git).parent) + os.pathsep + os.environ.get("PATH", "")

    if not (root / "run_pipeline.py").is_file():
        raise SystemExit(f"run_pipeline.py not found in {root}")

    if args.dry_run:
        print(f"Would init git in {root}")
        print(f"Would create https://github.com/{GITHUB_USER}/{args.repo_name}")
        return 0

    if not (root / ".git").exists():
        run([git, "init", "-b", "main"], root)

    run([git, "config", "user.name", GITHUB_USER], root, check=False)
    run([git, "config", "user.email", f"{GITHUB_USER}@users.noreply.github.com"], root, check=False)

    run([git, "add", "-A"], root)
    status = run([git, "status", "--porcelain"], root)
    if not status.stdout.strip():
        print("Nothing to commit (working tree clean).")
    else:
        run(
            [
                git,
                "commit",
                "-m",
                "Add aneurysm rupture genomics analysis pipeline\n\n"
                "End-to-end run_pipeline.py workflow: tissue DESeq2, PHASES-G shortlist,\n"
                "blood-bridge panel, external validation, supplementary tables/figures.\n"
                "Excludes manuscript Word generation and large data/results (see .gitignore).",
            ],
            root,
        )

    remote = f"https://github.com/{GITHUB_USER}/{args.repo_name}.git"
    view = subprocess.run(
        ["gh", "repo", "view", f"{GITHUB_USER}/{args.repo_name}"],
        cwd=str(root),
        capture_output=True,
        text=True,
    )
    if view.returncode != 0:
        vis = "private" if args.private else "public"
        run(
            [
                "gh",
                "repo",
                "create",
                args.repo_name,
                f"--{vis}",
                "--description",
                DESCRIPTION,
                "--source=.",
                "--remote=origin",
            ],
            root,
        )
    else:
        remotes = run([git, "remote", "-v"], root)
        if "origin" not in remotes.stdout:
            run([git, "remote", "add", "origin", remote], root)

    run([git, "push", "-u", "origin", "main"], root)

    url = f"https://github.com/{GITHUB_USER}/{args.repo_name}"
    print(f"\nSuccess: {url}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
