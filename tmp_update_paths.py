from pathlib import Path

FILES = [
    "comprehensive_analysis_pipeline.py",
    "comprehensive_multi_omics_integration.py",
    "execute_comprehensive_integration.py",
    "multi_omics_integration_pipeline.py",
    "phase2_population_prs.py",
    "phase3_validate_blood_signature.py",
    "complete_validation_prs_integration.py",
    "verify_performance_validity.py",
    "complete_validation_pathway_analysis.py",
    "comprehensive_validation_and_analysis.py",
    "data/complete_dataset_inventory.json",
    "setup_deseq2_analysis.py",
    "ml_classification_blood_genes.py",
    "deseq2_blood_analysis.py",
    "run_full_blood_de_analysis.py",
    "complete_blood_integration.py",
    "rerun_blood_integration_refined.py",
    "IMMEDIATE_ACTIONS_COMPLETE.md",
    "IMMEDIATE_ACTIONS_SUMMARY.md",
    "integrate_blood_transcriptomics.py",
]

OLD = "sample_metadata_refined.tsv"
NEW = "sample_metadata_with_labels.tsv"


def main() -> None:
    for rel in FILES:
        path = Path(rel)
        text = path.read_text(encoding="utf-8")
        if OLD not in text:
            print(f"[WARN] No match in {rel}")
            continue
        path.write_text(text.replace(OLD, NEW), encoding="utf-8")
        print(f"[OK] Updated {rel}")


if __name__ == "__main__":
    main()

