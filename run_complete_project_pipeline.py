"""
Complete Project Pipeline Execution
Runs all phases of the project as outlined in Cerebral_Aneurysm_Genomics_Project.md
"""
import os
import sys
import subprocess
from pathlib import Path
import json

def print_section(title):
    """Print a formatted section header"""
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)

def run_script(script_name, description, required=True):
    """Run a Python script and handle errors"""
    print_section(f"RUNNING: {description}")
    print(f"Script: {script_name}")
    
    if not os.path.exists(script_name):
        if required:
            print(f"[ERROR] Required script not found: {script_name}")
            return False
        else:
            print(f"[SKIP] Optional script not found: {script_name}")
            return None
    
    try:
        result = subprocess.run(
            [sys.executable, script_name],
            capture_output=True,
            text=True,
            timeout=1800  # 30 minute timeout
        )
        
        if result.returncode == 0:
            print(f"[OK] {description} completed successfully")
            # Show last few lines of output
            if result.stdout:
                lines = result.stdout.strip().split('\n')
                if len(lines) > 0:
                    print("\nLast output:")
                    for line in lines[-5:]:
                        if line.strip():
                            print(f"  {line}")
            return True
        else:
            print(f"[WARNING] {description} had issues (exit code: {result.returncode})")
            if result.stderr:
                error_lines = result.stderr.strip().split('\n')
                print("Error output (last 5 lines):")
                for line in error_lines[-5:]:
                    if line.strip():
                        print(f"  {line}")
            return False
    except subprocess.TimeoutExpired:
        print(f"[ERROR] {description} timed out (>30 minutes)")
        return False
    except Exception as e:
        print(f"[ERROR] {description} failed: {e}")
        return False

def check_data_availability():
    """Check if required data files are available"""
    print_section("CHECKING DATA AVAILABILITY")
    
    required_files = {
        "GSE122897 TPM data": "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz",
        "GSE122897 metadata": "data/sc_rnaseq/GSE122897/sample_metadata.tsv",
        "GSE159610 data": "data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz",
        "GWAS data": "data/gwas/bakker_2020/bakker_2020_summary_stats.txt.gz"
    }
    
    available = []
    missing = []
    
    for name, path in required_files.items():
        if os.path.exists(path):
            print(f"[OK] {name}: Available")
            available.append(name)
        else:
            print(f"[MISSING] {name}: {path}")
            missing.append(name)
    
    return available, missing

def process_supplementary_data():
    """Process supplementary ZIP, TAR, and PDF files"""
    print_section("PROCESSING SUPPLEMENTARY DATA")
    
    # Check for supplementary files
    supp_zip = "data/supplementary/cvae186/cvae186_supplementary_data.zip"
    tar_files = [
        "data/sc_rnaseq/GSE186927/GSE186927_RAW.tar",
        "data/microarray/GSE13353/GSE13353_RAW.tar",
        "data/microarray/GSE54083/GSE54083_RAW.tar",
        "data/microarray/GSE75436/GSE75436_RAW.tar"
    ]
    
    print("\nChecking supplementary files...")
    if os.path.exists(supp_zip):
        print(f"[OK] ZIP file found: {supp_zip}")
        # Check if already extracted
        if os.path.exists("results/supplementary_data_check"):
            print("  Already extracted")
        else:
            print("  [NOTE] Extraction may be needed")
    
    for tar_file in tar_files:
        if os.path.exists(tar_file):
            print(f"[OK] TAR file found: {tar_file}")
        else:
            print(f"[MISSING] TAR file: {tar_file}")
    
    # Run integration script
    return run_script("integrate_all_new_datasets.py", "Supplementary Data Integration", required=False)

def main():
    """Main pipeline execution"""
    
    print("=" * 80)
    print("COMPLETE PROJECT PIPELINE EXECUTION")
    print("Following Cerebral_Aneurysm_Genomics_Project.md phases")
    print("=" * 80)
    
    # Create results directory
    os.makedirs("results", exist_ok=True)
    
    results = {}
    
    # Phase 1: Data Acquisition & Integration
    print_section("PHASE 1: DATA ACQUISITION & INTEGRATION")
    available, missing = check_data_availability()
    results['data_check'] = len(available), len(missing)
    
    # Process supplementary data
    results['supplementary'] = process_supplementary_data()
    
    # Phase 2: Polygenic Risk Score Development
    print_section("PHASE 2: POLYGENIC RISK SCORE DEVELOPMENT")
    results['prs'] = run_script(
        "build_polygenic_risk_score.py",
        "Build Polygenic Risk Score from GWAS",
        required=False
    )
    
    # Phase 3: Blood Transcriptomics Integration
    print_section("PHASE 3: BLOOD TRANSCRIPTOMICS INTEGRATION")
    results['blood_integration'] = run_script(
        "complete_blood_integration.py",
        "Complete Blood Transcriptomics Integration",
        required=False
    )
    
    # Also run blood-tissue integration
    results['blood_tissue'] = run_script(
        "integrate_blood_tissue_comprehensive.py",
        "Blood-Tissue Cross-Dataset Integration",
        required=False
    )
    
    # Phase 4: Endothelial Dysfunction Profiling
    print_section("PHASE 4: ENDOTHELIAL DYSFUNCTION PROFILING")
    results['ec_dysfunction'] = run_script(
        "endothelial_dysfunction_analysis.py",
        "Endothelial Dysfunction Analysis",
        required=False
    )
    
    # Phase 5: Multi-Omics Integration
    print_section("PHASE 5: MULTI-OMICS INTEGRATION")
    results['multi_omics'] = run_script(
        "multi_omics_integration.py",
        "Multi-Omics Integration (PRS + Blood + EC + Clinical)",
        required=False
    )
    
    # Run comprehensive integration
    results['comprehensive'] = run_script(
        "integrate_all_findings_enhanced.py",
        "Comprehensive Integration of All Findings",
        required=False
    )
    
    # Pathway Enrichment
    print_section("PATHWAY ENRICHMENT ANALYSIS")
    results['pathway'] = run_script(
        "pathway_enrichment_analysis.py",
        "Pathway Enrichment Analysis",
        required=False
    )
    
    # Summary
    print_section("EXECUTION SUMMARY")
    
    print("\nPhase Results:")
    phase_names = {
        'prs': 'Phase 2: PRS Development',
        'blood_integration': 'Phase 3: Blood Integration',
        'blood_tissue': 'Phase 3: Blood-Tissue Integration',
        'ec_dysfunction': 'Phase 4: EC Dysfunction',
        'multi_omics': 'Phase 5: Multi-Omics Integration',
        'comprehensive': 'Comprehensive Integration',
        'pathway': 'Pathway Enrichment'
    }
    
    for key, name in phase_names.items():
        status = results.get(key)
        if status is True:
            print(f"  ✅ {name}: SUCCESS")
        elif status is False:
            print(f"  ❌ {name}: FAILED")
        else:
            print(f"  ⏭️  {name}: SKIPPED/NOT AVAILABLE")
    
    # Check output files
    print_section("OUTPUT FILES CHECK")
    
    expected_files = [
        "results/gse122897_integrated_de.tsv",
        "results/prs_weights.tsv",
        "results/blood_tissue_integration_fdr_0.1.tsv",
        "results/multi_omics_prioritized_genes_literature_annotated.tsv",
        "results/pathway_enrichment_fdr_0.1.tsv"
    ]
    
    for file_path in expected_files:
        if os.path.exists(file_path):
            size = os.path.getsize(file_path) / 1024  # KB
            print(f"  ✅ {file_path} ({size:.1f} KB)")
        else:
            print(f"  ⚠️  {file_path} (not found)")
    
    print("\n" + "=" * 80)
    print("PIPELINE EXECUTION COMPLETE")
    print("=" * 80)
    print("\nNext steps:")
    print("1. Review results in results/ directory")
    print("2. Check for any missing output files")
    print("3. Run additional analyses as needed")
    print("4. Generate final summary reports")

if __name__ == "__main__":
    main()






