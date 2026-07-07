"""
Identify Data Sources for Proper Multi-Omics Integration

Requirements for true multi-omics integration:
1. Same patients with multiple data types (blood + tissue + genotype)
2. Same phenotype (rupture prediction)
3. Sufficient sample size for validation
"""
import sys
import os

# Force UTF-8 encoding
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

import pandas as pd
import numpy as np

print("=" * 80)
print("MULTI-OMICS DATA SOURCES IDENTIFICATION")
print("=" * 80)

# ============================================================================
# CURRENT DATASETS ANALYSIS
# ============================================================================

def analyze_current_datasets():
    """Analyze currently available datasets"""
    print("\n" + "=" * 80)
    print("1. CURRENT DATASETS ANALYSIS")
    print("=" * 80)
    
    datasets = []
    
    # Blood dataset (GSE159610)
    print("\n  Blood Dataset (GSE159610):")
    print(f"    Type: Blood RNA-seq")
    print(f"    Task: Aneurysm vs Control")
    print(f"    Samples: ~47 (23 aneurysm, 24 control)")
    print(f"    Data: ✓ Available")
    print(f"    Genotype: ✗ Not available")
    print(f"    Tissue: ✗ Not available")
    
    datasets.append({
        'Name': 'GSE159610',
        'Type': 'Blood RNA-seq',
        'Task': 'Aneurysm vs Control',
        'Samples': 47,
        'Has_Blood': True,
        'Has_Tissue': False,
        'Has_Genotype': False,
        'Can_Integrate': False,
        'Reason': 'No tissue or genotype data'
    })
    
    # Tissue dataset (GSE122897)
    print("\n  Tissue Dataset (GSE122897):")
    print(f"    Type: Tissue RNA-seq")
    print(f"    Task: Ruptured vs Unruptured")
    print(f"    Samples: ~50 (ruptured vs unruptured)")
    print(f"    Data: ✓ Available")
    print(f"    Genotype: ✗ Not available")
    print(f"    Blood: ✗ Not available")
    
    datasets.append({
        'Name': 'GSE122897',
        'Type': 'Tissue RNA-seq',
        'Task': 'Ruptured vs Unruptured',
        'Samples': 50,
        'Has_Blood': False,
        'Has_Tissue': True,
        'Has_Genotype': False,
        'Can_Integrate': False,
        'Reason': 'No blood or genotype data'
    })
    
    # Tissue dataset (GSE66238)
    print("\n  Tissue Dataset (GSE66238):")
    print(f"    Type: Tissue RNA-seq")
    print(f"    Task: Aneurysm vs Control (superficial temporal artery)")
    print(f"    Samples: Unknown")
    print(f"    Data: ✓ Available")
    print(f"    Genotype: ✗ Not available")
    print(f"    Blood: ✗ Not available")
    
    datasets.append({
        'Name': 'GSE66238',
        'Type': 'Tissue RNA-seq',
        'Task': 'Aneurysm vs Control',
        'Samples': 'Unknown',
        'Has_Blood': False,
        'Has_Tissue': True,
        'Has_Genotype': False,
        'Can_Integrate': False,
        'Reason': 'No blood or genotype data'
    })
    
    print("\n  Summary:")
    print(f"    Current datasets: Different patients, different phenotypes")
    print(f"    Cannot directly integrate without same-patient data")
    
    return datasets

# ============================================================================
# POTENTIAL DATA SOURCES FOR MULTI-OMICS INTEGRATION
# ============================================================================

def identify_potential_sources():
    """Identify potential data sources for multi-omics integration"""
    print("\n" + "=" * 80)
    print("2. POTENTIAL DATA SOURCES FOR MULTI-OMICS INTEGRATION")
    print("=" * 80)
    
    sources = []
    
    # 1. dbGaP - Controlled Access
    print("\n  1. dbGaP (Database of Genotypes and Phenotypes)")
    print("     URL: https://www.ncbi.nlm.nih.gov/gap/")
    print("     Access: Controlled (requires institutional approval)")
    print("     Potential for:")
    print("       - Genotype data (SNP arrays, WGS)")
    print("       - Paired samples (blood + tissue)")
    print("       - Clinical metadata")
    print("     Requirements:")
    print("       - PI account")
    print("       - IRB approval")
    print("       - Data Use Agreement")
    print("     Status: Requires application")
    
    sources.append({
        'Source': 'dbGaP',
        'Type': 'Controlled Access',
        'Provides': ['Genotype', 'Blood', 'Tissue', 'Clinical'],
        'Access': 'Controlled (PI + IRB)',
        'URL': 'https://www.ncbi.nlm.nih.gov/gap/',
        'Status': 'Requires Application'
    })
    
    # 2. E-GEOD-15629 (from web search)
    print("\n  2. E-GEOD-15629 (ArrayExpress)")
    print("     Type: Gene Expression (Microarray)")
    print("     Task: Ruptured vs Unruptured vs Control")
    print("     Samples: 8 ruptured + 6 unruptured + 5 controls")
    print("     Provides: Tissue expression")
    print("     Limitations: Microarray (not RNA-seq), no blood/genotype")
    print("     URL: https://www.omicsdi.org/dataset/biostudies-arrayexpress/E-GEOD-15629")
    print("     Status: Available but limited for multi-omics")
    
    sources.append({
        'Source': 'E-GEOD-15629',
        'Type': 'Public (ArrayExpress)',
        'Provides': ['Tissue'],
        'Access': 'Open',
        'URL': 'https://www.omicsdi.org/dataset/biostudies-arrayexpress/E-GEOD-15629',
        'Status': 'Available (tissue only)',
        'Note': 'Microarray, not RNA-seq'
    })
    
    # 3. Contact Authors
    print("\n  3. Contact Study Authors")
    print("     Potential sources:")
    print("       - Poppenberg et al. (GSE159610): May have additional samples")
    print("       - Ruigrok et al. (GSE122897): May have blood samples")
    print("       - Bakker et al. (GWAS): May have genotype + phenotype data")
    print("     Benefits:")
    print("       - Paired samples (blood + tissue)")
    print("       - Genotype data")
    print("       - Clinical metadata")
    print("     Status: Requires collaboration")
    
    sources.append({
        'Source': 'Author Collaboration',
        'Type': 'Collaboration',
        'Provides': ['Blood', 'Tissue', 'Genotype', 'Clinical'],
        'Access': 'Collaboration Required',
        'Contacts': ['Poppenberg (GSE159610)', 'Ruigrok (GSE122897)', 'Bakker (GWAS)'],
        'Status': 'Requires Contact'
    })
    
    # 4. UK Biobank (if available)
    print("\n  4. UK Biobank")
    print("     Type: Population-based cohort")
    print("     Provides:")
    print("       - Genotype data (if approved)")
    print("       - Clinical data")
    print("       - May have blood samples")
    print("     Limitations:")
    print("       - Requires approval")
    print("       - May not have aneurysm-specific samples")
    print("       - Tissue data unlikely")
    print("     Status: Requires application")
    
    sources.append({
        'Source': 'UK Biobank',
        'Type': 'Controlled Access',
        'Provides': ['Genotype', 'Clinical', 'Blood (possible)'],
        'Access': 'Controlled (approved access)',
        'URL': 'https://www.ukbiobank.ac.uk/',
        'Status': 'Requires Application',
        'Note': 'May not have aneurysm-specific data'
    })
    
    # 5. Prospective Study
    print("\n  5. Prospective Clinical Study")
    print("     Type: Collect new data")
    print("     Requirements:")
    print("       - IRB approval")
    print("       - Patient recruitment")
    print("       - Sample collection (blood + tissue)")
    print("       - Genotyping (SNP array or WGS)")
    print("     Timeline: 2-5 years")
    print("     Benefits:")
    print("       - Paired samples")
    print("       - Complete multi-omics data")
    print("       - Controlled collection")
    print("     Status: Future option")
    
    sources.append({
        'Source': 'Prospective Study',
        'Type': 'New Data Collection',
        'Provides': ['Blood', 'Tissue', 'Genotype', 'Clinical'],
        'Access': 'IRB + Recruitment',
        'Timeline': '2-5 years',
        'Status': 'Future Option'
    })
    
    return sources

# ============================================================================
# RECOMMENDED APPROACH
# ============================================================================

def recommend_approach():
    """Recommend approach for multi-omics integration"""
    print("\n" + "=" * 80)
    print("3. RECOMMENDED APPROACH FOR MULTI-OMICS INTEGRATION")
    print("=" * 80)
    
    print("\n  IMMEDIATE (Use Current Data):")
    print("    1. Separate Models:")
    print("       a. Blood Model: Aneurysm vs Control (AUC 0.968) ✓")
    print("          - Task: Screening/detection")
    print("       b. Tissue Model: Ruptured vs Unruptured (build now)")
    print("          - Task: Risk stratification")
    print("    2. Feature-Level Integration:")
    print("       - Identify common signature genes (blood + tissue)")
    print("       - Cross-validate features across datasets")
    print("       - Use for biomarker panel development")
    
    print("\n  SHORT-TERM (1-6 months):")
    print("    1. Apply for dbGaP Access:")
    print("       - Search for aneurysm studies with genotype data")
    print("       - Request access for genotype + phenotype")
    print("       - May provide paired blood/tissue samples")
    print("    2. Contact Authors:")
    print("       - Poppenberg (blood): May have additional samples")
    print("       - Ruigrok (tissue): May have blood samples from same patients")
    print("       - Bakker (GWAS): May have genotype + phenotype data")
    print("    3. Integrate Available Data:")
    print("       - If genotype data obtained: Add PRS to blood model")
    print("       - If paired samples obtained: True multi-omics integration")
    
    print("\n  LONG-TERM (1-2 years):")
    print("    1. Prospective Study:")
    print("       - Collect blood + tissue from same patients")
    print("       - Genotype all patients (SNP array)")
    print("       - Follow for rupture outcomes")
    print("       - Complete multi-omics integration")
    print("    2. Multi-Center Collaboration:")
    print("       - Pool resources across institutions")
    print("       - Share samples and data")
    print("       - Larger sample size for validation")
    
    print("\n  INTERIM SOLUTION (Current Capability):")
    print("    1. Use Blood Model for Screening:")
    print("       - AUC 0.968 for aneurysm detection")
    print("       - Clinical utility for screening")
    print("    2. Use Tissue Model for Risk Stratification:")
    print("       - Separate model for ruptured vs unruptured")
    print("       - Apply after aneurysm detection")
    print("    3. Sequential Approach:")
    print("       - Step 1: Blood test (aneurysm present?)")
    print("       - Step 2: If positive, tissue analysis (rupture risk?)")
    print("       - Step 3: If high risk, intervention")
    
    print("\n  Key Insight:")
    print("    True multi-omics integration requires:")
    print("    ✓ Same patients across all omics layers")
    print("    ✓ Same phenotype (rupture prediction)")
    print("    ✓ Sufficient sample size (n > 100)")
    print("    ✗ Current datasets: Different patients, different phenotypes")

# ============================================================================
# DATA SOURCES CHECKLIST
# ============================================================================

def create_data_sources_checklist():
    """Create checklist for data sources"""
    print("\n" + "=" * 80)
    print("4. DATA SOURCES CHECKLIST")
    print("=" * 80)
    
    checklist = pd.DataFrame([
        {
            'Source': 'dbGaP',
            'Check': 'Search for aneurysm studies',
            'Action': 'Apply for access',
            'Provides': 'Genotype + phenotype',
            'Timeline': '2-4 weeks approval'
        },
        {
            'Source': 'Poppenberg et al.',
            'Check': 'Contact for additional samples',
            'Action': 'Request paired tissue samples',
            'Provides': 'Blood + tissue (if available)',
            'Timeline': 'Variable'
        },
        {
            'Source': 'Ruigrok et al.',
            'Check': 'Contact for blood samples',
            'Action': 'Request paired blood samples',
            'Provides': 'Blood + tissue (if available)',
            'Timeline': 'Variable'
        },
        {
            'Source': 'Bakker et al.',
            'Check': 'Check for genotype data',
            'Action': 'Request genotype + phenotype',
            'Provides': 'Genotype + clinical',
            'Timeline': 'Variable'
        },
        {
            'Source': 'E-GEOD-15629',
            'Check': 'Download microarray data',
            'Action': 'Use for tissue validation',
            'Provides': 'Tissue expression',
            'Timeline': 'Immediate'
        },
        {
            'Source': 'Prospective Study',
            'Check': 'IRB approval',
            'Action': 'Collect new samples',
            'Provides': 'Complete multi-omics',
            'Timeline': '2-5 years'
        }
    ])
    
    print("\n" + checklist.to_string(index=False))
    
    # Save checklist
    checklist_file = "results/multi_omics_data_sources_checklist.tsv"
    checklist.to_csv(checklist_file, sep='\t', index=False)
    print(f"\n  [OK] Saved checklist to {checklist_file}")
    
    return checklist

# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main pipeline"""
    
    os.makedirs("results", exist_ok=True)
    
    # Analyze current datasets
    current_datasets = analyze_current_datasets()
    
    # Identify potential sources
    potential_sources = identify_potential_sources()
    
    # Recommend approach
    recommend_approach()
    
    # Create checklist
    checklist = create_data_sources_checklist()
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    print("\n  Current Status:")
    print(f"    - Blood model: ✓ Available (AUC 0.968)")
    print(f"    - Tissue model: ⚠ Being built")
    print(f"    - Multi-omics: ✗ Not possible with current data")
    print(f"      Reason: Different patients, different phenotypes")
    
    print("\n  Next Steps:")
    print(f"    1. Build tissue model (ruptured vs unruptured) - IN PROGRESS")
    print(f"    2. Apply for dbGaP access (genotype data)")
    print(f"    3. Contact authors (paired samples)")
    print(f"    4. Consider prospective study (long-term)")
    
    print("\n  Output Files:")
    print(f"    - results/multi_omics_data_sources_checklist.tsv")
    
    # Save summary
    summary_file = "results/multi_omics_data_sources_summary.md"
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write("# Multi-Omics Data Sources Summary\n\n")
        f.write("## Current Status\n\n")
        f.write("- Blood model: Available (AUC 0.968)\n")
        f.write("- Tissue model: Being built\n")
        f.write("- Multi-omics: Not possible with current data\n\n")
        f.write("## Potential Sources\n\n")
        for source in potential_sources:
            f.write(f"### {source['Source']}\n")
            f.write(f"- Type: {source['Type']}\n")
            f.write(f"- Provides: {', '.join(source['Provides'])}\n")
            f.write(f"- Access: {source['Access']}\n")
            f.write(f"- Status: {source['Status']}\n\n")
    
    print(f"    - {summary_file}")
    
    return current_datasets, potential_sources, checklist

if __name__ == "__main__":
    main()

