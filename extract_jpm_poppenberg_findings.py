"""
Extract and Integrate Findings from Poppenberg et al. 2023 JPM Paper
Extracts: 39 gene list, PAT score formula, GEO accession info, and integrates with existing analyses
"""
import pandas as pd
import numpy as np
import re
from pathlib import Path
import json

def extract_39_genes():
    """Extract the 39 differentially expressed genes from the paper"""
    
    # Gene list from Table 2 in the paper
    genes_data = [
        # Genes with decreased expression in "growing" (negative Log2FC)
        {'Gene_ID': 'ENSG00000115155.16', 'Gene_Symbol': 'OTOF', 'Log2FC': -1.52, 'q_value': 1.00e-10, 'Direction': 'Decreased'},
        {'Gene_ID': 'ENSG00000206047.2', 'Gene_Symbol': 'DEFA1', 'Log2FC': -1.30, 'q_value': 2.06e-3, 'Direction': 'Decreased'},
        {'Gene_ID': 'ENSG00000206127.10', 'Gene_Symbol': 'GOLGA8O', 'Log2FC': -1.22, 'q_value': 2.19e-3, 'Direction': 'Decreased'},
        {'Gene_ID': 'ENSG00000175445.14', 'Gene_Symbol': 'LPL', 'Log2FC': -1.14, 'q_value': 7.92e-5, 'Direction': 'Decreased'},
        {'Gene_ID': 'ENSG00000076864.19', 'Gene_Symbol': 'RAP1GAP', 'Log2FC': -0.96, 'q_value': 2.59e-3, 'Direction': 'Decreased'},
        {'Gene_ID': 'ENSG00000137959.15', 'Gene_Symbol': 'IFI44L', 'Log2FC': -0.93, 'q_value': 7.39e-9, 'Direction': 'Decreased'},
        {'Gene_ID': 'ENSG00000165949.12', 'Gene_Symbol': 'IFI27', 'Log2FC': -0.87, 'q_value': 6.56e-5, 'Direction': 'Decreased'},
        {'Gene_ID': 'ENSG00000224940.8', 'Gene_Symbol': 'PRRT4', 'Log2FC': -0.86, 'q_value': 6.32e-4, 'Direction': 'Decreased'},
        {'Gene_ID': 'ENSG00000147689.16', 'Gene_Symbol': 'FAM83A', 'Log2FC': -0.75, 'q_value': 2.76e-3, 'Direction': 'Decreased'},
        {'Gene_ID': 'ENSG00000121933.17', 'Gene_Symbol': 'ADORA3', 'Log2FC': -0.74, 'q_value': 4.21e-2, 'Direction': 'Decreased'},
        {'Gene_ID': 'ENSG00000100320.22', 'Gene_Symbol': 'RBFOX2', 'Log2FC': -0.59, 'q_value': 7.27e-3, 'Direction': 'Decreased'},
        
        # Genes with increased expression in "growing" (positive Log2FC)
        {'Gene_ID': 'ENSG00000144648.14', 'Gene_Symbol': 'ACKR2', 'Log2FC': 0.59, 'q_value': 3.45e-2, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000165895.17', 'Gene_Symbol': 'ARHGAP42', 'Log2FC': 0.59, 'q_value': 4.09e-8, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000143479.15', 'Gene_Symbol': 'DYRK3', 'Log2FC': 0.59, 'q_value': 6.97e-9, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000183111.11', 'Gene_Symbol': 'ARHGEF37', 'Log2FC': 0.60, 'q_value': 6.90e-4, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000181856.14', 'Gene_Symbol': 'SLC2A4', 'Log2FC': 0.61, 'q_value': 6.72e-3, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000154114.12', 'Gene_Symbol': 'TBCEL', 'Log2FC': 0.62, 'q_value': 2.26e-11, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000188175.9', 'Gene_Symbol': 'HEPACAM2', 'Log2FC': 0.62, 'q_value': 3.50e-7, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000108309.12', 'Gene_Symbol': 'RUNDC3A', 'Log2FC': 0.63, 'q_value': 1.02e-3, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000162543.5', 'Gene_Symbol': 'UBXN10', 'Log2FC': 0.63, 'q_value': 3.89e-9, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000122644.12', 'Gene_Symbol': 'ARL4A', 'Log2FC': 0.66, 'q_value': 3.01e-9, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000152484.13', 'Gene_Symbol': 'USP12', 'Log2FC': 0.66, 'q_value': 7.20e-11, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000102904.14', 'Gene_Symbol': 'TSNAXIP1', 'Log2FC': 0.67, 'q_value': 9.05e-4, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000196407.11', 'Gene_Symbol': 'THEM5', 'Log2FC': 0.67, 'q_value': 2.28e-2, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000134184.12', 'Gene_Symbol': 'GSTM1', 'Log2FC': 0.69, 'q_value': 2.86e-4, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000165061.14', 'Gene_Symbol': 'ZMAT4', 'Log2FC': 0.72, 'q_value': 3.47e-2, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000111339.10', 'Gene_Symbol': 'ART4', 'Log2FC': 0.75, 'q_value': 2.49e-9, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000198892.6', 'Gene_Symbol': 'SHISA4', 'Log2FC': 0.78, 'q_value': 6.50e-3, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000163827.12', 'Gene_Symbol': 'LRRC2', 'Log2FC': 0.80, 'q_value': 9.63e-11, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000188582.8', 'Gene_Symbol': 'PAQR9', 'Log2FC': 0.80, 'q_value': 1.78e-5, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000066468.20', 'Gene_Symbol': 'FGFR2', 'Log2FC': 0.88, 'q_value': 2.91e-2, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000179914.4', 'Gene_Symbol': 'ITLN1', 'Log2FC': 1.09, 'q_value': 7.95e-7, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000176641.10', 'Gene_Symbol': 'RNF152', 'Log2FC': 1.18, 'q_value': 2.88e-13, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000145423.4', 'Gene_Symbol': 'SFRP2', 'Log2FC': 1.40, 'q_value': 2.44e-13, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000188536.12', 'Gene_Symbol': 'HBA2', 'Log2FC': 1.45, 'q_value': 5.37e-6, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000276345.1', 'Gene_Symbol': 'AC004556.1', 'Log2FC': 1.54, 'q_value': 2.16e-6, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000206172.8', 'Gene_Symbol': 'HBA1', 'Log2FC': 1.68, 'q_value': 2.23e-4, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000249992.1', 'Gene_Symbol': 'TMEM158', 'Log2FC': 1.69, 'q_value': 2.08e-14, 'Direction': 'Increased'},
        {'Gene_ID': 'ENSG00000180537.12', 'Gene_Symbol': 'RNF182', 'Log2FC': 1.78, 'q_value': 5.49e-17, 'Direction': 'Increased'},
    ]
    
    df = pd.DataFrame(genes_data)
    
    # Add significance flags
    df['Significant_01'] = df['q_value'] < 0.01
    df['Significant_001'] = df['q_value'] < 0.001
    df['Significant_0001'] = df['q_value'] < 0.0001
    
    # Add fold change category
    df['Abs_Log2FC'] = np.abs(df['Log2FC'])
    df['FC_Category'] = df['Abs_Log2FC'].apply(lambda x: 'High' if x > 1.0 else 'Moderate')
    
    return df

def extract_pat_score_info():
    """Extract PAT score formula information"""
    
    pat_info = {
        'name': 'Predicted Aneurysm Trajectory (PAT)',
        'description': 'Regression model that predicts the rate at which an IA will grow',
        'reference': 'Chien et al. 2019, J. Neurosurg. 2019, 132:1077-1087',
        'validation_dataset': '520 longitudinally followed IA cases',
        'factors': [
            'Initial IA size',
            'Age',
            'Sex',
            'Cigarette smoking',
            'Hypothyroidism',
            'Follow up time in months'
        ],
        'key_findings': [
            'Smoking and hypothyroidism had significant effects on growth trajectories of larger IAs (>=7 mm)',
            'Initial aneurysm size was significantly related to growth',
            'Larger IAs (>=7 mm) tended to grow significantly faster than small IAs (<3 mm)',
            'IAs in males showed faster growth rates than IAs in females',
            'PAT was markedly greater in ruptured IAs compared to unruptured IAs'
        ],
        'calculation_method': 'Multivariate adaptive splines for analysis of longitudinal data in R',
        'usage_in_study': 'Median PAT score used to dichotomize into "growing" (PAT >= 4.6) vs "stable" (< 4.6)',
        'note': 'Full coefficients available in Chien et al. 2019 manuscript Table',
    }
    
    return pat_info

def check_geo_accession():
    """Check for GEO accession numbers in the extracted text"""
    
    text_file = Path("results/jpm_analysis/jpm-13-00266/jpm-13-00266_extracted_text.txt")
    
    if text_file.exists():
        with open(text_file, 'r', encoding='utf-8') as f:
            text = f.read()
        
        # Search for GEO patterns
        geo_pattern = r'GSE\d+'
        geo_matches = re.findall(geo_pattern, text)
        
        # Also check for data availability statements
        data_availability = []
        if 'data available' in text.lower():
            # Find context around data availability
            lines = text.split('\n')
            for i, line in enumerate(lines):
                if 'data available' in line.lower() or 'data is available' in line.lower():
                    context_start = max(0, i-2)
                    context_end = min(len(lines), i+3)
                    context = '\n'.join(lines[context_start:context_end])
                    data_availability.append(context)
        
        return {
            'geo_accessions': list(set(geo_matches)),
            'data_availability_statement': data_availability,
            'note': 'Data is available upon reasonable request to corresponding author (vincentt@buffalo.edu)'
        }
    
    return {'geo_accessions': [], 'note': 'Could not find GEO accession numbers in paper'}

def integrate_with_existing_analyses(gene_df):
    """Integrate extracted genes with existing analyses"""
    
    print("\n" + "=" * 80)
    print("INTEGRATING WITH EXISTING ANALYSES")
    print("=" * 80)
    
    # Check for existing DE results
    existing_de_files = [
        "results/gse122897_differential_expression.tsv",
        "results/gse122897_top_genes.tsv",
    ]
    
    existing_results = {}
    
    for file_path in existing_de_files:
        if Path(file_path).exists():
            print(f"\n[INFO] Found existing results: {file_path}")
            df = pd.read_csv(file_path, sep='\t')
            existing_results[file_path] = df
            
            # Check overlap with Poppenberg genes
            if 'Gene' in df.columns:
                # Convert to string and uppercase for comparison
                existing_genes = set(df['Gene'].astype(str).str.upper())
                poppenberg_genes = set(gene_df['Gene_Symbol'].astype(str).str.upper())
                overlap = existing_genes & poppenberg_genes
                if overlap:
                    print(f"  [OK] Found {len(overlap)} overlapping genes: {', '.join(sorted(overlap))}")
    
    # Identify top priority genes (mentioned as common across multiple analyses)
    priority_genes = {
        'Common_Across_All_Metrics': ['DEFA1', 'HBA1', 'HBA2', 'OTOF', 'SFRP2'],
        'Common_PAT_Size_PHASES': ['DEFA1', 'HBA1', 'HBA2', 'OTOF', 'SFRP2'],
        'High_Prognostic_Potential': ['DEFA1', 'HBA1', 'HBA2', 'SFRP2', 'IFI27', 'LPL'],
        'Common_with_Other_Studies': ['FAM83A', 'GSTM1', 'RNF152', 'RNF182', 'RUNDC3A']
    }
    
    # Mark priority genes
    gene_df['Priority_Status'] = 'Standard'
    for category, genes in priority_genes.items():
        for gene in genes:
            mask = gene_df['Gene_Symbol'].str.upper() == gene.upper()
            if mask.any():
                gene_df.loc[mask, 'Priority_Status'] = category.replace('_', ' ')
    
    return gene_df, existing_results, priority_genes

def save_extracted_findings():
    """Save all extracted findings"""
    
    print("\n" + "=" * 80)
    print("EXTRACTING JPM POPPENBERG ET AL. 2023 FINDINGS")
    print("=" * 80)
    
    # Create output directory
    output_dir = Path("results/poppenberg_2023_integration")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract 39 genes
    print("\n[1] Extracting 39 differentially expressed genes...")
    gene_df = extract_39_genes()
    
    print(f"  [OK] Extracted {len(gene_df)} genes")
    print(f"  - Decreased in growing: {(gene_df['Direction'] == 'Decreased').sum()}")
    print(f"  - Increased in growing: {(gene_df['Direction'] == 'Increased').sum()}")
    print(f"  - High fold-change (>1.0): {(gene_df['FC_Category'] == 'High').sum()}")
    
    # Extract PAT score info
    print("\n[2] Extracting PAT score information...")
    pat_info = extract_pat_score_info()
    print(f"  [OK] PAT score factors: {', '.join(pat_info['factors'])}")
    
    # Check for GEO accession
    print("\n[3] Checking for GEO accession numbers...")
    geo_info = check_geo_accession()
    if geo_info['geo_accessions']:
        print(f"  [OK] Found GEO accessions: {', '.join(geo_info['geo_accessions'])}")
    else:
        print(f"  [INFO] {geo_info['note']}")
    
    # Integrate with existing analyses
    print("\n[4] Integrating with existing analyses...")
    gene_df_integrated, existing_results, priority_genes = integrate_with_existing_analyses(gene_df)
    
    # Save gene list
    gene_file = output_dir / "poppenberg_2023_39_genes.tsv"
    gene_df_integrated.to_csv(gene_file, sep='\t', index=False)
    print(f"\n[OK] Saved gene list to: {gene_file}")
    
    # Save PAT score info
    pat_file = output_dir / "pat_score_formula_info.json"
    with open(pat_file, 'w') as f:
        json.dump(pat_info, f, indent=2)
    print(f"[OK] Saved PAT score info to: {pat_file}")
    
    # Save GEO info
    geo_file = output_dir / "geo_accession_info.json"
    with open(geo_file, 'w') as f:
        json.dump(geo_info, f, indent=2)
    print(f"[OK] Saved GEO info to: {geo_file}")
    
    # Save priority genes
    priority_file = output_dir / "priority_genes.json"
    with open(priority_file, 'w') as f:
        json.dump(priority_genes, f, indent=2)
    print(f"[OK] Saved priority genes to: {priority_file}")
    
    # Create summary report
    summary = {
        'paper': {
            'title': 'RNA Expression Signatures of Intracranial Aneurysm Growth Trajectory Identified in Circulating Whole Blood',
            'authors': 'Poppenberg, K.E. et al.',
            'journal': 'Journal of Personalized Medicine',
            'year': 2023,
            'doi': '10.3390/jpm13020266'
        },
        'dataset': {
            'n_samples': 66,
            'training': 46,
            'testing': 20,
            'growing': 33,
            'stable': 33,
            'pat_threshold': 4.6
        },
        'findings': {
            'n_genes': len(gene_df),
            'n_decreased': (gene_df['Direction'] == 'Decreased').sum(),
            'n_increased': (gene_df['Direction'] == 'Increased').sum(),
            'model_auc_training': 0.85,
            'model_auc_testing': 0.86,
            'model_accuracy_testing': 0.80
        },
        'priority_genes': priority_genes,
        'pat_score': pat_info,
        'geo_accession': geo_info
    }
    
    summary_file = output_dir / "extraction_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[OK] Saved summary to: {summary_file}")
    
    # Print summary
    print("\n" + "=" * 80)
    print("EXTRACTION SUMMARY")
    print("=" * 80)
    print(f"\nPaper: {summary['paper']['title']}")
    print(f"Authors: {summary['paper']['authors']}")
    print(f"DOI: {summary['paper']['doi']}")
    print(f"\nDataset: {summary['dataset']['n_samples']} samples ({summary['dataset']['training']} training, {summary['dataset']['testing']} testing)")
    print(f"\nFindings:")
    print(f"  - {summary['findings']['n_genes']} differentially expressed genes")
    print(f"  - {summary['findings']['n_decreased']} decreased in growing aneurysms")
    print(f"  - {summary['findings']['n_increased']} increased in growing aneurysms")
    print(f"  - Model AUC: Training={summary['findings']['model_auc_training']}, Testing={summary['findings']['model_auc_testing']}")
    print(f"\nPriority Genes (common across analyses):")
    for category, genes in priority_genes.items():
        print(f"  - {category}: {', '.join(genes)}")
    
    print("\n" + "=" * 80)
    print("EXTRACTION COMPLETE!")
    print("=" * 80)
    print(f"\nAll results saved to: {output_dir}/")
    
    return gene_df_integrated, pat_info, geo_info, summary

if __name__ == "__main__":
    save_extracted_findings()

