"""
Verify extracted PDF data and check for overlaps with existing datasets
"""
import pandas as pd
import json
from pathlib import Path
from typing import Dict, List


def load_extracted_gwas_data() -> Dict[str, pd.DataFrame]:
    """Load all extracted GWAS data"""
    gwas_dir = Path("data/supplementary_pdfs/extracted_gwas")
    gwas_files = list(gwas_dir.glob("*_gwas_parsed.tsv"))
    
    gwas_data = {}
    for gwas_file in gwas_files:
        try:
            df = pd.read_csv(gwas_file, sep='\t')
            pdf_name = gwas_file.stem.replace('_gwas_parsed', '')
            gwas_data[pdf_name] = df
            print(f"  ✓ Loaded {len(df)} SNPs from {pdf_name}")
        except Exception as e:
            print(f"  ⚠ Error loading {gwas_file}: {e}")
    
    return gwas_data


def identify_paper_sources(gwas_data: Dict[str, pd.DataFrame]) -> Dict:
    """Identify and verify paper sources"""
    print("\n" + "=" * 80)
    print("PAPER IDENTIFICATION AND VERIFICATION")
    print("=" * 80)
    
    paper_info = {}
    
    for pdf_name, df in gwas_data.items():
        info = {
            'pdf_name': pdf_name,
            'num_snps': len(df),
            'study': df['study'].iloc[0] if 'study' in df.columns else 'Unknown',
            'pubmed_id': df['pubmed_id'].iloc[0] if 'pubmed_id' in df.columns else None,
        }
        
        # Get unique chromosomes
        if 'chr' in df.columns:
            info['chromosomes'] = sorted(df['chr'].dropna().unique().tolist())
        
        # Get p-value statistics
        if 'p' in df.columns:
            p_vals = df['p'].dropna()
            if len(p_vals) > 0:
                info['p_value_stats'] = {
                    'min': float(p_vals.min()),
                    'max': float(p_vals.max()),
                    'significant_5e8': int((p_vals < 5e-8).sum()),
                    'significant_5e7': int((p_vals < 5e-7).sum()),
                    'significant_1e5': int((p_vals < 1e-5).sum()),
                }
        
        paper_info[pdf_name] = info
        
        print(f"\n📄 {pdf_name}")
        print(f"  Study: {info['study']}")
        if info['pubmed_id']:
            print(f"  PubMed ID: {info['pubmed_id']}")
        print(f"  SNPs: {info['num_snps']}")
        if 'chromosomes' in info:
            print(f"  Chromosomes: {info['chromosomes']}")
        if 'p_value_stats' in info:
            stats = info['p_value_stats']
            print(f"  P-value range: {stats['min']:.2e} - {stats['max']:.2e}")
            print(f"  Genome-wide significant (p < 5e-8): {stats['significant_5e8']}")
    
    return paper_info


def check_overlaps(gwas_data: Dict[str, pd.DataFrame]) -> Dict:
    """Check for overlaps with existing datasets"""
    print("\n" + "=" * 80)
    print("OVERLAP ANALYSIS")
    print("=" * 80)
    
    overlap_results = {}
    
    # Get all SNP IDs from extracted data
    all_extracted_snps = set()
    for pdf_name, df in gwas_data.items():
        if 'snp' in df.columns:
            snps = set(df['snp'].dropna().astype(str))
            all_extracted_snps.update(snps)
            overlap_results[pdf_name] = {
                'snps': list(snps),
                'num_snps': len(snps)
            }
    
    print(f"\nTotal unique SNPs extracted: {len(all_extracted_snps)}")
    
    # Check for overlaps between PDFs
    if len(gwas_data) > 1:
        pdf_names = list(gwas_data.keys())
        for i, pdf1 in enumerate(pdf_names):
            for pdf2 in pdf_names[i+1:]:
                snps1 = set(gwas_data[pdf1]['snp'].dropna().astype(str)) if 'snp' in gwas_data[pdf1].columns else set()
                snps2 = set(gwas_data[pdf2]['snp'].dropna().astype(str)) if 'snp' in gwas_data[pdf2].columns else set()
                overlap = snps1.intersection(snps2)
                if overlap:
                    print(f"\n  Overlap between {pdf1} and {pdf2}: {len(overlap)} SNPs")
                    print(f"    Common SNPs: {', '.join(list(overlap)[:10])}")
    
    # Check against known GWAS datasets (if available)
    print("\n  Checking against known GWAS datasets...")
    
    # Known risk loci from literature
    known_risk_snps = {
        'rs1333040': {'chr': 9, 'study': 'Multiple', 'description': 'CDKN2B-AS1, 9p21.3'},
        'rs10958409': {'chr': 8, 'study': 'Bilguvar 2008', 'description': 'Chromosome 8q11.23'},
        'rs700651': {'chr': 2, 'study': 'Bilguvar 2008', 'description': 'Chromosome 2q'},
        'rs10116277': {'chr': 9, 'study': 'Bilguvar 2008', 'description': '9p21.3'},
        'rs2383207': {'chr': 9, 'study': 'Bilguvar 2008', 'description': '9p21.3'},
    }
    
    found_known = {}
    for pdf_name, df in gwas_data.items():
        if 'snp' in df.columns:
            snps = set(df['snp'].dropna().astype(str))
            found = {snp: known_risk_snps[snp] for snp in snps if snp in known_risk_snps}
            if found:
                found_known[pdf_name] = found
                print(f"\n  {pdf_name}: Found {len(found)} known risk SNPs")
                for snp, info in found.items():
                    print(f"    - {snp}: {info['description']}")
    
    overlap_results['known_risk_snps_found'] = found_known
    
    return overlap_results


def identify_new_data(gwas_data: Dict[str, pd.DataFrame], paper_info: Dict) -> Dict:
    """Identify any new data not in existing datasets"""
    print("\n" + "=" * 80)
    print("NEW DATA IDENTIFICATION")
    print("=" * 80)
    
    new_data_summary = {}
    
    for pdf_name, df in gwas_data.items():
        info = paper_info.get(pdf_name, {})
        
        new_data = {
            'pdf_name': pdf_name,
            'study': info.get('study', 'Unknown'),
            'total_snps': len(df),
            'novel_snps': [],
            'highly_significant_snps': []
        }
        
        # Identify highly significant SNPs
        if 'p' in df.columns:
            sig_snps = df[df['p'] < 5e-8].copy() if 'p' in df.columns else pd.DataFrame()
            if len(sig_snps) > 0:
                new_data['highly_significant_snps'] = sig_snps[['snp', 'chr', 'pos', 'p']].to_dict('records')
                print(f"\n📄 {pdf_name}")
                print(f"  Genome-wide significant SNPs (p < 5e-8): {len(sig_snps)}")
                for _, row in sig_snps.head(10).iterrows():
                    print(f"    - {row['snp']}: chr{row.get('chr', '?')}, p={row['p']:.2e}")
        
        new_data_summary[pdf_name] = new_data
    
    return new_data_summary


def generate_summary_report(gwas_data: Dict, paper_info: Dict, overlap_results: Dict, new_data: Dict):
    """Generate comprehensive summary report"""
    report = {
        'extraction_date': pd.Timestamp.now().isoformat(),
        'summary': {
            'total_pdfs_processed': len(gwas_data),
            'total_snps_extracted': sum(len(df) for df in gwas_data.values()),
            'unique_snps': len(set().union(*[set(df['snp'].dropna().astype(str)) 
                                           for df in gwas_data.values() 
                                           if 'snp' in df.columns]))
        },
        'papers': paper_info,
        'overlaps': overlap_results,
        'new_data': new_data
    }
    
    # Save report
    report_file = Path("data/supplementary_pdfs/extracted_gwas/verification_report.json")
    with open(report_file, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=2, ensure_ascii=False, default=str)
    
    print("\n" + "=" * 80)
    print("VERIFICATION COMPLETE")
    print("=" * 80)
    print(f"\n✓ Summary:")
    print(f"  - PDFs processed: {report['summary']['total_pdfs_processed']}")
    print(f"  - Total SNPs extracted: {report['summary']['total_snps_extracted']}")
    print(f"  - Unique SNPs: {report['summary']['unique_snps']}")
    print(f"\n✓ Report saved to: {report_file}")


def main():
    """Main verification function"""
    print("=" * 80)
    print("PDF EXTRACTION VERIFICATION")
    print("=" * 80)
    
    # Load extracted data
    print("\n📁 Loading extracted GWAS data...")
    gwas_data = load_extracted_gwas_data()
    
    if not gwas_data:
        print("\n❌ No GWAS data found to verify")
        return
    
    # Identify paper sources
    paper_info = identify_paper_sources(gwas_data)
    
    # Check overlaps
    overlap_results = check_overlaps(gwas_data)
    
    # Identify new data
    new_data = identify_new_data(gwas_data, paper_info)
    
    # Generate report
    generate_summary_report(gwas_data, paper_info, overlap_results, new_data)


if __name__ == "__main__":
    main()






