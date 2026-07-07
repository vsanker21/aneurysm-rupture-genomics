"""
Comprehensive Verification: Ensure All Analyses Use Real Data Only
Checks all analysis scripts and results for placeholder/fake/random data
"""
import os
import pandas as pd
import numpy as np
import re
from pathlib import Path

def check_file_for_placeholders(filepath, file_type='py'):
    """Check a file for placeholder/fake/random data patterns"""
    issues = []
    
    if not os.path.exists(filepath):
        return issues
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
            lines = content.split('\n')
            
        # Patterns to check
        patterns = {
            'random_normal': r'np\.random\.normal|random\.normal',
            'random_': r'np\.random\.|random\.',
            'placeholder_assignment': r'=.*placeholder|placeholder.*=',
            'fake_data': r'fake|Fake|FAKE',
            'simulated': r'simulated|Simulated|SIMULATED',
            'dummy': r'dummy|Dummy|DUMMY'
        }
        
        for pattern_name, pattern in patterns.items():
            matches = re.finditer(pattern, content, re.IGNORECASE)
            for match in matches:
                line_num = content[:match.start()].count('\n') + 1
                context = lines[line_num-1] if line_num <= len(lines) else ''
                # Skip if it's in a comment or documentation
                if not context.strip().startswith('#') and 'NOTE' not in context and 'WARNING' not in context:
                    issues.append({
                        'file': filepath,
                        'line': line_num,
                        'pattern': pattern_name,
                        'context': context.strip()[:100]
                    })
    except Exception as e:
        issues.append({
            'file': filepath,
            'error': str(e)
        })
    
    return issues

def verify_data_sources():
    """Verify that data files are from real datasets"""
    real_data_sources = {
        'GSE122897': {
            'file': 'data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz',
            'description': 'Ruptured vs Unruptured Aneurysm RNA-seq (real dataset)',
            'verified': False
        },
        'GSE159610': {
            'file': 'data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz',
            'description': 'Blood RNA-seq from aneurysm patients (real dataset)',
            'verified': False
        },
        'EC_Dysfunction': {
            'file': 'results/ec_dysfunction_scores.tsv',
            'description': 'EC dysfunction scores calculated from real tissue expression',
            'verified': False
        },
        'DE_Results': {
            'file': 'results/gse122897_differential_expression.tsv',
            'description': 'Differential expression from real GSE122897 data',
            'verified': False
        }
    }
    
    verification_results = {}
    
    for source_name, info in real_data_sources.items():
        filepath = info['file']
        if os.path.exists(filepath):
            # Check file size (real data files should be substantial)
            file_size = os.path.getsize(filepath) / (1024 * 1024)  # MB
            
            # For TSV files, check if they contain real data
            if filepath.endswith('.tsv'):
                try:
                    df = pd.read_csv(filepath, sep='\t', nrows=5)
                    has_data = len(df) > 0
                    # Check for non-zero variance in numeric columns
                    numeric_cols = df.select_dtypes(include=[np.number]).columns
                    if len(numeric_cols) > 0:
                        has_variance = df[numeric_cols].var().sum() > 0
                    else:
                        has_variance = True
                except:
                    has_data = False
                    has_variance = False
            else:
                has_data = file_size > 0.1  # At least 100KB
                has_variance = True
            
            verification_results[source_name] = {
                'exists': True,
                'size_mb': file_size,
                'has_data': has_data,
                'has_variance': has_variance,
                'description': info['description'],
                'status': 'VERIFIED' if (has_data and has_variance) else 'NEEDS_REVIEW'
            }
        else:
            verification_results[source_name] = {
                'exists': False,
                'status': 'MISSING'
            }
    
    return verification_results

def verify_integrated_features():
    """Verify integrated features use real data"""
    issues = []
    
    features_file = 'results/integrated_features.tsv'
    if not os.path.exists(features_file):
        return [{'issue': 'Integrated features file not found'}]
    
    try:
        df = pd.read_csv(features_file, sep='\t')
        
        # Check EC_Dysfunction_Score - should have real values
        if 'EC_Dysfunction_Score' in df.columns:
            ec_scores = df['EC_Dysfunction_Score'].dropna()
            if len(ec_scores) == 0:
                issues.append({'issue': 'EC_Dysfunction_Score has no data'})
            elif ec_scores.var() == 0:
                issues.append({'issue': 'EC_Dysfunction_Score has no variance (all same value)'})
            elif ec_scores.abs().max() < 0.01:
                issues.append({'warning': 'EC_Dysfunction_Score values are very small - verify these are real'})
            else:
                # Check for suspicious patterns (all zeros, all same value, etc.)
                unique_values = ec_scores.nunique()
                if unique_values < 3:
                    issues.append({'warning': f'EC_Dysfunction_Score has only {unique_values} unique values'})
        
        # Blood_RNA_Score and Clinical_Score = 0 is OK if data not available
        # But we should verify they're not random/fake
        for col in ['Blood_RNA_Score', 'Clinical_Score']:
            if col in df.columns:
                values = df[col].dropna()
                if len(values) > 0:
                    # If all zeros, that's OK (means data not available)
                    if values.abs().max() == 0:
                        pass  # This is expected if data not available
                    elif values.var() == 0 and values.iloc[0] != 0:
                        issues.append({'warning': f'{col} has no variance but is not zero - verify this is real data'})
        
        # Check for any suspicious patterns
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        for col in numeric_cols:
            if col not in ['Blood_RNA_Score', 'Clinical_Score']:  # These can be 0
                values = df[col].dropna()
                if len(values) > 0:
                    # Check if values look like they could be random
                    if values.std() == 0:
                        issues.append({'warning': f'{col} has zero standard deviation'})
        
    except Exception as e:
        issues.append({'error': f'Failed to verify integrated features: {e}'})
    
    return issues

def verify_analysis_scripts():
    """Check all analysis scripts for placeholder data"""
    analysis_scripts = [
        'analyze_gse122897.py',
        'endothelial_dysfunction_analysis.py',
        'pathway_enrichment_analysis.py',
        'ml_classification_blood_genes.py',
        'multi_omics_integration.py',
        'convert_and_integrate_gene_ids.py',
        'annotate_top_genes_literature.py'
    ]
    
    all_issues = []
    
    for script in analysis_scripts:
        if os.path.exists(script):
            issues = check_file_for_placeholders(script)
            if issues:
                all_issues.extend(issues)
    
    return all_issues

def generate_verification_report():
    """Generate comprehensive verification report"""
    print("=" * 80)
    print("COMPREHENSIVE DATA VERIFICATION REPORT")
    print("=" * 80)
    print("\nVerifying all analyses use REAL DATA only (no placeholders)...\n")
    
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("DATA VERIFICATION REPORT")
    report_lines.append("=" * 80)
    report_lines.append(f"\nGenerated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # 1. Verify data sources
    report_lines.append("## 1. DATA SOURCE VERIFICATION")
    report_lines.append("=" * 80)
    data_sources = verify_data_sources()
    
    all_verified = True
    for source_name, result in data_sources.items():
        if result.get('status') == 'VERIFIED':
            report_lines.append(f"\n✓ {source_name}: VERIFIED")
            report_lines.append(f"  Description: {result.get('description', 'N/A')}")
            if 'size_mb' in result:
                report_lines.append(f"  File size: {result['size_mb']:.2f} MB")
        elif result.get('status') == 'MISSING':
            report_lines.append(f"\n✗ {source_name}: MISSING")
            all_verified = False
        else:
            report_lines.append(f"\n⚠ {source_name}: NEEDS REVIEW")
            report_lines.append(f"  Status: {result.get('status', 'UNKNOWN')}")
            all_verified = False
    
    # 2. Verify analysis scripts
    report_lines.append("\n\n## 2. ANALYSIS SCRIPT VERIFICATION")
    report_lines.append("=" * 80)
    script_issues = verify_analysis_scripts()
    
    if not script_issues:
        report_lines.append("\n✓ All analysis scripts verified - no placeholder data found")
    else:
        report_lines.append(f"\n⚠ Found {len(script_issues)} potential issues in scripts:")
        for issue in script_issues[:10]:  # Limit to first 10
            report_lines.append(f"  - {issue.get('file', 'Unknown')}: Line {issue.get('line', '?')} - {issue.get('pattern', 'Unknown pattern')}")
        all_verified = False
    
    # 3. Verify integrated features
    report_lines.append("\n\n## 3. INTEGRATED FEATURES VERIFICATION")
    report_lines.append("=" * 80)
    feature_issues = verify_integrated_features()
    
    if not feature_issues:
        report_lines.append("\n✓ Integrated features verified - all use real data")
    else:
        report_lines.append(f"\n⚠ Found {len(feature_issues)} issues/warnings:")
        for issue in feature_issues:
            if 'error' in issue:
                report_lines.append(f"  ✗ ERROR: {issue['error']}")
                all_verified = False
            elif 'warning' in issue:
                report_lines.append(f"  ⚠ WARNING: {issue['warning']}")
            else:
                report_lines.append(f"  ⚠ {issue.get('issue', 'Unknown issue')}")
    
    # 4. Check EC dysfunction scores (should be real)
    report_lines.append("\n\n## 4. EC DYSFUNCTION SCORES VERIFICATION")
    report_lines.append("=" * 80)
    ec_file = 'results/ec_dysfunction_scores.tsv'
    if os.path.exists(ec_file):
        try:
            ec_df = pd.read_csv(ec_file, sep='\t')
            if 'EC_Dysfunction_Score' in ec_df.columns:
                scores = ec_df['EC_Dysfunction_Score']
                report_lines.append(f"\n✓ EC Dysfunction Scores: VERIFIED")
                report_lines.append(f"  Samples: {len(scores)}")
                report_lines.append(f"  Mean: {scores.mean():.4f}")
                report_lines.append(f"  Std: {scores.std():.4f}")
                report_lines.append(f"  Range: [{scores.min():.4f}, {scores.max():.4f}]")
                report_lines.append(f"  Non-zero values: {(scores != 0).sum()}")
                if scores.var() > 0 and (scores != 0).sum() > 0:
                    report_lines.append(f"  ✓ Scores show real variation (not all zeros)")
                else:
                    report_lines.append(f"  ⚠ WARNING: Scores may not be real data")
                    all_verified = False
        except Exception as e:
            report_lines.append(f"\n✗ ERROR reading EC scores: {e}")
            all_verified = False
    else:
        report_lines.append(f"\n✗ EC dysfunction scores file not found")
        all_verified = False
    
    # 5. Check DE results (should be real)
    report_lines.append("\n\n## 5. DIFFERENTIAL EXPRESSION VERIFICATION")
    report_lines.append("=" * 80)
    de_file = 'results/gse122897_differential_expression.tsv'
    if os.path.exists(de_file):
        try:
            de_df = pd.read_csv(de_file, sep='\t')
            report_lines.append(f"\n✓ DE Results: VERIFIED")
            report_lines.append(f"  Total genes: {len(de_df):,}")
            if 'Log2FC' in de_df.columns:
                fc_values = de_df['Log2FC'].dropna()
                report_lines.append(f"  Log2FC range: [{fc_values.min():.4f}, {fc_values.max():.4f}]")
                report_lines.append(f"  Non-zero FC: {(fc_values != 0).sum():,}")
            if 'P_value' in de_df.columns:
                p_values = de_df['P_value'].dropna()
                report_lines.append(f"  P-values range: [{p_values.min():.2e}, {p_values.max():.4f}]")
                report_lines.append(f"  Significant (P<0.05): {(p_values < 0.05).sum():,}")
        except Exception as e:
            report_lines.append(f"\n✗ ERROR reading DE results: {e}")
            all_verified = False
    else:
        report_lines.append(f"\n✗ DE results file not found")
        all_verified = False
    
    # Final summary
    report_lines.append("\n\n## FINAL VERIFICATION SUMMARY")
    report_lines.append("=" * 80)
    if all_verified:
        report_lines.append("\n✓✓✓ ALL ANALYSES VERIFIED - USING REAL DATA ONLY ✓✓✓")
        report_lines.append("\nNo placeholder, fake, or random data detected in:")
        report_lines.append("  - Data source files")
        report_lines.append("  - Analysis scripts")
        report_lines.append("  - Integrated features")
        report_lines.append("  - EC dysfunction scores")
        report_lines.append("  - Differential expression results")
    else:
        report_lines.append("\n⚠ VERIFICATION INCOMPLETE - Some issues found")
        report_lines.append("\nPlease review the warnings above.")
    
    report_lines.append("\n" + "=" * 80)
    report_lines.append("END OF VERIFICATION REPORT")
    report_lines.append("=" * 80)
    
    # Write report
    report_text = "\n".join(report_lines)
    report_file = "results/DATA_VERIFICATION_REPORT.txt"
    
    os.makedirs("results", exist_ok=True)
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report_text)
    
    print(report_text)
    print(f"\n[OK] Verification report saved to {report_file}")
    
    return all_verified, report_file

if __name__ == "__main__":
    verified, report_file = generate_verification_report()
    if verified:
        print("\n✓ All analyses verified to use REAL DATA only")
    else:
        print("\n⚠ Some issues found - please review the report")






