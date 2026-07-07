"""
Power Analysis for GSE122897 Differential Expression
Assesses detectable effect sizes given current sample size
"""
import pandas as pd
import numpy as np
import os
from scipy import stats
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

def load_de_results():
    """Load differential expression results"""
    de_file = "results/gse122897_differential_expression.tsv"
    if not os.path.exists(de_file):
        print(f"[ERROR] DE results file not found: {de_file}")
        return None
    
    df = pd.read_csv(de_file, sep='\t')
    return df

def calculate_power(n1, n2, effect_size, alpha=0.05):
    """
    Calculate statistical power for two-sample t-test
    
    Parameters:
    n1, n2: sample sizes for two groups
    effect_size: Cohen's d (standardized effect size)
    alpha: significance level
    
    Returns:
    power: statistical power (1 - beta)
    """
    # Degrees of freedom
    df = n1 + n2 - 2
    
    # Non-centrality parameter
    ncp = effect_size * np.sqrt((n1 * n2) / (n1 + n2))
    
    # Critical t-value
    t_crit = stats.t.ppf(1 - alpha/2, df)
    
    # Power calculation using non-central t-distribution
    power = 1 - stats.nct.cdf(t_crit, df, ncp) + stats.nct.cdf(-t_crit, df, ncp)
    
    return power

def power_analysis_for_sample_size(n1=22, n2=21):
    """
    Perform power analysis for current sample sizes
    """
    print("=" * 80)
    print("POWER ANALYSIS")
    print("=" * 80)
    
    print(f"\nCurrent sample sizes:")
    print(f"  Group 1 (Ruptured): {n1}")
    print(f"  Group 2 (Unruptured): {n2}")
    print(f"  Total: {n1 + n2}")
    
    # Calculate power for different effect sizes
    effect_sizes = np.arange(0.2, 3.1, 0.1)  # Cohen's d from 0.2 to 3.0
    powers = []
    
    for d in effect_sizes:
        power = calculate_power(n1, n2, d, alpha=0.05)
        powers.append(power)
    
    # Create power curve
    power_df = pd.DataFrame({
        'Effect_size_Cohens_d': effect_sizes,
        'Power': powers
    })
    
    # Find effect sizes at different power levels
    power_80_idx = np.argmin(np.abs(np.array(powers) - 0.80))
    power_90_idx = np.argmin(np.abs(np.array(powers) - 0.90))
    
    effect_80 = effect_sizes[power_80_idx]
    effect_90 = effect_sizes[power_90_idx]
    
    print(f"\nDetectable effect sizes:")
    print(f"  At 80% power: Cohen's d = {effect_80:.2f}")
    print(f"  At 90% power: Cohen's d = {effect_90:.2f}")
    
    # Convert Cohen's d to Log2FC (approximate)
    # For RNA-seq, typical SD of log2 expression is ~0.5-1.0
    # Log2FC ≈ Cohen's d × SD
    sd_estimate = 0.75  # Typical SD for log2 TPM
    log2fc_80 = effect_80 * sd_estimate
    log2fc_90 = effect_90 * sd_estimate
    
    print(f"\nApproximate detectable Log2FC (assuming SD ≈ 0.75):")
    print(f"  At 80% power: Log2FC ≈ {log2fc_80:.2f}")
    print(f"  At 90% power: Log2FC ≈ {log2fc_90:.2f}")
    
    # Plot power curve
    plt.figure(figsize=(10, 6))
    plt.plot(effect_sizes, powers, 'b-', linewidth=2, label='Power curve')
    plt.axhline(y=0.80, color='r', linestyle='--', label='80% power')
    plt.axhline(y=0.90, color='g', linestyle='--', label='90% power')
    plt.axvline(x=effect_80, color='r', linestyle=':', alpha=0.5)
    plt.axvline(x=effect_90, color='g', linestyle=':', alpha=0.5)
    plt.xlabel("Effect Size (Cohen's d)", fontsize=12)
    plt.ylabel('Statistical Power (1 - β)', fontsize=12)
    plt.title(f'Power Analysis: n1={n1}, n2={n2}, α=0.05', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    output_file = "results/power_analysis_curve.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n[OK] Saved power curve to {output_file}")
    plt.close()
    
    # Save power table
    power_df.to_csv("results/power_analysis_table.tsv", sep='\t', index=False)
    print(f"[OK] Saved power analysis table to results/power_analysis_table.tsv")
    
    return power_df, effect_80, effect_90

def analyze_observed_effects(df):
    """Analyze observed effect sizes in the data"""
    print("\n" + "=" * 80)
    print("OBSERVED EFFECT SIZE ANALYSIS")
    print("=" * 80)
    
    # Calculate Cohen's d for all genes (approximate)
    # Cohen's d ≈ mean difference / pooled SD
    # For our case, we can approximate using Log2FC and assume similar SDs
    
    # Get significant genes at different thresholds
    for fdr_threshold in [0.05, 0.1, 0.2]:
        sig_genes = df[df['P_adj_BH'] < fdr_threshold].copy()
        if len(sig_genes) == 0:
            continue
        
        print(f"\nFDR < {fdr_threshold}: {len(sig_genes)} genes")
        print(f"  Mean |Log2FC|: {sig_genes['Log2FC'].abs().mean():.3f}")
        print(f"  Median |Log2FC|: {sig_genes['Log2FC'].abs().median():.3f}")
        print(f"  Max |Log2FC|: {sig_genes['Log2FC'].abs().max():.3f}")
        print(f"  Min |Log2FC|: {sig_genes['Log2FC'].abs().min():.3f}")
        
        # Estimate Cohen's d (approximate: Log2FC / typical_SD)
        typical_sd = 0.75
        estimated_d = sig_genes['Log2FC'].abs() / typical_sd
        
        print(f"  Estimated mean Cohen's d: {estimated_d.mean():.3f}")
        print(f"  Estimated median Cohen's d: {estimated_d.median():.3f}")

def main():
    """Main power analysis pipeline"""
    import os
    
    # Load DE results
    df = load_de_results()
    
    # Perform power analysis
    power_df, effect_80, effect_90 = power_analysis_for_sample_size(n1=22, n2=21)
    
    # Analyze observed effects
    df = load_de_results()
    if df is not None:
        analyze_observed_effects(df)
    
    print("\n" + "=" * 80)
    print("POWER ANALYSIS COMPLETE")
    print("=" * 80)
    print("\nRecommendations:")
    print("1. Current sample size (22 vs 21) provides good power for moderate-large effects")
    print("2. For small effects (d < 0.5), consider increasing sample size")
    print("3. Focus on genes with |Log2FC| > 0.5-1.0 for robust findings")

if __name__ == "__main__":
    main()

