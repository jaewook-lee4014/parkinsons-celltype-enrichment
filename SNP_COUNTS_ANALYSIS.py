#!/usr/bin/env python3
"""
Analyze SNP counts from completed LDSC datasets for reliability assessment
"""

# Olig_unique SNP counts by chromosome
olig_unique = {
    1: 216, 2: 137, 3: 134, 4: 67, 5: 118, 6: 92, 7: 96, 8: 113, 9: 107, 10: 115,
    11: 127, 12: 96, 13: 67, 14: 80, 15: 84, 16: 69, 17: 59, 18: 44, 19: 11, 20: 73, 21: 17, 22: 15
}

# Neg_unique SNP counts by chromosome  
neg_unique = {
    1: 949, 2: 754, 3: 645, 4: 435, 5: 595, 6: 683, 7: 503, 8: 335, 9: 422, 10: 555,
    11: 783, 12: 450, 13: 205, 14: 396, 15: 307, 16: 290, 17: 477, 18: 163, 19: 175, 20: 192, 21: 72, 22: 110
}

# Chromosome lengths (Mb) for hg19
chr_lengths = {
    1: 249.25, 2: 242.19, 3: 198.30, 4: 191.15, 5: 180.92, 6: 171.12, 7: 159.14, 8: 146.36, 9: 141.21, 10: 135.53,
    11: 134.45, 12: 133.85, 13: 115.17, 14: 107.35, 15: 102.53, 16: 90.35, 17: 81.20, 18: 78.08, 19: 59.13, 20: 63.03, 21: 48.13, 22: 51.30
}

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import pandas as pd

def analyze_snp_distribution():
    """Analyze SNP distribution patterns for reliability assessment"""
    
    print("🔍 SNP Distribution Analysis for Reliability Assessment")
    print("=" * 60)
    
    # Calculate totals
    olig_total = sum(olig_unique.values())
    neg_total = sum(neg_unique.values())
    
    print(f"📊 Total SNP Counts:")
    print(f"  Olig_unique: {olig_total:,} SNPs")
    print(f"  Neg_unique:  {neg_total:,} SNPs")
    print(f"  Ratio (Neg/Olig): {neg_total/olig_total:.2f}x")
    
    # SNP density analysis
    print(f"\n📈 SNP Density Analysis:")
    
    # Calculate densities (SNPs per Mb)
    olig_densities = [olig_unique[chr] / chr_lengths[chr] for chr in range(1, 23)]
    neg_densities = [neg_unique[chr] / chr_lengths[chr] for chr in range(1, 23)]
    
    print(f"  Olig_unique density: {np.mean(olig_densities):.2f} ± {np.std(olig_densities):.2f} SNPs/Mb")
    print(f"  Neg_unique density:  {np.mean(neg_densities):.2f} ± {np.std(neg_densities):.2f} SNPs/Mb")
    
    # Correlation with chromosome length
    chr_nums = list(range(1, 23))
    lengths = [chr_lengths[chr] for chr in chr_nums]
    olig_counts = [olig_unique[chr] for chr in chr_nums]
    neg_counts = [neg_unique[chr] for chr in chr_nums]
    
    olig_corr, olig_p = pearsonr(lengths, olig_counts)
    neg_corr, neg_p = pearsonr(lengths, neg_counts)
    
    print(f"\n🔗 Correlation with Chromosome Length:")
    print(f"  Olig_unique: r = {olig_corr:.3f} (p = {olig_p:.2e})")
    print(f"  Neg_unique:  r = {neg_corr:.3f} (p = {neg_p:.2e})")
    
    # Inter-dataset correlation
    inter_corr, inter_p = pearsonr(olig_counts, neg_counts)
    print(f"  Inter-dataset: r = {inter_corr:.3f} (p = {inter_p:.2e})")
    
    # Outlier detection
    print(f"\n⚠️  Potential Outliers:")
    
    # Z-score analysis for Olig_unique
    olig_z = np.abs((np.array(olig_counts) - np.mean(olig_counts)) / np.std(olig_counts))
    olig_outliers = [chr_nums[i] for i, z in enumerate(olig_z) if z > 2]
    
    # Z-score analysis for Neg_unique  
    neg_z = np.abs((np.array(neg_counts) - np.mean(neg_counts)) / np.std(neg_counts))
    neg_outliers = [chr_nums[i] for i, z in enumerate(neg_z) if z > 2]
    
    if olig_outliers:
        print(f"  Olig_unique outliers (>2σ): Chr{olig_outliers}")
    else:
        print(f"  Olig_unique: No significant outliers")
        
    if neg_outliers:
        print(f"  Neg_unique outliers (>2σ): Chr{neg_outliers}")
    else:
        print(f"  Neg_unique: No significant outliers")
    
    # Expected vs observed patterns
    print(f"\n🎯 Expected vs Observed Patterns:")
    
    # Check if larger chromosomes have more SNPs (expected)
    large_chrs = [1, 2, 3, 4, 5]  # Top 5 largest
    small_chrs = [18, 19, 20, 21, 22]  # Bottom 5 smallest
    
    olig_large_mean = np.mean([olig_unique[chr] for chr in large_chrs])
    olig_small_mean = np.mean([olig_unique[chr] for chr in small_chrs])
    neg_large_mean = np.mean([neg_unique[chr] for chr in large_chrs])
    neg_small_mean = np.mean([neg_unique[chr] for chr in small_chrs])
    
    print(f"  Large chromosomes (1-5):")
    print(f"    Olig_unique: {olig_large_mean:.1f} SNPs/chr")
    print(f"    Neg_unique:  {neg_large_mean:.1f} SNPs/chr")
    print(f"  Small chromosomes (18-22):")
    print(f"    Olig_unique: {olig_small_mean:.1f} SNPs/chr")
    print(f"    Neg_unique:  {neg_small_mean:.1f} SNPs/chr")
    
    # Reliability assessment
    print(f"\n🏆 Reliability Assessment:")
    
    reliability_score = 0
    max_score = 6
    
    # 1. Positive correlation with chromosome length (expected)
    if olig_corr > 0.3 and olig_p < 0.05:
        reliability_score += 1
        print(f"  ✅ Olig_unique correlates with chr length")
    else:
        print(f"  ❌ Olig_unique weak correlation with chr length")
    
    if neg_corr > 0.3 and neg_p < 0.05:
        reliability_score += 1
        print(f"  ✅ Neg_unique correlates with chr length")
    else:
        print(f"  ❌ Neg_unique weak correlation with chr length")
    
    # 2. Reasonable SNP counts (not too extreme)
    if 500 < olig_total < 5000:
        reliability_score += 1
        print(f"  ✅ Olig_unique has reasonable total count")
    else:
        print(f"  ❌ Olig_unique extreme total count")
    
    if 5000 < neg_total < 50000:
        reliability_score += 1
        print(f"  ✅ Neg_unique has reasonable total count")
    else:
        print(f"  ❌ Neg_unique extreme total count")
    
    # 3. Consistent density patterns
    if np.std(olig_densities) / np.mean(olig_densities) < 1.0:  # CV < 100%
        reliability_score += 1
        print(f"  ✅ Olig_unique has consistent density")
    else:
        print(f"  ❌ Olig_unique highly variable density")
    
    if np.std(neg_densities) / np.mean(neg_densities) < 1.0:
        reliability_score += 1
        print(f"  ✅ Neg_unique has consistent density")
    else:
        print(f"  ❌ Neg_unique highly variable density")
    
    # Final assessment
    reliability_pct = (reliability_score / max_score) * 100
    
    print(f"\n📋 Final Reliability Score: {reliability_score}/{max_score} ({reliability_pct:.1f}%)")
    
    if reliability_pct >= 83:
        print("🟢 HIGH RELIABILITY - Data patterns are consistent with biological expectations")
    elif reliability_pct >= 67:
        print("🟡 MODERATE RELIABILITY - Some concerns but generally acceptable")
    elif reliability_pct >= 50:
        print("🟠 LOW RELIABILITY - Multiple issues detected, use with caution")
    else:
        print("🔴 UNRELIABLE - Data patterns inconsistent with expectations")
    
    return {
        'olig_total': olig_total,
        'neg_total': neg_total,
        'olig_correlation': (olig_corr, olig_p),
        'neg_correlation': (neg_corr, neg_p),
        'inter_correlation': (inter_corr, inter_p),
        'reliability_score': reliability_score,
        'reliability_pct': reliability_pct
    }

if __name__ == "__main__":
    results = analyze_snp_distribution()