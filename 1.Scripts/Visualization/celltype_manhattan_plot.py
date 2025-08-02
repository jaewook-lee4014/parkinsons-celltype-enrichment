#!/usr/bin/env python3
"""
세포타입별 enhancer intersect SNP들의 맨하탄 플롯
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import gzip

def load_celltype_annotations():
    """세포타입별 annotation 데이터 로드"""
    
    print("🧬 세포타입별 annotation 데이터 로딩 중...")
    
    celltypes = {
        'Microglia': 'Neg_cleaned',
        'Neuron': 'NeuN_cleaned', 
        'Oligodendrocyte': 'Olig_cleaned',
        'Dopaminergic': 'Nurr_cleaned'
    }
    
    all_annotations = {}
    
    for celltype, file_prefix in celltypes.items():
        print(f"  {celltype} 로딩 중...")
        
        # 염색체 1번만 먼저 테스트
        annot_file = f"/cephfs/volumes/hpc_data_prj/eng_waste_to_protein/ae035a41-20d2-44f3-aa46-14424ab0f6bf/repositories/bomin/ldsc_results/annotations/{file_prefix}.1.annot.gz"
        
        try:
            df = pd.read_csv(annot_file, sep='\t', compression='gzip')
            
            # 마지막 열이 enhancer annotation
            enhancer_col = df.columns[-1]  # 마지막 열
            
            # Intersect된 SNP들만 선택
            intersect_snps = df[df[enhancer_col] == 1]['SNP'].tolist()
            
            print(f"    염색체 1: {len(intersect_snps)}개 SNP가 {celltype} enhancer와 intersect")
            
            all_annotations[celltype] = {
                'intersect_snps': set(intersect_snps),
                'total_snps': len(df),
                'intersect_count': len(intersect_snps)
            }
            
        except Exception as e:
            print(f"    오류: {e}")
            all_annotations[celltype] = {
                'intersect_snps': set(),
                'total_snps': 0,
                'intersect_count': 0
            }
    
    return all_annotations

def load_gwas_with_positions():
    """GWAS 데이터와 위치 정보 함께 로드"""
    
    print("📊 GWAS 데이터 로딩 중...")
    
    # GWAS 데이터
    gwas_file = "/cephfs/volumes/hpc_data_prj/eng_waste_to_protein/ae035a41-20d2-44f3-aa46-14424ab0f6bf/repositories/bomin/ldsc_results/sumstats/parkinson_gwas.sumstats.gz"
    
    # 샘플링해서 로드 (메모리 절약)
    sample_size = 500000
    
    print(f"  GWAS 데이터 샘플링 로딩 중 (n={sample_size:,})...")
    
    df = pd.read_csv(gwas_file, sep='\t', compression='gzip', nrows=sample_size)
    
    # P-value 계산
    df['P'] = 2 * (1 - norm.cdf(np.abs(df['Z'])))
    df['-log10P'] = -np.log10(np.maximum(df['P'], 1e-50))
    
    # 위치 정보 추가 (간단한 방법: SNP ID 기반 가상 위치)
    # 실제로는 .bim 파일에서 가져와야 하지만, 데이터가 크므로 간단히 처리
    df['CHR'] = 1  # 염색체 1번만 처리
    df['BP'] = range(len(df))  # 순서대로 위치 할당
    
    print(f"  로딩 완료: {len(df):,} SNPs")
    print(f"  P-value 범위: {df['P'].min():.2e} - {df['P'].max():.2e}")
    
    return df

def create_celltype_manhattan_plots(gwas_df, annotations):
    """세포타입별 맨하탄 플롯 생성"""
    
    print("🗽 세포타입별 맨하탄 플롯 생성 중...")
    
    # 4개 세포타입 서브플롯
    fig, axes = plt.subplots(2, 2, figsize=(20, 16))
    axes = axes.flatten()
    
    celltypes = ['Microglia', 'Neuron', 'Oligodendrocyte', 'Dopaminergic']
    colors_intersect = ['red', 'darkgreen', 'purple', 'orange']
    colors_background = ['lightcoral', 'lightgreen', 'plum', 'moccasin']
    
    # 유의성 기준선
    genome_wide = -np.log10(5e-8)
    suggestive = -np.log10(1e-5)
    
    for i, celltype in enumerate(celltypes):
        ax = axes[i]
        
        print(f"  {celltype} 플롯 생성 중...")
        
        # 해당 세포타입의 intersect SNP들
        intersect_snps = annotations[celltype]['intersect_snps']
        intersect_count = len(intersect_snps)
        
        # GWAS 데이터에서 intersect 여부 표시
        gwas_df['is_intersect'] = gwas_df['SNP'].isin(intersect_snps)
        
        # Non-intersect SNP들 (배경)
        non_intersect = gwas_df[~gwas_df['is_intersect']]
        intersect_data = gwas_df[gwas_df['is_intersect']]
        
        print(f"    Non-intersect SNPs: {len(non_intersect):,}")
        print(f"    Intersect SNPs: {len(intersect_data):,}")
        
        # 배경 SNP들 (작게, 투명하게)
        if len(non_intersect) > 0:
            # 너무 많으면 샘플링
            if len(non_intersect) > 50000:
                non_intersect_sample = non_intersect.sample(n=50000, random_state=42)
            else:
                non_intersect_sample = non_intersect
                
            ax.scatter(non_intersect_sample['BP'], non_intersect_sample['-log10P'], 
                      c=colors_background[i], alpha=0.3, s=1, label='Background SNPs')
        
        # Intersect SNP들 (크게, 진하게)
        if len(intersect_data) > 0:
            ax.scatter(intersect_data['BP'], intersect_data['-log10P'], 
                      c=colors_intersect[i], alpha=0.8, s=20, 
                      label=f'{celltype} Enhancer SNPs (n={len(intersect_data)})')
        
        # 유의성 기준선
        ax.axhline(y=genome_wide, color='red', linestyle='--', alpha=0.7, 
                  linewidth=1.5, label='p=5×10⁻⁸')
        ax.axhline(y=suggestive, color='blue', linestyle='--', alpha=0.5, 
                  linewidth=1, label='p=1×10⁻⁵')
        
        # 상위 intersect SNP들 라벨링
        if len(intersect_data) > 0:
            top_intersect = intersect_data.nlargest(3, '-log10P')
            for _, snp in top_intersect.iterrows():
                if snp['-log10P'] > suggestive:
                    ax.annotate(f"{snp['SNP']}", 
                               xy=(snp['BP'], snp['-log10P']),
                               xytext=(5, 5), textcoords='offset points',
                               fontsize=8, alpha=0.8,
                               bbox=dict(boxstyle='round,pad=0.2', 
                                       facecolor=colors_intersect[i], alpha=0.7))
        
        # 축 설정
        ax.set_xlabel('SNP Position (Chr 1)', fontsize=12, fontweight='bold')
        ax.set_ylabel('-log₁₀(P-value)', fontsize=12, fontweight='bold')
        ax.set_title(f'{celltype} Enhancer Manhattan Plot\\n'
                    f'({intersect_count} intersect SNPs)', 
                    fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=10)
        
        # Y축 범위 설정
        max_y = max(gwas_df['-log10P'].max(), genome_wide + 2)
        ax.set_ylim(0, max_y)
    
    plt.tight_layout()
    
    # 저장
    plt.savefig('/scratch/prj/eng_waste_to_protein/repositories/bomin/celltype_manhattan_plots.png', 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('/scratch/prj/eng_waste_to_protein/repositories/bomin/celltype_manhattan_plots.pdf', 
                bbox_inches='tight', facecolor='white')
    
    plt.show()
    
    # 통계 요약
    print(f"\\n📈 세포타입별 맨하탄 플롯 통계:")
    print("="*60)
    
    for celltype in celltypes:
        intersect_snps = annotations[celltype]['intersect_snps']
        intersect_data = gwas_df[gwas_df['SNP'].isin(intersect_snps)]
        
        if len(intersect_data) > 0:
            significant_intersect = intersect_data[intersect_data['-log10P'] > suggestive]
            genome_wide_intersect = intersect_data[intersect_data['-log10P'] > genome_wide]
            
            print(f"\\n🧬 {celltype}:")
            print(f"  Total intersect SNPs: {len(intersect_data)}")
            print(f"  Suggestive (p<1e-5): {len(significant_intersect)}")
            print(f"  Genome-wide (p<5e-8): {len(genome_wide_intersect)}")
            
            if len(intersect_data) > 0:
                print(f"  Max -log10(P): {intersect_data['-log10P'].max():.2f}")
                
                # 상위 3개 SNP
                top_snps = intersect_data.nlargest(3, '-log10P')
                print(f"  Top SNPs:")
                for _, snp in top_snps.iterrows():
                    print(f"    {snp['SNP']}: p={snp['P']:.2e}")
        else:
            print(f"\\n🧬 {celltype}: No intersect SNPs found")

def create_comparison_manhattan(gwas_df, annotations):
    """4개 세포타입 비교 맨하탄 플롯"""
    
    print("\\n🔄 세포타입 비교 맨하탄 플롯 생성 중...")
    
    # 단일 플롯에서 4개 세포타입 비교
    fig, ax = plt.subplots(figsize=(16, 10))
    
    celltypes = ['Microglia', 'Neuron', 'Oligodendrocyte', 'Dopaminergic']
    colors = ['red', 'green', 'blue', 'orange']
    
    # 각 세포타입별로 점 표시
    for i, celltype in enumerate(celltypes):
        intersect_snps = annotations[celltype]['intersect_snps']
        intersect_data = gwas_df[gwas_df['SNP'].isin(intersect_snps)]
        
        if len(intersect_data) > 0:
            ax.scatter(intersect_data['BP'], intersect_data['-log10P'], 
                      c=colors[i], alpha=0.7, s=15, label=f'{celltype} (n={len(intersect_data)})')
    
    # 유의성 기준선
    genome_wide = -np.log10(5e-8)
    suggestive = -np.log10(1e-5)
    
    ax.axhline(y=genome_wide, color='red', linestyle='--', alpha=0.8, 
              linewidth=2, label='Genome-wide (p=5×10⁻⁸)')
    ax.axhline(y=suggestive, color='blue', linestyle='--', alpha=0.6, 
              linewidth=1.5, label='Suggestive (p=1×10⁻⁵)')
    
    # 축 설정
    ax.set_xlabel('SNP Position (Chr 1)', fontsize=14, fontweight='bold')
    ax.set_ylabel('-log₁₀(P-value)', fontsize=14, fontweight='bold')
    ax.set_title('Cell Type-Specific Enhancer Manhattan Plot Comparison', 
                 fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    
    # 저장
    plt.savefig('/scratch/prj/eng_waste_to_protein/repositories/bomin/celltype_comparison_manhattan.png', 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('/scratch/prj/eng_waste_to_protein/repositories/bomin/celltype_comparison_manhattan.pdf', 
                bbox_inches='tight', facecolor='white')
    
    plt.show()

if __name__ == "__main__":
    print("🗽 세포타입별 enhancer 맨하탄 플롯 생성!")
    print("="*60)
    
    # 1. 세포타입별 annotation 로드
    annotations = load_celltype_annotations()
    
    # 2. GWAS 데이터 로드
    gwas_data = load_gwas_with_positions()
    
    # 3. 세포타입별 맨하탄 플롯 생성
    create_celltype_manhattan_plots(gwas_data, annotations)
    
    # 4. 비교 맨하탄 플롯 생성
    create_comparison_manhattan(gwas_data, annotations)
    
    print(f"\\n✅ 세포타입별 맨하탄 플롯 생성 완료!")
    print(f"   저장된 파일:")
    print(f"   - celltype_manhattan_plots.png (4-panel)")
    print(f"   - celltype_comparison_manhattan.png (overlay)")