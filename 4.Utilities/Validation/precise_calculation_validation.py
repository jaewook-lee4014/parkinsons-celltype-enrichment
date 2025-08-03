#!/usr/bin/env python3
"""
Enrichment와 P-value 계산 과정의 정밀한 검증
"""
import numpy as np
from scipy import stats
import pandas as pd
from pathlib import Path

def validate_current_calculations():
    """현재 사용 중인 계산 과정 검증"""
    
    print("🔬 현재 Enrichment & P-value 계산 과정 정밀 검증")
    print("="*80)
    
    # 현재 quick_correct_visualization.py에서 사용하는 데이터
    realistic_enrichment = {
        'Microglia': {
            'unique': 68.3, 'all': 23.7,
            'unique_se': 12.5, 'all_se': 4.2
        },
        'Neuron': {
            'unique': 45.2, 'all': 31.5,
            'unique_se': 8.7, 'all_se': 5.3
        },
        'Oligodendrocyte': {
            'unique': 112.7, 'all': 18.9,
            'unique_se': 25.3, 'all_se': 3.8
        },
        'Dopaminergic': {
            'unique': 156.4, 'all': 12.3,
            'unique_se': 42.1, 'all_se': 2.9
        }
    }
    
    p_values = {
        'Microglia': {'unique_p': 0.0000038, 'all_p': 0.012},
        'Neuron': {'unique_p': 0.001, 'all_p': 0.008},
        'Oligodendrocyte': {'unique_p': 0.003, 'all_p': 0.015},
        'Dopaminergic': {'unique_p': 0.02, 'all_p': 0.035}
    }
    
    print("\n1️⃣ Enrichment 값 검증")
    print("-" * 60)
    print("❓ 문제점: 이 값들이 어디서 왔는가?")
    print("   - 실제 LDSC 결과 파일에서 추출? ❌")
    print("   - 임의로 설정한 값? ✅")
    print("   → 현재는 가상의 '현실적인' 값을 사용 중")
    
    print("\n2️⃣ P-value 계산 검증")
    print("-" * 60)
    
    for cell_type in realistic_enrichment.keys():
        enr = realistic_enrichment[cell_type]
        p_vals = p_values[cell_type]
        
        print(f"\n{cell_type}:")
        
        # Z-score 역계산
        unique_z = stats.norm.ppf(1 - p_vals['unique_p']/2)
        all_z = stats.norm.ppf(1 - p_vals['all_p']/2)
        
        # Enrichment에서 SE를 사용한 Z-score 계산
        unique_z_calc = (enr['unique'] - 1) / enr['unique_se']
        all_z_calc = (enr['all'] - 1) / enr['all_se']
        
        print(f"  Unique: Enrichment={enr['unique']:.1f}±{enr['unique_se']:.1f}")
        print(f"    주어진 p={p_vals['unique_p']:.6f} → Z={unique_z:.2f}")
        print(f"    Enrichment로 계산한 Z={(enr['unique']-1)/enr['unique_se']:.2f}")
        print(f"    불일치! ❌")
        
        print(f"  All: Enrichment={enr['all']:.1f}±{enr['all_se']:.1f}")
        print(f"    주어진 p={p_vals['all_p']:.3f} → Z={all_z:.2f}")
        print(f"    Enrichment로 계산한 Z={(enr['all']-1)/enr['all_se']:.2f}")
        print(f"    불일치! ❌")

def analyze_ldsc_calculation_process():
    """실제 LDSC 계산 과정 분석"""
    
    print("\n\n📚 실제 LDSC Enrichment 계산 과정")
    print("="*80)
    
    print("\n1. LDSC의 Enrichment 정의:")
    print("   Enrichment = τ_c / (M_c/M)")
    print("   여기서:")
    print("   - τ_c: category c의 per-SNP heritability")
    print("   - M_c: category c의 SNP 개수")
    print("   - M: 전체 SNP 개수")
    
    print("\n2. 통계적 검정:")
    print("   - τ_c는 regression coefficient로 추정")
    print("   - Standard error는 jackknife 방법으로 계산")
    print("   - P-value는 τ_c = 0 귀무가설 검정")
    
    print("\n3. Enrichment P-value:")
    print("   - Enrichment = 1 검정이 아님!")
    print("   - τ_c = 0 검정 (coefficient가 0인지)")
    print("   - 따라서 enrichment와 p-value가 직접 연결되지 않음")

def calculate_correct_values():
    """올바른 계산 방법 제시"""
    
    print("\n\n✅ 올바른 계산 방법")
    print("="*80)
    
    # 가상의 실제 LDSC 출력 예시
    print("\n예시: 실제 LDSC 결과 형식")
    print("-" * 60)
    
    example_results = {
        'Microglia_unique': {
            'Prop_SNPs': 0.00092,
            'Prop_h2': 0.0496,
            'Enrichment': 54.15,
            'Enrichment_std_error': 9.87,
            'Coefficient': 0.0215,
            'Coefficient_std_error': 0.0054,
            'Coefficient_z_score': 3.98,
            'Coefficient_p_value': 0.0000682
        }
    }
    
    res = example_results['Microglia_unique']
    
    print(f"Category: Microglia_unique")
    print(f"Proportion of SNPs: {res['Prop_SNPs']:.5f}")
    print(f"Proportion of h2: {res['Prop_h2']:.4f}")
    print(f"Enrichment: {res['Enrichment']:.2f} ({res['Enrichment_std_error']:.2f})")
    print(f"Coefficient: {res['Coefficient']:.4f} ({res['Coefficient_std_error']:.4f})")
    print(f"Coefficient z-score: {res['Coefficient_z_score']:.2f}")
    print(f"Coefficient p-value: {res['Coefficient_p_value']:.2e}")
    
    print("\n중요: P-value는 coefficient 검정에서 나옴!")
    print("      Enrichment 자체의 유의성이 아님!")

def propose_accurate_data():
    """정확한 데이터 생성 방법 제안"""
    
    print("\n\n🎯 정확한 시각화를 위한 제안")
    print("="*80)
    
    print("\n1. 실제 LDSC 실행:")
    print("   - 각 cell type의 unique/all annotation")
    print("   - Parkinson's GWAS summary statistics")
    print("   - Baseline model과 함께 분석")
    
    print("\n2. 결과 파일에서 추출:")
    print("   - .results 파일의 Enrichment, Enrichment_std_error")
    print("   - .results 파일의 Coefficient_p_value")
    
    print("\n3. 현재 사용 중인 가상 데이터의 문제:")
    print("   - Enrichment와 SE가 p-value와 수학적으로 일치하지 않음")
    print("   - P-value가 coefficient 검정이 아닌 임의 값")
    
    # 더 정확한 가상 데이터 생성
    print("\n4. 개선된 가상 데이터 (수학적으로 일관성 있게):")
    
    improved_data = {}
    
    for cell_type, params in [
        ('Microglia', {'unique_enr': 68.3, 'unique_coef': 0.025, 'unique_coef_se': 0.005}),
        ('Neuron', {'unique_enr': 45.2, 'unique_coef': 0.018, 'unique_coef_se': 0.006}),
        ('Oligodendrocyte', {'unique_enr': 112.7, 'unique_coef': 0.032, 'unique_coef_se': 0.008}),
        ('Dopaminergic', {'unique_enr': 156.4, 'unique_coef': 0.041, 'unique_coef_se': 0.012})
    ]:
        z_score = params['unique_coef'] / params['unique_coef_se']
        p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
        
        improved_data[cell_type] = {
            'enrichment': params['unique_enr'],
            'coefficient': params['unique_coef'],
            'coefficient_se': params['unique_coef_se'],
            'z_score': z_score,
            'p_value': p_value
        }
        
        print(f"\n{cell_type} Unique:")
        print(f"  Enrichment: {params['unique_enr']:.1f}")
        print(f"  Coefficient: {params['unique_coef']:.4f} ± {params['unique_coef_se']:.4f}")
        print(f"  Z-score: {z_score:.2f}")
        print(f"  P-value: {p_value:.2e}")

def identify_critical_issues():
    """중요한 문제점 식별"""
    
    print("\n\n⚠️ 현재 계산의 중요 문제점")
    print("="*80)
    
    print("\n1. Enrichment와 P-value의 불일치:")
    print("   - Enrichment SE로 계산한 Z-score와")
    print("   - P-value로 역계산한 Z-score가 다름")
    
    print("\n2. P-value의 의미 혼동:")
    print("   - LDSC에서 p-value는 coefficient = 0 검정")
    print("   - Enrichment = 1 검정이 아님")
    
    print("\n3. 가상 데이터의 한계:")
    print("   - 실제 LDSC 결과가 아닌 추정값")
    print("   - 수학적 일관성 부족")
    
    print("\n4. 해결 방안:")
    print("   - 실제 LDSC 실행 결과 사용")
    print("   - 또는 수학적으로 일관된 가상 데이터 생성")

if __name__ == "__main__":
    # 1. 현재 계산 검증
    validate_current_calculations()
    
    # 2. LDSC 계산 과정 분석
    analyze_ldsc_calculation_process()
    
    # 3. 올바른 계산 방법
    calculate_correct_values()
    
    # 4. 정확한 데이터 제안
    propose_accurate_data()
    
    # 5. 중요 문제점 식별
    identify_critical_issues()
    
    print("\n\n" + "="*80)
    print("📌 결론:")
    print("현재 시각화의 enrichment와 p-value는 수학적으로 일치하지 않습니다.")
    print("실제 LDSC 결과를 사용하거나, 수학적으로 일관된 데이터를 생성해야 합니다.")
    print("="*80)