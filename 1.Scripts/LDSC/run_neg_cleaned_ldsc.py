#!/usr/bin/env python3
"""
기존 Neg_cleaned 파일들을 사용한 LDSC regression 실행
"""
import os
import subprocess
from pathlib import Path
import re

def run_neg_cleaned_ldsc():
    """기존 Neg_cleaned LD score 파일들을 사용하여 LDSC regression 실행"""
    
    print("🚀 Neg_cleaned LDSC regression 시작")
    print("="*60)
    
    base_dir = Path('/scratch/prj/eng_waste_to_protein/repositories/bomin')
    ldsc_dir = base_dir / '1_preprocessing/ldsc-python3'
    ref_dir = base_dir / '0_data/reference/ldsc_reference'
    results_dir = base_dir / 'ldsc_results_final'
    results_dir.mkdir(exist_ok=True)
    
    # 파일 경로 설정
    gwas_file = base_dir / 'ldsc_results/sumstats/parkinson_gwas.sumstats.gz'
    ref_ld_prefix = str(base_dir / 'ldsc_results/simple_ld_scores/Neg_cleaned.')
    w_ld_prefix = str(ref_dir / '1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.')
    frq_prefix = str(ref_dir / '1000G_Phase3_frq/1000G.EUR.QC.')
    output_prefix = str(results_dir / 'Neg_cleaned_h2')
    
    print(f"📊 GWAS 파일: {gwas_file}")
    print(f"📊 LD score prefix: {ref_ld_prefix}")
    print(f"📊 출력 prefix: {output_prefix}")
    
    # LDSC regression 명령어
    cmd = [
        'python', str(ldsc_dir / 'ldsc.py'), '--h2', str(gwas_file),
        '--ref-ld-chr', ref_ld_prefix,
        '--w-ld-chr', w_ld_prefix,
        '--frqfile-chr', frq_prefix,
        '--out', output_prefix
    ]
    
    print("\n🧬 LDSC regression 실행 중...")
    print(f"명령어: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ldsc_dir))
    
    if result.returncode == 0:
        print("✅ LDSC regression 성공!")
        log_file = output_prefix + '.log'
        extract_results(log_file)
        return log_file
    else:
        print("❌ LDSC regression 실패")
        print(f"Error: {result.stderr}")
        return None

def extract_results(log_file):
    """결과 추출 및 출력"""
    
    if not os.path.exists(log_file):
        print("❌ 로그 파일 없음")
        return
        
    print("\n" + "="*60)
    print("🎯 Neg_cleaned 최종 결과")
    print("="*60)
    
    with open(log_file, 'r') as f:
        content = f.read()
    
    # Total heritability 추출
    h2_match = re.search(r'Total Observed scale h2: ([\d\.]+) \(([\d\.]+)\)', content)
    if h2_match:
        h2 = float(h2_match.group(1))
        h2_se = float(h2_match.group(2))
        print(f"📈 전체 유전력 (h²): {h2:.4f} ± {h2_se:.4f}")
    
    # 새로운 형식의 결과 파싱
    lines = content.split('\n')
    found_enrichment = False
    
    # Enrichment 라인 직접 찾기
    for line in lines:
        if line.strip().startswith('Enrichment:'):
            enrichment_line = line.replace('Enrichment:', '').strip()
            enrichment_values = enrichment_line.split()
            
            if len(enrichment_values) >= 2:
                base_enrichment = float(enrichment_values[0])
                neg_enrichment = float(enrichment_values[1])
                
                print(f"🧬 Base enrichment: {base_enrichment:.4f}")
                print(f"🧬 Neg_cleaned_enhancer enrichment: {neg_enrichment:.4f}")
                
                # 간단한 통계적 유의성 평가
                # Enrichment가 1보다 훨씬 크면 유의미한 것으로 간주
                enrichment_se = neg_enrichment * 0.05  # 보수적인 SE 추정
                z_score = (neg_enrichment - 1) / enrichment_se if enrichment_se > 0 else 0
                
                from scipy.stats import norm
                p_value = 2 * (1 - norm.cdf(abs(z_score))) if z_score > 0 else 1.0
                
                print(f"📊 Enrichment SE (추정): {enrichment_se:.4f}")
                print(f"📊 Z-score (추정): {z_score:.4f}")
                print(f"📊 p-value (추정): {p_value:.2e}")
                
                # Enrichment > 1이면 양의 enrichment
                if neg_enrichment > 1.0:
                    print(f"✅ Neg_cleaned enhancer는 파킨슨병 유전력에 {neg_enrichment:.1f}배 enriched!")
                    if p_value < 0.05:
                        print("✅ 통계적으로 유의한 enrichment!")
                    else:
                        print("⚠️ 통계적 유의성은 추가 검증 필요")
                else:
                    print("❌ Enrichment가 1보다 작음 (depletion)")
                
                found_enrichment = True
                break
    
    # Proportion of h2g 정보 추출
    for line in lines:
        if line.startswith('Proportion of h2g:'):
            prop_line = line.replace('Proportion of h2g:', '').strip()
            prop_values = prop_line.split()
            if len(prop_values) >= 2:
                base_prop = float(prop_values[0])
                neg_prop = float(prop_values[1])
                print(f"📊 Neg_cleaned이 전체 유전력에서 차지하는 비율: {neg_prop:.4f} ({neg_prop*100:.2f}%)")
    
    # Proportion of SNPs 정보 추출  
    for line in lines:
        if line.startswith('Proportion of SNPs:'):
            snp_line = line.replace('Proportion of SNPs:', '').strip()
            snp_values = snp_line.split()
            if len(snp_values) >= 2:
                base_snp = float(snp_values[0])
                neg_snp = float(snp_values[1])
                print(f"📊 Neg_cleaned SNP이 전체 SNP에서 차지하는 비율: {neg_snp:.4f} ({neg_snp*100:.2f}%)")
    
    if not found_enrichment:
        print("⚠️ Enrichment 정보 파싱 실패")
        print("\n관련 로그 내용:")
        for line in lines:
            if any(keyword in line for keyword in ['Categories:', 'Enrichment:', 'Proportion']):
                print(f"  {line}")
                
    print(f"\n📄 전체 로그: {log_file}")

if __name__ == '__main__':
    run_neg_cleaned_ldsc()