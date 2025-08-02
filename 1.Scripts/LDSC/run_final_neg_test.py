#!/usr/bin/env python3
"""
최종 Neg_cleaned LDSC regression 실행 및 결과 분석
"""
import subprocess
import os
import re

def fix_all_m_files():
    """모든 염색체의 M 파일을 정확히 수정"""
    
    print('🔧 모든 염색체 M 파일 수정 중...')
    
    for chr_num in range(1, 23):
        baseline_m_file = f'/scratch/prj/eng_waste_to_protein/repositories/bomin/0_data/reference/ldsc_reference/baselineLD.{chr_num}.l2.M'
        celltype_m_file = f'/scratch/prj/eng_waste_to_protein/repositories/bomin/ldsc_results/simple_ld_scores/Neg_cleaned.{chr_num}.l2.M'
        output_m_file = f'/scratch/prj/eng_waste_to_protein/repositories/bomin/ldsc_results/test_combined_ld_scores/Neg_cleaned.{chr_num}.l2.M'
        
        if not os.path.exists(baseline_m_file) or not os.path.exists(celltype_m_file):
            print(f'  ⚠️ Chr{chr_num}: 파일 없음')
            continue
            
        # BaselineLD M 읽기 (97개 값)
        with open(baseline_m_file, 'r') as f:
            baseline_values = f.read().strip().split('\t')
        
        # Neg_cleaned M에서 enhancer 값 추출
        with open(celltype_m_file, 'r') as f:
            lines = f.read().strip().split('\n')
            if len(lines) >= 2:
                enhancer_value = lines[1]  # 두 번째 줄
            else:
                enhancer_value = lines[0].split('\t')[1] if '\t' in lines[0] else lines[0]
        
        # 결합: 97개 + 1개 = 98개
        combined_values = baseline_values + [enhancer_value]
        
        # 한 줄로 저장
        with open(output_m_file, 'w') as f:
            f.write('\t'.join(combined_values) + '\n')
        
        print(f'{chr_num}', end=' ', flush=True)
    
    print('\n✅ 모든 M 파일 수정 완료')

def run_ldsc_regression():
    """LDSC regression 실행"""
    
    print('🧬 LDSC regression 실행 중...')
    
    ldsc_path = '/scratch/prj/eng_waste_to_protein/repositories/bomin/1_preprocessing/ldsc-python3/ldsc.py'
    gwas_file = '/scratch/prj/eng_waste_to_protein/repositories/bomin/ldsc_results/sumstats/parkinson_gwas_fixed.sumstats.gz'
    ref_ld_prefix = '/scratch/prj/eng_waste_to_protein/repositories/bomin/ldsc_results/test_combined_ld_scores/Neg_cleaned.'
    w_ld_prefix = '/scratch/prj/eng_waste_to_protein/repositories/bomin/0_data/reference/ldsc_reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.'
    frq_prefix = '/scratch/prj/eng_waste_to_protein/repositories/bomin/0_data/reference/ldsc_reference/1000G_Phase3_frq/1000G.EUR.QC.'
    output_prefix = '/scratch/prj/eng_waste_to_protein/repositories/bomin/ldsc_results/test_results/Neg_cleaned_FINAL'
    
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    
    cmd = [
        'python', ldsc_path, '--h2', gwas_file,
        '--ref-ld-chr', ref_ld_prefix,
        '--w-ld-chr', w_ld_prefix,
        '--overlap-annot',
        '--frqfile-chr', frq_prefix,
        '--out', output_prefix
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, 
                          cwd='/scratch/prj/eng_waste_to_protein/repositories/bomin/1_preprocessing/ldsc-python3')
    
    if result.returncode == 0:
        print('✅ LDSC regression 성공!')
        return output_prefix + '.log'
    else:
        print('❌ LDSC regression 실패')
        print(f'Error: {result.stderr[:1000]}')
        return None

def extract_results(log_file):
    """로그 파일에서 enrichment와 p-value 추출"""
    
    if not log_file or not os.path.exists(log_file):
        print('❌ 로그 파일이 없습니다.')
        return
    
    print('📊 결과 분석 중...')
    
    with open(log_file, 'r') as f:
        log_content = f.read()
    
    # 결과 추출
    print('\n' + '='*60)
    print('🧬 파킨슨병 GWAS - Neg_cleaned enhancer 분석 결과')
    print('='*60)
    
    # Total heritability
    h2_match = re.search(r'Total Observed scale h2: ([\d\.]+) \(([\d\.]+)\)', log_content)
    if h2_match:
        h2_value = h2_match.group(1)
        h2_se = h2_match.group(2)
        print(f'📈 전체 유전력 (h²): {h2_value} ± {h2_se}')
    
    # Neg_cleaned enhancer enrichment 찾기
    lines = log_content.split('\n')
    for i, line in enumerate(lines):
        if 'Neg_cleaned_enhancer' in line and ('Enrichment:' in line or 'Coefficient:' in line):
            print(f'🎯 {line.strip()}')
        elif 'Categories:' in line:
            # 카테고리별 결과 출력
            print(f'\n📋 카테고리별 결과:')
            for j in range(1, min(10, len(lines) - i)):
                next_line = lines[i + j].strip()
                if next_line and not next_line.startswith('Analysis') and 'L2' in next_line:
                    if 'Neg_cleaned_enhancer' in next_line:
                        print(f'🎯 **{next_line}**')
                    elif j <= 5:  # 처음 몇 개만 출력
                        print(f'   {next_line}')
    
    print(f'\n📄 전체 로그: {log_file}')
    return True

def main():
    """메인 실행"""
    
    print('🚀 최종 Neg_cleaned LDSC 분석 시작')
    
    # 1. M 파일 수정
    fix_all_m_files()
    
    # 2. LDSC regression 실행
    log_file = run_ldsc_regression()
    
    # 3. 결과 추출
    if log_file:
        extract_results(log_file)
    else:
        print('❌ 분석 실패')

if __name__ == '__main__':
    main()