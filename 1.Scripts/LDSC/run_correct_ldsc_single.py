#!/usr/bin/env python3
"""
올바른 방법으로 단일 데이터셋 LDSC 실행
"""
import subprocess
import re
from pathlib import Path
from scipy.stats import norm

def generate_correct_ld_scores(dataset_name):
    """올바른 LD score 생성"""
    print(f"🧬 {dataset_name} 올바른 LD score 생성...")
    
    base_dir = Path('/scratch/prj/eng_waste_to_protein/repositories/bomin')
    ldsc_dir = base_dir / '1_preprocessing/ldsc-python3'
    annot_dir = base_dir / 'ldsc_results/correct_annotations'
    output_dir = base_dir / 'ldsc_results/correct_ld_scores'
    output_dir.mkdir(exist_ok=True)
    
    completed = 0
    
    # 먼저 몇 개 염색체만 테스트
    test_chromosomes = [1, 2, 3]
    
    for chr_num in test_chromosomes:
        output_prefix = str(output_dir / f'{dataset_name}.{chr_num}')
        
        # 이미 완료된 경우 건너뛰기
        if Path(f"{output_prefix}.l2.ldscore.gz").exists():
            print(f"  Chr {chr_num}: 이미 완료")
            completed += 1
            continue
            
        print(f"  Chr {chr_num} 처리 중...")
        
        # 입력 파일 경로
        annot_file = annot_dir / f'{dataset_name}.{chr_num}.annot.gz'
        bfile_prefix = str(base_dir / f'0_data/reference/ldsc_reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr_num}')
        
        if not annot_file.exists():
            print(f"    ❌ Annotation 파일 없음")
            continue
            
        if not Path(f"{bfile_prefix}.bed").exists():
            print(f"    ❌ Plink 파일 없음")
            continue
        
        # LDSC 명령어
        cmd = [
            'python', 'ldsc.py', '--l2',
            '--bfile', bfile_prefix,
            '--ld-wind-cm', '1',
            '--annot', str(annot_file),
            '--out', output_prefix,
            '--print-snps', str(base_dir / '0_data/reference/ldsc_reference/listHM3.txt')
        ]
        
        try:
            result = subprocess.run(cmd, cwd=str(ldsc_dir), 
                                  capture_output=True, text=True, timeout=600)
            
            if result.returncode == 0:
                completed += 1
                print(f"    ✅ 완료")
            else:
                print(f"    ❌ 실패: {result.stderr[:100]}")
                
        except subprocess.TimeoutExpired:
            print(f"    ⏰ 시간 초과")
        except Exception as e:
            print(f"    ❌ 예외: {e}")
    
    print(f"완료된 염색체: {completed}/{len(test_chromosomes)}")
    return completed == len(test_chromosomes)

def run_correct_ldsc_regression(dataset_name):
    """올바른 LDSC regression 실행"""
    print(f"\n📊 {dataset_name} 올바른 LDSC regression...")
    
    base_dir = Path('/scratch/prj/eng_waste_to_protein/repositories/bomin')
    ldsc_dir = base_dir / '1_preprocessing/ldsc-python3'
    
    # 파일 경로 설정
    gwas_file = str(base_dir / 'ldsc_results/sumstats/parkinson_gwas.sumstats.gz')
    ref_ld_prefix = str(base_dir / f'ldsc_results/correct_ld_scores/{dataset_name}.')
    w_ld_prefix = str(base_dir / '0_data/reference/ldsc_reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.')
    frq_prefix = str(base_dir / '0_data/reference/ldsc_reference/1000G_Phase3_frq/1000G.EUR.QC.')
    output_prefix = str(base_dir / f'ldsc_results_correct/{dataset_name}_h2')
    
    # 출력 디렉토리 생성
    Path(output_prefix).parent.mkdir(exist_ok=True)
    
    # LD score 파일 확인 (테스트용 3개 염색체)
    ld_files = list(Path(base_dir / 'ldsc_results/correct_ld_scores').glob(f'{dataset_name}.*.l2.ldscore.gz'))
    print(f"사용 가능한 LD score 파일: {len(ld_files)}개")
    
    if len(ld_files) < 3:
        print(f"❌ LD score 파일 부족")
        return None
    
    # LDSC 명령어
    cmd = [
        'python', 'ldsc.py', '--h2', gwas_file,
        '--ref-ld-chr', ref_ld_prefix,
        '--w-ld-chr', w_ld_prefix,
        '--frqfile-chr', frq_prefix,
        '--out', output_prefix
    ]
    
    try:
        result = subprocess.run(cmd, cwd=str(ldsc_dir), 
                              capture_output=True, text=True, timeout=900)
        
        if result.returncode == 0:
            print(f"✅ LDSC regression 성공")
            return output_prefix + '.log'
        else:
            print(f"❌ LDSC regression 실패")
            print(f"Error: {result.stderr}")
            return None
            
    except Exception as e:
        print(f"❌ 예외 발생: {e}")
        return None

def extract_correct_results(log_file, dataset_name):
    """98번째 annotation (우리 enhancer) 결과 추출"""
    if not log_file or not Path(log_file).exists():
        return None
        
    with open(log_file, 'r') as f:
        content = f.read()
    
    results = {'dataset': dataset_name}
    lines = content.split('\n')
    
    # Total heritability
    h2_match = re.search(r'Total Observed scale h2: ([\d\.]+) \(([\d\.]+)\)', content)
    if h2_match:
        results['total_h2'] = float(h2_match.group(1))
        results['total_h2_se'] = float(h2_match.group(2))
    
    # annotation 개수 확인
    categories_line = None
    for line in lines:
        if line.startswith('Categories:'):
            categories_line = line
            categories = line.replace('Categories:', '').strip().split()
            results['n_annotations'] = len(categories)
            print(f"  총 annotation 수: {len(categories)}")
            break
    
    # 98번째 (마지막) annotation = 우리 enhancer
    for line in lines:
        if line.strip().startswith('Enrichment:'):
            enrichment_values = line.replace('Enrichment:', '').strip().split()
            if len(enrichment_values) >= 98:
                results['enhancer_enrichment'] = float(enrichment_values[-1])  # 마지막 값
                print(f"  Enhancer enrichment: {enrichment_values[-1]}")
    
    for line in lines:
        if line.startswith('Proportion of h2g:'):
            prop_values = line.replace('Proportion of h2g:', '').strip().split()
            if len(prop_values) >= 98:
                results['h2g_proportion'] = float(prop_values[-1])
                
        elif line.startswith('Proportion of SNPs:'):
            snp_values = line.replace('Proportion of SNPs:', '').strip().split()
            if len(snp_values) >= 98:
                results['snp_proportion'] = float(snp_values[-1])
    
    # Coefficient 및 SE (마지막 값)
    for line in lines:
        if line.startswith('Coefficients:'):
            coeff_values = line.replace('Coefficients:', '').strip().split()
            if len(coeff_values) >= 98:
                results['enhancer_coeff'] = float(coeff_values[-1])
                
        elif line.startswith('Coefficient SE:'):
            se_values = line.replace('Coefficient SE:', '').strip().split()
            if len(se_values) >= 98:
                results['enhancer_coeff_se'] = float(se_values[-1])
    
    # p-value 계산
    if 'enhancer_coeff' in results and 'enhancer_coeff_se' in results:
        z_score = results['enhancer_coeff'] / results['enhancer_coeff_se']
        results['z_score'] = z_score
        results['p_value'] = 2 * (1 - norm.cdf(abs(z_score)))
    
    return results

def main():
    """메인 실행"""
    dataset = 'Neg_cleaned'
    
    print("🔬 올바른 LDSC 방법 테스트")
    print("="*50)
    
    # 1단계: LD scores 생성
    ld_success = generate_correct_ld_scores(dataset)
    
    if not ld_success:
        print(f"❌ LD score 생성 실패")
        return
    
    # 2단계: LDSC regression
    log_file = run_correct_ldsc_regression(dataset)
    
    if not log_file:
        print(f"❌ LDSC regression 실패")
        return
    
    # 3단계: 결과 추출
    results = extract_correct_results(log_file, dataset)
    
    if results:
        print(f"\n🎯 {dataset} 수정된 결과:")
        print("-" * 40)
        if 'n_annotations' in results:
            print(f"사용된 annotation 수: {results['n_annotations']}")
        if 'enhancer_enrichment' in results:
            print(f"Enrichment: {results['enhancer_enrichment']:.1f}x")
        if 'p_value' in results:
            print(f"p-value: {results['p_value']:.2e}")
        if 'h2g_proportion' in results:
            print(f"유전력 기여도: {results['h2g_proportion']*100:.2f}%")
        if 'snp_proportion' in results:
            print(f"SNP 비율: {results['snp_proportion']*100:.3f}%")
            
        print(f"\n✅ 성공! 이제 전체 22개 염색체로 확장 가능")
    else:
        print(f"❌ 결과 추출 실패")

if __name__ == '__main__':
    main()