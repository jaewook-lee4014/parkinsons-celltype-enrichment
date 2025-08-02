# LDSC Cell-Type Specific Enhancer Enrichment Analysis Workflow

## 📋 Overview

이 문서는 LDSC (Linkage Disequilibrium Score Regression)를 사용하여 파킨슨병 GWAS 데이터에서 세포타입별 enhancer 영역의 heritability enrichment를 분석하는 전체 workflow를 정리합니다.

## 🎯 연구 목표

- **주요 목표**: 뇌 세포타입별 enhancer 영역이 파킨슨병 유전력에 기여하는 정도 정량화
- **가설**: 특정 뇌 세포타입의 enhancer 영역에 파킨슨병 관련 변이가 enriched되어 있을 것
- **기대 결과**: 세포타입별 enrichment 값과 통계적 유의성 확인

## 📊 데이터셋

### Input 데이터
1. **GWAS 요약통계**: 파킨슨병 GWAS 데이터 (`parkinson_gwas.sumstats.gz`)
2. **세포타입별 enhancer BED 파일**: 8개 데이터셋
   - Neg_cleaned, Neg_unique
   - NeuN_cleaned, NeuN_unique  
   - Nurr_cleaned, Nurr_unique
   - Olig_cleaned, Olig_unique
3. **Reference Panel**: BaselineLD v2.2 (97개 functional annotations)
4. **1000G EUR Phase 3**: LD 계산 및 weight files

### 데이터 구조
```
/scratch/prj/eng_waste_to_protein/repositories/bomin/
├── 0_data/
│   ├── raw/cleaned_data/           # 원본 BED 파일들
│   └── reference/ldsc_reference/   # BaselineLD, 1000G 파일들
├── ldsc_results/
│   ├── annotations/               # 생성된 annotation 파일들
│   ├── simple_ld_scores/         # LD score 파일들
│   └── sumstats/                 # 전처리된 GWAS 데이터
└── ldsc_results_final/           # 최종 결과 파일들
```

## 🔬 방법론

### 1. LDSC 이론적 배경

LDSC는 다음 선형 회귀 모델을 사용합니다:

```
E[χ²ⱼ] = 1 + N × Σₖ τₖ × lₖⱼ
```

**변수 설명:**
- `χ²ⱼ`: j번째 SNP의 GWAS 카이제곱 통계량
- `N`: GWAS 샘플 크기
- `τₖ`: k번째 annotation의 per-SNP heritability coefficient
- `lₖⱼ`: j번째 SNP의 k번째 annotation에 대한 LD score

### 2. Enrichment 계산 공식

```
Enrichment = (Proportion of h²g) / (Proportion of SNPs)
```

**의미:**
- 해당 annotation이 설명하는 유전력 비율을 SNP 비율로 나눈 값
- 1보다 크면 평균보다 중요한 영역 (enriched)
- 1보다 작으면 평균보다 덜 중요한 영역 (depleted)

### 3. 통계적 유의성 검정

**Z-score 계산:**
```
Z = coefficient / coefficient_SE
```

**p-value 계산:**
```
p-value = 2 × (1 - Φ(|Z|))
```
여기서 Φ는 표준정규분포 누적분포함수

## 🛠 단계별 Workflow

### Step 1: 환경 설정 및 데이터 준비

```bash
# 작업 디렉토리 설정
cd /scratch/prj/eng_waste_to_protein/repositories/bomin

# LDSC Python3 버전 확인
cd 1_preprocessing/ldsc-python3
python ldsc.py --help
```

### Step 2: BaselineLD 파일 검증

```python
def verify_baseline_files():
    """BaselineLD 파일 형식 및 개수 검증"""
    ref_dir = Path('0_data/reference/ldsc_reference')
    
    # 파일 존재 확인
    assert len(list(ref_dir.glob('baselineLD.*.l2.ldscore.gz'))) == 22
    assert len(list(ref_dir.glob('baselineLD.*.l2.M'))) == 22
    
    # 첫 번째 염색체로 형식 검증
    with gzip.open(ref_dir / 'baselineLD.1.l2.ldscore.gz', 'rt') as f:
        header = f.readline().strip().split('\t')
        n_annot = len(header) - 3  # CHR, SNP, BP 제외
    
    # 97개 BaselineLD annotation 확인
    assert n_annot == 97
    return n_annot
```

### Step 3: 세포타입별 Annotation 생성

**이미 생성된 annotation 파일 사용** (좌표계 문제 해결됨):

```python
def use_existing_annotations(dataset_name):
    """기존 생성된 annotation 파일 사용"""
    annot_dir = Path('ldsc_results/annotations')
    
    # 22개 염색체 annotation 파일 확인
    for chr_num in range(1, 23):
        annot_file = annot_dir / f'{dataset_name}.{chr_num}.annot.gz'
        assert annot_file.exists()
    
    return True
```

**Annotation 파일 구조:**
```
CHR  BP    SNP         BaselineLD_annotations...  CellType_enhancer
1    11008 rs575272151 [97 columns]              0
1    15274 rs62635286  [97 columns]              1
```

### Step 4: LD Score 계산

**기존 LD score 파일 사용** (이미 계산 완료):

```bash
# LD score 파일 확인
ls ldsc_results/simple_ld_scores/Neg_cleaned.*.l2.ldscore.gz
ls ldsc_results/simple_ld_scores/Neg_cleaned.*.l2.M
ls ldsc_results/simple_ld_scores/Neg_cleaned.*.l2.M_5_50
```

### Step 5: LDSC Regression 실행

```python
def run_ldsc_regression(dataset_name):
    """LDSC regression 실행"""
    
    # 파일 경로 설정
    gwas_file = 'ldsc_results/sumstats/parkinson_gwas.sumstats.gz'
    ref_ld_prefix = f'ldsc_results/simple_ld_scores/{dataset_name}.'
    w_ld_prefix = '0_data/reference/ldsc_reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.'
    frq_prefix = '0_data/reference/ldsc_reference/1000G_Phase3_frq/1000G.EUR.QC.'
    output_prefix = f'ldsc_results_final/{dataset_name}_h2'
    
    # LDSC 명령어 실행
    cmd = [
        'python', 'ldsc.py', '--h2', gwas_file,
        '--ref-ld-chr', ref_ld_prefix,
        '--w-ld-chr', w_ld_prefix,
        '--frqfile-chr', frq_prefix,
        '--out', output_prefix
    ]
    
    result = subprocess.run(cmd, cwd='1_preprocessing/ldsc-python3')
    return output_prefix + '.log' if result.returncode == 0 else None
```

### Step 6: 결과 추출 및 분석

```python
def extract_enrichment_results(log_file, dataset_name):
    """LDSC 결과에서 enrichment 정보 추출"""
    
    with open(log_file, 'r') as f:
        content = f.read()
    
    results = {}
    
    # Total heritability 추출
    h2_match = re.search(r'Total Observed scale h2: ([\d\.]+) \(([\d\.]+)\)', content)
    if h2_match:
        results['total_h2'] = float(h2_match.group(1))
        results['total_h2_se'] = float(h2_match.group(2))
    
    # Enrichment 값 추출
    lines = content.split('\n')
    for line in lines:
        if line.strip().startswith('Enrichment:'):
            enrichment_values = line.replace('Enrichment:', '').strip().split()
            if len(enrichment_values) >= 2:
                results['base_enrichment'] = float(enrichment_values[0])
                results['enhancer_enrichment'] = float(enrichment_values[1])
    
    # Proportion 정보 추출
    for line in lines:
        if line.startswith('Proportion of h2g:'):
            prop_values = line.replace('Proportion of h2g:', '').strip().split()
            if len(prop_values) >= 2:
                results['h2g_proportion'] = float(prop_values[1])
                
        elif line.startswith('Proportion of SNPs:'):
            snp_values = line.replace('Proportion of SNPs:', '').strip().split()
            if len(snp_values) >= 2:
                results['snp_proportion'] = float(snp_values[1])
    
    # Coefficient 및 SE 추출 (p-value 계산용)
    for line in lines:
        if line.startswith('Coefficients:'):
            coeff_values = line.replace('Coefficients:', '').strip().split()
            if len(coeff_values) >= 2:
                results['enhancer_coeff'] = float(coeff_values[1])
                
        elif line.startswith('Coefficient SE:'):
            se_values = line.replace('Coefficient SE:', '').strip().split()
            if len(se_values) >= 2:
                results['enhancer_coeff_se'] = float(se_values[1])
    
    # p-value 계산
    if 'enhancer_coeff' in results and 'enhancer_coeff_se' in results:
        z_score = results['enhancer_coeff'] / results['enhancer_coeff_se']
        from scipy.stats import norm
        results['p_value'] = 2 * (1 - norm.cdf(abs(z_score)))
    
    return results
```

## 📈 결과 해석

### 실제 분석 결과 (Neg_cleaned vs Neg_unique)

| 항목 | Neg_cleaned | Neg_unique | 해석 |
|------|-------------|------------|------|
| **Enrichment** | 53.3배 | 54.1배 | 둘 다 매우 강한 enrichment |
| **p-value** | 5.79e-06 | 2.05e-03 | 둘 다 통계적으로 유의함 |
| **유전력 기여도** | 14.86% | 4.96% | Neg_cleaned가 더 포괄적 |
| **SNP 비율** | 0.28% | 0.092% | Neg_unique가 더 선별적 |

### 생물학적 해석

1. **강한 Enrichment (50배 이상)**:
   - 뇌 enhancer 영역에 파킨슨병 관련 변이가 고도로 집중
   - 일반적인 enhancer enrichment (3-10배)보다 훨씬 강함

2. **Neg_cleaned vs Neg_unique**:
   - **Neg_cleaned**: 더 포괄적이고 안정적인 신호
   - **Neg_unique**: 더 선별적이고 특이적인 신호

3. **통계적 신뢰성**:
   - 두 데이터셋 모두 p < 0.01로 유의함
   - Neg_cleaned가 더 강한 통계적 증거 제공

## 🔧 주요 기술적 해결사항

### 1. 좌표계 문제 해결
- **문제**: BED 파일과 BaselineLD SNP 간 좌표 불일치
- **해결**: 기존 생성된 annotation 파일 사용 (좌표 매핑 완료)

### 2. 파일 형식 호환성
- **문제**: LDSC annotation 파일 형식 요구사항
- **해결**: BaselineLD (97개) + enhancer (1개) = 98개 annotation으로 통합

### 3. GWAS 데이터 호환성
- **문제**: SNP ID 형식 불일치 (chr:pos vs rsID)
- **해결**: rsID 형식의 GWAS 파일 사용

### 4. 통계적 검증
- **문제**: 50배 이상의 높은 enrichment 값 검증
- **해결**: BaselineLD 단독 분석과 비교하여 상대적 안정성 확인

## 📝 품질 관리 체크리스트

### 데이터 검증
- [ ] BaselineLD 파일 22개 염색체 모두 존재
- [ ] Annotation 파일 형식 및 컬럼 수 확인 (98개)
- [ ] LD score 파일 생성 완료 확인
- [ ] GWAS 파일 SNP ID 형식 확인

### 분석 검증
- [ ] LDSC regression 성공적 완료
- [ ] Total heritability 합리적 범위 (0.01-0.02)
- [ ] Enrichment 값 추출 성공
- [ ] p-value 계산 정확성 확인

### 결과 검증
- [ ] Enrichment > 1 (예상 방향)
- [ ] 통계적 유의성 확인 (p < 0.05)
- [ ] 다른 세포타입과 비교 가능
- [ ] 생물학적 해석 타당성

## 🔄 확장 가능성

### 추가 분석
1. **나머지 6개 세포타입 분석**: NeuN, Nurr, Olig 각각의 cleaned/unique
2. **세포타입별 비교**: 어떤 세포타입이 가장 강한 enrichment를 보이는가?
3. **Conditional 분석**: 여러 세포타입을 동시에 고려한 분석
4. **다른 질병과 비교**: 알츠하이머, 헌팅턴병 등과의 비교

### 방법론 개선
1. **더 정교한 annotation**: 세포타입별 특이성 향상
2. **Multi-trait 분석**: 여러 파킨슨병 관련 phenotype 동시 분석
3. **Functional validation**: 실험적 검증을 위한 후보 영역 제시

## 📚 참고문헌 및 도구

### 주요 도구
- **LDSC**: Bulik-Sullivan et al. (2015) Nature Genetics
- **BaselineLD v2.2**: Gazal et al. (2017) Nature Genetics
- **1000 Genomes Phase 3**: European ancestry reference panel

### 분석 환경
- **Python**: 3.9+
- **필수 패키지**: pandas, numpy, scipy, pathlib
- **시스템**: Linux HPC environment

---

**문서 작성일**: 2025-07-30
**분석 완료**: Neg_cleaned, Neg_unique
**다음 단계**: 나머지 6개 세포타입 분석 진행