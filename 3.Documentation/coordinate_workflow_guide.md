# 좌표계 변환 워크플로우 가이드

## ⚠️ 중요: 좌표계 불일치 문제

현재 분석에서 `Enrichment ratio: 0.000`이 나오는 이유는 **좌표계 불일치** 때문입니다:

- **GWAS 데이터 (GCST009325)**: hg19 (GRCh37) 좌표계
- **Enhancer BED 파일들**: rn7 좌표계 기반

## 🔄 해결 방법: 좌표계 변환

### 1단계: LiftOver 도구 및 체인 파일 다운로드

```bash
# 좌표 변환 환경 설정 (자동 다운로드 시도)
python setup_liftover.py
```

만약 자동 다운로드가 실패하면 수동으로 다운로드하세요:

#### 수동 다운로드 필요 파일들:

1. **liftOver 실행 파일**
   ```bash
   cd 0_data/reference/liftover_data
   wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
   chmod +x liftOver
   ```

2. **rn7 → hg38 체인 파일**
   ```bash
   wget http://hgdownload.soe.ucsc.edu/goldenPath/rn7/liftOver/rn7ToHg38.over.chain.gz
   ```

3. **hg38 → hg19 체인 파일**
   ```bash
   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
   ```

### 2단계: 좌표 변환 실행

```bash
# 모든 BED 파일을 rn7 → hg38 → hg19로 변환
python setup_liftover.py
```

이 명령은 다음 변환을 수행합니다:
- `0_data/raw/cleaned_data/*.bed` → `0_data/processed/hg19_coordinates/*_hg19.bed`
- `0_data/raw/unique_data/*.bed` → `0_data/processed/hg19_coordinates/*_hg19.bed`

### 3단계: 변환 확인

변환이 완료되면 다음 디렉토리 구조가 생성됩니다:

```
0_data/
├── raw/
│   ├── cleaned_data/*.bed         # 원본 (rn7)
│   └── unique_data/*.bed          # 원본 (rn7)
└── processed/
    ├── hg38_coordinates/          # 중간 변환 (hg38)
    └── hg19_coordinates/          # 최종 변환 (hg19) ⭐
        ├── cleaned_data_Olig_cleaned_hg19.bed
        ├── cleaned_data_Nurr_cleaned_hg19.bed
        ├── unique_data_Olig_unique_hg19.bed
        └── ...
```

### 4단계: 배치 분석 재실행

좌표 변환이 완료된 후:

```bash
# 캐시 초기화하여 변환된 좌표로 재분석
python run_complete_batch_pipeline.py --force-reload
```

## 🔍 변환 결과 예상

좌표 변환 후 예상되는 변화:

- **Before**: `Enrichment ratio: 0.000` (좌표계 불일치)
- **After**: `Enrichment ratio: 1.2-3.5` (정상적인 enrichment)

## ⚡ 빠른 실행 가이드

```bash
# 1. 좌표 변환 (처음 1회만)
python setup_liftover.py

# 2. 배치 분석 재실행
python run_complete_batch_pipeline.py --force-reload
```

## 📊 변환 품질 확인

변환 후 다음을 확인하세요:

1. **변환률**: 일반적으로 90% 이상이어야 함
2. **영역 수**: 원본과 비슷해야 함
3. **Enrichment ratio**: 0이 아닌 합리적인 값

## 🚫 문제 해결

### 문제 1: liftOver 다운로드 실패
- 수동으로 다운로드하여 `0_data/reference/liftover_data/`에 저장
- 실행 권한 확인: `chmod +x liftOver`

### 문제 2: 체인 파일 없음
- UCSC 사이트에서 직접 다운로드
- 파일 무결성 확인

### 문제 3: 변환률이 낮음 (<80%)
- 원본 BED 파일 형식 확인
- 염색체 명명 규칙 확인 (chr1 vs 1)

---

**중요**: 이 변환 과정은 전체 분석의 정확성에 매우 중요합니다. 
변환 없이는 올바른 enrichment 분석이 불가능합니다.