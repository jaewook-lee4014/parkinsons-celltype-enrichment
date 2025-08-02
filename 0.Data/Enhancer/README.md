# Enhancer Data Directory

이 폴더는 뇌 세포타입별 enhancer BED 파일들을 저장합니다.

## 파일 설명

### 세포타입별 Enhancer 파일

- **Olig_cleaned.bed / Olig_unique.bed**: Oligodendrocyte enhancer
- **Nurr_cleaned.bed / Nurr_unique.bed**: Dopaminergic neuron enhancer  
- **NeuN_cleaned.bed / NeuN_unique.bed**: General neuron enhancer
- **Neg_cleaned.bed / Neg_unique.bed**: Microglia enhancer

### 처리 방법

- **cleaned**: 전처리된 enhancer 영역
- **unique**: 고유한 enhancer 영역만 추출

## 좌표계

- 원본: rn7 (rat genome)
- 변환 필요: rn7 → hg38 → hg19 (human genome)
- 변환 도구: `1.Scripts/Utils/setup_liftover.py`