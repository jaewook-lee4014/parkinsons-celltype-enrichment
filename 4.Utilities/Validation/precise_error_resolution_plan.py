#!/usr/bin/env python3
"""
매우 정밀한 오류 해결 방안 제시
"""
import numpy as np
from pathlib import Path
import re

def analyze_current_issues():
    """현재 시각화의 모든 문제점 정밀 분석"""
    
    print("🔍 현재 시각화 코드의 정밀한 문제점 분석")
    print("="*80)
    
    issues = {
        "1. 데이터 소스 문제": {
            "문제점": [
                "가상의 enrichment 값 사용 (실제 LDSC 결과가 아님)",
                "Coefficient와 SE가 임의로 설정됨",
                "실제 annotation 파일의 SNP 수와 불일치"
            ],
            "영향": "과학적 신뢰성 부족",
            "심각도": "⚠️ 매우 높음"
        },
        
        "2. 계산 로직 문제": {
            "문제점": [
                "BED 파일 크기를 45 bytes/region으로 가정 (부정확)",
                "총 SNP 수를 1000만으로 임의 설정",
                "Enrichment factor가 순환 논리 (enrichment = factor)"
            ],
            "영향": "계산 결과의 정확성 결여",
            "심각도": "⚠️ 높음"
        },
        
        "3. P-value 계산 문제": {
            "문제점": [
                "Coefficient 기반 p-value 계산은 맞지만, coefficient 자체가 가상",
                "Enrichment SE와 coefficient SE의 관계가 불명확",
                "Multiple testing correction 없음"
            ],
            "영향": "통계적 유의성 판단 오류",
            "심각도": "⚠️ 중간"
        },
        
        "4. 시각화 문제": {
            "문제점": [
                "레전드가 너무 많아 복잡함 (8개 항목)",
                "Error bar가 실제 신뢰구간인지 SE인지 불명확",
                "색상 구분이 직관적이지 않음"
            ],
            "영향": "결과 해석의 어려움",
            "심각도": "⚠️ 낮음"
        }
    }
    
    for category, details in issues.items():
        print(f"\n{category}")
        print("-" * 60)
        print("문제점:")
        for issue in details["문제점"]:
            print(f"  • {issue}")
        print(f"영향: {details['영향']}")
        print(f"심각도: {details['심각도']}")
    
    return issues

def propose_precise_solutions():
    """정밀한 해결 방안 제시"""
    
    print("\n\n💡 정밀한 오류 해결 방안")
    print("="*80)
    
    solutions = {
        "해결방안 1: 실제 LDSC 실행": {
            "단계": [
                "1. LDSC 소프트웨어 설치 확인",
                "2. Annotation 파일 준비 (unique/all enhancer BED → annotation)",
                "3. GWAS summary statistics 준비",
                "4. Baseline model (53 annotations) 다운로드",
                "5. LDSC 실행 스크립트 작성",
                "6. 결과 파일에서 정확한 값 추출"
            ],
            "필요 명령어": """
# Annotation 생성
python make_annot.py \\
    --bed-file Microglia_unique.bed \\
    --bimfile 1000G.EUR.QC.{chr}.bim \\
    --annot-file Microglia_unique.{chr}.annot.gz

# LD score 계산  
python ldsc.py \\
    --l2 \\
    --bfile 1000G.EUR.QC.{chr} \\
    --ld-wind-cm 1 \\
    --annot Microglia_unique.{chr}.annot.gz \\
    --out Microglia_unique.{chr}

# Enrichment 분석
python ldsc.py \\
    --h2 PD_GWAS.sumstats.gz \\
    --ref-ld-chr baseline.,Microglia_unique. \\
    --w-ld-chr weights. \\
    --overlap-annot \\
    --frqfile-chr 1000G.EUR.QC.{chr}.frq \\
    --out PD_Microglia_unique
""",
            "결과 파일": "PD_Microglia_unique.results"
        },
        
        "해결방안 2: 기존 LDSC 결과 활용": {
            "단계": [
                "1. 이미 실행된 LDSC 결과 파일 찾기",
                "2. .results 파일 파싱",
                "3. Enrichment, Enrichment_std_error, Coefficient_p_value 추출",
                "4. 모든 셀타입과 조건에 대해 반복"
            ],
            "파일 패턴": [
                "*_unique.results",
                "*_all.results", 
                "*_cleaned.results"
            ],
            "파싱 코드": """
import pandas as pd

def parse_ldsc_results(results_file):
    df = pd.read_table(results_file)
    # Category, Prop._SNPs, Prop._h2, Enrichment, Enrichment_std_error, Coefficient_p_value
    return {
        'enrichment': df['Enrichment'].iloc[-1],
        'enrichment_se': df['Enrichment_std_error'].iloc[-1],
        'p_value': df['Coefficient_p_value'].iloc[-1],
        'prop_snps': df['Prop._SNPs'].iloc[-1],
        'prop_h2': df['Prop._h2'].iloc[-1]
    }
"""
        },
        
        "해결방안 3: Annotation 파일 기반 정확한 계산": {
            "단계": [
                "1. .annot.gz 파일에서 실제 SNP 수 계산",
                "2. 전체 baseline SNP 수 확인",
                "3. BED 파일의 실제 영역 수 계산",
                "4. 정확한 proportion 계산"
            ],
            "코드": """
import gzip
import numpy as np

def count_annotated_snps(annot_file):
    with gzip.open(annot_file, 'rt') as f:
        # Skip header
        next(f)
        # Count SNPs in annotation
        count = sum(1 for line in f if line.strip().split()[-1] == '1')
    return count

def get_total_snps(baseline_annot):
    with gzip.open(baseline_annot, 'rt') as f:
        next(f)
        total = sum(1 for _ in f)
    return total
"""
        },
        
        "해결방안 4: 수학적 일관성 보장": {
            "검증 항목": [
                "Enrichment = Prop_h2 / Prop_SNPs",
                "Z-score = Coefficient / Coefficient_SE", 
                "P-value = 2 * (1 - norm.cdf(|Z-score|))",
                "SE propagation이 올바른지 확인"
            ],
            "자동 검증 코드": """
def validate_consistency(data):
    # Enrichment 검증
    calc_enr = data['prop_h2'] / data['prop_snps']
    assert abs(calc_enr - data['enrichment']) < 0.1
    
    # P-value 검증 (coefficient 기반)
    z_score = data['coefficient'] / data['coefficient_se']
    calc_p = 2 * (1 - stats.norm.cdf(abs(z_score)))
    assert abs(np.log10(calc_p) - np.log10(data['p_value'])) < 0.5
"""
        }
    }
    
    for solution_name, details in solutions.items():
        print(f"\n📌 {solution_name}")
        print("-" * 70)
        
        if "단계" in details:
            print("실행 단계:")
            for step in details["단계"]:
                print(f"  {step}")
        
        if "필요 명령어" in details:
            print("\n필요 명령어:")
            print(details["필요 명령어"])
        
        if "코드" in details:
            print("\n구현 코드:")
            print(details["코드"])
    
    return solutions

def create_robust_visualization_framework():
    """견고한 시각화 프레임워크 제안"""
    
    print("\n\n🏗️ 견고한 시각화 프레임워크")
    print("="*80)
    
    framework = """
#!/usr/bin/env python3
'''
견고한 LDSC Enrichment 시각화 프레임워크
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
import logging

class LDSCEnrichmentVisualizer:
    def __init__(self, results_dir):
        self.results_dir = Path(results_dir)
        self.data = {}
        self.logger = logging.getLogger(__name__)
        
    def load_ldsc_results(self, cell_types, conditions=['unique', 'all']):
        '''실제 LDSC 결과 파일 로드'''
        for cell_type in cell_types:
            self.data[cell_type] = {}
            for condition in conditions:
                results_file = self.results_dir / f"{cell_type}_{condition}.results"
                if results_file.exists():
                    self.data[cell_type][condition] = self._parse_results(results_file)
                else:
                    self.logger.warning(f"Missing: {results_file}")
                    
    def _parse_results(self, results_file):
        '''LDSC results 파일 파싱'''
        df = pd.read_table(results_file)
        # 마지막 행이 우리가 추가한 annotation
        row = df.iloc[-1]
        
        data = {
            'enrichment': row['Enrichment'],
            'enrichment_se': row['Enrichment_std_error'],
            'prop_snps': row['Prop._SNPs'],
            'prop_h2': row['Prop._h2'],
            'coefficient': row.get('Coefficient', np.nan),
            'coefficient_se': row.get('Coefficient_std_error', np.nan),
            'p_value': row['Coefficient_p_value']
        }
        
        # 일관성 검증
        self._validate_data(data)
        return data
    
    def _validate_data(self, data):
        '''데이터 일관성 검증'''
        # Enrichment 계산 검증
        calc_enrichment = data['prop_h2'] / data['prop_snps']
        if abs(calc_enrichment - data['enrichment']) > 0.1:
            self.logger.warning("Enrichment calculation mismatch!")
            
        # 데이터 범위 검증
        if data['enrichment'] < 0 or data['enrichment'] > 1000:
            self.logger.warning(f"Unusual enrichment: {data['enrichment']}")
            
    def create_visualization(self, output_file):
        '''최종 시각화 생성'''
        # 데이터 준비
        cell_types = list(self.data.keys())
        
        # 그림 생성
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # ... 시각화 코드 ...
        
        # 저장
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        self.logger.info(f"Saved: {output_file}")

# 사용 예시
if __name__ == "__main__":
    visualizer = LDSCEnrichmentVisualizer("ldsc_results/")
    visualizer.load_ldsc_results(['Microglia', 'Neuron', 'Oligodendrocyte', 'Dopaminergic'])
    visualizer.create_visualization("final_enrichment_analysis.png")
"""
    
    print(framework)
    
    return framework

def provide_immediate_fixes():
    """즉시 적용 가능한 수정사항"""
    
    print("\n\n⚡ 즉시 적용 가능한 수정사항")
    print("="*80)
    
    fixes = {
        "1. Annotation 파일 기반 SNP 수 계산": {
            "문제": "BED 파일 크기로 추정하는 것은 부정확",
            "해결": "실제 annotation 파일에서 SNP 수 직접 계산",
            "코드": """
# check_intersect_ratios.py 수정
annot_file = 'Neg.1.annot.gz'
with gzip.open(annot_file, 'rt') as f:
    header = next(f).strip().split()
    col_idx = header.index('Neg')
    intersect_snps = sum(1 for line in f if line.strip().split()[col_idx] == '1')
"""
        },
        
        "2. 실제 LDSC 로그 파일 확인": {
            "문제": "가상의 enrichment 값 사용",
            "해결": "기존 LDSC 실행 로그에서 값 추출",
            "파일 패턴": [
                "*.log",
                "*_ldsc.log", 
                "PD_*.log"
            ]
        },
        
        "3. 전체 SNP 수 정확히 확인": {
            "문제": "1000만 SNP 가정",
            "해결": "Baseline annotation에서 실제 수 확인",
            "명령어": "zcat baseline.1.annot.gz | tail -n +2 | wc -l"
        }
    }
    
    for fix_name, details in fixes.items():
        print(f"\n{fix_name}")
        print("-" * 60)
        for key, value in details.items():
            print(f"{key}: {value}")
    
    return fixes

if __name__ == "__main__":
    print("🔧 매우 정밀한 오류 해결 방안\n")
    
    # 1. 현재 문제점 분석
    issues = analyze_current_issues()
    
    # 2. 정밀한 해결 방안
    solutions = propose_precise_solutions()
    
    # 3. 견고한 프레임워크
    framework = create_robust_visualization_framework()
    
    # 4. 즉시 적용 가능한 수정
    fixes = provide_immediate_fixes()
    
    print("\n\n" + "="*80)
    print("📋 우선순위별 실행 계획:")
    print("="*80)
    print("1. 즉시: Annotation 파일에서 실제 SNP 수 확인")
    print("2. 단기: 기존 LDSC 결과 파일 찾아서 파싱")
    print("3. 중기: 견고한 시각화 프레임워크 구현")
    print("4. 장기: 필요시 LDSC 재실행")