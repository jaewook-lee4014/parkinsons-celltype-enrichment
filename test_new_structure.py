#!/usr/bin/env python3
"""
새로운 데이터 구조 테스트 스크립트
"""

import sys
from pathlib import Path

# Add the script directory to path
sys.path.append(str(Path(__file__).parent / "1.Scripts" / "LDSC"))

from ldsc_analysis_system import LDSCConfig, AnnotationGenerator

def main():
    """새로운 LDSC 시스템 테스트"""
    print("🧪 새로운 LDSC 데이터 구조 테스트")
    
    # Initialize config
    config = LDSCConfig()
    
    # Test datasets definition
    print(f"\n📋 정의된 데이터셋 ({len(config.datasets)}개):")
    for i, dataset in enumerate(config.datasets, 1):
        print(f"  {i}. {dataset}")
    
    # Check enhancer files
    print(f"\n📁 Enhancer BED 파일 확인:")
    missing_files = []
    for dataset in config.datasets:
        bed_file = config.enhancer_dir / f"{dataset}.bed"
        if bed_file.exists():
            print(f"  ✅ {dataset}.bed")
        else:
            print(f"  ❌ {dataset}.bed")
            missing_files.append(dataset)
    
    # Check GWAS file
    print(f"\n🧬 GWAS 데이터 확인:")
    if config.gwas_file.exists():
        print(f"  ✅ {config.gwas_file.name}")
    else:
        print(f"  ❌ {config.gwas_file.name}")
    
    # Check directories
    print(f"\n📂 디렉토리 구조:")
    dirs_to_check = [
        ("Annotations", config.annotations_dir),
        ("LD Scores", config.ld_scores_dir), 
        ("Results", config.results_dir),
        ("Reference", config.reference_dir)
    ]
    
    for name, path in dirs_to_check:
        if path.exists():
            print(f"  ✅ {name}: {path}")
        else:
            print(f"  ❌ {name}: {path}")
    
    # Summary
    print(f"\n📊 요약:")
    print(f"  - 데이터셋: {len(config.datasets)}개")
    print(f"  - 누락된 BED 파일: {len(missing_files)}개")
    if missing_files:
        print(f"    누락: {', '.join(missing_files)}")
    
    # Test reference validation
    print(f"\n🔍 Reference 파일 검증:")
    try:
        is_valid = config.validate_reference_files()
        if is_valid:
            print(f"  ✅ 모든 reference 파일 확인됨")
        else:
            print(f"  ❌ 일부 reference 파일 누락")
    except Exception as e:
        print(f"  ⚠️ 검증 중 오류: {e}")
    
    print(f"\n🎯 다음 단계: LD score 계산 실행 준비됨")

if __name__ == "__main__":
    main()