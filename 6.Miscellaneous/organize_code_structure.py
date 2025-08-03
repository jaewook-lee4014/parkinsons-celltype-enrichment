#!/usr/bin/env python3
"""
코드 파일들을 속성에 따라 정리하는 스크립트
"""
import os
import shutil
from pathlib import Path
import re

def organize_code_files():
    """코드 파일들을 속성별로 정리"""
    
    base_dir = Path("/cephfs/volumes/hpc_data_prj/eng_waste_to_protein/ae035a41-20d2-44f3-aa46-14424ab0f6bf/repositories/bomin")
    
    # 새로운 디렉토리 구조
    directories = {
        "2.Analysis": {
            "LDSC": [],
            "Enrichment": [],
            "Annotation": [],
            "Overlap": []
        },
        "3.Visualization": {
            "Enrichment_Plots": [],
            "Manhattan_Plots": [],
            "Concept_Diagrams": []
        },
        "4.Utilities": {
            "Data_Processing": [],
            "Validation": [],
            "File_Search": []
        },
        "5.Documentation": {
            "Reports": [],
            "Summaries": []
        }
    }
    
    # 파일 분류 패턴
    patterns = {
        "LDSC": ["ldsc", "ld_score", "heritability", "h2"],
        "Enrichment": ["enrichment", "enrich"],
        "Annotation": ["annot", "annotation"],
        "Overlap": ["overlap", "intersect", "unique", "cleaned"],
        "Enrichment_Plots": ["visualization", "plot", "celltype.*vis", "bar.*plot"],
        "Manhattan_Plots": ["manhattan"],
        "Concept_Diagrams": ["concept", "diagram", "explain"],
        "Data_Processing": ["process", "convert", "calculate", "analyze"],
        "Validation": ["valid", "verify", "check", "test", "precise"],
        "File_Search": ["find", "search", "trace", "extensive"],
        "Reports": ["report", "summary"],
        "Summaries": ["\\.md$", "readme"]
    }
    
    # 현재 디렉토리의 Python 파일 목록
    python_files = list(base_dir.glob("*.py"))
    
    print(f"📁 총 {len(python_files)}개의 Python 파일 발견")
    print("="*80)
    
    # 파일 분류
    file_mapping = {}
    unclassified = []
    
    for file_path in python_files:
        filename = file_path.name.lower()
        classified = False
        
        # 각 카테고리별로 패턴 매칭
        for category, category_patterns in patterns.items():
            for pattern in category_patterns:
                if re.search(pattern, filename):
                    # 적절한 디렉토리 찾기
                    for main_dir, sub_dirs in directories.items():
                        if category in sub_dirs:
                            file_mapping[file_path] = (main_dir, category)
                            classified = True
                            break
                    if classified:
                        break
            if classified:
                break
        
        if not classified:
            unclassified.append(file_path)
    
    # 디렉토리 생성
    print("\n📂 디렉토리 구조 생성:")
    for main_dir, sub_dirs in directories.items():
        main_path = base_dir / main_dir
        main_path.mkdir(exist_ok=True)
        print(f"\n{main_dir}/")
        
        for sub_dir in sub_dirs:
            sub_path = main_path / sub_dir
            sub_path.mkdir(exist_ok=True)
            print(f"  └── {sub_dir}/")
    
    # 파일 이동
    print("\n\n📦 파일 이동:")
    print("-"*80)
    
    for file_path, (main_dir, sub_dir) in file_mapping.items():
        target_dir = base_dir / main_dir / sub_dir
        target_path = target_dir / file_path.name
        
        try:
            shutil.move(str(file_path), str(target_path))
            print(f"✅ {file_path.name} → {main_dir}/{sub_dir}/")
        except Exception as e:
            print(f"❌ {file_path.name}: {e}")
    
    # 분류되지 않은 파일
    if unclassified:
        print(f"\n\n⚠️ 분류되지 않은 파일 ({len(unclassified)}개):")
        misc_dir = base_dir / "6.Miscellaneous"
        misc_dir.mkdir(exist_ok=True)
        
        for file_path in unclassified:
            try:
                target_path = misc_dir / file_path.name
                shutil.move(str(file_path), str(target_path))
                print(f"  → {file_path.name} → 6.Miscellaneous/")
            except Exception as e:
                print(f"  ❌ {file_path.name}: {e}")
    
    # 특별 파일들 처리
    special_files = {
        "main.py": ".",  # 루트에 유지
        "README.md": ".",  # 루트에 유지
        "requirements.txt": ".",  # 루트에 유지
    }
    
    for filename, target in special_files.items():
        src = misc_dir / filename if (misc_dir / filename).exists() else None
        if src and src.exists():
            dst = base_dir / filename
            shutil.move(str(src), str(dst))
            print(f"\n♻️ {filename} → 루트 디렉토리로 복원")
    
    print("\n\n✨ 코드 정리 완료!")
    
    # 정리 결과 요약
    print("\n📊 정리 결과:")
    print("="*80)
    
    for main_dir in directories.keys():
        main_path = base_dir / main_dir
        if main_path.exists():
            total_files = sum(1 for _ in main_path.rglob("*.py"))
            print(f"{main_dir}: {total_files}개 파일")
            
            for sub_dir in directories[main_dir]:
                sub_path = main_path / sub_dir
                if sub_path.exists():
                    file_count = len(list(sub_path.glob("*.py")))
                    if file_count > 0:
                        print(f"  └── {sub_dir}: {file_count}개")

if __name__ == "__main__":
    organize_code_files()