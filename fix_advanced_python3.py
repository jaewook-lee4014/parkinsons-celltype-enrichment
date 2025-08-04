#!/usr/bin/env python3
"""
고급 Python 2 -> Python 3 호환성 수정 스크립트
"""

import os
import re
from pathlib import Path

def fix_advanced_python3_issues():
    """고급 Python 2/3 호환성 문제 수정"""
    
    ldsc_dir = Path("/scratch/prj/eng_waste_to_protein/repositories/bomin/1.Scripts/LDSC/ldsc")
    
    # Find all Python files
    python_files = list(ldsc_dir.rglob("*.py"))
    
    print(f"🔧 고급 Python 3 호환성 수정 시작 - {len(python_files)}개 파일")
    
    fixes_applied = 0
    files_modified = 0
    
    for py_file in python_files:
        try:
            # Read file
            with open(py_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            original_content = content
            
            # Fix 1: pandas .ix deprecated -> .loc
            content = re.sub(r'\.ix\[', '.loc[', content)
            
            # Fix 2: Line continuation issues (\)
            # Fix multiline format strings
            content = re.sub(
                r"'(.+?)\{.+?\}(.+?)'.format\(\\\)",
                r"'\1{}\2'.format(",
                content,
                flags=re.DOTALL
            )
            
            # Fix 3: Python 2 map() returns list -> Python 3 returns iterator
            # Wrap map() in list() where needed
            content = re.sub(r'map\(([^)]+)\)', r'list(map(\1))', content)
            
            # Fix 4: print >> syntax that wasn't caught
            content = re.sub(r'print\(>>', 'print(', content)
            
            # Fix 5: More complex print >> patterns
            content = re.sub(r'print\s*>>\s*open\(([^)]+)\)\s*,\s*(.+)', r'print(\2, file=open(\1))', content)
            
            # Fix 6: String methods that changed
            content = re.sub(r'\.rstrip\(\'\\n\'\)', '.rstrip()', content)
            
            # Fix 7: Fix import issues
            if 'import ConfigParser' in content:
                content = content.replace('import ConfigParser', 'import configparser as ConfigParser')
            
            # Write back if changed
            if content != original_content:
                with open(py_file, 'w', encoding='utf-8') as f:
                    f.write(content)
                
                files_modified += 1
                print(f"  ✅ 고급 수정: {py_file.name}")
                fixes_applied += 1
            
        except Exception as e:
            print(f"  ❌ 오류 {py_file.name}: {e}")
    
    print(f"\n📊 고급 수정 완료:")
    print(f"  - 수정된 파일: {files_modified}개")
    print(f"  - 적용된 수정사항: {fixes_applied}개")
    
    return files_modified > 0

if __name__ == "__main__":
    success = fix_advanced_python3_issues()
    if success:
        print("✅ 고급 Python 3 호환성 수정 완료!")
    else:
        print("⚠️ 수정할 항목이 없거나 오류 발생")