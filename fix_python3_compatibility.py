#!/usr/bin/env python3
"""
LDSC Python 2 -> Python 3 호환성 수정 스크립트
"""

import os
import re
from pathlib import Path

def fix_python3_compatibility():
    """LDSC 코드를 Python 3과 호환되도록 수정"""
    
    ldsc_dir = Path("/scratch/prj/eng_waste_to_protein/repositories/bomin/1.Scripts/LDSC/ldsc")
    
    # Find all Python files
    python_files = list(ldsc_dir.rglob("*.py"))
    
    print(f"🔧 Python 3 호환성 수정 시작 - {len(python_files)}개 파일")
    
    fixes_applied = 0
    files_modified = 0
    
    for py_file in python_files:
        try:
            # Read file
            with open(py_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            original_content = content
            
            # Fix 1: print statements -> print function
            # Match "print something" but not "print(...)"
            content = re.sub(r'\bprint\s+([^(][^\n]*)', r'print(\1)', content)
            
            # Fix 2: print >> file syntax
            content = re.sub(r'print\s*>>\s*([^,\n]+),\s*(.+)', r'print(\2, file=\1)', content)
            
            # Fix 3: xrange -> range
            content = re.sub(r'\bxrange\b', 'range', content)
            
            # Fix 4: Fix tabs/spaces (convert tabs to spaces)
            lines = content.split('\n')
            fixed_lines = []
            for line in lines:
                # Convert tabs to 4 spaces at the beginning of lines
                if line.startswith('\t'):
                    leading_tabs = len(line) - len(line.lstrip('\t'))
                    line = '    ' * leading_tabs + line.lstrip('\t')
                fixed_lines.append(line)
            content = '\n'.join(fixed_lines)
            
            # Write back if changed
            if content != original_content:
                with open(py_file, 'w', encoding='utf-8') as f:
                    f.write(content)
                
                files_modified += 1
                print(f"  ✅ 수정됨: {py_file.name}")
                
                # Count specific fixes
                if 'print(' in content and 'print(' not in original_content:
                    fixes_applied += 1
                if 'range(' in content and 'xrange(' in original_content:
                    fixes_applied += 1
            
        except Exception as e:
            print(f"  ❌ 오류 {py_file.name}: {e}")
    
    print(f"\n📊 수정 완료:")
    print(f"  - 수정된 파일: {files_modified}개")
    print(f"  - 적용된 수정사항: {fixes_applied}개")
    
    return files_modified > 0

if __name__ == "__main__":
    success = fix_python3_compatibility()
    if success:
        print("✅ Python 3 호환성 수정 완료!")
    else:
        print("⚠️ 수정할 항목이 없거나 오류 발생")