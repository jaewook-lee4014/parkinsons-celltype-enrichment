#!/usr/bin/env python3
"""
파킨슨병 GWAS - 세포타입별 Enhancer Enrichment 분석
Main execution script
"""

import os
import sys
import argparse
from pathlib import Path

# Add scripts directories to path
sys.path.append(str(Path(__file__).parent / "1.Scripts" / "LDSC"))
sys.path.append(str(Path(__file__).parent / "1.Scripts" / "Visualization"))
sys.path.append(str(Path(__file__).parent / "1.Scripts" / "Utils"))

def main():
    parser = argparse.ArgumentParser(
        description="파킨슨병 GWAS 세포타입별 Enhancer Enrichment 분석",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
사용 예시:
  # 1. 좌표 변환
  python main.py --step coordinate
  
  # 2. LDSC 분석
  python main.py --step ldsc
  
  # 3. 시각화
  python main.py --step visualize
  
  # 전체 파이프라인 실행
  python main.py --all
        """
    )
    
    parser.add_argument(
        '--step', 
        choices=['coordinate', 'ldsc', 'visualize'],
        help='실행할 단계 선택'
    )
    parser.add_argument(
        '--all', 
        action='store_true',
        help='전체 파이프라인 실행'
    )
    
    args = parser.parse_args()
    
    if args.all:
        print("=" * 60)
        print("파킨슨병 GWAS 세포타입별 Enhancer Enrichment 분석")
        print("전체 파이프라인 실행")
        print("=" * 60)
        
        # Step 1: Coordinate conversion
        print("\n[1/3] 좌표계 변환...")
        from setup_liftover import main as setup_liftover_main
        setup_liftover_main()
        
        # Step 2: LDSC analysis
        print("\n[2/3] LDSC 분석...")
        from ldsc_analysis_system import main as ldsc_main
        ldsc_main()
        
        # Step 3: Visualization
        print("\n[3/3] 시각화...")
        from celltype_manhattan_plot import main as manhattan_main
        manhattan_main()
        
        print("\n✅ 전체 파이프라인 완료!")
        
    elif args.step == 'coordinate':
        print("좌표계 변환 실행...")
        from setup_liftover import main as setup_liftover_main
        setup_liftover_main()
        
    elif args.step == 'ldsc':
        print("LDSC 분석 실행...")
        from ldsc_analysis_system import main as ldsc_main
        ldsc_main()
        
    elif args.step == 'visualize':
        print("시각화 실행...")
        from celltype_manhattan_plot import main as manhattan_main
        manhattan_main()
        
    else:
        parser.print_help()

if __name__ == "__main__":
    main()