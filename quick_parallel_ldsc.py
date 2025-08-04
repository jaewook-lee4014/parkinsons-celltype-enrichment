#!/usr/bin/env python3
"""
현재 로그인 노드에서 완료된 annotation에 대해 LD score 계산
"""

import os
import sys
import subprocess
from pathlib import Path
import time
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
import glob

# 로깅 설정
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class QuickLDSCCalculator:
    def __init__(self):
        self.base_dir = Path("/cephfs/volumes/hpc_data_prj/eng_waste_to_protein/ae035a41-20d2-44f3-aa46-14424ab0f6bf/repositories/bomin")
        self.ldsc_dir = self.base_dir / "1.Scripts/LDSC/ldsc"
        self.annotations_dir = self.base_dir / "ldsc_results/annotations"
        self.reference_dir = Path("/scratch/prj/eng_waste_to_protein/repositories/bomin/0.Data/Reference/ldsc_reference")
        self.results_dir = Path("/scratch/prj/eng_waste_to_protein/repositories/bomin/combined_ld_scores")
        
        # 완료된 데이터셋들 (현재 완료 확인된 것들)
        self.completed_datasets = [
            "Olig_unique",    
            "Neg_unique",     
            "Neg_cleaned",    
            "Olig_cleaned",   
            "NeuN_cleaned",   
            "NeuN_unique"     
        ]
        
        # 로그인 노드 최적화 (너무 많은 병렬 처리는 피함)
        self.max_workers = 2
        
    def check_completed_annotations(self):
        """완료된 annotation 파일들 확인"""
        completed_annots = {}
        
        for dataset in self.completed_datasets:
            logger.info(f"🔍 {dataset} annotation 파일 확인 중...")
            
            annot_files = list(self.annotations_dir.glob(f"{dataset}.*.annot.gz"))
            if len(annot_files) >= 22:
                completed_annots[dataset] = len(annot_files)
                logger.info(f"  ✅ {dataset}: {len(annot_files)}개 파일 발견")
            elif len(annot_files) > 0:
                completed_annots[dataset] = len(annot_files)
                logger.info(f"  🔄 {dataset}: {len(annot_files)}개 파일 발견 (부분 완료)")
            else:
                logger.warning(f"  ❌ {dataset}: annotation 파일 없음")
        
        return completed_annots
    
    def calculate_ld_score_single(self, dataset_name, chromosome):
        """단일 염색체에 대한 LD score 계산"""
        try:
            # annotation 파일 확인
            annot_file = self.annotations_dir / f"{dataset_name}.{chromosome}.annot.gz"
            if not annot_file.exists():
                logger.warning(f"    ⚠️ {dataset_name} Chr{chromosome}: annotation 파일 없음")
                return False
            
            # 결과 파일이 이미 존재하면 스킵
            ldscore_file = self.results_dir / f"{dataset_name}.{chromosome}.l2.ldscore.gz"
            if ldscore_file.exists():
                logger.info(f"    ✅ {dataset_name} Chr{chromosome}: 이미 완료됨")
                return True
                
            # 결과 디렉토리 생성
            self.results_dir.mkdir(parents=True, exist_ok=True)
            
            # LD score 계산 명령어
            ldscore_cmd = [
                "python", str(self.ldsc_dir / "ldsc.py"),
                "--l2",
                "--bfile", f"{self.reference_dir}/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chromosome}",
                "--ld-wind-cm", "1",
                "--annot", str(annot_file),
                "--out", str(self.results_dir / f"{dataset_name}.{chromosome}"),
                "--print-snps", str(self.reference_dir / "w_hm3.snplist")
            ]
            
            logger.info(f"    🔗 {dataset_name} Chr{chromosome}: LD score 계산 시작...")
            
            # 실행
            result = subprocess.run(
                ldscore_cmd, 
                capture_output=True, 
                text=True, 
                cwd=str(self.ldsc_dir),
                timeout=1200  # 20분 타임아웃
            )
            
            if result.returncode == 0:
                logger.info(f"    ✅ {dataset_name} Chr{chromosome}: 완료")
                return True
            else:
                logger.error(f"    ❌ {dataset_name} Chr{chromosome}: 실패")
                logger.error(f"    stderr: {result.stderr[:300]}")
                return False
                
        except subprocess.TimeoutExpired:
            logger.error(f"    ⏰ {dataset_name} Chr{chromosome}: 타임아웃 (20분)")
            return False
        except Exception as e:
            logger.error(f"    ❌ {dataset_name} Chr{chromosome}: 오류 - {e}")
            return False
    
    def process_single_dataset(self, dataset_name, available_chromosomes):
        """단일 데이터셋 처리 (사용 가능한 염색체만)"""
        logger.info(f"\n🧬 {dataset_name} LD scores 계산 시작")
        
        # 사용 가능한 염색체 찾기
        chromosomes = []
        for chr_num in range(1, 23):
            annot_file = self.annotations_dir / f"{dataset_name}.{chr_num}.annot.gz"
            if annot_file.exists():
                chromosomes.append(chr_num)
        
        logger.info(f"  📋 {dataset_name}: {len(chromosomes)}개 염색체 처리 예정")
        
        success_count = 0
        
        # 순차 처리 (로그인 노드에서 너무 많은 병렬 처리 방지)
        for chromosome in chromosomes:
            success = self.calculate_ld_score_single(dataset_name, chromosome)
            if success:
                success_count += 1
                
        logger.info(f"  📊 {dataset_name} 완료: {success_count}/{len(chromosomes)} 염색체")
        return success_count, len(chromosomes)
    
    def run_calculation(self):
        """LD score 계산 실행"""
        logger.info("🚀 LD Score 계산 시작")
        
        # 환경 확인
        logger.info("🔧 환경 확인 중...")
        logger.info(f"  작업 디렉토리: {os.getcwd()}")
        logger.info(f"  LDSC 디렉토리: {self.ldsc_dir}")
        logger.info(f"  Annotations 디렉토리: {self.annotations_dir}")
        logger.info(f"  Reference 디렉토리: {self.reference_dir}")
        logger.info(f"  결과 디렉토리: {self.results_dir}")
        
        # 필요한 디렉토리 존재 확인
        if not self.ldsc_dir.exists():
            logger.error(f"❌ LDSC 디렉토리가 없습니다: {self.ldsc_dir}")
            return False
        
        if not self.annotations_dir.exists():
            logger.error(f"❌ Annotations 디렉토리가 없습니다: {self.annotations_dir}")
            return False
            
        if not self.reference_dir.exists():
            logger.error(f"❌ Reference 디렉토리가 없습니다: {self.reference_dir}")
            return False
        
        # 완료된 annotation 확인
        completed_annots = self.check_completed_annotations()
        
        if not completed_annots:
            logger.error("❌ 완료된 annotation 파일이 없습니다")
            return False
        
        logger.info(f"📋 처리할 데이터셋: {list(completed_annots.keys())}")
        
        # 각 데이터셋 처리
        total_success = 0
        total_chromosomes = 0
        
        for dataset_name in completed_annots.keys():
            success, total = self.process_single_dataset(dataset_name, completed_annots[dataset_name])
            total_success += success
            total_chromosomes += total
        
        logger.info(f"\n🎉 LD Score 계산 완료!")
        logger.info(f"  ✅ 성공: {total_success}/{total_chromosomes} 염색체")
        logger.info(f"  📊 처리된 데이터셋: {len(completed_annots)}")
        
        return total_success > 0

def main():
    """메인 실행 함수"""
    calculator = QuickLDSCCalculator()
    
    success = calculator.run_calculation()
    
    if success:
        logger.info("🎯 LD Score 계산이 성공적으로 시작되었습니다!")
        return True
    else:
        logger.error("❌ LD Score 계산이 실패했습니다")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)