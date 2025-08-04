#!/usr/bin/env python3
"""
HPC 노드에서 완료된 annotation에 대해 병렬로 LD score 계산
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

class ParallelLDSCCalculator:
    def __init__(self):
        self.base_dir = Path("/cephfs/volumes/hpc_data_prj/eng_waste_to_protein/ae035a41-20d2-44f3-aa46-14424ab0f6bf/repositories/bomin")
        self.ldsc_dir = self.base_dir / "1.Scripts/LDSC/ldsc"
        self.reference_dir = Path("/scratch/prj/eng_waste_to_protein/repositories/bomin/0.Data/Reference/ldsc_reference")
        self.results_dir = Path("/scratch/prj/eng_waste_to_protein/repositories/bomin/combined_ld_scores")
        
        # 완료된 데이터셋들
        self.completed_datasets = [
            "Olig_unique",    # 1,937 SNPs
            "Neg_unique",     # 9,496 SNPs  
            "Neg_cleaned",    # 30,885 SNPs
            "Olig_cleaned",   # 17,888 SNPs
            "NeuN_cleaned",   # 42,295 SNPs
            "NeuN_unique"     # 20,184 SNPs (방금 완료)
        ]
        
        # HPC 노드 최적화
        self.max_workers = 4  # 2 CPU cores -> 4 threads
        
    def find_completed_annotations(self):
        """완료된 annotation 파일들 찾기"""
        completed_annots = {}
        
        for dataset in self.completed_datasets:
            logger.info(f"🔍 {dataset} annotation 파일 확인 중...")
            
            annot_files = list(self.results_dir.glob(f"{dataset}.*.annot.gz"))
            if len(annot_files) >= 22:  # 22개 염색체 모두 있는지 확인
                completed_annots[dataset] = sorted(annot_files)
                logger.info(f"  ✅ {dataset}: {len(annot_files)}개 파일 발견")
            else:
                logger.warning(f"  ⚠️ {dataset}: {len(annot_files)}개 파일만 발견 (22개 필요)")
        
        return completed_annots
    
    def calculate_ld_score_single(self, dataset_name, chromosome):
        """단일 염색체에 대한 LD score 계산"""
        try:
            # 결과 파일이 이미 존재하면 스킵
            ldscore_file = self.results_dir / f"{dataset_name}.{chromosome}.l2.ldscore.gz"
            if ldscore_file.exists():
                logger.info(f"    ✅ {dataset_name} Chr{chromosome}: 이미 완료됨")
                return True
            
            # annotation 파일 확인
            annot_file = self.results_dir / f"{dataset_name}.{chromosome}.annot.gz"
            if not annot_file.exists():
                logger.error(f"    ❌ {dataset_name} Chr{chromosome}: annotation 파일 없음")
                return False
            
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
                timeout=600  # 10분 타임아웃
            )
            
            if result.returncode == 0:
                logger.info(f"    ✅ {dataset_name} Chr{chromosome}: 완료")
                return True
            else:
                logger.error(f"    ❌ {dataset_name} Chr{chromosome}: 실패")
                logger.error(f"    stderr: {result.stderr[:200]}")
                return False
                
        except subprocess.TimeoutExpired:
            logger.error(f"    ⏰ {dataset_name} Chr{chromosome}: 타임아웃")
            return False
        except Exception as e:
            logger.error(f"    ❌ {dataset_name} Chr{chromosome}: 오류 - {e}")
            return False
    
    def calculate_dataset_parallel(self, dataset_name):
        """특정 데이터셋의 모든 염색체를 병렬로 계산"""
        logger.info(f"\n🧬 {dataset_name} LD scores 병렬 계산 시작")
        
        chromosomes = list(range(1, 23))  # 1-22번 염색체
        success_count = 0
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # 모든 염색체 작업 제출
            future_to_chr = {
                executor.submit(self.calculate_ld_score_single, dataset_name, chr): chr 
                for chr in chromosomes
            }
            
            # 결과 수집
            for future in as_completed(future_to_chr):
                chromosome = future_to_chr[future]
                try:
                    success = future.result()
                    if success:
                        success_count += 1
                except Exception as e:
                    logger.error(f"    ❌ Chr{chromosome} 처리 중 오류: {e}")
        
        logger.info(f"  📊 {dataset_name} 완료: {success_count}/22 염색체")
        return success_count == 22
    
    def run_parallel_calculation(self):
        """모든 완료된 데이터셋에 대해 병렬 LD score 계산"""
        logger.info("🚀 병렬 LD Score 계산 시작")
        logger.info(f"  💻 HPC 노드: erc-hpc-comp048")
        logger.info(f"  🧵 최대 동시 작업: {self.max_workers}")
        
        # 완료된 annotation 확인
        completed_annots = self.find_completed_annotations()
        
        if not completed_annots:
            logger.error("❌ 완료된 annotation 파일이 없습니다")
            return False
        
        logger.info(f"📋 처리할 데이터셋: {list(completed_annots.keys())}")
        
        # 각 데이터셋별로 순차 처리 (데이터셋 내에서는 병렬)
        success_datasets = []
        
        for dataset_name in completed_annots.keys():
            success = self.calculate_dataset_parallel(dataset_name)
            if success:
                success_datasets.append(dataset_name)
        
        logger.info(f"\n🎉 병렬 계산 완료!")
        logger.info(f"  ✅ 성공: {len(success_datasets)}/{len(completed_annots)} 데이터셋")
        logger.info(f"  📊 성공한 데이터셋: {success_datasets}")
        
        return len(success_datasets) > 0

def main():
    """메인 실행 함수"""
    calculator = ParallelLDSCCalculator()
    
    # 환경 확인
    logger.info("🔧 환경 확인 중...")
    logger.info(f"  작업 디렉토리: {os.getcwd()}")
    logger.info(f"  LDSC 디렉토리: {calculator.ldsc_dir}")
    logger.info(f"  Reference 디렉토리: {calculator.reference_dir}")
    logger.info(f"  결과 디렉토리: {calculator.results_dir}")
    
    # 필요한 디렉토리 존재 확인
    if not calculator.ldsc_dir.exists():
        logger.error(f"❌ LDSC 디렉토리가 없습니다: {calculator.ldsc_dir}")
        return False
    
    if not calculator.reference_dir.exists():
        logger.error(f"❌ Reference 디렉토리가 없습니다: {calculator.reference_dir}")
        return False
    
    # 결과 디렉토리 생성
    calculator.results_dir.mkdir(parents=True, exist_ok=True)
    
    # 병렬 계산 실행
    success = calculator.run_parallel_calculation()
    
    if success:
        logger.info("🎯 병렬 LD Score 계산이 성공적으로 완료되었습니다!")
        return True
    else:
        logger.error("❌ 병렬 LD Score 계산이 실패했습니다")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)