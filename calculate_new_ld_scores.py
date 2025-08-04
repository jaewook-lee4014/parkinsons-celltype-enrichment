#!/usr/bin/env python3
"""
새로운 annotation에 대한 LD score 계산
8개 데이터셋 × 22개 염색체 = 176개 LD score 파일 생성
"""

import sys
import os
import subprocess
from pathlib import Path
import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class NewLDScoreCalculator:
    def __init__(self):
        # Paths
        self.base_cephfs = Path("/cephfs/volumes/hpc_data_prj/eng_waste_to_protein/ae035a41-20d2-44f3-aa46-14424ab0f6bf/repositories/bomin")
        self.base_scratch = Path("/scratch/prj/eng_waste_to_protein/repositories/bomin")
        
        # LDSC software
        self.ldsc_dir = self.base_scratch / "1.Scripts" / "LDSC" / "ldsc"
        
        # Input annotations (already generated)
        self.annotations_dir = self.base_cephfs / "ldsc_results" / "annotations"
        
        # Output directory for NEW LD scores
        self.output_dir = self.base_scratch / "new_ld_scores"
        self.output_dir.mkdir(exist_ok=True)
        
        # Reference data
        self.reference_dir = self.base_scratch / "0.Data" / "Reference" / "ldsc_reference"
        self.plink_files = self.reference_dir / "1000G_EUR_Phase3_plink" / "1000G.EUR.QC"
        self.snp_list = self.reference_dir / "hm3_no_MHC.list.txt"
        
        # 8 datasets
        self.datasets = [
            "Olig_cleaned", "Olig_unique",
            "Neg_cleaned", "Neg_unique", 
            "NeuN_cleaned", "NeuN_unique",
            "Nurr_cleaned", "Nurr_unique"
        ]
        
        logger.info(f"🧬 LD Score Calculator 초기화")
        logger.info(f"📋 처리할 데이터셋: {len(self.datasets)}개")
        logger.info(f"📂 출력 디렉토리: {self.output_dir}")
    
    def check_annotations(self):
        """annotation 파일들이 모두 존재하는지 확인"""
        logger.info("🔍 Annotation 파일 확인 중...")
        
        total_annotations = 0
        missing_annotations = []
        
        for dataset in self.datasets:
            dataset_count = 0
            for chr_num in range(1, 23):
                annot_file = self.annotations_dir / f"{dataset}.{chr_num}.annot.gz"
                if annot_file.exists():
                    dataset_count += 1
                else:
                    missing_annotations.append(f"{dataset}.{chr_num}")
            
            logger.info(f"  {dataset}: {dataset_count}/22 chromosomes")
            total_annotations += dataset_count
        
        logger.info(f"📊 총 annotation 파일: {total_annotations}/176")
        
        if missing_annotations:
            logger.warning(f"⚠️ 누락된 annotation: {len(missing_annotations)}개")
            for missing in missing_annotations[:10]:  # Show first 10
                logger.warning(f"    {missing}")
            if len(missing_annotations) > 10:
                logger.warning(f"    ... and {len(missing_annotations)-10} more")
        
        return total_annotations, missing_annotations
    
    def calculate_ld_score_single(self, dataset: str, chromosome: int):
        """단일 염색체의 LD score 계산"""
        try:
            # Check if annotation exists
            annot_file = self.annotations_dir / f"{dataset}.{chromosome}.annot.gz"
            if not annot_file.exists():
                return f"{dataset} Chr{chromosome}: annotation 파일 없음"
            
            # Check if output already exists
            output_file = self.output_dir / f"{dataset}.{chromosome}.l2.ldscore.gz"
            if output_file.exists():
                return f"{dataset} Chr{chromosome}: 이미 완료됨"
            
            # LD score calculation command
            cmd = [
                "python3", "/users/k23070952/.local/bin/ldsc.py",
                "--l2",
                "--bfile", f"{self.plink_files}.{chromosome}",
                "--ld-wind-cm", "1",
                "--annot", str(annot_file),
                "--out", str(self.output_dir / f"{dataset}.{chromosome}"),
                "--print-snps", str(self.snp_list)
            ]
            
            # Execute command
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=1800  # 30 minutes timeout
            )
            
            if result.returncode == 0:
                return f"{dataset} Chr{chromosome}: 성공"
            else:
                error_msg = result.stderr[:200] if result.stderr else "Unknown error"
                return f"{dataset} Chr{chromosome}: 실패 - {error_msg}"
                
        except subprocess.TimeoutExpired:
            return f"{dataset} Chr{chromosome}: 타임아웃 (30분)"
        except Exception as e:
            return f"{dataset} Chr{chromosome}: 오류 - {str(e)}"
    
    def calculate_all_ld_scores(self, max_workers=2):
        """모든 LD score 계산 (병렬 처리)"""
        logger.info(f"🚀 LD Score 계산 시작 - {max_workers}개 워커 사용")
        
        # Create task list
        tasks = []
        for dataset in self.datasets:
            for chromosome in range(1, 23):
                tasks.append((dataset, chromosome))
        
        logger.info(f"📋 총 작업: {len(tasks)}개 (8 datasets × 22 chromosomes)")
        
        completed = 0
        successful = 0
        failed = 0
        
        # Progress tracking
        start_time = time.time()
        
        # Execute tasks with thread pool
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_task = {
                executor.submit(self.calculate_ld_score_single, dataset, chr_num): (dataset, chr_num)
                for dataset, chr_num in tasks
            }
            
            # Process results as they complete
            for future in as_completed(future_to_task):
                dataset, chr_num = future_to_task[future]
                result = future.result()
                
                completed += 1
                
                if "성공" in result:
                    successful += 1
                    logger.info(f"✅ {result} [{completed}/{len(tasks)}]")
                elif "이미 완료됨" in result:
                    successful += 1
                    logger.info(f"⏭️ {result} [{completed}/{len(tasks)}]")
                else:
                    failed += 1
                    logger.error(f"❌ {result} [{completed}/{len(tasks)}]")
                
                # Progress update every 10 tasks
                if completed % 10 == 0 or completed in [1, 5]:
                    elapsed = time.time() - start_time
                    rate = completed / elapsed if elapsed > 0 else 0
                    eta = (len(tasks) - completed) / rate if rate > 0 else 0
                    
                    logger.info(f"📊 진행률: {completed}/{len(tasks)} ({completed/len(tasks)*100:.1f}%)")
                    logger.info(f"⏱️ 처리 속도: {rate:.2f} tasks/sec, 예상 완료: {eta/60:.1f}분")
        
        # Final summary
        elapsed = time.time() - start_time
        logger.info(f"\n🎉 LD Score 계산 완료!")
        logger.info(f"  ✅ 성공: {successful}/{len(tasks)}")
        logger.info(f"  ❌ 실패: {failed}/{len(tasks)}")
        logger.info(f"  ⏱️ 총 시간: {elapsed/60:.1f}분")
        
        return successful, failed
    
    def verify_output(self):
        """생성된 LD score 파일들 검증"""
        logger.info("🔍 생성된 LD score 파일 검증 중...")
        
        total_expected = len(self.datasets) * 22  # 8 * 22 = 176
        total_found = 0
        
        for dataset in self.datasets:
            dataset_count = 0
            for chr_num in range(1, 23):
                ldscore_file = self.output_dir / f"{dataset}.{chr_num}.l2.ldscore.gz"
                if ldscore_file.exists():
                    dataset_count += 1
                    total_found += 1
            
            logger.info(f"  {dataset}: {dataset_count}/22 파일")
        
        logger.info(f"📊 총 LD score 파일: {total_found}/{total_expected}")
        
        if total_found == total_expected:
            logger.info("✅ 모든 LD score 파일이 성공적으로 생성되었습니다!")
        else:
            logger.warning(f"⚠️ {total_expected - total_found}개 파일이 누락되었습니다.")
        
        return total_found, total_expected

def main():
    """메인 실행 함수"""
    calculator = NewLDScoreCalculator()
    
    # Step 1: Check annotations
    total_annots, missing_annots = calculator.check_annotations()
    
    if len(missing_annots) > 20:  # Too many missing
        logger.error("❌ 너무 많은 annotation 파일이 누락되었습니다. 먼저 annotation을 생성하세요.")
        return False
    
    # Step 2: Calculate LD scores
    successful, failed = calculator.calculate_all_ld_scores(max_workers=2)
    
    # Step 3: Verify output
    found, expected = calculator.verify_output()
    
    # Final status
    if found >= expected * 0.9:  # 90% success rate
        logger.info("🎯 LD Score 계산이 성공적으로 완료되었습니다!")
        return True
    else:
        logger.error("❌ LD Score 계산이 불완전합니다. 오류를 확인하세요.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)