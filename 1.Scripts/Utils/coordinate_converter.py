#!/usr/bin/env python3
"""
염색체 좌표계 변환 유틸리티
========================
rn7 → hg38 → hg19 좌표계 변환
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import subprocess
import tempfile
import os
from typing import Optional, Dict, Any, Tuple
import warnings
warnings.filterwarnings('ignore')

logger = logging.getLogger(__name__)


class CoordinateConverter:
    """염색체 좌표계 변환 클래스"""
    
    def __init__(self):
        self.base_dir = Path(".")
        self.reference_dir = self.base_dir / "0_data" / "reference"
        self.liftover_dir = self.reference_dir / "liftover_data"
        
        # LiftOver 체인 파일 경로
        self.chain_files = {
            'rn7_to_hg38': self.liftover_dir / "rn7ToHg38.over.chain.gz",
            'hg38_to_hg19': self.liftover_dir / "hg38ToHg19.over.chain.gz"
        }
        
        # LiftOver 실행 파일 경로
        self.liftover_cmd = self.liftover_dir / "liftOver"
        
        logger.info("좌표계 변환기 초기화")
    
    def check_liftover_availability(self) -> Dict[str, bool]:
        """LiftOver 도구 및 체인 파일 확인"""
        status = {
            'liftover_executable': self.liftover_cmd.exists() and os.access(self.liftover_cmd, os.X_OK),
            'rn7_to_hg38_chain': self.chain_files['rn7_to_hg38'].exists(),
            'hg38_to_hg19_chain': self.chain_files['hg38_to_hg19'].exists()
        }
        
        logger.info(f"LiftOver 가용성 체크:")
        for component, available in status.items():
            logger.info(f"  {component}: {'✅' if available else '❌'}")
        
        return status
    
    def detect_coordinate_system(self, bed_file: Path) -> str:
        """BED 파일의 좌표계 추정"""
        logger.info(f"좌표계 추정 중: {bed_file}")
        
        # 샘플 영역 읽기
        sample_df = pd.read_csv(bed_file, sep='\t', header=None, nrows=100,
                               names=['chr', 'start', 'end', 'name'])
        
        # 염색체 1번의 좌표 범위 확인
        chr1_data = sample_df[sample_df['chr'] == 'chr1']
        if len(chr1_data) == 0:
            return "unknown"
        
        max_pos = chr1_data['end'].max()
        min_pos = chr1_data['start'].min()
        
        logger.info(f"  염색체 1번 좌표 범위: {min_pos:,} - {max_pos:,}")
        
        # 휴리스틱 기반 좌표계 추정
        if max_pos > 240000000:  # hg19/hg38 (약 249Mb)
            if max_pos > 248000000:
                return "hg19_or_hg38"
            else:
                return "hg19_or_hg38"
        elif max_pos > 260000000:  # rn7 (약 285Mb)
            return "rn7"
        else:
            return "unknown"
    
    def convert_bed_coordinates(self, input_bed: Path, output_bed: Path,
                               from_assembly: str, to_assembly: str) -> bool:
        """BED 파일 좌표 변환"""
        logger.info(f"좌표 변환: {from_assembly} → {to_assembly}")
        
        # 체인 파일 매핑
        chain_mapping = {
            ('rn7', 'hg38'): 'rn7_to_hg38',
            ('hg38', 'hg19'): 'hg38_to_hg19'
        }
        
        chain_key = (from_assembly, to_assembly)
        if chain_key not in chain_mapping:
            logger.error(f"지원하지 않는 변환: {from_assembly} → {to_assembly}")
            return False
        
        chain_file = self.chain_files[chain_mapping[chain_key]]
        if not chain_file.exists():
            logger.error(f"체인 파일을 찾을 수 없습니다: {chain_file}")
            return False
        
        # 임시 파일 생성
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_input:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_unmapped:
                try:
                    # LiftOver 실행
                    cmd = [
                        str(self.liftover_cmd),
                        str(input_bed),
                        str(chain_file),
                        str(output_bed),
                        tmp_unmapped.name
                    ]
                    
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    
                    if result.returncode == 0:
                        logger.info(f"좌표 변환 성공: {output_bed}")
                        
                        # 변환 통계
                        if output_bed.exists():
                            original_count = sum(1 for _ in open(input_bed))
                            converted_count = sum(1 for _ in open(output_bed))
                            success_rate = converted_count / original_count * 100
                            
                            logger.info(f"  변환 성공률: {converted_count}/{original_count} ({success_rate:.1f}%)")
                        
                        return True
                    else:
                        logger.error(f"LiftOver 실행 실패: {result.stderr}")
                        return False
                        
                finally:
                    # 임시 파일 정리
                    try:
                        os.unlink(tmp_input.name)
                        os.unlink(tmp_unmapped.name)
                    except:
                        pass
    
    def setup_liftover_environment(self):
        """LiftOver 환경 설정"""
        logger.info("LiftOver 환경 설정 중...")
        
        # 디렉토리 생성
        self.liftover_dir.mkdir(parents=True, exist_ok=True)
        
        # 필요한 파일들 다운로드 URL
        download_urls = {
            'liftOver': 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver',
            'rn7ToHg38.over.chain.gz': 'http://hgdownload.soe.ucsc.edu/goldenPath/rn7/liftOver/rn7ToHg38.over.chain.gz',
            'hg38ToHg19.over.chain.gz': 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz'
        }
        
        logger.info("자동 다운로드는 지원하지 않습니다.")
        logger.info("다음 파일들을 수동으로 다운로드하세요:")
        for filename, url in download_urls.items():
            target_path = self.liftover_dir / filename
            logger.info(f"  {filename}: {url}")
            logger.info(f"    → {target_path}")
        
        return False


class EnhancedDataManager:
    """좌표계 변환이 포함된 데이터 매니저"""
    
    def __init__(self, enhancer_file: Path):
        self.enhancer_file = enhancer_file
        self.converter = CoordinateConverter()
        self.cache_dir = Path("coordinate_conversion_cache")
        self.cache_dir.mkdir(exist_ok=True)
        
        logger.info(f"향상된 데이터 매니저 초기화: {enhancer_file}")
    
    def get_converted_enhancer_data(self, target_assembly: str = "hg19") -> Optional[pd.DataFrame]:
        """좌표계 변환된 enhancer 데이터 반환"""
        
        # 캐시 파일 확인
        cache_file = self.cache_dir / f"{self.enhancer_file.stem}_{target_assembly}.pkl"
        if cache_file.exists():
            logger.info(f"캐시된 변환 데이터 로딩: {cache_file}")
            return pd.read_pickle(cache_file)
        
        # 원본 좌표계 추정
        original_assembly = self.converter.detect_coordinate_system(self.enhancer_file)
        logger.info(f"추정된 원본 좌표계: {original_assembly}")
        
        if original_assembly == target_assembly:
            logger.info("좌표계 변환이 필요하지 않습니다.")
            enhancer_df = pd.read_csv(self.enhancer_file, sep='\t', header=None,
                                    names=['CHR', 'START', 'END', 'NAME'])
            # 데이터 정리
            enhancer_df['CHR'] = enhancer_df['CHR'].str.replace('chr', '')
            numeric_mask = enhancer_df['CHR'].str.isnumeric()
            enhancer_df = enhancer_df[numeric_mask].copy()
            enhancer_df['CHR'] = enhancer_df['CHR'].astype(int)
            enhancer_df = enhancer_df[enhancer_df['CHR'].isin(range(1, 23))]
            
            # 캐시 저장
            enhancer_df.to_pickle(cache_file)
            return enhancer_df
        
        # 좌표계 변환 필요
        logger.warning("좌표계 변환이 필요하지만 LiftOver 도구가 설정되지 않았습니다.")
        logger.info("현재는 원본 데이터를 그대로 사용합니다.")
        logger.info("정확한 분석을 위해서는 좌표계 변환을 수행하세요.")
        
        # 임시로 원본 데이터 반환
        enhancer_df = pd.read_csv(self.enhancer_file, sep='\t', header=None,
                                names=['CHR', 'START', 'END', 'NAME'])
        # 데이터 정리
        enhancer_df['CHR'] = enhancer_df['CHR'].str.replace('chr', '')
        numeric_mask = enhancer_df['CHR'].str.isnumeric()
        enhancer_df = enhancer_df[numeric_mask].copy()
        enhancer_df['CHR'] = enhancer_df['CHR'].astype(int)
        enhancer_df = enhancer_df[enhancer_df['CHR'].isin(range(1, 23))]
        
        return enhancer_df


def main():
    """테스트 및 검증"""
    converter = CoordinateConverter()
    
    # 환경 체크
    status = converter.check_liftover_availability()
    
    if not all(status.values()):
        logger.warning("LiftOver 환경이 완전하지 않습니다.")
        converter.setup_liftover_environment()
    
    # 샘플 파일 좌표계 추정
    sample_files = [
        Path("0_data/raw/cleaned_data/Olig_cleaned.bed"),
        Path("0_data/raw/unique_data/Olig_unique.bed")
    ]
    
    for sample_file in sample_files:
        if sample_file.exists():
            assembly = converter.detect_coordinate_system(sample_file)
            logger.info(f"{sample_file}: {assembly}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()