#!/usr/bin/env python3
"""
LiftOver 도구 설정 스크립트
=========================
rn7 → hg38 → hg19 변환을 위한 환경 구축
"""

import os
import subprocess
from pathlib import Path
import urllib.request
import gzip
import shutil
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def setup_liftover_environment():
    """LiftOver 환경 설정"""
    base_dir = Path(".")
    liftover_dir = base_dir / "0_data" / "reference" / "liftover_data"
    liftover_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"LiftOver 환경 설정: {liftover_dir}")
    
    # 필요한 파일들
    files_to_download = {
        'liftOver': {
            'url': 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver',
            'executable': True
        },
        'rn7ToHg38.over.chain.gz': {
            'url': 'http://hgdownload.soe.ucsc.edu/goldenPath/rn7/liftOver/rn7ToHg38.over.chain.gz',
            'executable': False
        },
        'hg38ToHg19.over.chain.gz': {
            'url': 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz',
            'executable': False
        }
    }
    
    # 다운로드 및 설정
    for filename, info in files_to_download.items():
        file_path = liftover_dir / filename
        
        if file_path.exists():
            logger.info(f"✅ {filename} 이미 존재")
            continue
        
        logger.info(f"📥 {filename} 다운로드 중...")
        try:
            urllib.request.urlretrieve(info['url'], file_path)
            
            if info['executable']:
                os.chmod(file_path, 0o755)
                logger.info(f"✅ {filename} 다운로드 및 실행권한 설정 완료")
            else:
                logger.info(f"✅ {filename} 다운로드 완료")
                
        except Exception as e:
            logger.error(f"❌ {filename} 다운로드 실패: {e}")
            logger.info(f"   수동 다운로드 필요: {info['url']}")
    
    # 설치 확인
    liftover_cmd = liftover_dir / "liftOver"
    if liftover_cmd.exists():
        try:
            result = subprocess.run([str(liftover_cmd)], capture_output=True, text=True)
            logger.info("✅ liftOver 도구 설치 확인됨")
        except Exception as e:
            logger.error(f"❌ liftOver 실행 실패: {e}")
    
    return liftover_dir


def convert_bed_file(input_bed: Path, output_bed: Path, 
                    chain_file: Path, liftover_cmd: Path):
    """BED 파일 좌표 변환"""
    logger.info(f"좌표 변환: {input_bed} → {output_bed}")
    
    # 변환되지 않은 영역을 저장할 파일
    unmapped_file = output_bed.with_suffix('.unmapped')
    
    # liftOver 실행
    cmd = [
        str(liftover_cmd),
        str(input_bed),
        str(chain_file),
        str(output_bed),
        str(unmapped_file)
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # 변환 결과 통계
        if input_bed.exists() and output_bed.exists():
            original_count = sum(1 for _ in open(input_bed))
            converted_count = sum(1 for _ in open(output_bed))
            success_rate = converted_count / original_count * 100
            
            logger.info(f"변환 완료: {converted_count}/{original_count} ({success_rate:.1f}%)")
            return True
        else:
            logger.error("변환 결과 파일이 생성되지 않았습니다.")
            return False
            
    except subprocess.CalledProcessError as e:
        logger.error(f"liftOver 실행 실패: {e}")
        logger.error(f"stderr: {e.stderr}")
        return False


def convert_all_enhancer_files():
    """모든 enhancer 파일 좌표 변환"""
    logger.info("=" * 60)
    logger.info("🔄 모든 enhancer 파일 좌표 변환 시작")
    logger.info("=" * 60)
    
    # liftOver 환경 설정
    liftover_dir = setup_liftover_environment()
    liftover_cmd = liftover_dir / "liftOver"
    
    if not liftover_cmd.exists():
        logger.error("liftOver 도구가 설치되지 않았습니다.")
        return False
    
    # 체인 파일들
    rn7_to_hg38_chain = liftover_dir / "rn7ToHg38.over.chain.gz"
    hg38_to_hg19_chain = liftover_dir / "hg38ToHg19.over.chain.gz"
    
    if not rn7_to_hg38_chain.exists() or not hg38_to_hg19_chain.exists():
        logger.error("필요한 체인 파일이 없습니다.")
        return False
    
    # 변환할 파일들
    base_dir = Path(".")
    raw_dir = base_dir / "0_data" / "raw"
    
    # 변환된 파일을 저장할 디렉토리 생성
    converted_dir = base_dir / "0_data" / "processed"
    converted_dir.mkdir(exist_ok=True)
    
    hg38_dir = converted_dir / "hg38_coordinates"
    hg19_dir = converted_dir / "hg19_coordinates"
    hg38_dir.mkdir(exist_ok=True)
    hg19_dir.mkdir(exist_ok=True)
    
    # 변환할 파일 목록
    bed_files = []
    for subdir in ['cleaned_data', 'unique_data']:
        subdir_path = raw_dir / subdir
        if subdir_path.exists():
            bed_files.extend(subdir_path.glob("*.bed"))
    
    logger.info(f"변환할 파일 수: {len(bed_files)}")
    
    success_count = 0
    
    for bed_file in bed_files:
        file_stem = f"{bed_file.parent.name}_{bed_file.stem}"
        
        logger.info(f"\n🔄 {bed_file.name} 변환 중...")
        
        # 1단계: rn7 → hg38
        hg38_file = hg38_dir / f"{file_stem}_hg38.bed"
        logger.info("  1단계: rn7 → hg38")
        
        if convert_bed_file(bed_file, hg38_file, rn7_to_hg38_chain, liftover_cmd):
            
            # 2단계: hg38 → hg19
            hg19_file = hg19_dir / f"{file_stem}_hg19.bed"
            logger.info("  2단계: hg38 → hg19")
            
            if convert_bed_file(hg38_file, hg19_file, hg38_to_hg19_chain, liftover_cmd):
                logger.info(f"✅ {bed_file.name} 변환 완료: {hg19_file}")
                success_count += 1
            else:
                logger.error(f"❌ {bed_file.name} 2단계 변환 실패")
        else:
            logger.error(f"❌ {bed_file.name} 1단계 변환 실패")
    
    logger.info(f"\n🎉 좌표 변환 완료: {success_count}/{len(bed_files)} 파일")
    logger.info(f"변환된 파일 위치: {hg19_dir}")
    
    return success_count == len(bed_files)


if __name__ == "__main__":
    convert_all_enhancer_files()