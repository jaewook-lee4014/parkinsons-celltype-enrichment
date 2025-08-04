#!/bin/bash
#SBATCH --job-name=parallel_ldsc
#SBATCH --account=kcl
#SBATCH --partition=nmes_gpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --output=parallel_ldsc_%j.out
#SBATCH --error=parallel_ldsc_%j.err

# HPC 노드에서 병렬 LD Score 계산 실행

echo "🚀 병렬 LD Score 계산 시작"
echo "작업 ID: $SLURM_JOB_ID"
echo "노드: $SLURM_NODELIST"
echo "CPU 코어: $SLURM_CPUS_PER_TASK"
echo "메모리: 4GB"
echo "시간: $(date)"
echo "=============================================="

# 작업 디렉토리로 이동
cd /cephfs/volumes/hpc_data_prj/eng_waste_to_protein/ae035a41-20d2-44f3-aa46-14424ab0f6bf/repositories/bomin

# Python 환경 활성화 (필요시)
# source activate ldsc_env

# 병렬 계산 실행
echo "🧬 병렬 LD Score 계산 실행 중..."
python3 parallel_ldsc_calculation.py

# 실행 결과 확인
if [ $? -eq 0 ]; then
    echo "✅ 병렬 계산 성공!"
else
    echo "❌ 병렬 계산 실패!"
fi

echo "=============================================="
echo "완료 시간: $(date)"