#!/bin/bash
#SBATCH --job-name=energetics_Diffusion
#SBATCH --output=%j.out
#SBATCH --partition=debug
#SBATCH --ntasks=8
#SBATCH --exclusive

export QT_QPA_PLATFORM=offscreen

# ===== Clean old files =====
rm -rf ./runs
rm -rf *.out
rm -rf log

# ===== Write start time =====
echo "===============================" >> log
echo "Job start time: $(date '+%Y-%m-%d %H:%M:%S')" >> log
echo "Running on: $(hostname)" >> log
echo "===============================" >> log

start_time=$(date +%s)

# ===== Run program =====
python Diffusion_v1.2.py >> log

# ===== Write end time =====
end_time=$(date +%s)
runtime=$((end_time - start_time))

echo "===============================" >> log
echo "Job end time: $(date '+%Y-%m-%d %H:%M:%S')" >> log
echo "Total runtime: ${runtime} seconds" >> log
echo "===============================" >> log