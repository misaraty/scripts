#!/bin/bash
#SBATCH --job-name=GoalDFT.jl_v5
#SBATCH --output=%j.out
#SBATCH --partition=debug
#SBATCH --ntasks=8
#SBATCH --exclusive

julia GoalDFT.jl && python plot_v3.py
