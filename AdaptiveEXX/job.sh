#!/bin/bash
#SBATCH --job-name=AdaptiveEXX
#SBATCH --output=%j.out
#SBATCH --partition=debug
#SBATCH --ntasks=8
#SBATCH --exclusive

julia AdaptiveEXX.jl && python plot_v4.py
