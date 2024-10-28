#!/bin/bash
#SBATCH --job-name=qe_test
#SBATCH --partition=cpu3_q
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
							  
qe='mpirun /opt/ohpc/pub/apps/qe-7.3/bin/pw.x'

# # $qe <relax.in> relax.out
# $qe <vc-relax.in> vc-relax.out
# $qe <scf.in> scf.out
$qe <nscf.in> nscf.out