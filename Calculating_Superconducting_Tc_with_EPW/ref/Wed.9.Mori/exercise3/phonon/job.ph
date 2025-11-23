#!/bin/bash
#SBATCH -J job.ph             # Job name
#SBATCH -N 1                  # Total # of nodes
#SBATCH --ntasks-per-node 24
#SBATCH -t 00:30:00           # Run time (hh:mm:ss)
#SBATCH -A DMR23030
#SBATCH -p skx
#SBATCH --reservation=NSF_Summer_School_Wed
 	
# Launch MPI code...
export PATHQE=/work2/05193/sabyadk/stampede3/EPWSchool2024/q-e
	
ibrun $PATHQE/bin/pw.x -nk 6 -in scf.in > scf.out
ibrun $PATHQE/bin/ph.x -nk 6 -in ph.in > ph.out

exit
