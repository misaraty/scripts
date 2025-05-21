#!/bin/bash
#SBATCH --job-name=epw_test
#SBATCH --partition=cpu3_q
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

cd ./phonon &&
mpirun pw.x -in scf.in > scf.out &&
mpirun ph.x -in ph.in > ph.out &&
python /opt/ohpc/pub/apps/q-e-qe-7.3/EPW/bin/pp.py pb &&
cd ../epw &&
mkdir pb.save &&
cp ../phonon/pb.save/charge-density.dat pb.save/ &&
cp ../phonon/pb.save/data-file-schema.xml pb.save/ &&
mpirun -np 1 pw.x -in nscf.in > nscf.out &&
mpirun -np 1 epw.x -in epw1.in > epw1.out &&
mpirun -np 1 epw.x -in epw2.in > epw2.out &&
python plot.py &&
cd ../
