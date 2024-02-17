#!/bin/bash
#SBATCH -N 4
#SBATCH --ntasks-per-node=64
#SBATCH -o slurm.txt
#SBATCH -e slurm.txt

module load QuantumESPRESSO/7.2-intel-2022a

mpirun -np 256 pw.x -inp defect_relaxation.in > defect_relaxation.out
