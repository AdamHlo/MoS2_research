#!/bin/bash
#SBATCH -N 4
#SBATCH --ntasks-per-node=64
#SBATCH -o slurm_electron_minimization.txt
#SBATCH -e slurm_electron_minimization.txt

module load QuantumESPRESSO/7.2-intel-2022a

mpirun -np 256 pw.x -inp relaxation.in > relaxation.out
