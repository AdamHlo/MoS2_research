#!/bin/bash
#SBATCH -N 8
#SBATCH --ntasks-per-node=64
#SBATCH -o slurm.txt
#SBATCH -e slurm.txt
#SBATCH --account=p376-23-1

module load QuantumESPRESSO/7.2-intel-2022a

mpirun -np 512 pw.x -inp md_defect.in > out/md_defect.out
