#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH -o slurm.txt
#SBATCH -e slurm.txt
#SBATCH --account=p376-23-1

module load QuantumESPRESSO/7.2-intel-2022a


JOB_NAME='441_unitc_relaxation'


DT=`date +"%d-%b-%Y_%H-%M-%S"`
INP_FILENAME="input.in"
OUT_FILENAME="output.out"
OUT_DIR="out/${JOB_NAME}_${DT}"

mkdir -p "$OUT_DIR"


cat > $INP_FILENAME << EOF
&CONTROL
  calculation = 'vc-relax'
  outdir = '${OUT_DIR}'
  prefix = 'mos2'
  pseudo_dir = './pseudo/'
  restart_mode = 'from_scratch'
  verbosity = 'high'
  wf_collect = .true.
/
&SYSTEM
  ecutwfc = 70
  ecutrho = 560
  ibrav = 0
  nat = 3
  ntyp = 2
/
&ELECTRONS
  conv_thr =   1.0000000000d-08
/
&IONS
/
&CELL
  cell_dofree='2Dxy'
/
ATOMIC_SPECIES
Mo    95.95    Mo_ONCV_PBE-1.0.oncvpsp.upf
S     32.06    s_pbe_v1.4.uspp.F.UPF

ATOMIC_POSITIONS angstrom
Mo    1.58    0.4561067126598042    25.0
S    1.58    2.280533563299022    26.58
S    1.58    2.280533563299022    23.42

CELL_PARAMETERS angstrom
3.16    0    0
-1.58    2.736640275958826    0
0    0    50

K_POINTS automatic
4 4 1 0 0 0
EOF


mpirun -np 64 pw.x -inp $INP_FILENAME > "${OUT_DIR}/${OUT_FILENAME}"
