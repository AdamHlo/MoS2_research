<%text>#!/bin/bash
#SBATCH -N 3
#SBATCH --ntasks-per-node=64
#SBATCH -o slurm.txt
#SBATCH -e slurm.txt
#SBATCH --account=p376-23-1

module load QuantumESPRESSO/7.2-intel-2022a</%text>

JOB_NAME='${job_name}'

<%text>
DT=`date +"%d-%b-%Y_%H-%M-%S"`
OUT_DIR="out/${JOB_NAME}_${DT}"
INP_FILE="${OUT_DIR}/input.in"
OUTPUT_FILE="${OUT_DIR}/output.out"

mkdir -p "$OUT_DIR"


cat > $INP_FILE << EOF
</%text>
&CONTROL
  calculation = 'scf'
  iprint = 1
  isave = 1
  ndr = ${ndr}
  ndw = ${ndw}
  nstep = 1
  <%text>outdir = '${OUT_DIR}'</%text>
  prefix = 'mos2'
  pseudo_dir = '../pw/pseudo/'
  restart_mode = 'from_scratch'
  verbosity = 'high'
/
&SYSTEM
  ecutwfc = ${ecutwfc}
  ecutrho = ${ecutrho}
  ibrav = ${ibrav}
  nat = ${nat}
  ntyp = ${ntyp}
  nr1b = 24,
  nr2b = 24,
  nr3b = 24
/
&ELECTRONS
  orthogonalization = 'ortho'
  electron_dynamics = 'cg'
  electron_maxstep = 100
  conv_thr = ${conv_thr}
/
&CELL
/
ATOMIC_SPECIES
Mo    95.95    Mo_ONCV_PBE-1.0.oncvpsp.upf
S     32.06    s_pbe_v1.4.uspp.F.UPF

${atomic_configuration}
EOF

mpirun -np 192 cp.x -inp $INP_FILE > $OUTPUT_FILE
