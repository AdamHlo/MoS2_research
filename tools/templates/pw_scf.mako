<%text>#!/bin/bash
#SBATCH -N 4
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
  <%text>outdir = '${OUT_DIR}'</%text>
  prefix = 'mos2'
  pseudo_dir = './pseudo/'
  restart_mode = 'from_scratch'
  verbosity = 'high'
  wf_collect = .true.
/
&SYSTEM
  ecutwfc = ${ecutwfc}
  ecutrho = ${ecutrho}
  ibrav = ${ibrav}
  nat = ${nat}
  ntyp = ${ntyp}
  nosym = ${nosym}
/
&ELECTRONS
  conv_thr = ${conv_thr}
  diagonalization = 'cg'
/
&CELL
/
ATOMIC_SPECIES
Mo    95.95    Mo_ONCV_PBE-1.0.oncvpsp.upf
S     32.06    s_pbe_v1.4.uspp.F.UPF

${atomic_configuration}

K_POINTS gamma
EOF

mpirun -np 256 pw.x -inp $INP_FILE > $OUTPUT_FILE
