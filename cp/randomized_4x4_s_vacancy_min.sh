#!/bin/bash
#SBATCH -N 3
#SBATCH --ntasks-per-node=64
#SBATCH -o slurm.txt
#SBATCH -e slurm.txt
#SBATCH --account=p376-23-1

module load QuantumESPRESSO/7.2-intel-2022a

JOB_NAME='randomized_4x4_cp'


DT=`date +"%d-%b-%Y_%H-%M-%S"`
OUT_DIR="out/${JOB_NAME}_${DT}"
INP_FILE="${OUT_DIR}/input.in"
OUTPUT_FILE="${OUT_DIR}/output.out"

mkdir -p "$OUT_DIR"


cat > $INP_FILE << EOF

&CONTROL
  calculation = 'cp'
  iprint = 1
  isave = 1
  ndr = 50
  ndw = 50
  nstep = 1
  outdir = '${OUT_DIR}'
  prefix = 'mos2'
  pseudo_dir = '../pw/pseudo/'
  restart_mode = 'from_scratch'
  verbosity = 'high'
/
&SYSTEM
  ecutwfc = 50
  ecutrho = 400
  ibrav = 0
  nat = 47
  ntyp = 2
  nr1b = 24,
  nr2b = 24,
  nr3b = 24
/
&ELECTRONS
  electron_dynamics = 'cg'
  electron_maxstep = 100
  conv_thr = 0.000001
/
&CELL
/
ATOMIC_SPECIES
Mo    95.95    Mo_ONCV_PBE-1.0.oncvpsp.upf
S     32.06    s_pbe_v1.4.uspp.F.UPF

ATOMIC_POSITIONS angstrom
Mo    1.622386100708352    0.434345914381836    25.005028843977584
S    1.5232743761662175    2.2928016633877863    26.55956933874013
S    1.5869430537058205    2.2938201508045983    23.570061026200413
Mo    0.05010733408664491    3.2268044517661316    25.041452115973623
S    0.06788706148364636    5.049095360825358    26.610515658193314
S    -0.0611902271467819    5.058780087209127    23.516141184184516
Mo    -1.627643364046756    5.937411139464601    24.98111498613477
S    -1.6346169177298044    7.77833682653642    26.60721473683568
S    -1.5814288850939366    7.688908846008571    23.494508758981286
Mo    -3.155433315423902    8.568951217733558    24.99390598330971
S    -3.1516834540901097    10.57948561342733    26.631065357971796
S    -3.093970731357702    10.488407706630964    23.42824290121949
Mo    4.720010591167794    0.533108605249183    24.98348605769457
S    4.833075380367041    2.276899861840346    26.57746031552407
S    4.763832400944368    2.3900323103726864    23.469297581421888
Mo    3.121199927118948    3.1991546535868847    25.045180434990566
S    3.212960241588171    5.0322520424805735    26.560378172430045
S    3.2046779750690617    5.000098772714196    23.504480872234492
Mo    1.576730285767859    5.9002047856208195    24.985546426002653
S    1.5904572183627188    7.730531866517308    26.47454497470717
S    1.6360959981893908    7.762171545752723    23.466582852964933
Mo    -0.028481862687732475    8.677082333221822    25.04621268204037
S    0.013122935483061912    10.564259460465442    26.662651080544435
S    0.02181718079280175    10.54832561436618    23.40268421327704
Mo    7.879137835806739    0.534227943874585    24.958859010638342
S    7.887965166415952    2.2316491350169216    23.402252779890773
Mo    6.393871653142948    3.0901743400476867    24.88531027491434
S    6.325478831744737    4.98724960489446    26.682196865460593
S    6.257500887329686    5.017056139622887    23.399163141288266
Mo    4.82329442964417    5.919132371938357    25.08869547667515
S    4.727546518404341    7.717439760542513    26.590669311843648
S    4.7447169285488044    7.847595667541017    23.506674150757657
Mo    3.145492441798787    8.678251241204704    25.054092004722257
S    3.2174149972018053    10.566415061590142    26.509908202302277
S    3.114623145118107    10.503395343823607    23.372629913693203
Mo    11.029328623649318    0.5750441603812542    24.969365572958274
S    10.89366528078753    2.3381988360292905    26.57667550347642
S    10.953402479651688    2.280982019385105    23.357428265919026
Mo    9.42047987911675    3.1627871760594957    24.96262283465893
S    9.424259606100403    4.937401175193002    26.56778697355968
S    9.555367525035278    4.997960647899881    23.442587943753885
Mo    7.969440399729769    5.937251206626313    25.001625772900496
S    7.891155432695152    7.719782566971478    26.642134202197273
S    7.843647128027467    7.713592215823791    23.40781916277521
Mo    6.3540770731980345    8.693617118860244    24.9701857493089
S    6.278979806658644    10.55077440430077    26.626090809176045
S    6.321972066014423    10.55187576411829    23.394138074445493

CELL_PARAMETERS angstrom
12.651922535508914    -4.3663250128404126e-07    0.0
-6.3259616458894525    10.956886648195539    0.0
0.0    0.0    50.00000719928987
EOF

mpirun -np 192 cp.x -inp $INP_FILE > $OUTPUT_FILE
