#!/bin/bash
#SBATCH -N 4
#SBATCH --ntasks-per-node=64
#SBATCH -o slurm.txt
#SBATCH -e slurm.txt
#SBATCH --account=p376-23-1

module load QuantumESPRESSO/7.2-intel-2022a


JOB_NAME='4x4_s_vacancy_cp_minimization'


DT=`date +"%d-%b-%Y_%H-%M-%S"`
OUT_DIR="out/${JOB_NAME}_${DT}"
INP_FILE="${OUT_DIR}/input.in"
OUTPUT_FILE="${OUT_DIR}/output.out"

mkdir -p "$OUT_DIR"


cat > $INP_FILE << EOF
&CONTROL
  calculation = 'cp'
  iprint = 1
  isave = 5
  max_seconds = 43200
  ndr = 50
  ndw = 50
  nstep = 100
  outdir = '$OUT_DIR'
  prefix = 'MoS2'
  pseudo_dir = '../pw/pseudo/'
  restart_mode = 'from_scratch'
  verbosity = 'high'
/
&SYSTEM
    ibrav = 0,
    nat = 47,   ! number of atoms in the unit cell
    ntyp = 2,   ! number of atom types in the unit cell
    ecutwfc = 40,  ! kin. E cutoff for wavefunctions
    ecutrho = 320,
    smearing = 'gaussian',
/
&ELECTRONS
  electron_dynamics = 'cg'
  electron_maxstep = 100
/
&IONS
ion_dynamics = 'none'
/
ATOMIC_SPECIES
Mo    95.95    Mo_ONCV_PBE-1.0.oncvpsp.upf
S     32.06    s_pbe_v1.4.uspp.F.UPF

ATOMIC_POSITIONS angstrom
Mo    1.5814904846075442    0.4624789663310789    25.00980583549507
S    1.581489750369667    2.2784360124190757    26.587520909543763
S    1.5814906496809087    2.292850466351797    23.451996055193053
Mo    -0.01135604387842685    3.194555873535763    25.01812506453629
S    -0.009903989027045009    5.016188627050577    26.591681632711392
S    0.00020123135506283494    5.022023820936632    23.459322892998724
Mo    -1.5882098750070923    5.925749070495656    25.01812592152558
S    -1.5851693787794294    7.763252989024683    26.587523506334435
S    -1.572685890649262    7.756045863392389    23.451993951483306
Mo    -3.1578326703464468    8.671229882391595    25.009805980231395
S    -3.144480948937866    10.51102991258774    26.56353882260719
S    -3.1611504564531177    10.501405367427175    23.441072388482894
Mo    4.770178603336136    0.47137699468734995    24.986135845350944
S    4.824285327639041    2.2803607527969882    26.552977541072504
S    4.753120228407201    2.3035441820067715    23.416833020843484
Mo    3.1743365627331457    3.194555196278774    25.01812627382427
S    3.1728830251550475    5.016189400012678    26.591680972739134
S    3.162777838590276    5.022023773128949    23.459324742730495
Mo    1.581489070298558    5.934981464464826    25.031104022111553
S    1.5814902532904027    7.772564281607689    26.591682614714216
S    1.5814896474919868    7.7608959448119395    23.45932319803891
Mo    0.004635868585480184    8.684638007603581    25.018126737503888
S    0.03789578873961902    10.570631292494859    26.552979363689285
S    0.02238970476948166    10.497407841954457    23.41683039293036
Mo    7.907451307531323    0.5093172200550129    24.946586319655978
S    7.9074525009804395    2.28268423375103    23.366968393421967
Mo    6.371672740256277    3.1693659848662477    24.94658759733222
S    6.367879981717273    4.953946064492534    26.552977208811843
S    6.312220136064212    5.003985180541051    23.416833817539633
Mo    4.751190065552052    5.925749046611082    25.018126778953153
S    4.748148970110832    7.763253575501975    26.587522387135696
S    4.73566588232878    7.756045815471187    23.451995954272316
Mo    3.158343662502187    8.684638371923304    25.018126599641672
S    3.1250832110930444    10.570631504759998    26.55297922514066
S    3.14058825408825    10.497407528735566    23.416829220289312
Mo    11.044724733763504    0.4713768290773028    24.986134820408818
S    10.990618287074891    2.280361237582921    26.552977276306574
S    11.061783172072413    2.30354448779322    23.41683012061734
Mo    9.443231585009864    3.169366562922714    24.946586346408747
S    9.447023601948494    4.953945495355315    26.552977171492785
S    9.502682609015956    5.003985633433151    23.416832238715944
Mo    7.907451823440337    5.905296820235694    24.98613632793092
S    7.907452009272781    7.739766711928713    26.563539462629613
S    7.907451502541246    7.759014924365625    23.441073159814696
Mo    6.3208136236444314    8.67123060317483    25.0098067236859
S    6.307460941276415    10.511030635948972    26.563539388163104
S    6.324131564043984    10.501404821803373    23.441071924695617

CELL_PARAMETERS angstrom
12.651922535508914    -4.3663250128404126e-07    0.0
-6.3259616458894525    10.956886648195539    0.0
0.0    0.0    50.00000719928987
EOF


mpirun -np 256 cp.x -inp $INP_FILE > $OUTPUT_FILE
