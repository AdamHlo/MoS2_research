#!/bin/bash
#SBATCH -N 4
#SBATCH --ntasks-per-node=64
#SBATCH -o slurm.txt
#SBATCH -e slurm.txt
#SBATCH --account=p376-23-1

module load QuantumESPRESSO/7.2-intel-2022a


JOB_NAME='4x4_s_vacancy_relaxation'


DT=`date +"%d-%b-%Y_%H-%M-%S"`
OUT_DIR="out/${JOB_NAME}_${DT}"
INP_FILE="${OUT_DIR}/input.in"
OUTPUT_FILE="${OUT_DIR}/output.out"

mkdir -p "$OUT_DIR"


cat > $INP_FILE << EOF
&CONTROL
  calculation = 'vc-relax'
  outdir = './out/'
  prefix = 'MoS2'
  pseudo_dir = './pseudo/'
  restart_mode = 'from_scratch'
  verbosity = 'high'
  wf_collect = .true.
/
&SYSTEM
  ecutwfc = 70
  ecutrho = 560
  ibrav = 0
  nat = 47
  ntyp = 2
  nosym = .true.
/
&ELECTRONS
  conv_thr =   1.0000000000d-07
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
Mo    1.5904639853567222    0.45912740504103866    25.0000017998224
S    1.5904639853567106    2.2956370252052567    26.561494720962884
S    1.5904639853567106    2.2956370252052567    23.4385088786824
Mo    4.440892098500626e-14    3.213891835287337    25.0000017998224
S    3.2862601528904634e-14    5.050401455451555    26.561494720962884
S    3.2862601528904634e-14    5.050401455451555    23.4385088786824
Mo    -1.5904639853566334    5.968656265533635    25.0000017998224
S    -1.590463985356645    7.805165885697853    26.561494720962884
S    -1.590463985356645    7.805165885697853    23.4385088786824
Mo    -3.180927970713311    8.723420695779934    25.0000017998224
S    -3.1809279707133227    10.55993031594415    26.561494720962884
S    -3.1809279707133227    10.55993031594415    23.4385088786824
Mo    -4.771391956069989    11.478185126026233    25.0000017998224
S    -4.77139195607    13.314694746190451    26.561494720962884
S    -4.77139195607    13.314694746190451    23.4385088786824
Mo    -6.361855941426667    14.232949556272532    25.0000017998224
S    -6.361855941426678    16.069459176436748    26.561494720962884
S    -6.361855941426678    16.069459176436748    23.4385088786824
Mo    4.77139195607012    0.4591274050410386    25.0000017998224
S    4.77139195607011    2.2956370252052567    26.561494720962884
S    4.77139195607011    2.2956370252052567    23.4385088786824
Mo    3.180927970713443    3.213891835287337    25.0000017998224
S    3.1809279707134315    5.050401455451555    26.561494720962884
S    3.1809279707134315    5.050401455451555    23.4385088786824
Mo    1.5904639853567653    5.968656265533635    25.0000017998224
S    1.5904639853567537    7.805165885697853    26.561494720962884
Mo    8.748557434046234e-14    8.723420695779934    25.0000017998224
S    7.593925488436071e-14    10.55993031594415    26.561494720962884
S    7.593925488436071e-14    10.55993031594415    23.4385088786824
Mo    -1.5904639853565903    11.478185126026233    25.0000017998224
S    -1.5904639853566018    13.314694746190451    26.561494720962884
S    -1.5904639853566018    13.314694746190451    23.4385088786824
Mo    -3.1809279707132685    14.232949556272532    25.0000017998224
S    -3.18092797071328    16.069459176436748    26.561494720962884
S    -3.18092797071328    16.069459176436748    23.4385088786824
Mo    7.952319926783519    0.45912740504103855    25.0000017998224
S    7.952319926783508    2.2956370252052567    26.561494720962884
S    7.952319926783508    2.2956370252052567    23.4385088786824
Mo    6.361855941426842    3.213891835287337    25.0000017998224
S    6.36185594142683    5.050401455451555    26.561494720962884
S    6.36185594142683    5.050401455451555    23.4385088786824
Mo    4.771391956070164    5.968656265533635    25.0000017998224
S    4.771391956070152    7.805165885697853    26.561494720962884
S    4.771391956070152    7.805165885697853    23.4385088786824
Mo    3.180927970713486    8.723420695779934    25.0000017998224
S    3.1809279707134746    10.55993031594415    26.561494720962884
S    3.1809279707134746    10.55993031594415    23.4385088786824
Mo    1.5904639853568083    11.478185126026233    25.0000017998224
S    1.5904639853567968    13.314694746190451    26.561494720962884
S    1.5904639853567968    13.314694746190451    23.4385088786824
Mo    1.305622276959184e-13    14.232949556272532    25.0000017998224
S    1.1901590823981678e-13    16.069459176436748    26.561494720962884
S    1.1901590823981678e-13    16.069459176436748    23.4385088786824
Mo    11.133247897496918    0.45912740504103844    25.0000017998224
S    11.133247897496908    2.2956370252052567    26.561494720962884
S    11.133247897496908    2.2956370252052567    23.4385088786824
Mo    9.54278391214024    3.213891835287337    25.0000017998224
S    9.54278391214023    5.050401455451555    26.561494720962884
S    9.54278391214023    5.050401455451555    23.4385088786824
Mo    7.952319926783563    5.968656265533635    25.0000017998224
S    7.952319926783551    7.805165885697853    26.561494720962884
S    7.952319926783551    7.805165885697853    23.4385088786824
Mo    6.361855941426885    8.723420695779934    25.0000017998224
S    6.361855941426874    10.55993031594415    26.561494720962884
S    6.361855941426874    10.55993031594415    23.4385088786824
Mo    4.771391956070207    11.478185126026233    25.0000017998224
S    4.771391956070196    13.314694746190451    26.561494720962884
S    4.771391956070196    13.314694746190451    23.4385088786824
Mo    3.1809279707135296    14.232949556272532    25.0000017998224
S    3.180927970713518    16.069459176436748    26.561494720962884
S    3.180927970713518    16.069459176436748    23.4385088786824
Mo    14.314175868210317    0.4591274050410384    25.0000017998224
S    14.314175868210306    2.2956370252052563    26.561494720962884
S    14.314175868210306    2.2956370252052563    23.4385088786824
Mo    12.723711882853639    3.2138918352873365    25.0000017998224
S    12.723711882853628    5.050401455451555    26.561494720962884
S    12.723711882853628    5.050401455451555    23.4385088786824
Mo    11.133247897496961    5.968656265533635    25.0000017998224
S    11.13324789749695    7.805165885697853    26.561494720962884
S    11.13324789749695    7.805165885697853    23.4385088786824
Mo    9.542783912140283    8.723420695779934    25.0000017998224
S    9.542783912140273    10.55993031594415    26.561494720962884
S    9.542783912140273    10.55993031594415    23.4385088786824
Mo    7.9523199267836056    11.478185126026233    25.0000017998224
S    7.952319926783594    13.314694746190451    26.561494720962884
S    7.952319926783594    13.314694746190451    23.4385088786824
Mo    6.361855941426928    14.232949556272532    25.0000017998224
S    6.361855941426916    16.069459176436748    26.561494720962884
S    6.361855941426916    16.069459176436748    23.4385088786824
Mo    17.495103838923715    0.45912740504103833    25.0000017998224
S    17.495103838923704    2.2956370252052563    26.561494720962884
S    17.495103838923704    2.2956370252052563    23.4385088786824
Mo    15.904639853567037    3.2138918352873365    25.0000017998224
S    15.904639853567026    5.050401455451555    26.561494720962884
S    15.904639853567026    5.050401455451555    23.4385088786824
Mo    14.31417586821036    5.968656265533635    25.0000017998224
S    14.314175868210349    7.805165885697853    26.561494720962884
S    14.314175868210349    7.805165885697853    23.4385088786824
Mo    12.723711882853681    8.723420695779934    25.0000017998224
S    12.72371188285367    10.55993031594415    26.561494720962884
S    12.72371188285367    10.55993031594415    23.4385088786824
Mo    11.133247897497004    11.478185126026233    25.0000017998224
S    11.133247897496993    13.314694746190451    26.561494720962884
S    11.133247897496993    13.314694746190451    23.4385088786824
Mo    9.542783912140326    14.232949556272532    25.0000017998224
S    9.542783912140315    16.069459176436748    26.561494720962884
S    9.542783912140315    16.069459176436748    23.4385088786824

CELL_PARAMETERS angstrom
19.085567824280393    -4.003125588536473e-16    0.0
-9.542783912140067    16.52858658147779    0.0
0.0    0.0    50.000003599644806
EOF


mpirun -np 256 pw.x -inp $INP_FILE > $OUTPUT_FILE
