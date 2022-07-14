#!/bin/sh
# reminder: from now on, what follows the character # is a comment
####################################################################
#
# define the following variables according to your needs
#
outdir=out/
pseudo_dir=./pseudo
# the following is not actually used:
# espresso_dir=top_directory_of_espresso_package
####################################################################

# make directory for processed data
mkdir -p _data 

# remove files to prevent overwrite
rm -f $processed_output $raw_output

# set up input and output files
input=pwscf.in
raw_output=_data/pwscf_kpoints.out
processed_output=_data/mos2_kpoints.csv

# modify input file for each set of kpoints and send output to si_kpoints.out 
for k_points in 8 9 10 11 12; do

echo "kpoints: ${k_points}"
# self-consistent calculation
cat > $input << EOF
 &CONTROL
  calculation = 'scf'
  etot_conv_thr =   6.0000000000d-05
  forc_conv_thr =   1.0000000000d-04
  outdir = './out/'
  prefix = 'MoS2'
  pseudo_dir = './pseudo/'
  tprnfor = .true.
  tstress = .true.
  verbosity = 'high'
/
&SYSTEM
  degauss =   1.4699723600d-02
  ecutrho =   2.8000000000d+02
  ecutwfc =   50.00000
  ibrav =     4
  celldm(1) = 3.122
  celldm(3) = 3.839
  nat = 6
  nosym = .false.
  ntyp = 2
  occupations = 'smearing'
  smearing = 'cold'
/
&ELECTRONS
  conv_thr =   1.2000000000d-09
  electron_maxstep = 80
  mixing_beta =   4.0000000000d-01
/
ATOMIC_SPECIES
Mo     95.96 Mo_ONCV_PBE-1.0.oncvpsp.upf
S      32.065 s_pbe_v1.4.uspp.F.UPF
ATOMIC_POSITIONS crystal
Mo           0.3333333300       0.6666666700       0.2500000000 
Mo           0.6666666700       0.3333333300       0.7500000000 
S            0.3333333300       0.6666666700       0.8551740000 
S            0.6666666700       0.3333333300       0.1448260000 
S            0.6666666600       0.3333333300       0.3551740000 
S            0.3333333400       0.6666666700       0.6448260000 
K_POINTS automatic
$k_points $k_points $k_points 0 0 0
EOF

   # If pw.x is not found, specify the correct value for $espresso_dir,
   # use $espresso_dir/bin/pw.x instead of pw.x

   pw.x -in $input > $raw_output

   # grep -e 'lattice parameter' -e ! $raw_output | \
      # awk '/lattice/{alat=$(NF-1)}/!/{print alat, $(NF-1)}' >> $processed_output
   # set each grep to variable
   cutoff="$(grep -e 'kinetic-energy cutoff' ${raw_output} | awk '{print $(NF-1)}')"
   total_energy="$(grep ! ${raw_output} | awk '{print $(NF-1)}')"
   wall_time="$(grep 'PWSCF.*WALL' ${raw_output} | awk '{printf $(NF-1)}' | head -c -1)" 

   # write energy cutoff, total energy, wall time and kpoints to csv file
   echo "${cutoff}, ${total_energy}, ${wall_time}, ${k_points}" >> $processed_output
done
