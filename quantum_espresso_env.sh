#!/bin/bash

module load slurm git gcc/10.3.0 intel-mkl/2020.4.304 openmpi/4.0.7 python/3.9 cmake
export ESPRESSO_PSEUDO=/mnt/home/landerson1/apps/q-e/pseudo
export ASE_ESPRESSO_COMMAND="mpirun -np 20 /mnt/home/landerson1/apps/q-e/bin/pw.x -in PREFIX.pwi > PREFIX.pwo" 
export PATH=~/apps/q-e/build/bin/:$PATH

export PYTHONPATH=~/apps/ase:$PYTHONPATH
export PATH=~/apps/ase/ase/build/cli:$PATH
