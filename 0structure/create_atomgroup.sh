#!/bin/sh -xeu
#PJM -L rscgrp=fx-debug
#PJM -L elapse=01:00:00
#PJM --mpi proc=1
#PJM -j

prefix=$TEMPLATE_PROJECT_PREFIX

source $prefix/util/load_amber
module load python/3.7.7

cd $prefix/0structure

cpptraj -i input/atominfo.cpptraj

python3 generate_atomgroup_sidechain.py \
    --sidechain-filename output/atomgroup_side.dat \
    --residue-filename output/atomgroup_residue.dat \
    --atominfo-filename output/atominfo.dat \
    --no-separated-residues 1 2
