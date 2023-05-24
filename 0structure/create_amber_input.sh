#!/bin/sh -xeu
#PJM -L rscgrp=fx-debug
#PJM -L elapse=01:00:00
#PJM --mpi proc=1
#PJM -j

# Note: This script has no compatibility with other environment like YtvA(light)
# If you want to reuse this script, please check not only settings but also antechamber options.

# ========
# Settings
# ========

# ------
# Ligand
# ------
ligand_name=FMN
charge=0
multiplicity=1

# ---------------------
# Files and directories
# ---------------------
prefix=$TEMPLATE_PROJECT_PREFIX

# input file
gaussian_esp_file=$prefix/FMN/output/B3LYP_dp_5/${ligand_name}.esp

# output file
output_dir=$prefix/0structure/output
mol2_file=$output_dir/${ligand_name}.mol2
frcmod_file=$output_dir/${ligand_name}.frcmod

# Create output directory and log directory
if [ ! -e $output_dir ]; then mkdir -p $output_dir; fi
if [ ! -e $prefix/0structure/log ]; then mkdir -p $prefix/0structure/log; fi

# =======================================================
# Converting gaussian output into force modification file
# =======================================================
source $prefix/util/load_amber
cd $prefix/0structure

# -----------
# esp -> mol2
# -----------
antechamber -i $gaussian_esp_file \
           -fi gesp \
           -o $mol2_file \
           -fo mol2 \
           -c resp \
           -at gaff2 \
           -nc $charge \
           -m $multiplicity \
           -rn $ligand_name \
           -pf y &> log/convert_gesp_to_mol2.log

# --------------
# mol2 -> frcmod
# --------------
parmchk2 -i $mol2_file \
         -f mol2 \
         -o $frcmod_file \
         -s gaff2 &> log/convert_mol2_to_frcmod.log

# ==========================
# Change atom name of flavin
# ==========================
if [ ! -e input/2PR5.pdb ]; then
  wget "https://files.rcsb.org/download/2PR5.pdb" -P input/
fi

python3 convert_fmn.py input/2PR5.pdb output/FMN_taken_out_of_original.pdb

# ====
# Leap
# ====
# leapin require
#  - output/FMN.mol2
#  - output/FMN.frcmod
#  - input/protein.pdb
#  - output/FMN_taken_out_of_original.pdb
# tleap output
#  - output/system.dry.pdb
#  - output/system.pdb
#  - output/system.prmtop
#  - output/system.crd
#  - output/leap.log
tleap -f input/leapin &> log/leap_stdout.log
