#!/bin/sh -xeu
#PJM -L rscgrp=cl-share
#PJM -L elapse=03:00:00
#PJM --mpi proc=1
#PJM -S
#PJM -j

# Get run and sampling number
# specify_job_id should be defined by submitter
run=${specify_job_id:0:2}
smp=${specify_job_id:2:1}

# Set appropriate environment
prefix=$TEMPLATE_PROJECT_PREFIX

module purge
source $prefix/util/load_curp

output_residue=$prefix/5curp/output/residue/${run}_${smp}
output_side=$prefix/5curp/output/side/${run}_${smp}
mkdir -p $output_residue
mkdir -p $output_side
cd $output_residue

echo "job start time" `date +'%Y%m%d %H:%M:%S'` >> $output_residue/time.log

curp analyze pickup-respairs \
  --input-prmtop-file $prefix/0structure/output/system.dry.prmtop \
  --input-prmtop-format amber \
  --interval 1000 \
  --cutoff 6.0 \
  --input-format netcdf \
  $prefix/4nve/output/${run}_${smp}/center.crd.nc \
  > $output_residue/group_pair_original.dat

echo "pickup res end" `date +'%Y%m%d %H:%M:%S'` >> $output_residue/time.log

python3 $prefix/5curp/generate_grouppair_sidechain.py \
  $output_residue/group_pair_original.dat 1 2 \
  > $output_side/group_pair_original.dat

echo "job end time  " `date +'%Y%m%d %H:%M:%S'` >> $output_residue/time.log
