#!/bin/sh -eu
#PJM -L rscgrp=cl-share
#PJM -L elapse=2:00:00
#PJM --mpi proc=1
#PJM -S

# Get run and sampling number
# specify_job_id should be defined by submitter
run=${specify_job_id:0:2}
smp=${specify_job_id:2:1}

# Set appropriate environment
prefix=$TEMPLATE_PROJECT_PREFIX

module purge
source $prefix/util/load_curp

set -xeu

output=$prefix/4nve/output/${run}_${smp}
# faster directory
tmp=/dev/shm/$PJM_JOBID

mkdir -p $output
mkdir -p $tmp
cd $tmp

echo "job start time" `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

cp $prefix/4nve/output/${run}_${smp}/md.vel.nc .

echo "copy end time " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

curp conv-trj -vel \
    -p $prefix/0structure/output/system.dry.prmtop \
    -pf amber \
    -i md.vel.nc \
    -if netcdf \
    --irange 1 -1 1 \
    -o $output/adjust.vel.nc \
    -of netcdf \
    --orange 2 -1 2 \
    adjust-vel > /dev/null

echo "adjust end time" `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

if [ -e $output/adjust.vel.nc ]; then rm $output/md.vel.nc; fi

echo "job end time  " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log
