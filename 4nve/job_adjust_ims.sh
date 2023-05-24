#!/bin/sh
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1:jobtype=core
#PBS -l walltime=167:00:00

prefix=$TEMPLATE_PROJECT_PREFIX

source $prefix/util/load_curp_ims

set -xeu

output=$prefix/4nve/output/${run}_${smp}
tmp=/work/users/$USER/$PBS_JOBID

mkdir -p $tmp
cd $tmp

echo "adjust job start time" `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

cp $output/md.vel.nc .

echo "copy end time " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

curp conv-trj -vel \
    -p $prefix/0structure/output/system.dry.prmtop \
    -pf amber \
    -i md.vel.nc \
    -if netcdf \
    --irange 1 -1 1 \
    -o $tmp/adjust.vel.nc \
    -of netcdf \
    --orange 1 -1 2 \
    adjust-vel > /dev/null

echo "adjust end time" `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

cp $tmp/adjust.vel.nc $output
if [ -e $output/adjust.vel.nc ]; then rm $output/md.vel.nc; fi
rm -rf $tmp

echo "adjust job end time" `date +'%Y%m%d %H:%M:%S'` >> $output/time.log
