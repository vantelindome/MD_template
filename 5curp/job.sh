#!/bin/sh
#PJM -L rscunit=lm
#PJM -L rscgrp=lm-middle
#PJM -L socket=1
#PJM --mpi proc=28
#PJM -L elapse=48:00:00
#PJM -S
#PJM -j

prefix=$TEMPLATE_PROJECT_PREFIX

module purge
source $prefix/util/load_curp

set -xeu

specify_job_id=$(printf %03d $PJM_BULKNUM)
run=${specify_job_id:0:2}
smp=${specify_job_id:2:1}

output=$prefix/5curp/output/$group_pair/${run}_${smp}
# faster directory
tmp=/dev/shm/$PJM_JOBID/$PJM_BULKNUM
# tmp=$prefix/5curp/output/vd/$PJM_JOBID/$PJM_BULKNUM

mkdir -p $output
mkdir -p $tmp
cd $tmp

echo "[${group_pair}] job start time" `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

# Copy atom group file if need
if [ $group_pair == "inter_side" ];then
    cp $prefix/5curp/input/atomgroup/side.dat atom_group.dat
elif [ $group_pair == "intra_dimer_even" ];then
    cp $prefix/5curp/input/atomgroup/dimer_even.dat atom_group.dat
elif [ $group_pair == "intra_dimer_odd" ];then
    cp $prefix/5curp/input/atomgroup/dimer_odd.dat atom_group.dat
elif [ $group_pair == "intra_whole" ];then
    cp $prefix/5curp/input/atomgroup/whole.dat atom_group.dat
fi

cp $prefix/0structure/output/system.dry.prmtop .
# ln -s $prefix/4nve/output/${run}_${smp}/center.crd.nc md.crd.nc
# ln -s $prefix/4nve/output/${run}_${smp}/adjust.vel.nc md.vel.nc
cp $prefix/4nve/output/${run}_${smp}/center.crd.nc md.crd.nc
cp $prefix/4nve/output/${run}_${smp}/adjust.vel.nc md.vel.nc
cp $prefix/5curp/input/config/generated_${group_pair}.cfg flux.cfg
cp $prefix/5curp/output/group_pair/${group_pair}.dat group_pair.dat
cp $prefix/5curp/output/group_pair/${group_pair}.dat $output/group_pair_used.dat

echo "[${group_pair}] copy end time " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

OMP_NUM_THREADS=1
mpiexec -np $PJM_MPI_PROC curp compute flux.cfg > $output/flux.log

echo "[${group_pair}] flux end time " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

cp flux_grp.nc $output/flux.nc

echo "[${group_pair}] copy end time " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

OMP_NUM_THREADS=4 # MPI process * OpenMP threads = Core number
mpiexec -np 7 \
  curp cal-tc \
    --frame-range 1 50000 1 --average-shift 2 \
    --dt 0.001 \
    -a acf.nc \
    -o $output/conductivity.dat \
    flux_grp.nc > $output/conductivity.log

echo "[${group_pair}] cond end time " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

cp acf.nc $output/

echo "[${group_pair}] copy end time " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

rm $tmp -rf

echo "[${group_pair}] rm tmp directory end time " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log
