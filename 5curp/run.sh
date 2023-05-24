#!/bin/bash
prefix=$TEMPLATE_PROJECT_PREFIX

module purge
source $prefix/util/load_curp

set -xeu

run=$1
smp=$2
group_pair_type=$3
output=$prefix/5curp/output/$group_pair_type/${run}_${smp}
# faster directory
tmp=/dev/shm/$PJM_JOBID/$PJM_BULKNUM

mkdir -p $output
mkdir -p $tmp
cd $tmp

echo "[${group_pair_type}] job start time" `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

# Copy atom group file if need
if [ $group_pair_type == "inter_side" ];then
    cp $prefix/5curp/input/atomgroup/side.dat atom_group.dat
elif [ $group_pair_type == "intra_dimer_even" ];then
    cp $prefix/5curp/input/atomgroup/dimer_even.dat atom_group.dat
elif [ $group_pair_type == "intra_dimer_odd" ];then
    cp $prefix/5curp/input/atomgroup/dimer_odd.dat atom_group.dat
elif [ $group_pair_type == "intra_whole" ];then
    cp $prefix/5curp/input/atomgroup/whole.dat atom_group.dat
fi

cp $prefix/0structure/output/system.dry.prmtop .
cp $prefix/4nve/output/${run}_${smp}/center.crd.nc md.crd.nc
cp $prefix/4nve/output/${run}_${smp}/adjust.vel.nc md.vel.nc
cp $prefix/5curp/input/config/${group_pair_type}.cfg flux.cfg
cp $prefix/5curp/output/group_pair/${group_pair_type}.dat group_pair.dat
cp $prefix/5curp/output/group_pair/${group_pair_type}.dat $output/group_pair_used.dat

echo "[${group_pair_type}] copy end time " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

mpirun -np 28 curp compute flux.cfg > $output/flux.log

echo "[${group_pair_type}] flux end time " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

cp flux_grp.nc $output/flux.nc

echo "[${group_pair_type}] copy end time " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

OMP_NUM_THREADS=4 # MPI process * OpenMP threads = Core number
mpirun -np 7 \
  curp cal-tc \
    --frame-range 1 50000 1 --average-shift 2 \
    --dt 0.001 \
    -a acf.nc \
    -o $output/conductivity.dat \
    flux_grp.nc > $output/conductivity.log

echo "[${group_pair_type}] cond end time " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

cp acf.nc $output/

echo "[${group_pair_type}] copy end time " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log
