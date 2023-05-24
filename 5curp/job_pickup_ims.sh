#!/bin/sh
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1:jobtype=core
#PBS -l walltime=01:00:00
prefix=$TEMPLATE_PROJECT_PREFIX

if [ "$PBS_O_WORKDIR" ]; then
  cd ${PBS_O_WORKDIR}
fi

module purge
source $prefix/util/load_curp_ims

set -xeu

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
  > $output_residue/group_pair.dat

echo "pickup res end" `date +'%Y%m%d %H:%M:%S'` >> $output_residue/time.log

python3 $prefix/5curp/generate_grouppair_sidechain.py \
  $output_residue/group_pair.dat 1 2 \
  > $output_side/group_pair.dat

echo "job end time  " `date +'%Y%m%d %H:%M:%S'` >> $output_residue/time.log
