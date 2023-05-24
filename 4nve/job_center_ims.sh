#!/bin/sh
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1:jobtype=core
#PBS -l walltime=24:00:00
prefix=$TEMPLATE_PROJECT_PREFIX

if [ "$PBS_O_WORKDIR" ]; then
  cd ${PBS_O_WORKDIR}
fi

module purge
module load amber/18/bugfix16

set -xeu

output=$prefix/4nve/output/${run}_${smp}
cd $output # Because output dir is equal input, you don't need to create it.

echo "job start time" `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

ln -s $prefix/0structure/output/system.dry.prmtop system.dry.prmtop
ln -s $prefix/0structure/output/system.dry.crd system.dry.crd
cp $prefix/4nve/center.in .

echo "copy end time " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

cpptraj -i center.in

echo "cpptraj end   " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

md5sum *.nc > checksum

echo "checksum end  " `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

if [ -e $output/center.crd.nc ]; then rm $output/md.crd.nc; fi

echo "remove crd end" `date +'%Y%m%d %H:%M:%S'` >> $output/time.log
