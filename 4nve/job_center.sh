#!/bin/sh -eu
#PJM -L rscgrp=cl-share
#PJM -L elapse=03:00:00
#PJM --mpi proc=1
#PJM -S

# Get run and sampling number
# specify_job_id should be defined by submitter
run=${specify_job_id:0:2}
smp=${specify_job_id:2:1}

# Set appropriate environment
prefix=$TEMPLATE_PROJECT_PREFIX

module purge
source $prefix/util/load_amber_on_cloud

output=$prefix/4nve/output/${run}_${smp}

cd $output

# Job start
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
rm system.dry.prmtop \
   system.dry.crd

echo "remove crd end" `date +'%Y%m%d %H:%M:%S'` >> $output/time.log
