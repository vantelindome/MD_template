#!/bin/sh -eu
#PJM -L rscgrp=cx-share
#PJM -L gpu=1
#PJM -L elapse=00:30:00
#PJM --mpi proc=1
#PJM -S

prefix=$TEMPLATE_PROJECT_PREFIX

module purge
source $prefix/util/load_amber_gpu

PMEMD=pmemd.cuda
output=$prefix/1minimization/output

mkdir -p $output

if [ -e $output/time.log ];then rm $output/time.log; fi
echo "job start time:" `date +'%Y%m%d %H:%M:%S'` >> $output/time.log

# Optimize system in multiple steps
coodinates=${prefix}/0structure/output/system.crd 
for step in optimize_hydrogen optimize_sidechain optimize_whole_system;do
  mkdir -p $output/$step
  cd $output/$step
  ${PMEMD} \
    -O \
    -i $prefix/1minimization/input/${step}.input \
    -p $prefix/0structure/output/system.prmtop \
    -c $coodinates \
    -ref $prefix/0structure/output/system.crd \
    -r $output/$step/md.rst \
    -o $output/$step/md.out \
    -inf $output/$step/md.info \
    -l $output/$step/md.log
  coodinates=$output/$step/md.rst # Update initial structure for the next step
done

echo "job end time  :" `date +'%Y%m%d %H:%M:%S'` >> $output/time.log
