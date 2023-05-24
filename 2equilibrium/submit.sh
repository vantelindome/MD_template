#!/bin/bash

run_start=$1
run_end=$2

pjsub -N yd_2_${run_start}_${run_end} \
  --bulk --sparam ${run_start}-${run_end} \
  $TEMPLATE_PROJECT_PREFIX/2equilibrium/job.sh
