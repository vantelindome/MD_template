#!/bin/bash

prefix=$TEMPLATE_PROJECT_PREFIX
run_start=$1
run_end=$2

for i in $( seq $run_start $run_end ); do
    specify_job_id=$(printf %03d $i)
    run_number=${specify_job_id:0:2}
    sampling_number=${specify_job_id:2:1}
    job_name="ytva_util_${run_number}_${sampling_number}"

    echo "${job_name} will be submitted"

    if [[ "$(hostname)" =~ flow ]];then
        pjsub -N $job_name \
            --step \
            -x specify_job_id=$specify_job_id \
            $prefix/4nve/job_center.sh \
            $prefix/5curp/job_pickup.sh
        
        pjsub -N $job_name \
            -x specify_job_id=$specify_job_id \
            $prefix/4nve/job_adjust.sh
    else
        variables="run=${run_number},smp=${sampling_number}"
        
        # Process coordinate and velocity severally
        jsub --step -N $job_name -v $variables \
            $prefix/4nve/job_center_ims.sh \
            $prefix/5curp/job_pickup_ims.sh
        jsub --step -N $job_name -v $variables \
            $prefix/4nve/job_adjust_ims.sh
    fi
done
