#!/bin/bash

prefix=$TEMPLATE_PROJECT_PREFIX

run_start=$1
run_end=$2

if [ -z "$2" ];then echo "check arguments";exit 1; fi


if [[ "$(hostname)" =~ flow ]];then
    for group_pair in "inter_residue" "inter_side" "intra_dimer_even" "intra_dimer_odd" "intra_residue" "intra_whole" ; do
        pjsub -N ytva_curp_${group_pair}_${run_start}_${run_end} \
            -x group_pair=$group_pair \
            --bulk --sparam ${run_start}-${run_end} \
            $prefix/5curp/job.sh
    done
else
    for i in $( seq $run_start $run_end ); do
        specify_job_id=$(printf %03d $i)
        run_number=${specify_job_id:0:2}
        sampling_number=${specify_job_id:2:1}

        variables="run=${run_number},smp=${sampling_number},group_type=${group_type}"
        job_name="ytva_${run_number}_${sampling_number}"

        jsub --step \
            -N $job_name \
            -v $variables \
            $prefix/5curp/job_flux_ims.sh \
            $prefix/5curp/job_cond_ims.sh
    done
fi
