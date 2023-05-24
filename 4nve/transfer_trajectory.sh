#/bin/bash

run_start=$1
run_end=$2

eval `ssh-agent`
ssh-add $HOME/.ssh/ims_rsa

target_output_prefix=ims:/home/users/$USER/flat/project/ytva_dark/4nve/output

for i in $( seq $run_start $run_end ); do
    specify_job_id=$(printf %03d $i)
    run_number=${specify_job_id:0:2}
    sampling_number=${specify_job_id:2:1}
    job_id=${run_number}_${sampling_number}
    scp output/$job_id/md.crd.nc \
        output/$job_id/md.vel.nc \
        $target_output_prefix/$job_id
    md5sum output/$job_id/*.nc > output/$job_id/checksum
done
