#!/bin/bash

module purge

if [[ "$(hostname)" =~ flow ]];then
    source $TEMPLATE_PROJECT_PREFIX/util/load_pinter
else
    source $TEMPLATE_PROJECT_PREFIX/util/load_pinter_ims
fi

set -eu

cd $TEMPLATE_PROJECT_PREFIX/5curp

if [ -e group_pair/ ]; then mkdir group_pair/; fi

for group_type in "residue" "side"; do
    pinter -i $(find output/$group_type -type f -name group_pair_original.dat) -o output/group_pair/inter_${group_type}.dat
done
