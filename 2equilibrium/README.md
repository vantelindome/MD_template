# Heating system

After optimization, we heat our system.

## Require other steps

- $prefix/0structure/output/system.crd
- $prefix/0structure/output/system.prmtop
- $prefix/1minimization/output/optimize_whole_system/md.rst

## Required from other steps

- $prefix/2equilibrium/output/md.rst
  Used by 3sampling step.

## Usage

Run `./submit.sh $run`. $run is integer.

## Execution time

Each run has 50k + 350k + 100k = 500k takes 5mins.
