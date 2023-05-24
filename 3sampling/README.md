# Sampling structures

After Heating, we obtain structures as initial structures of 4nve

## Require other steps

- `$prefix/0structure/output/system.crd`
- `$prefix/0structure/output/system.prmtop`
- `$prefix/2equilibrium/output/md.rst`

## Required from other steps

- `$prefix/3sampling/output/$run/${step}/md.rst`
  Used by 4nve step.

## Usage

Run `./submit.sh $run`. $run is integer.

## Execution time

25M + 10M + 50M = 85M steps
