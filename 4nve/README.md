# Sampling trajectories

Finally, we can start production run to get trajectories for CURP.

## Require other steps

- `$prefix/0structure/output/system.crd`
- `$prefix/0structure/output/system.prmtop`
- `$prefix/3sampling/output/$run/sampling/md.rst`

## Required from other steps

- `$prefix/4nve/output/${run}_${step}/md.crd.nc`
- `$prefix/4nve/output/${run}_${step}/md.vel.nc`
  Used by 5hflux step.

## Usage

Run `./submit.sh $bulk_start_number $bulk_end_number`. `$bulk_start_number` and `$bulk_end_number` are integer. See below corresponding bulk number and run and sample numbers.

|bulk number|run|sample|
|--|--|--|
|0 |0 |0 |
|5 |0 |5 |
|10|1 |0 |
|23|2 |3 |

## Execution time

1M steps take 1 hour and 10 mins to copy data from `/local` disk of Type 2 node to `/data` disk.
