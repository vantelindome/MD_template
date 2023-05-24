# Calcualte Heat flux and conductivity

Using curp, we calculate heat flux and its conductivity

## Require other steps

- `$prefix/0structure/output/system.prmtop`
- `$prefix/4nve/output/${run}_${sample}/md.nc`

## Required from other steps

- `$prefix/5heat_flux/output/${run}_${sample}/flux.nc`
  Used by 6heat_conductivity step.

## Usage

Run `./submit.sh $run $sample`. $run and $sample are integer.

## Process

### Create whole group file

I wrote `input/atomgroup/whole.dat` by my hand.

```dat
[00001_whole]
1-4144
```

### Create dimer group file

```shell
cd $TEMPLATE_PROJECT_PREFIX/5curp/input/atomgroup
python3 ../../utilize_script/get_dimer_residue_group.py $TEMPLATE_PROJECT_PREFIX/0structure/output/atominfo.dat dimer
```

Note: `utilize_script/get_dimer_residue_group.py` is specific script instead universal. It has many magic number.(assuming FMN and the positions)

### Create side chain group file

I did it at `0structure` and copy it `input/atomgroup/side.dat`

### Create intra group pair file

```shell
cd $TEMPLATE_PROJECT_PREFIX/5curp
python3 utilize_script/get_self_pair.py input/atomgroup/residue.dat output/group_pair/intra_residue.dat
python3 utilize_script/get_self_pair.py input/atomgroup/whole.dat output/group_pair/intra_whole.dat
python3 utilize_script/get_self_pair.py input/atomgroup/dimer_even.dat output/group_pair/intra_dimer_even.dat
python3 utilize_script/get_self_pair.py input/atomgroup/dimer_odd.dat output/group_pair/intra_dimer_odd.dat
```

## LOG

2022/02/26 15:24

lm-middleに30-34のinter-residue投入。このとき0-2のinter-sideがまだ走っていたので、10-5-3 = 2個走らせられるため、 30-31のinter-sideを投入。
