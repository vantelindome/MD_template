# TEMPLATE_SYSTEM_NAME

<!-- Write details of your system -->

## Requirements

- `TEMPLATE_PROJECT_PREFIX`
- conda environment

### Environment variable to detect project directory path

`TEMPLATE_PROJECT_PREFIX` to detect project directory.

If this repository located in `/home/flat35hd99/ytva/` and the directory name is `dark` (See below tree), `TEMPLATE_PROJECT_PREFIX=/home/flat35hd99/ytva/dark`. Don't add slash ("/") in the last of the variable.

I recommend you to add `export TEMPLATE_PROJECT_PREFIX=/path/to/directory` into `.bashrc` or other script that is loaded in initialization of CLI.

```
ytva/
├── dark <- This repository
│   ├── 0structure
│   ├── 1minimization
│   ├── 2equilibrium
│   ├── 3sampling
│   ├── 4nve
│   ├── 5heat_flux
│   ├── 6heat_conductivity
│   ├── FMN
│   └── README.md
└── light
```

### conda environment

```bash
# Install miniconda in $HOME/miniconda3_x86_64

# If you want to install other place,
# please change util/load_curp

# conda environment name to intall curp
env_name=py2curp

# If possible, please use newer versions without python
conda create -n $env_name python=2.7 netcdf4=1.4.2 pygraphviz=1.3
conda activate $env_name
pip install curp==1.3.1
```

## Process

1. 0structure
   1. `create_amber_input.sh`
   2. `create_atomgroup_sidechain.sh`
2. 1minimization
   1. `pjsub job.sh`
3. 2equilibrium
   1. `submit.sh 0 4`
4. 3sampling
   1. `submit.sh 0 4`
5. 4nve
   1. `submit.sh 0 49`
   2. `submit_utility.sh 0 49` on IMS
6. 5curp
   1. `./generate_representation_of_grouppair.sh` on IMS.
   2. `submit.sh 0 49 residue` and `submit.sh 0 49 side` on IMS
