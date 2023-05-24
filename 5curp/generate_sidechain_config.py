import sys
from string import Template

# Usage
# python generate_sidechain_config.py \
#            input/config/template_inter_side.cfg \
#            input/atomgroup/side.dat \
#            input/config/generated_inter_side.cfg


def main(template, atomgroup, output):
    atoms = ""
    with open(atomgroup, mode="r") as f:
        for l in f.readlines():
            if l.startswith("["):
                continue
            atoms += l.replace("\n", " ")

    with open(template, mode="r") as f:
        config_template = Template(f.read())

    with open(output, mode="w") as f:
        f.write(config_template.substitute(target_atoms=atoms))


if __name__ == "__main__":
    template_file = sys.argv[1]
    atomgroup_file = sys.argv[2]
    output_file = sys.argv[3]
    main(template_file, atomgroup_file, output_file)

