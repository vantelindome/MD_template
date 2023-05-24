
import sys, re

# Usage
# python3 get_self_pair.py atomgroup.dat selfpair.dat

def main(atomgroup_file, output_file):
    is_comment = re.compile('^\s*#')
    is_section = re.compile('^\[.*\]')
    groups = []
    with open(atomgroup_file, mode="r") as f:
        for l in f.readlines():
            if is_comment.match(l): continue
            if is_section.match(l): groups.append(is_section.search(l).group())
    
    result = [f"{g}\n{g[1:-1]}" for g in groups]
    with open(output_file, mode="w") as f:
        f.write("\n".join(result))


if __name__ == '__main__':
    atomgroup_file = sys.argv[1]
    output_file = sys.argv[2]
    main(atomgroup_file, output_file)
