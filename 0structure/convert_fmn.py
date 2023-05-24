import sys, re

pdb2mol2 = {
    "N1 ": "N1 ",
    "C2 ": "C1 ",
    "O2 ": "O1 ",
    "N3 ": "N2 ",
    "C4 ": "C2 ",
    "O4 ": "O2 ",
    "C4A": "C3 ",
    "N5 ": "N3 ",
    "C5A": "C4 ",
    "C6 ": "C5 ",
    "C7 ": "C6 ",
    "C7M": "C7 ",
    "C8 ": "C8 ",
    "C8M": "C9 ",
    "C9 ": "C10",
    "C9A": "C11",
    "N10": "N4 ",
    "C10": "C12",
    "C1'": "C13",
    "C2'": "C14",
    "O2'": "O3 ",
    "C3'": "C15",
    "O3'": "O4 ",
    "C4'": "C16",
    "O4'": "O5 ",
    "C5'": "C17",
    "O5'": "O6 ",
    "P  ": "P1 ",
    "O1P": "O7 ",
    "O2P": "O8 ",
    "O3P": "O9 ",
}

def cli(original, output):
    # Get FMN lines
    fmn_lines = []
    with open(original) as f:
        for l in f.readlines():
            if re.search("^HETATM.*FMN", l):
                for k, v in pdb2mol2.items():
                    if l[13:16] == k:
                        fmn_lines.append(l.replace(k,v))
                        break

    # Write FMN
    with open(output, "w") as f_output:
        f_output.writelines(fmn_lines)

if __name__ == "__main__":
    original_pdb = sys.argv[1]
    output = sys.argv[2]
    cli(original_pdb, output)

