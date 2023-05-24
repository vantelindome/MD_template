import matplotlib.pyplot as plt

# Get structure name from residue number of chain B
number2structure_at_chain_B = (
    ["unfold"] * (10 - 6 + 1)
    + ["beta"] * (15 - 11 + 1)
    + ["loop"] * (23 - 16 + 1)
    + ["beta"] * (27 - 24 + 1)
    + ["alpha"] * (36 - 28 + 1)
    + ["unfold"] * (37 - 37 + 1)
    + ["alpha"] * (42 - 38 + 1)
    + ["unfold"] * (45 - 43 + 1)
    + ["alpha"] * (51 - 46 + 1)
    + ["unfold"] * (55 - 52 + 1)
    + ["alpha"] * (70 - 56 + 1)
    + ["unfold"] * (72 - 71 + 1)
    + ["beta"] * (80 - 73 + 1)
    + ["loop"] * (85 - 81 + 1)
    + ["beta"] * (98 - 86 + 1)
    + ["loop"] * (100 - 99 + 1)
    + ["beta"] * (110 - 101 + 1)
    + ["alpha"] * (131 - 111 + 1)
)

# Get structure name from residue number of chain A
number2structure_at_chain_A = number2structure_at_chain_B + ["unfold"] * (132 - 132 + 1)

number2structure = (
    ["FMN"] * 2 + number2structure_at_chain_A + number2structure_at_chain_B
)

number2chain = ["FMN"] * (2 - 1 + 1) + ["A"] * (129 - 3 + 1) + ["B"] * (255 - 130 + 1)

chain2color = {"FMN": "#FF6699", "A": "#99FF66", "B": "#6699FF"}

structure2color = {
    "unfold": "white",
    "loop": "white",
    "beta": "#4099D4",
    "alpha": "#FC9A9E",
    "FMN": "white",
}

# + FMN x2
# + chain A x2
# - one extra residue of chain B
_chainB = (
    ["unfold"] * (10 - 6 + 1)
    + ["A_beta"] * (15 - 11 + 1)
    + ["loop"] * (23 - 16 + 1)
    + ["B_beta"] * (27 - 24 + 1)
    + ["C_alpha"] * (36 - 28 + 1)
    + ["unfold"] * (37 - 37 + 1)
    + ["D_alpha"] * (42 - 38 + 1)
    + ["unfold"] * (45 - 43 + 1)
    + ["E_alpha"] * (51 - 46 + 1)
    + ["unfold"] * (55 - 52 + 1)
    + ["F_alpha"] * (70 - 56 + 1)
    + ["unfold"] * (72 - 71 + 1)
    + ["G_beta"] * (80 - 73 + 1)
    + ["loop"] * (85 - 81 + 1)
    + ["H_beta"] * (98 - 86 + 1)
    + ["loop"] * (100 - 99 + 1)
    + ["I_beta"] * (110 - 101 + 1)
    + ["J_alpha"] * (131 - 111 + 1)
)
number2named_structure = ["FMN"] * 2 + _chainB + ["unfold"] + _chainB

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
named_structure2color = {
    "unfold": "#ddd",
    "loop": "#ddd",
    "FMN": "#ddd",
    "A_beta": colors[0],
    "B_beta": colors[1],
    "C_alpha": colors[2],
    "D_alpha": colors[3],
    "E_alpha": colors[4],
    "F_alpha": colors[5],
    "G_beta": colors[6],
    "H_beta": colors[7],
    "I_beta": colors[8],
    "J_alpha": colors[9],
}

short2single = {
    "FMN": "FMN",
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "HIP": "H",  # http://ambermd.org/Questions/HIS.html
    "HIE": "H",
    "HID": "H",
}

conserved_residue_number = [
    28,  # F at C-alpha
    32,  # T at C-alpha
    34,  # Y at loop between C-alpha and D-alpha
    44,  # C at loop between D-alpha and E-alpha
    47,  # L at loop between E-alpha and F-alpha
    76,  # N at G-beta
    81,  # G at loop between G-beta and H-beta
    86,  # N at H-beta
    105,  # Q at I-beta
]

if __name__ == "__main__":
    print(len(number2named_structure))
