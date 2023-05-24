import sys
import pandas as pd

# Usage
# python3 get_dimer_residue_group.py $TEMPLATE_PROJECT_PREFIX/0structure/output/atominfo.dat dimer

def get_dictionaries(df):
    # return {name: "residue_name", atoms: [2,3,4,5]}
    result = []
    for residue_id in df.residue_id.unique():
        residue = df[df["residue_id"] == residue_id]
        result.append({
            "name": residue.iloc[0]["residue_name"],
            "id": residue_id,
            "atoms": list(residue["atom_id"])
        })
    return result

def create_atomgroup(atominfo_file, output_file, except_chain_ids):
    atominfo = pd.read_table(atominfo_file, delim_whitespace=True, header=0, names=["atom_id", "atom_name", "residue_id", "residue_name", "chain_id"], usecols=[0,1,2,3,4])
    # Remove chains not used
    for chain_id in except_chain_ids:
        atominfo = atominfo[atominfo["chain_id"] != chain_id]
    
    result_strings = []
    
    for chain_id in atominfo.chain_id.unique():
        chain = atominfo[atominfo.chain_id == chain_id]
        even = chain[chain["residue_id"] % 2 == 0]
        odd = chain[chain["residue_id"] % 2 == 1]

        if len(even.residue_id.unique()) == len(odd.residue_id.unique()):
            longer = get_dictionaries(odd)
            shorter = get_dictionaries(even)
            last_of_longer = None
        elif len(even.residue_id.unique()) >= len(odd.residue_id.unique()):
            longer = get_dictionaries(even)
            shorter = get_dictionaries(odd)
            last_of_longer = longer.pop(-1)
        elif len(even.residue_id.unique()) <= len(odd.residue_id.unique()):
            longer = get_dictionaries(odd)
            shorter = get_dictionaries(even)
            last_of_longer = longer.pop(-1)
        
        # Start from begin
        result = ""
        for a, b in zip(longer, shorter):
            if a["id"] < b["id"]:
                left = a
                right = b
            else:
                left = b
                right = a
            group_name = f"{left['id']:0>5}_{left['name']}-{right['id']:0>5}_{right['name']}"
            atom_list = list(map(str, left['atoms'] + right['atoms']))
            result += f"[{group_name}]\n{' '.join(atom_list)}\n"
        
        result_strings.append(result)
        
        # Start from begin + 1
        if last_of_longer is not None:
            longer.append(last_of_longer)
            longer.pop(0)
        else:
            shorter.pop(0)
            longer.pop(-1)
        result = ""
        for a, b in zip(longer, shorter):
            if a["id"] < b["id"]:
                left = a
                right = b
            else:
                left = b
                right = a
            group_name = f"{left['id']:0>5}_{left['name']}-{right['id']:0>5}_{right['name']}"
            atom_list = list(map(str, left['atoms'] + right['atoms']))
            result += f"[{group_name}]\n{' '.join(atom_list)}\n"
        
        result_strings.append(result)

    with open(output_file + "_odd.dat", mode="w") as f:
        f.write(result_strings[0] + result_strings[3])

    with open(output_file + "_even.dat", mode="w") as f:
        f.write(result_strings[1] + result_strings[2])
        

if __name__ == '__main__':
    except_chain_ids = [1,2]
    try:
        atominfo_file = sys.argv[1]
        output_file = sys.argv[2]
    except Exception:
        print("Two arguments required.")
        exit()
    
    create_atomgroup(atominfo_file, output_file, except_chain_ids)
