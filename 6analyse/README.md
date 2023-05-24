# Analysis

## Require other steps

- Conductivity data
  - `$prefix/5curp/output/{inter_residue | inter_side | intra_dimer_even | intra_dimer_odd | intra_residue | intra_whole}/{data_id}/conductivity.dat`

## Discussion

### Error convergence

<img src="https://github.com/flat35hd99/ytva_dark/blob/output/graphs/residue_and_side/inter_side-sem.png" width="500" height="500" alt="Error bar" >

From this picture, error size looks small.(I have not known how to evaluate error). I think I can go next model(YtvA in light) and compare with it.

### Compare inter residue and inter side

### Are there any cooperativity between monomers ?

A question(発問) from Yamato.

<img src="https://github.com/flat35hd99/ytva_dark/blob/output/graphs/residue_and_side/inter_residue-comparison_between_monomers.png" width="500" height="500" alt="Comparison between each monomer" >

Other quesitons about this topic are

 - Where large differences between monomers ?
   - specific structure such as loop, J-alpha
   - conserved residue or not
   - dense/sparse(疎密) place in network

### Are there any interaction between monomers ? If so, where ?

We can see interactions of inter sidechains effect is larger than that of inter residue. In other words, the effects of backbone is huge.

<img src="https://github.com/flat35hd99/ytva_dark/raw/output/graphs/residue_and_side/inter_side-network_between_monomers.png" width="500" height="500" alt="Interactions of inter sidechain between monomers" >

<img src="https://github.com/flat35hd99/ytva_dark/raw/output/graphs/residue_and_side/inter_residue-network_between_monomers.png" width="500" height="500" alt="Interactions of inter residue between monomers" >

### How each conserved residue plays its role ?

- `Q503`
  - According to 
  - > A Conserved Glutamine Plays a Central Role in LOV Domain Signal Transmission
  - , Q503 have a role.
  - In `2PR5` structure,
    - `108` at original pdb file
    - `105` at `atominfo.dat`
    - `103` at network graph
  - It is a part of pathway inter-residue and inter-sidechain from FMN to J-alpha.
- Cystein that make covalent bond with FMN is not included in the pathway.

### Where is pathway ?

Here.

<img src="https://github.com/flat35hd99/ytva_dark/blob/output/graphs/residue_and_side/inter_side-network_pathway.png" width="500" height="500" alt="Signal pathway between sidechain" >

### Does pathway pass through conserved residues ? If so, where ?
