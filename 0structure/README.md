# Create amber input files

In this directory, I created amber input files(system.crd and system.prmtop) from Gaussian output and 2PR5 pdb data.

You can reproduce my works following process:

1. Use H++ server
   1. Access H++ server.
   2. Input pdb code, 2PR5.
   3. Change pH to `7.0` and click PROCESS...
   4. Wait a few moments
   5. Click VIEW RESULTS.
   6. Download `0.15_80_10_pH7.0_2pr5.pqr.from_ambpdb` as `input/protein.pdb`
2. Run `create_amber_input.sh`
3. Run `parmed -p output/system.dry.prmtop -c output/system.dry.crd -i input/parmed.input --no-splash > log/summary_of_system_dry.log` to check the mass of system + ligand
4. Check the number of Ions using [SLTCAP](https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html)
   1. Protein mass (kDa) is `29.9633760`
   2. Solution salt concentration (mM/l) is `150`
   3. Net charge of solutes (proton charge units) is `-21`
   4. Number of water molecules is `12365`. This value is from `Added 12373 residues.` on `log/leap_stdout.log`
   5. I left other blanks as blanks
   6. The result show `Your system requires 21.43 anions and 42.43 cations.`
5. Change ion number
```diff
addIons2 system Na+ 0
- addIonsRand system Na+ 20 Cl- 20
+ addIonsRand system Na+ 21 Cl- 21
```
6. Run `create_amber_input.sh` again
7. Run `md5sum output/* > checksum` in `0structure` directory
8. To calculate side chain heat conductivity, we want atomgroup file. Run `create_atomgroup_sidechain.sh`
9. Re-run `md5sum output/* > checksum` in `0structure` directory

`parmed` step output

```
Amino Acid Residues:   253
Nucleic Acid Residues: 0
Number of cations:     0
Number of anions:      0
Num. of solvent mols:  0
Num. of unknown res:   2
Total charge (e-):     -21.0000
Total mass (amu):      29963.3760
Number of atoms:       4148
Number of residues:    255
Residue set:           ALA, ARG, ASN, ASP, CYS, FMN, GLN
                       GLU, GLY, HID, HIE, HIP, ILE, LEU
                       LYS, MET, PHE, PRO, SER, THR, TRP
                       TYR, VAL
Residue count:         ALA: 8, ARG: 6, ASN: 18, ASP: 20, CYS: 2, FMN: 2, GLN: 16
                       GLU: 26, GLY: 14, HID: 2, HIE: 1, HIP: 1, ILE: 20, LEU: 18
                       LYS: 18, MET: 6, PHE: 8, PRO: 10, SER: 3, THR: 22, TRP: 2
                       TYR: 10, VAL: 22
```

## Directory and files

- `create_amber_input.sh` is the script to do all you have to do in this step. 
- `log/` hold logs that record proceess in `create_amber_input.sh`
- `input/` hold input files of `create_amber_input.sh`
- `output/` hold output files of `create_amber_input.sh`.

## Requirements

- ESP charge information, calculated by Gaussian(FMN.esp)
- PDB file
  - Modify original 2PR5 using H++ server
    - pH = 7.0 

## Outputs

### antechamber outputs

- FMN.mol2
- FMN.frcmod

### tleap outputs

- system.dry.pdb
- system.pdb

The former is not solvated structure. The latter is solvated structure.

- system.crd
- system.prmtop

These files are created by saveAmberParam in tleap. `system.crd` is initial *input coordinates* of the system. `system.prmtop` records system *topology and parameters*. These files are also known as `inpcrd` and `prmtop`, respectively.

You can check details how these files are created at `input/leapin` file.

- leap.log

This is another log of tleap. I think readability of this file is not better than `log/leap_stdout.log`.

## Next step

1minimization is the next step. The step use

- `$prefix/0structure/output/system.prmtop`
- `$prefix/0structure/output/system.crd`

