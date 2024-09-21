# example_usage.py

from chem_formater import (
    smiles_to_pdbqt,
    pdbqt_to_pdb,
    pdb_to_sdf,
    pdbqt_to_smiles,
    pdb_to_smiles,
    OUTPUT_DIR
)
import os

# Example SMILES string
smiles = "CC(=O)NC1=CC=C(C=C1)O"  # Ethanol

# Convert SMILES to PDBQT
pdbqt_filename = "ethanol.pdbqt"
pdbqt_path = smiles_to_pdbqt(smiles, pdbqt_filename)

# Convert PDBQT to PDB
pdb_filename = "ethanol.pdb"
pdb_path = pdbqt_to_pdb(pdbqt_filename, pdb_filename)

# Convert PDB to SDF
sdf_filename = "ethanol.sdf"
sdf_path = pdb_to_sdf(pdb_filename, sdf_filename)

# Convert PDBQT to SMILES
smiles_from_pdbqt = pdbqt_to_smiles(pdbqt_filename)
print("SMILES from PDBQT:", smiles_from_pdbqt)

# Convert PDB to SMILES
smiles_from_pdb = pdb_to_smiles(pdb_filename)
print("SMILES from PDB:", smiles_from_pdb)

# Print the OUTPUT_DIR where files are saved
print("All files are saved in:", OUTPUT_DIR)
