# chem_formater.py

import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

# Define the output directory in the current working directory
OUTPUT_DIR = os.path.join(os.getcwd(), "output_files")
os.makedirs(OUTPUT_DIR, exist_ok=True)

def smiles_to_pdbqt(smiles, output_filename):
    """
    Convert SMILES to PDBQT with minimum energy conformation.
    Output file is saved in OUTPUT_DIR with the specified output_filename.
    Requires Open Babel installed and accessible via command line.
    """
    # Generate RDKit molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates
    embed_status = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if embed_status != 0:
        raise ValueError("Embedding molecule failed.")

    # Optimize geometry
    optimize_status = AllChem.UFFOptimizeMolecule(mol)
    if optimize_status != 0:
        raise ValueError("Geometry optimization failed.")

    # Write to temporary PDB file in OUTPUT_DIR
    pdb_temp_path = os.path.join(OUTPUT_DIR, "temp.pdb")
    Chem.MolToPDBFile(mol, pdb_temp_path)

    # Output PDBQT path in OUTPUT_DIR
    pdbqt_output_path = os.path.join(OUTPUT_DIR, output_filename)

    # Convert PDB to PDBQT using Open Babel
    obabel_cmd = ["obabel", pdb_temp_path, "-O", pdbqt_output_path, "--partialcharge", "gasteiger"]
    result = subprocess.run(obabel_cmd, capture_output=True, text=True)

    # Remove temporary PDB file
    os.remove(pdb_temp_path)

    if result.returncode != 0:
        raise RuntimeError(f"Open Babel conversion failed: {result.stderr}")

    print(f"PDBQT file saved to {pdbqt_output_path}")
    return pdbqt_output_path

def pdbqt_to_pdb(input_filename, output_filename):
    """
    Convert PDBQT to PDB using Open Babel.
    Input and output files are in OUTPUT_DIR.
    """
    pdbqt_input_path = os.path.join(OUTPUT_DIR, input_filename)
    pdb_output_path = os.path.join(OUTPUT_DIR, output_filename)

    obabel_cmd = ["obabel", pdbqt_input_path, "-O", pdb_output_path]
    result = subprocess.run(obabel_cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"Open Babel conversion failed: {result.stderr}")

    print(f"PDB file saved to {pdb_output_path}")
    return pdb_output_path

def pdb_to_sdf(input_filename, output_filename):
    """
    Convert PDB to SDF using RDKit.
    Input and output files are in OUTPUT_DIR.
    """
    pdb_input_path = os.path.join(OUTPUT_DIR, input_filename)
    sdf_output_path = os.path.join(OUTPUT_DIR, output_filename)

    mol = Chem.MolFromPDBFile(pdb_input_path, removeHs=False)
    if mol is None:
        raise ValueError("Failed to read PDB file.")

    w = Chem.SDWriter(sdf_output_path)
    if w is None:
        raise IOError("Failed to create SDF writer.")

    w.write(mol)
    w.close()
    print(f"SDF file saved to {sdf_output_path}")
    return sdf_output_path

def pdbqt_to_smiles(input_filename):
    """
    Convert PDBQT to SMILES using Open Babel.
    Input file is in OUTPUT_DIR.
    """
    pdbqt_input_path = os.path.join(OUTPUT_DIR, input_filename)

    obabel_cmd = ["obabel", pdbqt_input_path, "-osmi"]
    result = subprocess.run(obabel_cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"Open Babel conversion failed: {result.stderr}")

    smiles = result.stdout.strip()
    if not smiles:
        raise ValueError("Failed to convert PDBQT to SMILES.")

    return smiles

def pdb_to_smiles(input_filename):
    """
    Convert PDB to SMILES using RDKit.
    Input file is in OUTPUT_DIR.
    """
    pdb_input_path = os.path.join(OUTPUT_DIR, input_filename)

    mol = Chem.MolFromPDBFile(pdb_input_path, removeHs=False)
    if mol is None:
        raise ValueError("Failed to read PDB file.")

    smiles = Chem.MolToSmiles(mol)
    if not smiles:
        raise ValueError("Failed to convert PDB to SMILES.")

    return smiles
