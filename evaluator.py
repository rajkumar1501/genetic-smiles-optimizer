# evaluator.py

from rdkit import Chem
from posebusters import PoseBusters
import pandas as pd
from pathlib import Path
import os

def evaluate_molecules_with_posebusters(sdf_file_path):
    """
    Evaluates molecules in an SDF file using PoseBusters in 'mol' mode.

    Args:
        sdf_file_path (str): Path to the SDF file containing the molecules.

    Returns:
        pd.DataFrame: DataFrame containing SMILES strings and overall PoseBusters result.
    """
    # Load molecules from SDF file using RDKit
    suppl = Chem.SDMolSupplier(sdf_file_path, removeHs=False)
    molecules = [mol for mol in suppl if mol is not None]

    if not molecules:
        print("No valid molecules found in the SDF file.")
        return pd.DataFrame(columns=['SMILES', 'All_PoseBusters_Parameters_Passed'])

    # Prepare PoseBusters in 'mol' mode
    buster = PoseBusters(config="mol")

    # Temporary directory to save individual molecule files
    temp_dir = Path("posebusters_temp")
    temp_dir.mkdir(exist_ok=True)

    results = []

    for idx, mol in enumerate(molecules):
        # Generate SMILES string
        smiles = Chem.MolToSmiles(mol)

        # Save molecule to temporary SDF file for PoseBusters
        mol_filename = temp_dir / f"molecule_{idx+1}.sdf"
        writer = Chem.SDWriter(str(mol_filename))
        writer.write(mol)
        writer.close()

        # Run PoseBusters on the molecule
        try:
            df = buster.bust([mol_filename], None, None, full_report=False)
            # Check if all PoseBusters parameters are True
            # Exclude the 'mol_pred_loaded' column from the check
            posebusters_checks = df.drop(columns=['mol_pred_loaded'])
            all_passed = posebusters_checks.all(axis=1).iloc[0]
        except Exception as e:
            print(f"An error occurred while processing molecule {idx+1}: {e}")
            all_passed = False

        # Append results
        results.append({
            'SMILES': smiles,
            'All_PoseBusters_Parameters_Passed': all_passed
        })

    # Clean up temporary files
    for file in temp_dir.glob("*"):
        file.unlink()
    temp_dir.rmdir()

    # Create DataFrame
    result_df = pd.DataFrame(results)

    return result_df
