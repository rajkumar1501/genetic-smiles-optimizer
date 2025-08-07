# docking_script.py

import os
import subprocess
from chem_formater import smiles_to_pdbqt

def extract_best_docking_score(autodock_output):
    """
    Extracts the best docking score from the AutoDock-GPU output.

    Args:
        autodock_output (str): The stdout output from AutoDock-GPU.

    Returns:
        float or None: The best docking score if found, else None.
    """
    import re
    best_score = None
    lines = autodock_output.strip().split('\n')
    for line in lines:
        line = line.strip()
        # Look for the line that contains "best energy"
        if 'best energy' in line:
            # Example line: "159 samples, best energy    -6.19 kcal/mol."
            match = re.search(r'best energy\s+([-\d\.]+)\s+kcal/mol', line)
            if match:
                best_score = float(match.group(1))
                break
    return best_score

def run_autodock_gpu(smiles, receptor_fld_path, output_dir, nrun=50, ligand_pdbqt_filename='ligand.pdbqt', best_pdbqt_filename='best.pdbqt'):
    """
    Runs AutoDock-GPU docking and returns the best docking score and the path to the best PDBQT file.

    Args:
        smiles (str): SMILES string of the ligand.
        receptor_fld_path (str): Path to the receptor grid map (.fld) file.
        output_dir (str): Directory to save output files.
        nrun (int): Number of docking runs.
        ligand_pdbqt_filename (str): Filename for the ligand PDBQT file.
        best_pdbqt_filename (str): Filename for the best pose PDBQT file.

    Returns:
        tuple: (best_pose_filepath, best_docking_score)
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Convert SMILES to PDBQT
    print("Converting SMILES to PDBQT...")
    ligand_pdbqt_full_path = os.path.join(output_dir, ligand_pdbqt_filename)
    smiles_to_pdbqt(smiles, ligand_pdbqt_full_path)

    # Verify that the ligand PDBQT file exists
    if not os.path.isfile(ligand_pdbqt_full_path):
        print(f"Ligand PDBQT file not found at {ligand_pdbqt_full_path}")
        return None, None

    # Build the AutoDock-GPU command
    autodock_gpu_executable = './autodock_linux/adgpu'  # Update path if necessary
    cmd = [
        autodock_gpu_executable,
        '--ffile', receptor_fld_path,
        '--lfile', ligand_pdbqt_full_path,
        '--nrun', str(nrun),
        '--gbest', '1',  # Output the best pose as a PDBQT file
    ]
    print(cmd)

    print("Running AutoDock-GPU docking...")
    # Run AutoDock-GPU
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check for errors
    if result.returncode != 0:
        print("AutoDock-GPU failed to run:")
        print(result.stderr)
        return None, None

    # Optional: Print AutoDock-GPU output
    print(result.stdout)

    # Extract the best docking score from the stdout
    best_docking_score = extract_best_docking_score(result.stdout)
    if best_docking_score is not None:
        print(f"Best docking score: {best_docking_score} kcal/mol")
    else:
        print("Failed to extract the best docking score.")

    # The best pose is saved as 'best.pdbqt' as specified by --gbest
    best_pose_filename = 'best.pdbqt'
    best_pose_filepath = os.path.join(os.getcwd(), best_pose_filename)
    renamed_best_pose_filepath = os.path.join(output_dir, best_pdbqt_filename)

    if os.path.exists(best_pose_filepath):
        os.rename(best_pose_filepath, renamed_best_pose_filepath)
        print(f"The best docking pose has been saved as {renamed_best_pose_filepath}")
    else:
        print("The best docking pose file was not found.")
        return None, None

    return renamed_best_pose_filepath, best_docking_score
