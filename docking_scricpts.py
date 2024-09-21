# docking_script.py

import os
import subprocess
import shutil
from chem_formater import (
    smiles_to_pdbqt,
    OUTPUT_DIR  # Ensure OUTPUT_DIR is accessible
)

def run_autodock_gpu(smiles, receptor_fld_path, output_dir, nrun=5, top_n=5):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Convert SMILES to PDBQT
    ligand_pdbqt_filename = 'ligand.pdbqt'
    print("Converting SMILES to PDBQT...")
    ligand_pdbqt_path = os.path.join(output_dir, ligand_pdbqt_filename)
    smiles_to_pdbqt(smiles, ligand_pdbqt_filename)
    ligand_pdbqt_full_path = os.path.join(output_dir, ligand_pdbqt_filename)

    # Build the AutoDock-GPU command
    autodock_gpu_executable = '/mnt/d/Projects/genetic-smiles-optimizer/autodock_linux/adgpu'  # Update path if necessary
    cmd = [
        autodock_gpu_executable,
        '--ffile', receptor_fld_path,
        '--lfile', ligand_pdbqt_full_path,
        '--nrun', str(nrun),
        '--xmloutput', '1',  # Enable XML output for easier parsing
        '--npdb','1',
    ]

    print("Running AutoDock-GPU docking...")
    # Run AutoDock-GP
    print(cmd)
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=output_dir)

    # Check for errors
    if result.returncode != 0:
        print("AutoDock-GPU failed to run:")
        print(result.stderr)
        return

    # Optional: Print AutoDock-GPU output
    print(result.stdout)

    # Parse the XML output to extract top N poses
    xml_output_file = os.path.join(output_dir, 'ligand.xml')
    if not os.path.isfile(xml_output_file):
        print("AutoDock-GPU did not generate the expected XML output.")
        return

    print("Parsing docking results to extract top poses...")
    top_poses = parse_autodock_gpu_xml(xml_output_file, top_n)

    # Save top poses as individual PDBQT files
    for idx, pose in enumerate(top_poses, 1):
        pose_filename = f'top_pose_{idx}.pdbqt'
        pose_filepath = os.path.join(output_dir, pose_filename)
        with open(pose_filepath, 'w') as f:
            f.write(pose)
        print(f"Saved pose {idx} to {pose_filepath}")

    print(f"Top {top_n} docking poses have been saved in {output_dir}")

def parse_autodock_gpu_xml(xml_file_path, top_n):
    """
    Parses the AutoDock-GPU XML output file to extract the top N docking poses.

    Args:
        xml_file_path (str): Path to the XML output file.
        top_n (int): Number of top poses to extract.

    Returns:
        List[str]: List of PDBQT strings for the top N poses.
    """
    import xml.etree.ElementTree as ET

    tree = ET.parse(xml_file_path)
    root = tree.getroot()

    # Namespace handling
    ns = {'ad': 'http://autosolvatedocking.org/schema'}

    # Find all docking poses
    poses = []
    for prediction in root.findall('.//ad:prediction', ns):
        pdbqt_lines = []
        for atom in prediction.findall('.//ad:atom', ns):
            pdbqt_line = atom.text
            pdbqt_lines.append(pdbqt_line)
        # Get the MODEL number and energy
        model_number = prediction.get('id')
        energy = float(prediction.find('.//ad:affinity', ns).text)
        # Assemble the PDBQT content
        pdbqt_content = f"MODEL {model_number}\n" + "\n".join(pdbqt_lines) + f"\nENDMDL\n"
        poses.append({'energy': energy, 'pdbqt': pdbqt_content})

    # Sort poses by energy (ascending)
    poses.sort(key=lambda x: x['energy'])

    # Extract top N poses
    top_poses = [pose['pdbqt'] for pose in poses[:top_n]]

    return top_poses

# Example usage
if __name__ == '__main__':
    # SMILES string of the ligand
    smiles = 'CC(=O)Oc1ccccc1C(=O)O'  # Aspirin, for example

    # Path to receptor grid maps (.fld file)
    receptor_fld_path = os.path.abspath('autodock_linux/Docking-4ieh-main/4ieh_protein.maps.fld')  # Update with your receptor .fld file path

    # Output directory (ensure it matches OUTPUT_DIR in chem_formater.py)
    output_dir = OUTPUT_DIR

    # Number of docking runs
    nrun = 5  # Increase for better sampling

    # Number of top poses to save
    top_n = 5

    run_autodock_gpu(smiles, receptor_fld_path, output_dir, nrun, top_n)
