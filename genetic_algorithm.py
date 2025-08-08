# genetic_algorithm.py

import os
import torch
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit import RDLogger
import subprocess

# Import the updated docking function and evaluator
from docking_script import run_autodock_gpu
from evaluator import evaluate_molecules_with_posebusters

from mutate_smiles import load_model, generate_samples
from molecule_generation import load_model_from_directory
from pathlib import Path

# Disable RDKit warnings
RDLogger.DisableLog('rdApp.*')

# Paths to models and configurations
SMILES_GENERATION_MODEL_PATH = 'molgen_weights'  # Update as needed
MUTATION_MODEL_CONFIG_PATH = 'paper_checkpoints/ecfp4_with_counts_with_rank/config.yml'
MUTATION_MODEL_CHECKPOINT_PATH = 'paper_checkpoints/ecfp4_with_counts_with_rank/weights.ckpt'
MUTATION_MODEL_VOCAB_PATH = 'paper_checkpoints/vocabulary.pkl'

# Path to receptor grid maps (.fld file)
RECEPTOR_FLD_PATH = os.path.abspath('autodock_linux/Docking-4ieh-main/4ieh_protein.maps.fld')  # Update path

# Output directory
OUTPUT_DIR = 'output_files'
os.makedirs(OUTPUT_DIR, exist_ok=True)

def is_valid_molecule(smiles, allowed_elements={'H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'I'}):
    """
    Checks if a SMILES string can be converted to a valid 3D molecule and contains only allowed elements.
    Returns True if successful, False otherwise.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        elements_in_mol = {atom.GetSymbol() for atom in mol.GetAtoms()}
        if not elements_in_mol.issubset(allowed_elements):
            return False
        mol = Chem.AddHs(mol)
        params = Chem.AllChem.ETKDG()
        params.randomSeed = 0xf00d
        embed_status = Chem.AllChem.EmbedMolecule(mol, params)
        if embed_status != 0:
            return False
        optimize_status = Chem.AllChem.UFFOptimizeMolecule(mol)
        if optimize_status != 0:
            return False
        return True
    except Exception:
        return False

def generate_initial_population(desired_size, batch_size=100):
    """
    Generate an initial population of SMILES strings.
    Continues generating SMILES in batches until the desired number of valid molecules is obtained.
    """
    model_kwargs = {}  # Add any model-specific keyword arguments if needed
    valid_smiles = []
    total_generated = 0
    try:
        with load_model_from_directory(SMILES_GENERATION_MODEL_PATH, **model_kwargs) as model:
            while len(valid_smiles) < desired_size:
                samples = model.sample(batch_size)
                total_generated += batch_size
                # Filter out invalid SMILES
                for smiles in samples:
                    if is_valid_molecule(smiles):
                        valid_smiles.append(smiles)
                        if len(valid_smiles) >= desired_size:
                            break
        print(f"Generated {len(valid_smiles)} valid SMILES out of {total_generated} generated.")
        return valid_smiles[:desired_size]
    except Exception as e:
        print(f"Error generating initial SMILES: {e}")
        return []

def dock_and_evaluate(smiles_list, iteration, output_dir):
    """
    Dock SMILES strings, evaluate poses, and return results.
    """
    results = []
    for idx, smiles in enumerate(smiles_list):
        print(f"\n=== Generation {iteration}, Molecule {idx+1}/{len(smiles_list)} ===")
        try:
            # Define unique filenames for this molecule
            ligand_pdbqt_filename = f'ligand_{iteration}_{idx}.pdbqt'
            best_pdbqt_filename = f'best_{iteration}_{idx}.pdbqt'
            ligand_pdbqt_path = os.path.join(output_dir, ligand_pdbqt_filename)
            best_pdbqt_path = os.path.join(output_dir, best_pdbqt_filename)

            # Run AutoDock-GPU and get the best docking score and best PDBQT file path
            best_pdbqt_path, docking_score = run_autodock_gpu(
                smiles=smiles,
                receptor_fld_path=RECEPTOR_FLD_PATH,
                output_dir=output_dir,
                nrun=10,
                ligand_pdbqt_filename=ligand_pdbqt_filename,
                best_pdbqt_filename=best_pdbqt_filename
            )
            if not best_pdbqt_path or docking_score is None:
                print(f"Docking failed or docking score not found for molecule: {smiles}")
                continue

            # Convert best pose PDBQT to SDF for PoseBusters
            best_sdf_filename = f'best_{iteration}_{idx}.sdf'
            best_sdf_path = os.path.join(output_dir, best_sdf_filename)

            # Convert PDBQT to SDF using Open Babel
            obabel_cmd = [
                "obabel",
                "-ipdbqt", best_pdbqt_path,
                "-osdf",
                "-O", best_sdf_path
            ]
            subprocess.run(obabel_cmd, capture_output=True, text=True)

            # Evaluate with PoseBusters
            evaluation_df = evaluate_molecules_with_posebusters(best_sdf_path)
            passed_posebusters = evaluation_df['All_PoseBusters_Parameters_Passed'].iloc[0]

            # Store results
            results.append({
                'SMILES': smiles,
                'Docking_Score': docking_score,
                'Passed_PoseBusters': passed_posebusters,
                'PDBQT_File': best_pdbqt_filename
            })
        except Exception as e:
            print(f"Error processing molecule {smiles}: {e}")
            continue
    return pd.DataFrame(results)

def select_top_molecules(results_df, top_n):
    """
    Select molecules that passed PoseBusters evaluation and have top docking scores.
    """
    if results_df.empty:
        print("No results to select from.")
        return []

    filtered_df = results_df[results_df['Passed_PoseBusters'] == True]
    if filtered_df.empty:
        print("No molecules passed PoseBusters evaluation.")
        return []

    top_molecules = filtered_df.nsmallest(top_n, 'Docking_Score')
    selected_smiles = top_molecules['SMILES'].tolist()
    return selected_smiles
def perform_smiles_crossover(smiles1, smiles2):
    """
    Perform crossover between two SMILES strings at random points.

    Args:
        smiles1 (str): SMILES string of the first parent.
        smiles2 (str): SMILES string of the second parent.

    Returns:
        str or None: Crossover SMILES string, or None if crossover failed.
    """
    import random
    # Convert SMILES strings to lists of tokens
    tokens1 = list(smiles1)
    tokens2 = list(smiles2)

    len1 = len(tokens1)
    len2 = len(tokens2)

    if len1 < 2 or len2 < 2:
        return None

    # Identify possible crossover points avoiding special characters
    valid_indices1 = [i for i, c in enumerate(tokens1) if c.isalnum()]
    valid_indices2 = [i for i, c in enumerate(tokens2) if c.isalnum()]

    if not valid_indices1 or not valid_indices2:
        return None

    cross_point1 = random.choice(valid_indices1)
    cross_point2 = random.choice(valid_indices2)

    # Create new SMILES by combining parts of the parents
    new_smiles = ''.join(tokens1[:cross_point1] + tokens2[cross_point2:])

    return new_smiles


# New Function: Crossover Operation
def crossover_smiles(parent_smiles_list, num_offspring, model, device='cpu'):
    """
    Perform crossover between pairs of parent SMILES strings to produce offspring.
    Uses the mutation model to generate valid molecules from crossover products.

    Args:
        parent_smiles_list (list): List of parent SMILES strings.
        num_offspring (int): Number of offspring to generate.
        model: Mutation model to generate valid molecules.
        device (str): 'cpu' or 'cuda'

    Returns:
        list: List of valid offspring SMILES strings.
    """
    import random
    offspring_smiles = []
    num_parents = len(parent_smiles_list)
    if num_parents < 2:
        print("Not enough parents for crossover.")
        return offspring_smiles

    max_attempts = num_offspring * 10  # To avoid infinite loops
    attempts = 0

    while len(offspring_smiles) < num_offspring and attempts < max_attempts:
        attempts += 1
        # Randomly select two parents
        parent1_smiles = random.choice(parent_smiles_list)
        parent2_smiles = random.choice(parent_smiles_list)
        if parent1_smiles == parent2_smiles:
            continue  # Skip if both parents are the same

        # Perform crossover on SMILES strings
        crossover_smiles_candidate = perform_smiles_crossover(parent1_smiles, parent2_smiles)
        if crossover_smiles_candidate is None:
            continue  # Crossover failed, try again

        # Use the mutation model to generate valid molecules from the crossover SMILES
        try:
            mutated_smiles_list = generate_samples(model, crossover_smiles_candidate, beam_size=5, device=device)
            # Filter valid molecules
            for mutated_smiles in mutated_smiles_list:
                if is_valid_molecule(mutated_smiles):
                    offspring_smiles.append(mutated_smiles)
                    if len(offspring_smiles) >= num_offspring:
                        break
        except Exception as e:
            print(f"Error generating samples for crossover SMILES '{crossover_smiles_candidate}': {e}")
            continue

    return offspring_smiles[:num_offspring]


def mutate_smiles_list(smiles_list, model, mutations_per_parent=5, device='cpu'):
    """
    Generate new SMILES strings by mutating the given list.
    Generates 'mutations_per_parent' mutated SMILES per input SMILES.
    """
    new_smiles = []
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Cannot understand SMILES: {smiles}")
                continue
        except Exception as e:
            print(f"Error parsing SMILES '{smiles}': {e}")
            continue
        try:
            mutated_smiles = generate_samples(model, smiles, beam_size=mutations_per_parent, device=device)
            new_smiles.extend(mutated_smiles)
        except Exception as e:
            print(f"Error generating samples for SMILES '{smiles}': {e}")
            continue
    return new_smiles


class GeneticAlgorithm:
    def __init__(self, population_size, generations, mutations_per_parent=5, crossover_offspring_per_generation=10):
        self.population_size = population_size
        self.generations = generations
        self.mutations_per_parent = mutations_per_parent
        self.crossover_offspring_per_generation = crossover_offspring_per_generation
        self.current_population = []
        self.device = 'cuda' if torch.cuda.is_available() else 'cpu'
        self.mutation_model = None

    def initialize_population(self):
        print("\n=== Generating Initial Population ===")
        self.current_population = generate_initial_population(self.population_size)
        print(f"Initial Population ({len(self.current_population)} molecules): {self.current_population}")

    def evaluate_population(self, generation, generation_output_dir):
        # Evaluate current population
        results_df = dock_and_evaluate(self.current_population, generation, generation_output_dir)
        # Save results
        results_csv = os.path.join(generation_output_dir, f'docking_results_generation_{generation}.csv')
        results_df.to_csv(results_csv, index=False)
        print(f"Results saved to {results_csv}")
        return results_df

    def select_population(self, results_df):
        # Check if 'Passed_PoseBusters' column exists
        if 'Passed_PoseBusters' not in results_df.columns:
            print("No valid results to select from.")
            return []

        # Select molecules passing PoseBusters and with best docking scores
        selected_smiles = select_top_molecules(results_df, top_n=self.population_size)
        print(f"Selected {len(selected_smiles)} molecules for next generation.")
        return selected_smiles

    def mutate_population(self, selected_smiles):
        # Generate new population by mutating selected molecules
        mutated_smiles = mutate_smiles_list(selected_smiles, self.mutation_model, mutations_per_parent=self.mutations_per_parent, device=self.device)
        print(f"Generated {len(mutated_smiles)} mutated molecules.")
        return mutated_smiles

    def crossover_population(self, selected_smiles):
    # Generate new population by performing crossover
        offspring_smiles = crossover_smiles(
            parent_smiles_list=selected_smiles,
            num_offspring=self.crossover_offspring_per_generation,
            model=self.mutation_model,
            device=self.device
        )
        print(f"Generated {len(offspring_smiles)} crossover offspring.")
        return offspring_smiles


    def run(self):
        # Load mutation model
        self.mutation_model = load_model(MUTATION_MODEL_CONFIG_PATH, MUTATION_MODEL_CHECKPOINT_PATH, MUTATION_MODEL_VOCAB_PATH, self.device)
        # Initialize population
        self.initialize_population()
        for generation in range(1, self.generations + 1):
            print(f"\n=== Generation {generation} ===")
            # Create output directory for this generation
            generation_output_dir = os.path.join(OUTPUT_DIR, f'generation_{generation}')
            os.makedirs(generation_output_dir, exist_ok=True)
            # Evaluate population
            results_df = self.evaluate_population(generation, generation_output_dir)
            # Selection
            selected_smiles = self.select_population(results_df)
            if not selected_smiles:
                print("No molecules passed PoseBusters evaluation. Stopping generations.")
                break
            # Mutation
            mutated_smiles = self.mutate_population(selected_smiles)
            # Crossover
            crossover_smiles_list = self.crossover_population(selected_smiles)
            # Combine new population
            self.current_population = mutated_smiles + crossover_smiles_list
            # Ensure population size does not exceed the desired size
            if len(self.current_population) > self.population_size:
                self.current_population = self.current_population[:self.population_size]
            print(f"New population size: {len(self.current_population)}")

def main():
    # Parameters
    population_size = 5  # Number of molecules in the population
    generations = 3  # Number of generations
    mutations_per_parent = 2  # Number of mutations per parent molecule
    crossover_offspring_per_generation = 3  # Number of crossover offspring to generate per generation

    ga = GeneticAlgorithm(
        population_size=population_size,
        generations=generations,
        mutations_per_parent=mutations_per_parent,
        crossover_offspring_per_generation=crossover_offspring_per_generation
    )
    ga.run()

if __name__ == '__main__':
    main()
