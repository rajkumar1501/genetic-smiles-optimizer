import random
from typing import List, Tuple
from smiles_generation import generate_initial_population, generate_similar_smiles
from docking import convert_smiles_to_pdbqt, perform_docking
from pose_evaluation import evaluate_pose

# Define constants for the genetic algorithm
POPULATION_SIZE = 100
NUM_GENERATIONS = 50
MUTATION_RATE = 0.1
CROSSOVER_RATE = 0.5
SELECTION_SIZE = 20
NUM_SAMPLES = 10

def genetic_algorithm(target_protein: str, reference_ligand: str) -> str:
    """
    Main genetic algorithm loop for optimizing SMILES strings.

    Args:
        target_protein (str): Path to the target protein PDB file.
        reference_ligand (str): Path to the reference ligand SDF file.

    Returns:
        str: The best SMILES string found after optimization.
    """
    # Step 1: Generate initial population of SMILES strings
    population = generate_initial_population(POPULATION_SIZE)
    
    for generation in range(NUM_GENERATIONS):
        print(f"Generation {generation + 1}/{NUM_GENERATIONS}")

        scored_population = []
        
        # Step 2: Evaluate fitness of each SMILES string in the population
        for smiles in population:
            pdbqt_file = f'{smiles}.pdbqt'
            docked_file = f'{smiles}_docked.pdbqt'
            
            convert_smiles_to_pdbqt(smiles, pdbqt_file)
            perform_docking(pdbqt_file, target_protein, docked_file)
            score = evaluate_pose(docked_file, reference_ligand, target_protein)
            
            scored_population.append((smiles, score))
        
        # Step 3: Select the top individuals based on fitness scores
        top_individuals = select_top_individuals(scored_population, SELECTION_SIZE)
        
        # Step 4: Generate the next generation
        next_generation = []
        
        while len(next_generation) < POPULATION_SIZE:
            if random.random() < CROSSOVER_RATE:
                parent1, parent2 = random.sample(top_individuals, 2)
                offspring = crossover(parent1, parent2)
            else:
                offspring = random.choice(top_individuals)
            
            if random.random() < MUTATION_RATE:
                offspring = mutate(offspring)
            
            next_generation.append(offspring)
        
        population = next_generation
    
    # Step 5: Return the best individual found
    best_individual = max(scored_population, key=lambda x: x[1])
    print(f'Best SMILES: {best_individual[0]}, Score: {best_individual[1]}')
    return best_individual[0]

def select_top_individuals(population: List[Tuple[str, float]], top_n: int) -> List[str]:
    """
    Select the top N individuals based on fitness scores.

    Args:
        population (List[Tuple[str, float]]): The population with their fitness scores.
        top_n (int): Number of top individuals to select.

    Returns:
        List[str]: List of top N SMILES strings.
    """
    population.sort(key=lambda x: x[1], reverse=True)
    return [individual[0] for individual in population[:top_n]]

def mutate(smiles: str) -> str:
    """
    Mutate a SMILES string by generating a similar SMILES string.

    Args:
        smiles (str): The original SMILES string.

    Returns:
        str: A mutated SMILES string.
    """
    similar_smiles_list = generate_similar_smiles([smiles], 1)
    return similar_smiles_list[0] if similar_smiles_list else smiles

def crossover(parent1: str, parent2: str) -> str:
    """
    Perform crossover between two SMILES strings to produce an offspring.

    Args:
        parent1 (str): The first parent SMILES string.
        parent2 (str): The second parent SMILES string.

    Returns:
        str: The offspring SMILES string resulting from the crossover.
    """
    # Simple crossover logic: split and combine
    split_point = random.randint(1, len(parent1) - 1)
    return parent1[:split_point] + parent2[split_point:]

# Example usage
if __name__ == "__main__":
    target_protein_file = 'protein.pdb'  # Replace with actual target protein file
    reference_ligand_file = 'crystal_ligand.sdf'  # Replace with actual reference ligand file
    best_smiles = genetic_algorithm(target_protein_file, reference_ligand_file)
    print(f"The optimized SMILES string is: {best_smiles}")
