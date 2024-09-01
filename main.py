from genetic_smiles_optimizer import (
    generate_initial_population, 
    generate_similar_smiles,
    smiles_to_pdbqt, 
    perform_docking, 
    pdbqt_to_sdf, 
    evaluate_pose, 
    calculate_fitness
)

def main():
    # Step 1: Generate Initial Population of SMILES
    print("Generating initial population of SMILES...")
    initial_population = generate_initial_population(10)
    print("Initial Population:", initial_population)

    # Step 2: Process Each SMILES String
    protein_maps_file = 'data/protein/protein.maps.fld'
    protein_file = 'data/protein/protein.pdbqt'

    for smiles in initial_population:
        print(f"\nProcessing SMILES: {smiles}")
        
        # Convert SMILES to PDBQT
        ligand_pdbqt_file = "ligand.pdbqt"
        print(f"Converting SMILES to PDBQT: {ligand_pdbqt_file}")
        smiles_to_pdbqt(smiles, ligand_pdbqt_file)
        
        # Perform docking
        print(f"Performing docking for {ligand_pdbqt_file}...")
        perform_docking(protein_maps_file, ligand_pdbqt_file)
        
        # Convert PDBQT to SDF
        sdf_file = "ligand.sdf"
        print(f"Converting PDBQT to SDF: {sdf_file}")
        pdbqt_to_sdf(ligand_pdbqt_file, sdf_file)
        
        # Evaluate the docked pose and calculate fitness
        print(f"Evaluating pose for {sdf_file}...")
        results = evaluate_pose(sdf_file, protein_file)
        fitness_scores = calculate_fitness(results)
        
        # Print the fitness scores
        for molecule, score in fitness_scores.items():
            print(f"Molecule: {molecule}, Fitness Score: {score}")

    # Step 3: Generate Similar SMILES and Handle Results
    print("\nGenerating similar SMILES...")
    similar_smiles_df = generate_similar_smiles(initial_population, 1000)
    print("Similar SMILES Generated:")
    print(similar_smiles_df)

    # The pipeline can be expanded further with more iterations, optimizations,
    # or even integrating user feedback based on the results.

if __name__ == "__main__":
    main()
