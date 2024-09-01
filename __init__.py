# __init__.py
# ----------------------------------------
# This file marks the directory as a Python package and is used to control what is 
# accessible when the package is imported. By selectively importing functions 
# from different modules within the package, we make the code cleaner and easier to use.
# ----------------------------------------

# Step 1: Import SMILES Generation Functions
# ------------------------------------------
# We start by importing functions related to SMILES generation. These functions are 
# essential for creating the initial pool of candidate molecules and for generating 
# similar molecules during the optimization process.
from .smiles_generation import generate_initial_population, generate_similar_smiles

# Step 2: Import Docking Functions
# --------------------------------
# Next, we import functions that handle the conversion of chemical formats and 
# the docking process itself. These functions are critical for preparing the molecules 
# for docking simulations and converting the results back into formats that can be analyzed.
from .docking import smiles_to_pdbqt, pdbqt_to_pdb, pdbqt_to_sdf, pdbqt_to_smiles, perform_docking

# Step 3: Import Pose Evaluation Functions
# ----------------------------------------
# Finally, we import functions for evaluating the docked poses. These functions 
# use PoseBusters to assess how well the docked molecules fit the target protein, 
# and to calculate fitness scores that help guide the optimization process.
from .pose_evaluation import evaluate_pose, calculate_fitness

# Summary
# ----------------------------------------
# By organizing the imports this way, we provide a clear and logical structure for 
# the package. Users of the package can easily access the necessary functions for 
# generating SMILES, performing docking, and evaluating poses without needing to know 
# the internal structure of the package. This makes the package more user-friendly 
# and helps to maintain clean and maintainable code.
