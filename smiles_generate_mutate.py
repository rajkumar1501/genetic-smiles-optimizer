from molecule_generation import load_model_from_directory
from functools import partial
import subprocess
import pandas as pd
from typing import List


from time import time
import argparse
import os
import yaml


from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity
from scipy.stats import pearsonr
import numpy as np
import torch
from mutate_smiles import load_model, generate_samples

from molecule_generation import load_model_from_directory

SMILES_GENERATION_MODEL_PATH = '/mnt/d/Projects/molgen_weights'

def generate_initial_population(size: int) -> List[str]:
    """
    Generate an initial population of SMILES strings.

    Args:
        size (int): The number of SMILES strings to generate.

    Returns:
        List[str]: A list of generated SMILES strings.
    """
    model_kwargs = {}  # Add any model-specific keyword arguments if needed
    
    try:
        with load_model_from_directory(SMILES_GENERATION_MODEL_PATH, **model_kwargs) as model:
            samples = model.sample(size)
        return samples
    except Exception as e:
        print(f"Error generating initial SMILES: {e}")
        return []
    
config_path = os.path.join('paper_checkpoints/ecfp4_with_counts_with_rank', "config.yml")
checkpoint_path = os.path.join('paper_checkpoints/ecfp4_with_counts_with_rank', "weights.ckpt")
vocabulary_path = "paper_checkpoints/vocabulary.pkl"

initial_population = generate_initial_population(2)
print("Initial Population:", initial_population)
device = "cpu"
if torch.cuda.is_available():
     device = "cuda"

model = load_model(config_path, checkpoint_path, vocabulary_path, device)

samples = []
for smi in initial_population:
    smi = smi.strip()
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"Cannot understand SMILES: {smi}")
            continue
    except BaseException:
        print(f"Cannot understand SMILES: {smi}")
        continue
        
    _samples = generate_samples(model, smi, beam_size=10, device=device)
    for new_smi in _samples:
                try:
                    mol = Chem.MolFromSmiles(new_smi)
                    if mol is None:
                        print(f"Cannot understand SMILES: {new_smi}")
                        continue
                    else:
                        samples.append(new_smi)
                except BaseException:
                        print(f"Cannot understand SMILES: {new_smi}")

                        continue

print(samples)


