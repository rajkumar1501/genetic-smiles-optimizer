# genetic_algorithm.py

import os
import random
import subprocess
import argparse
import shutil
from pathlib import Path
from typing import Optional, List

import torch
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit import DataStructs


# External modules you provide
from docking_script import run_autodock_gpu
from evaluator import evaluate_molecules_with_posebusters
from mutate_smiles import load_model, generate_samples
from molecule_generation import load_model_from_directory

# Silence RDKit
RDLogger.DisableLog("rdApp.*")

# ----------------------- Config (defaults; can be overridden by CLI) -----------------------

SMILES_GENERATION_MODEL_PATH = Path("molgen_weights")  # update as needed
MUTATION_MODEL_CONFIG_PATH = Path("paper_checkpoints/ecfp4_with_counts_with_rank/config.yml")
MUTATION_MODEL_CHECKPOINT_PATH = Path("paper_checkpoints/ecfp4_with_counts_with_rank/weights.ckpt")
MUTATION_MODEL_VOCAB_PATH = Path("paper_checkpoints/vocabulary.pkl")

RECEPTOR_FLD_PATH = Path("example/4uxl.maps.fld").resolve()

OUTPUT_DIR = Path("output_files")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

ALLOWED_ELEMENTS = {"H", "C", "N", "O", "F", "P", "S", "Cl", "I"}

# ----------------------- Utils -----------------------

def canonicalize_smiles(smiles: str) -> Optional[str]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        return None

def dedupe_smiles(smiles_list: List[str]) -> List[str]:
    seen = set()
    uniq = []
    for s in smiles_list:
        c = canonicalize_smiles(s)
        if c and c not in seen:
            seen.add(c)
            uniq.append(c)
    return uniq

def is_valid_molecule(smiles: str, allowed_elements=ALLOWED_ELEMENTS) -> bool:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        elements_in_mol = {atom.GetSymbol() for atom in mol.GetAtoms()}
        if not elements_in_mol.issubset(allowed_elements):
            return False
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDG()
        params.randomSeed = 0xF00D
        if AllChem.EmbedMolecule(mol, params) != 0:
            return False
        # UFF returns 0 on success
        if AllChem.UFFOptimizeMolecule(mol) != 0:
            return False
        return True
    except Exception:
        return False

# ----------------------- Generation -----------------------

def generate_initial_population(desired_size: int, batch_size: int = 8) -> List[str]:
    valid_smiles = []
    total_generated = 0
    try:
        with load_model_from_directory(SMILES_GENERATION_MODEL_PATH) as model:
            while len(valid_smiles) < desired_size:
                samples = model.sample(batch_size)
                total_generated += batch_size
                for s in samples:
                    if is_valid_molecule(s):
                        valid_smiles.append(s)
                        if len(valid_smiles) >= desired_size:
                            break
        print(f"Generated {len(valid_smiles)} valid SMILES out of {total_generated} generated.")
        return dedupe_smiles(valid_smiles)[:desired_size]
    except Exception as e:
        print(f"Error generating initial SMILES: {e}")
        return []

# ----------------------- Dock + Evaluate -----------------------

def dock_and_evaluate(smiles_list: List[str], iteration: int, generations: int, output_dir: Path, nrun: int = 100) -> pd.DataFrame:
    results = []
    output_dir.mkdir(parents=True, exist_ok=True)

    for idx, smiles in enumerate(smiles_list):
        print(f"\n=== Generation {iteration}/{generations}, Molecule {idx+1}/{len(smiles_list)} ===")
        try:
            ligand_pdbqt = output_dir / f"ligand_{iteration}_{idx}.pdbqt"
            best_pdbqt = output_dir / f"best_{iteration}_{idx}.pdbqt"

            best_path, docking_score = run_autodock_gpu(
                smiles=smiles,
                receptor_fld_path=str(RECEPTOR_FLD_PATH),
                output_dir=str(output_dir),
                nrun=nrun,
                ligand_pdbqt_filename=str(ligand_pdbqt.name),
                best_pdbqt_filename=str(best_pdbqt.name),
            )
            if not best_path or docking_score is None:
                print(f"Docking failed or score missing for: {smiles}")
                continue

            # Convert PDBQT -> SDF via Open Babel
            best_sdf = output_dir / f"best_{iteration}_{idx}.sdf"
            obabel_cmd = [
                "obabel",
                "-ipdbqt",
                str(best_path),
                "-osdf",
                "-O",
                str(best_sdf),
            ]
            conv = subprocess.run(obabel_cmd, capture_output=True, text=True)
            if conv.returncode != 0 or not best_sdf.exists():
                print(f"Open Babel failed for {best_path}: {conv.stderr.strip()}")
                continue

            evaluation_df = evaluate_molecules_with_posebusters(str(best_sdf))
            col = "All_PoseBusters_Parameters_Passed"
            passed = bool(evaluation_df[col].iloc[0]) if col in evaluation_df.columns else False

            results.append(
                {
                    "SMILES": canonicalize_smiles(smiles) or smiles,
                    "Docking_Score": docking_score,
                    "Passed_PoseBusters": passed,
                    "PDBQT_File": str(Path(best_path).resolve()),
                    "SDF_File": str(best_sdf.resolve()),
                }
            )
        except Exception as e:
            print(f"Error processing molecule {smiles}: {e}")
            continue

    return pd.DataFrame(results)


def _smiles_to_fp(smiles: str, radius: int = 2, n_bits: int = 2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)


def _select_top_diverse_strict(
    df: pd.DataFrame,
    top_n: int,
    sim_threshold: float = 0.6,
    radius: int = 2,
    n_bits: int = 2048,
) -> List[str]:
    """
    Greedy + STRICT diversity:
    - Sort by Docking_Score ascending (best first).
    - Add a candidate only if its max Tanimoto to all selected ≤ sim_threshold.
    - If no remaining candidate fits, STOP (even if < top_n selected).
    """
    if df.empty:
        return []

    pool = df[df["Passed_PoseBusters"] == True].copy()
    if pool.empty:
        return []

    pool.sort_values("Docking_Score", inplace=True, ascending=True)
    pool["SMILES"] = pool["SMILES"].apply(lambda s: canonicalize_smiles(s) or s)

    # precompute fingerprints
    rows, fps = [], []
    for _, row in pool.iterrows():
        fp = AllChem.GetMorganFingerprintAsBitVect(
            Chem.MolFromSmiles(row["SMILES"]), radius, nBits=n_bits
        )
        if fp is not None:
            rows.append(row)
            fps.append(fp)
    if not rows:
        return []

    selected_idx = []
    remaining_idx = list(range(len(rows)))

    # seed with best scoring candidate
    selected_idx.append(0)
    remaining_idx.remove(0)

    while len(selected_idx) < top_n and remaining_idx:
        added = False
        for i in list(remaining_idx):
            cand_fp = fps[i]
            sims = [DataStructs.TanimotoSimilarity(cand_fp, fps[j]) for j in selected_idx]
            if max(sims) <= sim_threshold:
                selected_idx.append(i)
                remaining_idx.remove(i)
                added = True
                if len(selected_idx) >= top_n:
                    break
        if not added:
            # nothing else fits the threshold → strict stop
            break

    return [rows[i]["SMILES"] for i in selected_idx]


def select_top_molecules(
    results_df: pd.DataFrame,
    top_n: int,
    use_diversity: bool = True,
    sim_threshold: float = 0.6,
    radius: int = 2,
    n_bits: int = 2048,
    strict_diversity: bool = True,   # NEW: default strict
) -> List[str]:
    if results_df is None or results_df.empty:
        print("No results to select from.")
        return []
    required = {"Passed_PoseBusters", "Docking_Score", "SMILES"}
    if not required.issubset(results_df.columns):
        print("Results missing required columns.")
        return []

    if use_diversity:
        if strict_diversity:
            return _select_top_diverse_strict(results_df, top_n, sim_threshold, radius, n_bits)
        else:
            # your previous non-strict version (which backfilled) can be kept as _select_top_diverse(...)
            return _select_top_diverse(results_df, top_n, sim_threshold, radius, n_bits)

    # score-only fallback
    filtered = results_df[results_df["Passed_PoseBusters"] == True]
    if filtered.empty:
        print("No molecules passed PoseBusters.")
        return []
    return filtered.nsmallest(top_n, "Docking_Score")["SMILES"].tolist()

# ----------------------- Crossover + Mutation -----------------------

def perform_smiles_crossover(smiles1: str, smiles2: str) -> Optional[str]:
    # very naive splice; used only as a seed for the mutation model
    t1, t2 = list(smiles1), list(smiles2)
    if len(t1) < 2 or len(t2) < 2:
        return None

    # choose alnum splice points
    idxs1 = [i for i, c in enumerate(t1) if c.isalnum()]
    idxs2 = [i for i, c in enumerate(t2) if c.isalnum()]
    if not idxs1 or not idxs2:
        return None

    i1 = random.choice(idxs1)
    i2 = random.choice(idxs2)
    child = "".join(t1[:i1] + t2[i2:])

    # quick sanity check before spending model time
    if Chem.MolFromSmiles(child) is None:
        return None
    return child

def crossover_smiles(parent_smiles_list: List[str], num_offspring: int, model, device: str = "cpu") -> List[str]:
    offspring = []
    parents = list(set(parent_smiles_list))
    if len(parents) < 2:
        print("Not enough parents for crossover.")
        return offspring

    attempts = 0
    max_attempts = max(50, num_offspring * 10)

    while len(offspring) < num_offspring and attempts < max_attempts:
        attempts += 1
        p1, p2 = random.sample(parents, 2)
        seed = perform_smiles_crossover(p1, p2)
        if not seed:
            continue
        try:
            # small beam; mutation model repairs invalid bits
            mutated = generate_samples(model, seed, beam_size=5, device=device)
            for m in mutated:
                if is_valid_molecule(m):
                    offspring.append(canonicalize_smiles(m) or m)
                    if len(offspring) >= num_offspring:
                        break
        except Exception as e:
            print(f"Error mutating crossover seed '{seed}': {e}")
            continue

    return dedupe_smiles(offspring)[:num_offspring]

def mutate_smiles_list(smiles_list: List[str], model, mutations_per_parent: int = 5, device: str = "cpu") -> List[str]:
    new_smiles = []
    for s in smiles_list:
        try:
            if Chem.MolFromSmiles(s) is None:
                print(f"Cannot parse SMILES: {s}")
                continue
            mutated = generate_samples(model, s, beam_size=mutations_per_parent, device=device)
            new_smiles.extend(mutated)
        except Exception as e:
            print(f"Error generating samples for SMILES '{s}': {e}")
            continue
    # validate + canonicalize once here to prune junk early
    return dedupe_smiles([m for m in new_smiles if is_valid_molecule(m)])

def _mutations_from_seed(seed: str, want: int, model, device: str) -> List[str]:
    """Generate up to `want` valid, deduped SMILES by mutating a seed."""
    if Chem.MolFromSmiles(seed) is None:
        print(f"Provided seed SMILES is invalid: {seed}")
        return []
    # oversample a bit to survive validation/dedupe
    beam = max(want * 5, 32)
    try:
        pool = generate_samples(model, seed, beam_size=beam, device=device)
    except Exception as e:
        print(f"Error mutating seed '{seed}': {e}")
        return []
    pool = [m for m in pool if is_valid_molecule(m)]
    pool = dedupe_smiles(pool)
    return pool[:want]

# ----------------------- GA -----------------------

class GeneticAlgorithm:
    def __init__(
        self,
        population_size: int,
        generations: int,
        mutations_per_parent: int = 5,
        crossover_offspring_per_generation: int = 10,
        seed_smiles: Optional[str] = None,
        nrun: int = 100,
    ):
        self.population_size = population_size
        self.generations = generations
        self.mutations_per_parent = mutations_per_parent
        self.crossover_offspring_per_generation = crossover_offspring_per_generation
        self.current_population: List[str] = []
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.mutation_model = None
        self.seed_smiles = seed_smiles
        self.nrun = nrun

    def initialize_population(self):
        print("\n=== Generating Initial Population ===")
        if self.seed_smiles:
            half_a = self.population_size // 2
            half_b = self.population_size - half_a

            # half from generative model
            gen_half = generate_initial_population(half_a)

            # half from mutations of the provided seed
            if self.mutation_model is None:
                raise RuntimeError("Mutation model must be loaded before seeding from SMILES.")
            mut_half = _mutations_from_seed(self.seed_smiles, half_b, self.mutation_model, self.device)

            # combine, dedupe, and top-up if needed
            combined = dedupe_smiles(gen_half + mut_half)

            # top up if dedupe shrank the list
            attempts = 0
            while len(combined) < self.population_size and attempts < 5:
                attempts += 1
                need = self.population_size - len(combined)
                extra_mut = _mutations_from_seed(self.seed_smiles, need, self.mutation_model, self.device)
                combined = dedupe_smiles(combined + extra_mut)
                if len(combined) >= self.population_size:
                    break
                extra_gen = generate_initial_population(need)
                new_combined = dedupe_smiles(combined + extra_gen)
                if len(new_combined) == len(combined):
                    # cannot top-up further
                    break
                combined = new_combined

            self.current_population = combined[:self.population_size]
            if len(self.current_population) < self.population_size:
                print(f"Warning: initial population shortfall {len(self.current_population)}/{self.population_size}.")
        else:
            # original behavior
            self.current_population = generate_initial_population(self.population_size)

        self.current_population = dedupe_smiles(self.current_population)[:self.population_size]
        print(f"Initial Population ({len(self.current_population)}): {self.current_population}")

    def evaluate_population(self, generation: int, generations: int, generation_output_dir: Path) -> pd.DataFrame:
        df = dock_and_evaluate(self.current_population, generation, generations, generation_output_dir, nrun=self.nrun)
        results_csv = generation_output_dir / f"docking_results_generation_{generation}.csv"
        df.to_csv(results_csv, index=False)
        print(f"Results saved to {results_csv}")
        return df

    def select_population(self, results_df: pd.DataFrame) -> List[str]:
        selected = select_top_molecules(
            results_df,
            top_n=max(1, int(self.population_size * 0.1)),
            use_diversity=getattr(self, "use_diversity", True),
            sim_threshold=getattr(self, "sim_thresh", 0.6),
            radius=getattr(self, "fp_radius", 2),
            n_bits=getattr(self, "fp_bits", 2048),
        )
        print(f"Selected {len(selected)} molecules for next generation. \n Candiadtes {selected}\n")
        return selected


    def mutate_population(self, selected_smiles: List[str]) -> List[str]:
        mutated = mutate_smiles_list(
            selected_smiles,
            self.mutation_model,
            mutations_per_parent=self.mutations_per_parent,
            device=self.device,
        )
        print(f"Generated {len(mutated)} mutated molecules.")
        return mutated

    def crossover_population(self, selected_smiles: List[str]) -> List[str]:
        offspring = crossover_smiles(
            parent_smiles_list=selected_smiles,
            num_offspring=self.crossover_offspring_per_generation,
            model=self.mutation_model,
            device=self.device,
        )
        print(f"Generated {len(offspring)} crossover offspring.")
        return offspring

    def run(self, seed: int = 1337):
        # deterministic-ish behavior
        random.seed(seed)
        try:
            import numpy as np
            np.random.seed(seed)
        except Exception:
            pass
        if self.device == "cuda":
            try:
                torch.manual_seed(seed)
                torch.cuda.manual_seed_all(seed)
            except Exception:
                pass

        # Load mutation model
        self.mutation_model = load_model(
            str(MUTATION_MODEL_CONFIG_PATH),
            str(MUTATION_MODEL_CHECKPOINT_PATH),
            str(MUTATION_MODEL_VOCAB_PATH),
            self.device,
        )

        # Initialize
        self.initialize_population()

        for generation in range(1, self.generations + 1):
            print(f"\n=== Generation {generation} / {self.generations} ===")
            gen_dir = OUTPUT_DIR / f"generation_{generation}"
            gen_dir.mkdir(parents=True, exist_ok=True)

            # Always evaluate current population
            results_df = self.evaluate_population(generation, self.generations, gen_dir)

            # If this is the final generation, stop: no selection/mutation/crossover
            if generation == self.generations:
                print("Final generation evaluated. Stopping without creating a new population.")
                break

            # Selection
            selected = self.select_population(results_df)
            if not selected:
                print("No molecules passed selection. Stopping.")
                break

            # Mutation + crossover
            mutated = self.mutate_population(selected)
            crossed = self.crossover_population(selected)

            # Combine, dedupe, shuffle BEFORE next round
            elite = dedupe_smiles(selected)
            candidates = dedupe_smiles(mutated + crossed)
            random.shuffle(candidates)  # shuffle after dedupe
            print(f"Candidates after dedupe: {len(candidates)}")

            if not candidates:
                self.current_population = elite[:self.population_size]
                print("No new candidates after dedupe; carrying elites forward.")
                continue

            elite_set = set(elite)
            offspring_only = [s for s in candidates if s not in elite_set]

            # Commit next population
            self.current_population = (elite + offspring_only)[:self.population_size]
            print(f"New population size: {len(self.current_population)}")

# ----------------------- Entry -----------------------

def main():
    global RECEPTOR_FLD_PATH, OUTPUT_DIR
    parser = argparse.ArgumentParser(description="Run Genetic Algorithm for molecule generation.")
    parser.add_argument("--population_size", type=int, default=10, help="Number of molecules in the population")
    parser.add_argument("--generations", type=int, default=5, help="Number of generations")
    parser.add_argument("--mutations_per_parent", type=int, default=5, help="Mutations per parent molecule")
    parser.add_argument("--crossover_offspring_per_generation", type=int, default=5,
                        help="Number of crossover offspring to generate per generation")
    parser.add_argument("--seed_smiles", type=str, default=None,
                        help="Optional seed SMILES. If provided, initial population is 50% seed-mutations and 50% model-generated.")
    parser.add_argument("--nrun", type=int, default=100, help="AutoDock-GPU nrun")
    parser.add_argument("--receptor_fld", type=str, default=str(RECEPTOR_FLD_PATH),
                        help="Path to receptor .fld file")
    parser.add_argument("--output_dir", type=str, default=str(OUTPUT_DIR),
                        help="Directory for outputs")
    parser.add_argument("--seed", type=int, default=1337, help="Random seed")
    parser.add_argument("--diversity_filter", action="store_true",
                    help="Enable diversity-aware selection using Tanimoto")
    parser.add_argument("--similarity_threshold", type=float, default=0.40,
                        help="Max Tanimoto similarity allowed among selected (default 0.4)")
    parser.add_argument("--fp_radius", type=int, default=2, help="ECFP radius (default 2)")
    parser.add_argument("--fp_bits", type=int, default=2048, help="ECFP nBits (default 2048)")


    args = parser.parse_args()

    # Validate external deps and inputs
    if shutil.which("obabel") is None:
        raise SystemExit("Open Babel 'obabel' not found in PATH.")
    fld = Path(args.receptor_fld)
    if not fld.exists():
        raise SystemExit(f"Receptor FLD not found: {fld}")

    # Apply runtime paths
    
    RECEPTOR_FLD_PATH = fld.resolve()
    OUTPUT_DIR = Path(args.output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.seed_smiles:
        control_dir = OUTPUT_DIR / "control"
        control_dir.mkdir(parents=True, exist_ok=True)

        print("\n=== Running control docking for provided seed SMILES ===")
        control_results_df = dock_and_evaluate(
            smiles_list=[args.seed_smiles],
            iteration=0,
            generations=0,  # use 0 to indicate control
            output_dir=control_dir,
            nrun=args.nrun
        )
        control_csv = control_dir / "docking_results_control.csv"
        control_results_df.to_csv(control_csv, index=False)
        print(f"Control docking results saved to {control_csv}")

    ga = GeneticAlgorithm(
        population_size=args.population_size,
        generations=args.generations,
        mutations_per_parent=args.mutations_per_parent,
        crossover_offspring_per_generation=args.crossover_offspring_per_generation,
        seed_smiles=args.seed_smiles,
        nrun=args.nrun,
    )
    ga.use_diversity = args.diversity_filter
    ga.sim_thresh = args.similarity_threshold
    ga.fp_radius = args.fp_radius
    ga.fp_bits = args.fp_bits
    ga.run(seed=args.seed)
    

if __name__ == "__main__":
    main()
