# Genetic SMILES Optimizer

This repository implements a **genetic algorithm pipeline for de novo ligand discovery**, combining deep generative models with docking-based evaluation and diversity control.

The workflow integrates:

* **RDKit** for SMILES parsing, conformer generation, and fingerprinting
* **PyTorch** for pretrained generative/mutation model inference
* **AutoDock-GPU** for docking
* **Open Babel** and **Gypsum-DL** for 3-D structure generation and format conversions
* **PoseBusters** for stereochemical validation

---

## Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/rajkumar1501/genetic-smiles-optimizer.git
   cd genetic-smiles-optimizer
   ```

2. **Create and activate the conda environment**
   The environment is defined in `environment.yml`:

   ```bash
   conda env create -f environment.yml
   conda activate mol2mol_env
   ```

3. **Prepare receptor grids (.fld)** for AutoDock-GPU using your target protein structure.

---

## Workflow

### 1. Initial Population Generation

* **Generative sampling**: molecules are drawn from a pretrained SMILES generative model.
* **Seed-based initialization**: if a `--seed_smiles` is provided, the initial population is split evenly between seed mutations (via MolFormer) and generative samples.

### 2. Generation Loop

Each generation proceeds through:

1. **Molecular preparation & docking**

   * SMILES → 3-D coordinates with RDKit & Open Babel
   * Docking via AutoDock-GPU with multiple runs (`--nrun`)
   * Best pose → SDF, validated by PoseBusters
   * Molecules failing PoseBusters checks are discarded

2. **Selection with diversity control**

   * Candidates ranked by docking score
   * Morgan fingerprints computed
   * Greedy selection ensures Tanimoto similarity ≤ `--similarity_threshold`
   * Top fraction retained as elite pool

3. **Reproduction (mutation & crossover)**

   * Mutation: MolFormer generates SMILES variants per parent
   * Crossover: two parents spliced → repaired/expanded by mutation model
   * Elites always carried forward; duplicates removed; population shuffled

### 3. Termination

* Loop continues for `--generations` or until no viable offspring survive
* Outputs include docking scores, PoseBusters results, diversity analysis, and final candidate set

---

## Usage

Run the pipeline with:

```bash
python genetic_algorithm.py [options]
```

### Arguments

* `--population_size` (int, default: **10**) – Number of molecules in the population
* `--generations` (int, default: **5**) – Number of generations
* `--mutations_per_parent` (int, default: **5**) – Mutations per parent
* `--crossover_offspring_per_generation` (int, default: **5**) – Number of crossover offspring
* `--seed_smiles` (str) – Optional seed SMILES string
* `--nrun` (int, default: **100**) – AutoDock-GPU runs per docking
* `--receptor_fld` (str, default: `example/4uxl.maps.fld`) – Path to receptor `.fld` file
* `--output_dir` (str, default: `output_files`) – Directory for results
* `--seed` (int, default: **1337**) – Random seed
* `--diversity_filter` – Enable diversity control in selection
* `--similarity_threshold` (float, default: **0.40**) – Max Tanimoto similarity allowed
* `--fp_radius` (int, default: **2**) – Morgan fingerprint radius
* `--fp_bits` (int, default: **2048**) – Number of fingerprint bits

---

## Examples

* **Default run**

  ```bash
  python genetic_algorithm.py
  ```

* **Larger population and more generations**

  ```bash
  python genetic_algorithm.py --population_size 50 --generations 20
  ```

* **Seed-guided initialization**

  ```bash
  python genetic_algorithm.py --seed_smiles "CCO"
  ```

* **Docking with custom receptor**

  ```bash
  python genetic_algorithm.py --receptor_fld my_receptor.fld --nrun 200
  ```

* **Enable diversity filter**

  ```bash
  python genetic_algorithm.py --diversity_filter --similarity_threshold 0.3
  ```

---

## Output

Results are stored in `output_files/` (or `--output_dir`):

* Docking scores per molecule
* PoseBusters validation flags
* Best poses in `.pdbqt` and `.sdf` formats
* Final population statistics
