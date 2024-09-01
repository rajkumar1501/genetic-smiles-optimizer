<<<<<<< HEAD
# Genetic SMILES Optimizer
=======
# README.md
[![DOI](https://zenodo.org/badge/833483470.svg)](https://zenodo.org/doi/10.5281/zenodo.12958255)
>>>>>>> 95854f02b79c56b8bf3a5a86f6e09e4b79557368

## Project Overview

The Genetic SMILES Optimizer is a computational tool designed to generate and optimize drug-like molecules represented by SMILES (Simplified Molecular Input Line Entry System) strings. The project leverages the power of genetic algorithms to evolve chemical structures towards desired properties, making it a valuable tool in the early stages of drug discovery.

### Key Components:

### 1. Drug Molecule Generation:
- **Objective**: Create an initial pool of diverse, chemically valid SMILES strings that represent potential drug-like molecules.
- **Approach**: Use a chemical generative model to produce a wide range of SMILES strings. These generated molecules serve as the starting point for the optimization process.

### 2. Drug Molecule Optimization:
- **Objective**: Optimize the generated SMILES strings to enhance their drug-like properties, such as binding affinity to a target protein and overall stability.
- **Approach**: Implement a genetic algorithm that iteratively refines the SMILES strings. The algorithm applies evolutionary techniques, including selection, mutation, and crossover, to evolve the molecules toward improved characteristics.

### 3. Scoring Function:
- **Docking Energy**: Assess the binding affinity of the molecules to a target protein using molecular docking simulations. This metric helps identify molecules with high potential as drug candidates.
- **PoseBusters**: Evaluate the structural validity and quality of the docked poses using PoseBusters. This tool ensures that the generated molecules are not only chemically feasible but also biologically relevant.

## Project Goals

The ultimate goal of this project is to develop a robust pipeline for the discovery and optimization of drug-like molecules. By combining SMILES generation with a genetic algorithm-driven optimization process, the tool aims to efficiently explore chemical space and identify promising candidates for further development in drug discovery.

This project represents a significant step toward automating the early stages of drug design, enabling researchers to rapidly generate and refine potential therapeutic compounds.
