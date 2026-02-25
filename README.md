# Neuromodulation and homeostasis: complementary mechanisms for robust neural function

## About this repository

This repository contains all code and data involved in **Neuromodulation and homeostasis: complementary mechanisms for robust neural function**.

## Getting started

This work uses the [Julia programming language](https://julialang.org/). To set up the environment, download Julia, then open a Julia REPL, navigate to this repository, and run:

```jl
include("dependencies.jl")
```

This will install all required packages.

## Repository structure

Each figure with associated simulations has its own folder. Within each folder, you will find:

- `.jl` files implementing model equations and utility functions
- A `data/` folder containing pre-generated data needed to reproduce the figures
- `.ipynb` notebooks that run the experiments and produce the figures

Figures 1 and 4 are conceptual schematics and do not have associated code.

| Folder | Description |
|--------|-------------|
| `figure 2 controlled/` | Controlled neuromodulation with homeostasis |
| `figure 2 sharp/` | Sharp neuromodulation with homeostasis |
| `figure 3 controlled/` | Controlled neuromodulation with homeostasis |
| `figure 3 sharp/` | Sharp neuromodulation with homeostasis |
| `figure 5/` | Extended analysis with blockades |
| `figure 6/` | Network-level simulations |

## Launching the notebooks

To open JupyterNotebook or JupyterLab using Julia, execute:

```jl
using IJulia
notebook() # or jupyterlab()
```

Then browse to the `.ipynb` files to run them. If this is your first time launching Jupyter, Julia will prompt you to install it through Conda.

## Author

Code written by Arthur Fyon.