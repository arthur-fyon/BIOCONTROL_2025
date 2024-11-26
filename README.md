# Neuromodulation and homeostasis: complementary mechanisms for robust neural function

## What do you find in this repository?

In this repository, you will be able to find all codes and data that were involved in **Neuromodulation and homeostasis: complementary mechanisms for robust neural function**.

## How to use the codes?
In this work, the Julia programming language was used. To use the codes, first download the latest version of Julia [here](https://julialang.org/) if it is not already the case.
Once this is done, you can run *dependencies.jl* to download and install all the packages used in this work. To do so, open Julia, *cd* to the folder where you can find *dependencies.jl* and run
```jl
include("dependencies.jl")
```

## How are the codes organized?
In this repository, you have one folder per image. In all folders, you will find *.jl* files in which the model equations are implemented, as well as other useful functions. You will also find a *data* folder in which all the data necessary to reproduce the figures of the article are already generated. Finally, notebooks *.ipynb* are also available. These notebooks contain the code for running the experiments as well as producing the figures.

## Special note to launch Jupyter with Julia
To open JupyterNotebook/JupyterLab using Julia, execute in Julia
```jl
using IJulia
notebook() # Or jupyterlab()
```

Then simply browse to the *.ipynb* files to run them. If it is the first time you ever launch Jupyter, Julia will ask you if you want to install it through Conda. Accept to start the installation of Jupyter.

**Codes written by Arthur Fyon**
