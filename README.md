# PIXL Martian Mineral Chemistry Analysis

Python-based analysis tools for investigating the mineral chemistry and spatial distribution of Martian surface materials measured by the Planetary Instrument for X-ray Lithochemistry (PIXL) aboard NASA's Perseverance rover.

The workflow processes PIXL elemental/oxide data to identify chemically distinct domains, characterize their composition and geometry, and investigate their spatial distribution. Particular emphasis is placed on identifying probable olivine-bearing domains and calculating their forsterite–fayalite composition.

### Main capabilities

- Processing and filtering of PIXL chemical datasets
- Separation of observations by target and Sol
- Chemical domain and edge detection
- Calculation of domain area, characteristic dimensions, and compositional statistics
- Identification and characterization of probable olivine domains
- Calculation of Fo–Fa compositions
- Harker and ternary compositional diagrams
- RGB chemical mapping using multiple oxides
- Grain/domain size distributions
- Geographic visualization of PIXL observations on Mars
- Mapping of the spatial distribution of olivine-bearing material

The project is developed in Python using scientific libraries including `pandas`, `NumPy`, `SciPy`, `Matplotlib`, `mpltern`, and `adjustText`.


### Setting up

When running the code for the first time, you will have to install the necessary packages from the environment.yml file into an environment in python.

The first thing you will want to do is:
conda activate PlanetaryEnv
to activate the environment in the IDE