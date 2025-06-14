[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15594168.svg)](https://doi.org/10.5281/zenodo.15594168)

# TexturedPerovskiteSolarCells.jl -- Numerical examples to analyse the opto-electronic behaviour of textured perovskite solar cells

**TexturedPerovskiteSolarCells.jl** contains all scripts and data needed to reproduce the numerical results and figures from the manuscript:
> *Unravelling the mystery of enhanced open-circuit voltages in nanotextured perovskite solar cells*
> 
> Dilara Abdel, Jacob Relle, Thomas Kirchartz, Patrick Jaap, Jürgen Fuhrmann, Sven Burger, Christiane Becker, Klaus Jäger, and Patricio Farrell.

## Overview

The simulations and postprocessing scripts in this repository rely on the Julia package [**ChargeTransport.jl**](https://github.com/WIAS-PDELib/ChargeTransport.jl), which solves the drift-diffusion equations using the Voronoi finite volume method. This numerical method is implemented via [**VoronoiFVM.jl**](https://github.com/WIAS-PDELib/VoronoiFVM.jl).

## Directory Structure

The repository is organised into the following main folders:

- **`data/`**
  Contains all optical photogeneration files, generated with [**JCMSuite**](https://jcmwave.com/jcmsuite). These files serve as input for the drift-diffusion simulations. This folder also includes additional data required to reproduce the figures in the manuscript, stored via Git LFS. These files can also be regenerated using the scripts provided in the `scripts/` folder.

- **`PostProcess/`**
  Includes all scripts used to postprocess simulation results and generate the figures presented in the manuscript.

- **`scripts/`**
  Contains the main simulation scripts, including definitions of physical parameters and configurations for generating solutions and J–V curves.

- **`src/`**
  Includes auxiliary scripts such as grid definitions and the photogeneration data reader.

## Usage

To reproduce the results from the paper, you will need a Julia installation (version ≥ 1.10.0).

1. **Clone this repository**
```
git clone https://github.com/NUMSEMIC/TexturedPerovskiteSolarCells.jl
cd TexturedPerovskiteSolarCells.jl
```

2. **Start julia and instantiate the project dependencies (first time only; it takes a while and it will download all necessary dependencies and the pyplot backend)**
```
julia --project
julia> using Pkg; Pkg.instantiate()
```

3. **Run simulations**
The following command generates a basic one-dimensional set-up and shows solution and J–V curve plots:

```
julia> include("scripts/SingleJunction.jl"); SingleJunction.main(gridDim = 1, plotting = true)
```

4. **Reproduce figures from the opto-electronic simulations in the manuscript**
Each figure can be generated by running the corresponding script from the `PostProcess/` folder. For example, to generate and save the figures from Figure 3:
```
julia> include("PostProcess/Fig3CharacteristicsStudy.jl"); Fig3CharacteristicsStudy.main(saveFig = true)
```
