# A Positive Feedback between Iceberg Calving and Local Sea Level through Hydrostatic Pressure
Code repository for:

**Ambelorun A. A., Robel, A. A., & Seroussi, H.** "A Positive Feedback between Iceberg Calving and Local Sea Level through Hydrostatic Pressure."

## Overview
This repository contains MATLAB scripts to reproduce all the figures in the manuscript. These scripts require the Ice-sheet and Sea-level System Model (ISSM) and can only be run in an environment where ISSM is installed.

## Repository structure

|  Folder  |  Description 
|  ---  |  ---  |
|  `initialization/`  |  Model initialization and supporting files.  |
|  `idealized_simulations/`  |  Idealized sea level perturbation experiments.  | 
|  `fullycoupled_simulations/`  |  Contains the fully coupled experiments in which sea-level changes from ISSM-SESAW are applied to different spatial extent of the ice-sheet domain.  |
|  `functions/`  |  Contains functions required by the simulation scripts.  |

Each simulation folder also contains the supporting model files, parameter files, scripts, and subdirectories required to run the corresponding experiments.

### Figure scripts

|  File  |  Description  |
|  ---  |  ---  |
|  `Figure2.m`  |  Reproduces Figure 2 from the completed idealized simulations.  |
|  `Figure3.m`  |  Reproduces Figure 3 from the completed fully coupled simulations.  |

## Requirements

- MATLAB
- Ice-sheet and Sea-level System Model (ISSM)
- A functioning ISSM-SESAW setup, including all external packages and dependencies required by ISSM-SESAW.
- [`cbrewer2`](https://github.com/scottclowe/cbrewer2) for figure colormaps
## Reproducing figures

1. Run `initialization/initializemodel.m` before running any simulations.
2. Run `idealized_simulations/idealized_sim.m`. Set `cf` and `gl` to the perturbation values used in the manuscript. For example, `cf = 0` and `gl = 5` applies a 5m sea-level perturbation at the grounding line.
3. Run `Figure2.m` to reproduce Figure 2 in the manuscript.
4. Run `fullycoupled_simulations/cfsimulations.m` for perturbations applied at the calving front and floating ice, `fullycoupled_simulations/glsimulations.m` for perturbations applied at the grounding line only, and `fullycoupled_simulations/fullsimulations.m` for perturbations applied across the full domain.
5. Run `Figure3.m` to reproduce Figure 3(panels a-e) in the manuscript.
