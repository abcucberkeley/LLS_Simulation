# LLS_Simulation
Lattice light sheet simulation and stripe pattern deconvolution simulation code for paper

**Characterization, Comparison, and Optimization of Lattice Light Sheets**, Gaoxiang Liu*, Xiongtao Ruan*, Daniel E. Milkie, Frederik Görlitz, Matthew Mueller, Wilmene Hercule, Alison Kililea, Eric Betzig#, Srigokul Upadhyayula# bioRxiv 2022.07.30.502108; doi: https://doi.org/10.1101/2022.07.30.502108

## Usage
The code is tested with Matlab R2022b. In matlab, run **simulation.m** and select a light sheet for simulation. 

The simulation of lattice light sheet generation is self-contained within the repo. To run the stripe pattern deconvolution simulation, the theoretical 3D detection PSF file can be downloaded from [this link](https://www.dropbox.com/s/a5gaz1tdj7g6ozm/Det_PSF_OTF_3D_510_NA1p0_px_38p346nm_RichardsWolf.mat?dl=0). After downloading the file, create a folder with name "RW_3d" under the simulation result root directory and put the file in this folder. Then the stripe pattern deconvolution simulation should work smoothly. 

If there are questions or issues, please feel free to contact Eric Betzig (betzige@janelia.hhmi.org) or Srigokul Upadhyayula (sup@berkeley.edu). 
