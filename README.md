# LLS_Simulation
Lattice light sheet simulation and stripe pattern deconvolution simulation code for paper

**Characterization, Comparison, and Optimization of Lattice Light Sheets**, Gaoxiang Liu*, Xiongtao Ruan*, Daniel E. Milkie, Frederik Görlitz, Matthew Mueller, Wilmene Hercule, Alison Kililea, Eric Betzig#, Srigokul Upadhyayula# bioRxiv 2022.07.30.502108; doi: https://doi.org/10.1101/2022.07.30.502108

## Usage

The code can be cloned from the source code with command:

```bash
git clone --recurse-submodules https://github.com/abcucberkeley/LLS_Simulation.git
```

or downloaded from https://github.com/abcucberkeley/LLS_Simulation/releases/download/v1.0.0/LLS_Simulation.tgz


The code is tested with Matlab R2022b. In matlab, run **simulation.m** and select a light sheet for simulation. 

The simulation of lattice light sheet generation is self-contained within the repo. To run the stripe pattern deconvolution simulation, the theoretical 3D detection PSF file can be downloaded from [this link](https://www.dropbox.com/s/a5gaz1tdj7g6ozm/Det_PSF_OTF_3D_510_NA1p0_px_38p346nm_RichardsWolf.mat?dl=0). After downloading the file, put it to the folder "RW_PSFs" under the code directory. Then the stripe pattern deconvolution simulation should work smoothly. 

If there are questions or issues, please feel free to contact Eric Betzig (betzige@janelia.hhmi.org) or Srigokul Upadhyayula (sup@berkeley.edu). 
