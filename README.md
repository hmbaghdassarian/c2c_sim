# c2c_sim
c2c_sim is a pipeline for simulating cell-cell interaction dynamics across contexts (e.g., time). See the Tutorial to get started

## Install 
1. clone the repo and create the environment using the .yml file

```
git clone git@github.com:hmbaghdassarian/c2c_sim.git
cd c2c_sim
conda create -n c2c_sim --file=c2c_sim.yml
conda activate c2c_sim
```

2. with your environment activated, install [StabEco](https://rdrr.io/github/keepsimpler/StabEco/) following the installation instructions. If install with the "remotes" package does not work, try to install with devtools in R:

>Option 1:
```
remotes::install_github("keepsimpler/StabEco")
```

>Option 2: 
```
devtools::install_github("YosefLab/SymSim@a9e24981ba3f68cc35b93ae34402e25bfa1971f9")
```

>Option 3:
```
devtools::install_github("YosefLab/SymSim@v0.0.0.9000")
```
3. To explore tensor outputs of the simulation, also install [cell2cell](https://github.com/earmingol/cell2cell)