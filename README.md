# c2c_sim
c2c_sim is a pipeline for simulating cell-cell interaction dynamics across contexts (e.g., time). See the Tutorial to get started

## Install 
1. clone the repo and create the environment using the .yml file

```
git clone git@github.com:hmbaghdassarian/c2c_sim.git
cd c2c_sim
conda env create -n c2c_sim --file=c2c_sim.yml
conda activate c2c_sim
```

2. with your environment activated, install [StabEco](https://rdrr.io/github/keepsimpler/StabEco/) following the installation instructions. First try installing the specific version package was written from using devtools, otherwise can try installing with either remotes or devtools:

>Option 1:
```
devtools::install_github("keepsimpler/StabEco@fe71129dda26db1a0578557808960460e09201f6")

```

>Option 2: 
```
remotes::install_github("keepsimpler/StabEco")
```

>Option 3:
```
devtools::install_github("keepsimpler/StabEco")
```
3. Explore directory shows analyses related to https://doi.org/10.1101/2021.09.20.461129