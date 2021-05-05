source('http://renozao.github.io/repotools/install.R') #1.10.1
library(repotools)
repotools::install.pkgs('NMF', force = TRUE, quiet = TRUE) #>=0.23.0

library(devtools)
devtools::install_github("jokergoo/circlize@feccef4f5f25f0ad389ccb6d3e36a2e6813ad7a8")
