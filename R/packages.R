suppressPackageStartupMessages({
  
  library(targets)
  library(tarchetypes)
  library(quarto)
  library(patchwork)
  library(rempsyc)
  library(ggplot2)
  library(data.table)
  library(crew)
  
  library(cmdstanr)
  library(ape)
  library(rethinking)
  library(bayesplot)
  library(phangorn)
  library(tarchetypes)
  library(ggplot2)
  library(ggtree)
  library(brms)
  library(stringr)
  library(igraph)
  library(tidybayes)
  library(ggrepel)
  library(ggdist)
  library(readxl)
  library(viridis)
  library(testthat)
  library(mice)
  library(phytools)
  library(sjPlot)
  library(stringi)
  library(tidyr)

  library(TreeTools)
  library(Rphylopars)
  library(rphylopic)
  library(lavaan)
  library(coevolve)
  library(phylobase)
  library(bayestestR)
  library(ggraph)
  library(ggridges)
  library(blavaan)
  library(dplyr)
  library(tibble)
  library(ggplotify)
  library(cowplot)
  
  library(marinelifehistdata)
  
  
  # Comment these out for supercomputer
  # library(rnaturalearth)
  # library(rnaturalearthdata)
  # library(sf)
  # library(stantargets)
  # library(rbbt)
  
  # If packages need to be installed 
  
  # install.packages(c(
  #   "targets", 
  #   "tarchetypes", 
  #   "rbbt", 
  #   "quarto", 
  #   "patchwork", 
  #   "rempsyc", 
  #   "ggplot2", 
  #   "data.table", 
  #   "crew", 
  #   "cmdstanr", 
  #   "ape", 
  #   "rethinking", 
  #   "bayesplot", 
  #   # "stantargets", 
  #   "phangorn", 
  #   "ggtree", 
  #   "brms", 
  #   "stringr", 
  #   "igraph", 
  #   "tidybayes", 
  #   "ggrepel", 
  #   "ggdist", 
  #   "readxl", 
  #   "viridis", 
  #   "testthat", 
  #   "mice", 
  #   "phytools", 
  #   "sjPlot", 
  #   "stringi", 
  #   "tidyr", 
  #   "rnaturalearth", 
  #   "rnaturalearthdata", 
  #   "sf", 
  #   "TreeTools", 
  #   "Rphylopars", 
  #   "rphylopic", 
  #   "lavaan", 
  #   "coevolve", 
  #   "phylobase", 
  #   "bayestestR", 
  #   "ggraph", 
  #   "ggridges", 
  #   "blavaan", 
  #   "dplyr", 
  #   "tibble", 
  #   "marinelifehistdata"
  # ))
  # install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
  # also devtools for rethinking
  # devtools::install_github("rmcelreath/rethinking")
  
  # if (!requireNamespace("BiocManager", quietly = TRUE))
  #   install.packages("BiocManager")
  # 
  # BiocManager::install("ggtree")
  # devtools::install_github("ropensci/rnaturalearth")
  # devtools::install_github("ropensci/rnaturalearthdata")
  # devtools::install_github("ScottClaessens/coevolve")
  # devtools::install_github("samellisq/marinelifehistdata")
  


})

print('Cleared packages')

