

library(targets)

###############
# Run analysis
# tar_make() 
###############

######################################
# Run analysis with parallel computing
tar_make_clustermq(workers=16)


# Visualize analysis
# tar_visnetwork(targets_only = TRUE, label = 'time') # visualizes analysis pipeline

# Examine object from analysis
# tar_read(printsByPopulation, branches=1)

# Diagnostics 
# View(tar_meta()) # useful tool for diagnostics
# View(tar_meta(targets_only = TRUE)) # simplified



