

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


# Check with Erik: are these messages an issue? Seems fairly sporadic...

# m$fit$output # check with Erik if these warnings are an issue? 
# [74] "Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:"                                                     
# [75] "Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/tmp/RtmplKeyyW/model-3a4a063d19797d.stan', line 204, column 2 to column 36)"
# [76] "If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,"                             
# [77] "but if this warning occurs often then your model may be either severely ill-conditioned or misspecified." 


# Additional coevolve figures

m <- tar_read(coev_Q_noTransform)

coev_plot_flowfield(m, var1='length.mean_F',
                    var2='Q',
                    limits=c(-5,5))

# coev_plot_pred_series(m)

coev_plot_pred_series(m,
                      eta_anc = list(length.mean_F=-0.5,
                                     Q=0))

coev_plot_pred_series(m,
                      eta_anc = list(length.mean_F=2,
                                     Q=-2))

coev_plot_pred_series(m,
                      eta_anc = list(length.mean_F=-0.5,
                                     Q=0.5))

coev_plot_pred_series(m,
                      eta_anc = list(length.mean_F=1,
                                     Q=1))

coev_plot_pred_series(m,
                      eta_anc = list(length.mean_F=-1,
                                     Q=-1))


coev_plot_pred_series(m,
                      eta_anc = list(length.mean_F=0,
                                     Q=0))




coev_plot_flowfield(tar_read(coev_Q), var1='log_length_F',
                    var2='Q',
                    limits=c(-2,2))





