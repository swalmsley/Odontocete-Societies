
library(targets)

# Source
tar_source('R') # will do all in 'R' folder

# Seed
tar_option_set(seed = 1234)

# Configuration - reduce uploads to cloud for improved efficiency
tar_config_set(seconds_meta_append = 15,
               seconds_meta_upload = 15,
               seconds_reporter = 0.5)

# Variables
suppressMessages(set_cmdstan_path(path='C:/Users/sjfwa/AppData/Local/R/cmdstan-2.33.1')) # cmdstan path for local machine
# cmdstanr::set_cmdstan_path('/home/sjfw/.cmdstan/cmdstan-2.36.0') # path for supercomputer


# ############################
# # Parallel computing - Canada Supercomputer (Comment out to run locally)
# options(clustermq.scheduler = 'slurm', clustermq.template='SLURM.tmpl')
# # Set targets options with SLURM resources
# tar_option_set(
#   resources = tar_resources(
#     clustermq = tar_resources_clustermq(
#       template = list(
#         cores = 4,            # Number of cores per task
#         memory = 8192,        # Memory in MB (e.g., 8192 for 8 GB)
#         time = '7-00:00:00'   # Time in D-HH:MM:SS
#       )
#     )
#   )
# )


# Set color theme
proj_color <- 'A'

# Targets
list(
  
  ####################################
  # Prepare data
  ###################################
  
  # load species list
  tar_target(species_file, './Input/speciesList.csv', format='file'),
  tar_target(species, read_data(species_file)),
  
  # load database of social network traits
  tar_target(database_file, 'Input/Literature_Review/5-Database.xlsx', format = 'file'),
  tar_target(database, read_data_excel(database_file)),  
  
  # load lifespan data from Ellis et al. 2023
  # https://www.biorxiv.org/content/10.1101/2023.02.22.529527v1.abstract
  tar_target(ellisLifespan_file, 'Input/Odontocete_Traits/Ellis_Lifespan/Ellis_Lifespan.xlsx', format = 'file'),
  tar_target(ellisLifespan_raw, read_data_excel(ellisLifespan_file)),
  tar_target(ellisLifespan, prep_ellisLifespanData(species, ellisLifespan_raw)),
  
  # prepare life history data from marinelifehistdata package
  # https://github.com/samellisq/marinelifehistdata
  tar_target(ellisLifeHistory, extract_Ellis_data(species)),
  # identify species to exclude from SSD comparison
  tar_target(exclude_from_SSD, SSD_exclusions(species, ellisLifeHistory)),

  # merge data sources
  tar_target(data_a, merge(species, database, all.x=TRUE)),
  tar_target(data_b, merge(data_a, ellisLifespan, by='Species', all.x=TRUE)),
  tar_target(data_c, merge(data_b, ellisLifeHistory, by='Species', all.x=TRUE)),
  
  # clean up final database and add data transformations
  tar_target(data, finalize_data(data_c, exclude_from_SSD)),
  
  # make versions of data incorporating exclusions specific to Q and S
  tar_target(data_Q, data[AssociationIndex %in% c('HWI','SRI','SAI','GAI (with HWI)') & !(Exclusions %in% c('Both','Q')) & SexFocus=='Assume mixed',,]),
  tar_target(data_S, data[S_method=='Likelihood' & AssociationIndex %in% c('HWI','SRI','SAI','GAI (with HWI)') & !(Exclusions %in% c('Both','S')) & SexFocus=='Assume mixed',,]),
  
  # create species-specific dataframe with IVSO at species level
  tar_target(data_IVSO, add_IVSO(data)),

  # tar_target(A_list, rep(list(A=A_pruned), 5),)
  tar_target(A_list, lapply(1:10, function(x) list(A = A))),  # Corrected list structure
  

  ####################################
  # Prepare phylogeny
  ####################################
  
  # load gene key
  tar_target(gene_key_file, 'Input/McGowen_Trees/Gene_Key.xlsx', format = 'file'),
  tar_target(gene_key, read_data_excel(gene_key_file)),  
  
  # load tree
  tar_target(tree, prepare_mcgowen_trees_noSister(species, gene_key, 'Input/McGowen_Trees/FigTree_parts_6_mcmctree_AR.tre')),
  tar_target(tree_with_additions, prepare_mcgowen_trees(species, gene_key, 'Input/McGowen_Trees/FigTree_parts_6_mcmctree_AR.tre')), # can use for visualization if desired, not necessary

  # check for species name mismatch
  tar_target(missing_from_tree, identify_mismatch(species, tree, database)),
  # Currently 19 species missing from tree
  
  # create variance-covariance matrix 
  tar_target(A, vcv(tree, corr=TRUE)),
  ###### Do I need to subset the variance-covariance matrix to included species or does it not matter?
  
  
  ####################################
  # Phylogenetic signal
  ####################################

  # Fit models to estimate phylogenetic signal

  tar_target(psQ, ps_model('Q', data_Q, A, repeatObs=FALSE)),
  tar_target(psQ_repeat, ps_model('Q', data_Q, A, repeatObs=TRUE)),

  tar_target(psS, ps_model('S', data_S, A, repeatObs=FALSE)),
  tar_target(psS_repeat, ps_model('S', data_S, A, repeatObs=TRUE)),

  tar_target(psLifespan, ps_model('lifespan_Post.Mean_F', data, A, repeatObs=FALSE)),

  tar_target(psAgeM, ps_model('age.mat_F', data, A, repeatObs=FALSE)),

  tar_target(psLength, ps_model('log_length_F', data, A, repeatObs=FALSE)),
  tar_target(psLength_noTransform, ps_model('length.mean_F', data, A, repeatObs=FALSE)),

  tar_target(psSSD, ps_model('SSD', data, A, repeatObs=FALSE)),


  ####################################
  # Single-trait phylogenetic models
  ####################################

  ## Lifespan
  tar_target(m1q, brm(Q ~ scale(lifespan_Post.Mean_F) + (1|AssociationIndex) + (1|Species) + (1|Population) + (1|gr(phylo, cov = A)),
                     family = 'Beta', data = data_Q, data2 = list(A=A),
                     prior = c(prior("normal(0,2)", class = "Intercept"),
                               prior("normal(0,1)", class = "b")),
                     control=list(adapt_delta=0.99, max_treedepth=15),
                     iter = 4000, chains = 4, cores = 4)),
  tar_target(m1s, brm(S ~ scale(lifespan_Post.Mean_F) + (1|AssociationIndex) + (1|Species) + (1|Population) + (1|gr(phylo, cov = A)),
                     family = 'Gamma', data = data_S, data2 = list(A=A),
                     prior = c(prior("normal(0,1)", class = "Intercept"),
                               prior("normal(0,1)", class = "b")),
                     control=list(adapt_delta=0.99, max_treedepth=15),
                     iter = 4000, chains = 4, cores = 4)),
  tar_target(m1q_IVSO, brm(Q_IVSO ~ scale(lifespan_Post.Mean_F) + (1|gr(phylo, cov = A)),
                      family = 'Gamma', data = data_IVSO[numObsQ>1,,], data2 = list(A=A),
                      prior = c(prior("normal(-2,1.5)", class = "Intercept"),
                                prior("normal(0,1)", class = "b"),
                                prior("exponential(1)", class='sd'),
                                prior("gamma(3,0.1)", class = "shape")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),
  tar_target(m1s_IVSO, brm(S_IVSO ~ scale(lifespan_Post.Mean_F) + (1|gr(phylo, cov = A)),
                           family = 'Gamma', data = data_IVSO[numObsS>1,,], data2 = list(A=A),
                           prior = c(prior("normal(-2,1.5)", class = "Intercept"),
                                     prior("normal(0,1)", class = "b"),
                                     prior("exponential(1)", class='sd')),
                           control=list(adapt_delta=0.99, max_treedepth=15),
                           iter = 4000, chains = 4, cores = 4)),



  ## Age at maturity
  tar_target(m2q, brm(Q ~ scale(age.mat_F) + (1|AssociationIndex) + (1|Species) + (1|Population) + (1|gr(phylo, cov = A)),
                      family = 'Beta', data = data_Q, data2 = list(A=A),
                      prior = c(prior("normal(0,2)", class = "Intercept"),
                                prior("normal(0,1)", class = "b")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),
  tar_target(m2s, brm(S ~ scale(age.mat_F) + (1|AssociationIndex) + (1|Species) + (1|Population) + (1|gr(phylo, cov = A)),
                      family = 'Gamma', data = data_S, data2 = list(A=A),
                      prior = c(prior("normal(0,1)", class = "Intercept"),
                                prior("normal(0,1)", class = "b")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),
  tar_target(m2q_IVSO, brm(Q_IVSO ~ scale(age.mat_F) + (1|gr(phylo, cov = A)),
                      family = 'Gamma', data = data_IVSO[numObsQ>1], data2 = list(A=A),
                      prior = c(prior("normal(-2,1.5)", class = "Intercept"),
                                prior("normal(0,1)", class = "b"),
                                prior("exponential(1)", class='sd'),
                                prior("gamma(3,0.1)", class = "shape")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),
  tar_target(m2s_IVSO, brm(S_IVSO ~ scale(age.mat_F) + (1|gr(phylo, cov = A)),
                      family = 'Gamma', data = data_IVSO[numObsS>1], data2 = list(A=A),
                      prior = c(prior("normal(-2,1.5)", class = "Intercept"),
                                prior("normal(0,1)", class = "b"),
                                prior("exponential(1)", class='sd'),
                                prior("gamma(3,0.1)", class = "shape")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),


  # Body length
  tar_target(m3q, brm(Q ~ scale(log_length_F) + (1|AssociationIndex) + (1|Species) + (1|Population) + (1|gr(phylo, cov = A)),
                      family = 'Beta', data = data_Q, data2 = list(A=A),
                      prior = c(prior("normal(0,2)", class = "Intercept"),
                                prior("normal(0,1)", class = "b")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),
  tar_target(m3s, brm(S ~ scale(log_length_F) + (1|AssociationIndex) + (1|Species) + (1|Population) + (1|gr(phylo, cov = A)),
                      family = 'Gamma', data = data_S, data2 = list(A=A),
                      prior = c(prior("normal(0,1)", class = "Intercept"),
                                prior("normal(0,1)", class = "b")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),
  tar_target(m3q_noTransform, brm(Q ~ scale(length.mean_F) + (1|AssociationIndex) + (1|Species) + (1|Population) + (1|gr(phylo, cov = A)),
                      family = 'Beta', data = data_Q, data2 = list(A=A),
                      prior = c(prior("normal(0,2)", class = "Intercept"),
                                prior("normal(0,1)", class = "b")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),
  tar_target(m3s_noTransform, brm(S ~ scale(length.mean_F) + (1|AssociationIndex) + (1|Species) + (1|Population) + (1|gr(phylo, cov = A)),
                      family = 'Gamma', data = data_S, data2 = list(A=A),
                      prior = c(prior("normal(0,1)", class = "Intercept"),
                                prior("normal(0,1)", class = "b")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),
  tar_target(m3q_IVSO, brm(Q_IVSO ~ scale(log_length_F) + (1|gr(phylo, cov = A)),
                      family = 'Gamma', data = data_IVSO[numObsQ>1,,], data2 = list(A=A),
                      prior = c(prior("normal(-2,1.5)", class = "Intercept"),
                                prior("normal(0,1)", class = "b"),
                                prior("exponential(1)", class='sd'),
                                prior("gamma(3,0.1)", class = "shape")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),
  tar_target(m3s_IVSO, brm(S_IVSO ~ scale(log_length_F) + (1|gr(phylo, cov = A)),
                      family = 'Gamma', data = data_IVSO[numObsS>1], data2 = list(A=A),
                      prior = c(prior("normal(-2,1.5)", class = "Intercept"),
                                prior("normal(0,1)", class = "b"),
                                prior("exponential(1)", class='sd'),
                                prior("gamma(3,0.1)", class = "shape")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),



  # Sexual size dimorphism
  tar_target(m4q, brm(Q ~ scale(SSD) + (1|AssociationIndex) + (1|Species) + (1|Population) + (1|gr(phylo, cov = A)),
                      family = 'Beta', data = data_Q, data2 = list(A=A),
                      prior = c(prior("normal(0,2)", class = "Intercept"),
                                prior("normal(0,1)", class = "b")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),
  tar_target(m4s, brm(S ~ scale(SSD) + (1|AssociationIndex) + (1|Species) + (1|Population) + (1|gr(phylo, cov = A)),
                      family = 'Gamma', data = data_S, data2 = list(A=A),
                      prior = c(prior("normal(0,1)", class = "Intercept"),
                                prior("normal(0,1)", class = "b")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),
  tar_target(m4q_IVSO, brm(Q_IVSO ~ scale(SSD) + (1|gr(phylo, cov = A)),
                      family = 'Gamma', data = data_IVSO[numObsQ>1], data2 = list(A=A),
                      prior = c(prior("normal(-2,1.5)", class = "Intercept"),
                                prior("normal(0,1)", class = "b"),
                                prior("exponential(1)", class='sd'),
                                prior("gamma(3,0.1)", class = "shape")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),
  tar_target(m4s_IVSO, brm(S_IVSO ~ scale(SSD) + (1|gr(phylo, cov = A)),
                      family = 'Gamma', data = data_IVSO[numObsS>1], data2 = list(A=A),
                      prior = c(prior("normal(-2,1.5)", class = "Intercept"),
                                prior("normal(0,1)", class = "b"),
                                prior("exponential(1)", class='sd'),
                                prior("gamma(3,0.1)", class = "shape")),
                      control=list(adapt_delta=0.99, max_treedepth=15),
                      iter = 4000, chains = 4, cores = 4)),


  ####################################
  # Coevolutionary models
  ####################################

  tar_target(data_coev_Q, create_coev_data(data, 'Q')),
  tar_target(data_coev_S, create_coev_data(data, 'S')),

  tar_target(coev_Q, fit_and_save_coev_model(data = data_coev_Q, var1 = 'Q', var2 = 'log_length_F', tree)),
  tar_target(coev_S, fit_and_save_coev_model(data = data_coev_S, var1 = 'S', var2 = 'log_length_F', tree)),
  
  tar_target(coev_Q_noTransform, fit_and_save_coev_model(data = data_coev_Q, var1 = 'Q', var2 = 'length.mean_F', tree)),
  tar_target(coev_S_noTransform, fit_and_save_coev_model(data = data_coev_S, var1 = 'S', var2 = 'length.mean_F', tree)),

  tar_target(coev_Q_lifespan, fit_and_save_coev_model(data = data_coev_Q, var1 = 'Q', var2 = 'lifespan_Post.Mean_F', tree)),
  tar_target(coev_S_lifespan, fit_and_save_coev_model(data = data_coev_S, var1 = 'S', var2 = 'lifespan_Post.Mean_F', tree)),

  tar_target(coev_Q_ageMat, fit_and_save_coev_model(data = data_coev_Q, var1 = 'Q', var2 = 'age.mat_F', tree)),
  tar_target(coev_S_ageMat, fit_and_save_coev_model(data = data_coev_S, var1 = 'S', var2 = 'age.mat_F', tree)),

  tar_target(coev_Q_SSD, fit_and_save_coev_model(data = data_coev_Q, var1 = 'Q', var2 = 'SSD', tree)),
  tar_target(coev_S_SSD, fit_and_save_coev_model(data = data_coev_S, var1 = 'S', var2 = 'SSD', tree)),
  


  ####################################
  # Figures
  ####################################
  
  # F1 - Map of study sites
  tar_target(Figure1, save_figure('./Manuscript/Figures/Figure1.png',w=12,h=7,
                                  create_map(data))),


  # F2 - Relationship between Q and S
  tar_target(Figure2, save_figure('./Manuscript/Figures/Figure2.png',w=5.5,h=5.5,
                                    (plot_qs(data)))),



  # F3 - Version 2: Univariate effects and ancestral states
  tar_target(fg_layout, layout <- "
             ABFFFFGH
             CDFFFFIJ
             EEFFFFKK
             "),
  tar_target(Figure3, save_figure('./Manuscript/Figures/Figure3.png',w=24,h=12,
                                  (plot_univariate(m1q, 'lifespan_Post.Mean_F', 'Q', 'Lifespan', 'Modularity', 1, 2) +
                                     plot_univariate(m2q, 'age.mat_F', 'Q', 'Age at maturity', 'Modularity', 1, 3) +
                                     plot_univariate(m3q, 'log_length_F', 'Q', expression('Body length ' * (log)), 'Modularity', 1, 4) +
                                     plot_univariate(m4q, 'SSD', 'Q', 'Sexual size dimorphism', 'Modularity', 1, 5) +
                                     ancestral_state(m3q, ancestral_length_raw=250, colnum=6, var='Q', label_position=c(0.75, 2), c(0,1)) +
                                     plot_tree_two_adjacent_traits(tree, data, 'Q', 'S', flip=FALSE, just_labels=FALSE) +
                                     plot_univariate(m1s, 'lifespan_Post.Mean_F', 'S', 'Lifespan', 'Social differentiation', 2.27, 2) +
                                     plot_univariate(m2s, 'age.mat_F', 'S', 'Age at maturity', 'Social differentiation', 2.27, 3) +
                                     plot_univariate(m3s, 'log_length_F', 'S', expression('Body length ' * (log)), 'Social differentiation', 2.27, 4) +
                                     plot_univariate(m4s, 'SSD', 'S', 'Sexual size dimorphism', 'Social differentiation', 2.27, 5) +
                                     ancestral_state(m3s, ancestral_length_raw=250, colnum=4, var='S', label_position=c(2, 1), c(0,3)) +
                                     plot_layout(design = fg_layout)) + plot_annotation(tag_levels='A'))),

  # F4 - Dynamic effects
  tar_target(Figure4, save_figure('./Manuscript/Figures/Figure4.png',w=14,h=7,
                                  (trait_change_plot(coev_Q_noTransform) |
                                     custom_coev_plot_flowfield(coev_Q_noTransform, 'length.mean_F', 'Q', nullclines=FALSE, limits=c(-5,5), var1_lab='Body length', var2_lab='Modularity')) +
                                     plot_annotation(tag_levels='A'))),
  
  # S1 - Phylogenetic signal
  tar_target(FigureS1, save_figure('./Manuscript/Figures/FigureS1.png',w=3,h=5,
                                   plot_phylogenetic_signal_ridges(psQ, psS, psLifespan, psAgeM, psLength, psSSD))),

  # S2 - Dynamic effects (all)
  tar_target(FigureS2, save_figure('./Manuscript/Figures/FigureS2.png',w=25,h=12,
                                   (trait_change_plot(coev_Q) / trait_change_plot(coev_S) |
                                      trait_change_plot(coev_Q_lifespan) / trait_change_plot(coev_S_lifespan) |
                                      trait_change_plot(coev_Q_ageMat) / trait_change_plot(coev_S_ageMat) |
                                      trait_change_plot(coev_Q_SSD) / trait_change_plot(coev_S_SSD))
                                      + plot_annotation(tag_levels='A'))),

  # Visual summaries of GDPMs
  tar_target(Figure_GDPM_1, save_figure('./Manuscript/Figures/GDPM_summaries/Figure_GDPM_1.png',w=11,h=11, (GDPM_visual_summary(coev_Q, 'Q', 'log_length_F')))),
  tar_target(Figure_GDPM_2, save_figure('./Manuscript/Figures/GDPM_summaries/Figure_GDPM_2.png',w=11,h=11, (GDPM_visual_summary(coev_S, 'S', 'log_length_F')))),

  tar_target(Figure_GDPM_1a, save_figure('./Manuscript/Figures/GDPM_summaries/Figure_GDPM_1a.png',w=11,h=11, (GDPM_visual_summary(coev_Q_noTransform, 'Q', 'length.mean_F')))),
  tar_target(Figure_GDPM_2a, save_figure('./Manuscript/Figures/GDPM_summaries/Figure_GDPM_2a.png',w=11,h=11, (GDPM_visual_summary(coev_S_noTransform, 'S', 'length.mean_F')))),

  tar_target(Figure_GDPM_3, save_figure('./Manuscript/Figures/GDPM_summaries/Figure_GDPM_3.png',w=11,h=11, (GDPM_visual_summary(coev_Q_lifespan, 'Q', 'lifespan_Post.Mean_F')))),
  tar_target(Figure_GDPM_4, save_figure('./Manuscript/Figures/GDPM_summaries/Figure_GDPM_4.png',w=11,h=11, (GDPM_visual_summary(coev_S_lifespan, 'S', 'lifespan_Post.Mean_F')))),

  tar_target(Figure_GDPM_5, save_figure('./Manuscript/Figures/GDPM_summaries/Figure_GDPM_5.png',w=11,h=11, (GDPM_visual_summary(coev_Q_ageMat, 'Q', 'age.mat_F')))),
  tar_target(Figure_GDPM_6, save_figure('./Manuscript/Figures/GDPM_summaries/Figure_GDPM_6.png',w=11,h=11, (GDPM_visual_summary(coev_S_ageMat, 'S', 'age.mat_F')))),

  tar_target(Figure_GDPM_7, save_figure('./Manuscript/Figures/GDPM_summaries/Figure_GDPM_7.png',w=11,h=11, (GDPM_visual_summary(coev_Q_SSD, 'Q', 'SSD')))),
  tar_target(Figure_GDPM_8, save_figure('./Manuscript/Figures/GDPM_summaries/Figure_GDPM_8.png',w=11,h=11, (GDPM_visual_summary(coev_S_SSD, 'S', 'SSD')))),


  ####################################
  # Write manuscript
  ####################################

  tar_quarto(
    supplement,
    file.path('Manuscript','Supplement_OdontoceteSocieties.qmd')),

  tar_quarto(
    paper,
    file.path('Manuscript','MS_OdontoceteSocieties.qmd'))

)











