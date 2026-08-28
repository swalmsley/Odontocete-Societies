
# Custom functions to be called in '_targets.R'


###########################################################################
# Basic functions
###########################################################################


# Read data ---------------------------------------------------------------
read_data <- function(path) {
  d <- data.table(read.csv(path))
  return(d)
}


# Read data (excel) -------------------------------------------------------
read_data_excel <- function(path) {
  d <- data.table(read_xlsx(path))
  return(d)
}



# Save plot ---------------------------------------------------------------
save_figure <- function(path, w, h, call) {
  png(path, width=w, height=h, units='in', res=1000)
  print(call)
  dev.off()
}



# In-line effect ----------------------------------------------------------
inline <- function(fit, var, p, CI) {
  
  est <- summary(fit,prob=p)$fixed[paste(var),'Estimate']
  low <- summary(fit,prob=p)$fixed[paste(var),'l-90% CI']
  high <- summary(fit,prob=p)$fixed[paste(var),'u-90% CI']
  
  # report credible interval
  if (CI) output <- paste(paste(format_number(est), ', ', sep = ''), 'CI ', paste(format_number(low), format_number(high), sep = ' -- '), sep = '')
  
  # report probability of directional effect
  if (!CI) output <- paste(paste(format_number(est), ', ', sep = ''), p_dir <- paste('pd = ', round(as.numeric(p_direction(fit, parameters = paste(var)))*100,1), '%',sep=''), sep = '')
  
  return(output)
  
}
# fit <- tar_read(m1q)
# var <- 'scalelifespan_Post.Mean_F'
# p = 0.9
# CI = FALSE



# In-line effect given samples --------------------------------------------
inline_simple <- function(samples, p, CI) {
  
  est <- mean(samples)
  low <- PI(samples, prob=p)[1]
  high <- PI(samples, prob=p)[2]
  
  # report credible interval
  if (CI) output <- paste(paste(format_number(est), ', ', sep = ''), 'CI ', paste(format_number(low), format_number(high), sep = ' -- '), sep = '')
  
  # report probability of directional effect
  # if (!CI) output <- paste(paste(format_number(est), ', ', sep = ''), p_dir <- paste('pd = ', round(as.numeric(p_direction(fit, parameters = paste(var)))*100,1), '%',sep=''), sep = '')
  
  return(output)
  
}
# samples <- as_draws_df(tar_read(psQ_repeat))$sd_phylo__Intercept
# var <- 'numSamplingPeriodsByYear'
# p = 0.9
# CI = FALSE




# Helper for inline effect ------------------------------------------------
format_number <- function(num) {
  if (abs(num) < 0.01) {
    
    # Convert to scientific notation with custom formatting
    formatted_number <- format(num, scientific = TRUE, trim = TRUE, digits = 2)
    # Extract mantissa and exponent
    parts <- strsplit(formatted_number, "e", fixed = TRUE)[[1]]
    # Print in the desired format
    out <- capture.output(cat(sprintf("%.2f x 10^%+d\n", as.numeric(parts[1]), as.numeric(parts[2]))))
    return(out)
    
  } else {
    return(round(num, 2))
  }
}



###########################################################################
# Project-specific functions
###########################################################################


# # Custom helper to deal with cmdstanr issues on supercomputer -------------
# fit_and_save_coev_model <- function(data, vars, tree) {
#   
#   model_name <- paste(vars, collapse = '-')
#   
#   save_coevfit(coev_fit(data = data,
#                         variables = setNames(as.list(rep("normal", length(vars))), vars), # currently assuming normal
#                         id = 'phylo',
#                         tree = tree,
#                         complete_cases = FALSE,
#                         estimate_residual = TRUE,
#                         chains=4,
#                         parallel_chains = 4,
#                         iter_warmup = 2000,
#                         iter_sampling = 2000,
#                         adapt_delta = 0.99,
#                         max_treedepth = 15,
#                         refresh = 0,
#                         seed = 1,
#                         output_dir = file.path(getwd(), "Output")), file = paste(file.path(getwd(), "Output"), '/', model_name, '.RDS', sep=''))
#   
#   readRDS(paste(file.path(getwd(), "Output"), '/', model_name, '.RDS', sep=''))
#   
# }
# # data <- tar_read(data_coev_Q)
# # vars <- c('Q', 'log_length_F')
# # tree <- tar_read(tree)


# Custom helper to deal with cmdstanr issues on supercomputer -------------
fit_and_save_coev_model <- function(data, vars, tree) {
  
  model_name <- paste(vars, collapse = '-')
  
  m <- coev_fit(data = data,
                        variables = setNames(as.list(rep("normal", length(vars))), vars), # currently assuming normal
                        id = 'phylo',
                        tree = tree,
                        complete_cases = FALSE,
                        estimate_residual = TRUE,
                        chains=4,
                        parallel_chains = 4,
                        iter_warmup = 2000,
                        iter_sampling = 2000,
                        adapt_delta = 0.99,
                        max_treedepth = 15,
                        refresh = 0,
                        seed = 1)
  
  return(m)

}
# data <- tar_read(data_coev_Q)
# vars <- c('Q', 'log_length_F')
# tree <- tar_read(tree)






# Create data for coevolutionary model ------------------------------------
create_coev_data <- function(data, measure) {
  
  # We want to retain all species for coevolutionary models
  # Instead of excluding observations, set to NA
  data[!((SexFocus)=='Assume mixed'), (measure):=NA ,] 
  data[!(AssociationIndex %in% c('HWI','SRI','SAI','GAI (with HWI)')), (measure):=NA,]
  data[Exclusions %in% c('Both',paste(measure)), (measure):=NA,]
  
  if (measure=='S') (data[!(S_method=='Likelihood'),(measure):=NA,])
  
  data <- data[!is.na(get(measure)) | !duplicated(Species), ] # for species with only NAs, just retain one observation of life history data
  
  data <- data[!(phylo %in% c('Platanista gangetica', 'Platanista minor')),,] # cut out species not located in tree
  
  return(data)
  
}
# data <- tar_read(data)
# measure <- 'S'



# Fit model for estimating phylogenetic signal ----------------------------
ps_model <- function(focal_trait, data, A, repeatObs) {
  
  # If response variable is a life history trait, restrict to one measure per species
  if (!(focal_trait %in% c('Q','S'))) (data <- unique(data[,.(phylo,Species,get(focal_trait))])) 
  if (!(focal_trait %in% c('Q','S'))) (colnames(data)[3] <- focal_trait) # rename focal trait appropriately 
   
  # Dynamically build the formula with focal_trait
  if (repeatObs==FALSE) (formula <- bf(paste0("scale(", focal_trait, ") ~ (1|gr(phylo, cov = A))")))
  if (repeatObs==TRUE) (formula <- bf(paste0("scale(", focal_trait, ") ~ (1|Species) + (1|gr(phylo, cov = A))")))
  
  # Slice out species not located on tree
  data <- data[!(phylo %in% c('Platanista gangetica', 'Platanista minor')),,]
  
  # fit model
  fit <- brm(formula,
             family = 'Gaussian',
             data = data, 
             data2 = list(A=A),
             prior = c(prior("normal(0,1)", class = "Intercept"),
                       prior("exponential(1)", class = "sigma"),
                       prior("exponential(1)", class = "sd")),
             control=list(adapt_delta=0.99),
             iter = 4000, warmup = 2000, chains = 4, cores = 4)
  
  return(fit)
  
}
# focal_trait <- 'SSD'
# data <- tar_read(data)
# A <- tar_read(A)
# repeatObs <- FALSE

# focal_trait <- 'log_length_F'
# data <- tar_read(data)
# A <- tar_read(A)
# repeatObs <- FALSE



# Calculate phylogenetic signal -------------------------------------------
phylogenetic_signal <- function (fit, repeat_obs, species_effect) {
  
  # Note: will depend on random effect structure of model in question 
  
  # With species-level random effect (in addition to phylogenetic one)
  if (repeat_obs & (!(species_effect))) (hyp <- paste("sd_phylo__Intercept^2 /", "(sd_phylo__Intercept^2 + sd_Species__Intercept^2 + sigma^2) = 0")) # if interested in phylo over and above species
  if (repeat_obs & species_effect) (hyp <- paste("sd_Species__Intercept^2 /", "(sd_phylo__Intercept^2 + sd_Species__Intercept^2 + sigma^2) = 0")) # if interested in species over and above phylo

  # With only phylogenetic random effect
  if (!(repeat_obs)) (hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0")
  
  # Run hypothesis test
  (hyp <- hypothesis(fit, hyp, class = NULL))
  
  return(data.table(hyp$samples))
  
}
# fit <- tar_read(psQ)
# repeat_obs <- FALSE
# species_effect <- FALSE

# fit <- tar_read(psQ_repeat)
# repeat_obs <- TRUE
# species_effect <- TRUE



# Plot posteriors of phylogenetic signal ----------------------------------
plot_phylogenetic_signal_ridges <- function(Q, S, Lifespan, AgeM, Length, SSD) {
  
  # calculate proportion of variance explained by phylogeny 
  Q_samples <- phylogenetic_signal(Q, repeat_obs=FALSE, species_effect=FALSE)
  Q_samples[, trait:='Q', ]
  
  S_samples <- phylogenetic_signal(S, repeat_obs=FALSE, species_effect=FALSE)
  S_samples[, trait:='S', ]
  
  Lifespan_samples <- phylogenetic_signal(Lifespan, repeat_obs=FALSE, species_effect=FALSE)
  Lifespan_samples[, trait:='Lifespan', ]
  
  AgeM_samples <- phylogenetic_signal(AgeM, repeat_obs=FALSE, species_effect=FALSE)
  AgeM_samples[, trait:='AgeM', ]
  
  Length_samples <- phylogenetic_signal(Length, repeat_obs=FALSE, species_effect=FALSE)
  Length_samples[, trait:='Length', ]
  
  SSD_samples <- phylogenetic_signal(SSD, repeat_obs=FALSE, species_effect=FALSE)
  SSD_samples[, trait:='SSD', ]
  
  # combine all samples into single dataframe
  combined_samples <- rbindlist(list(Q_samples, S_samples, Lifespan_samples, AgeM_samples, Length_samples, SSD_samples), use.names = TRUE, fill = TRUE)
  
  # calculate mean signal by trait
  combined_samples[, Rbar:=round(mean(H1),2), by=trait]
  Rbar_labels <- combined_samples[, .(Rbar = unique(Rbar)), by = trait]
  
  # change order of traits
  combined_samples$trait <- factor(combined_samples$trait, levels=c('SSD', 'Length', 'AgeM', 'Lifespan', 'S', 'Q')) 
  
  # plot
  g <- ggplot(combined_samples, aes(x=H1, y=trait))+
    geom_density_ridges(alpha=0.5, fill=viridis(10, option=proj_color)[3], scale=2) + 
    geom_text(
      data = Rbar_labels,
      aes(x = 0.2, y = trait, label = sapply(Rbar, function(r) bquote(bar(R^2) == .(r)))),
      parse = TRUE, # Necessary for math expressions
      color = viridis(10, option=proj_color)[3], size = 3.5, hjust = 0, vjust=-2
    ) +    
    labs(x=expression(Phylogenetic*' '*R^2), y='') +
    coord_cartesian(clip = "off") +
    theme_classic() + 
    theme(axis.text.y=element_text(size=12))
  
  return(g)
  
}
# Q <- tar_read(psQ)
# S <- tar_read(psS)
# Lifespan <- tar_read(psLifespan)
# AgeM <- tar_read(psAgeM)
# Length <- tar_read(psLength)
# SSD <- tar_read(psSSD)



# Coefficients table (brms) -----------------------------------------------
brms_table <- function(fit) {
  
  fixed <- data.table(as.data.frame(summary(fit, prob=0.9)$fixed)) 
  r1 <- data.table(as.data.frame(summary(fit, prob=0.9)$random$phylo)) 
  r2 <- data.table(as.data.frame(summary(fit, prob=0.9)$random$Species)) 
  raw.table <- rbindlist(list(fixed,r1,r2))
  
  stats.table <- cbind(c(row.names(summary(fit)$fixed),
                         'Phylogeny',
                         'Species'), 
                       raw.table)
  
  # may want to include phi or population-level effects for some models if we end up including these
  
  #stats.table <- cbind(row.names(summary(fit)$fixed), stats.table)
  stats.table[,c('Rhat','Bulk_ESS','Tail_ESS'):=NULL,]
  names(stats.table) <- c("Term", "Estimate", "SE", "PI-Lower", "PI-Upper")
  
  nice_table(stats.table)
  
}
# fit <- tar_read(m1q)



# Prep Ellis Lifespan data ------------------------------------------------
prep_ellisLifespanData <- function(species, d) {
  
  # replace species names due to wacky encoding issue
  for (i in 1:nrow(d)) {
    for (j in 1:length(species$Species)) {
      if ((stri_compare(d[i,Species,], species$Species[j]))==1) (d[i,SpeciesCorrected:=species$Species[j]])
    }
  }
  d[,Species:=SpeciesCorrected,]

  # make more informative column names
  d[,lifespan_median:=Median,]
  d[,lifespan_IC195:=IC195,]
  d[,lifespan_uC195:=uC195,]
  d[,lifespan_Post.Mean:=Post.Mean,]
  d[,lifespan_Post.sd:=Post.sd,]
  
  # make wide for M and F trait values
  wide <- dcast(d, Species ~ Sex, value.var = c("lifespan_median", "lifespan_IC195", "lifespan_uC195", "lifespan_Post.Mean", "lifespan_Post.sd"))
  
  return(wide)
  
}
# d <- tar_read(ellisLifespan_raw)
# species <- tar_read(species)



# Finalize dataframe for analysis -----------------------------------------
finalize_data <- function(d, exclude_from_SSD) {
  
  # add column for phylogenetic effects
  d[,phylo:=Species,]
  
  # modify 'phylo' to sister species that are actually in tree, for modelling
  d[Species=='Berardius minimus', phylo:='Berardius bairdii',]
  d[Species %in% c('Cephalorhynchus hectori','Cephalorhynchus eutropia'), phylo:='Cephalorhynchus commersonii',]
  d[Species=='Lagenorhynchus cruciger', phylo:='Lagenorhynchus australis',]
  d[Species=='Indopacetus pacificus', phylo:='Tasmacetus shepherdi',]
  d[Species=='Mesoplodon eueu', phylo:='Mesoplodon mirus',]
  d[Species=='Mesoplodon hotaula', phylo:='Mesoplodon ginkgodens',]
  d[Species=='Mesoplodon traversii', phylo:='Mesoplodon hectori',]
  d[Species=='Neophocaena asiaeorientalis', phylo:='Neophocaena phocaenoides',]
  d[Species=='Phocoena sinus', phylo:='Phocoena spinipinnis',]
  d[Species=='Sotalia fluviatilis', phylo:='Sotalia guianensis',]
  d[Species %in% c('Sousa plumbea', 'Sousa sahulensis', 'Sousa teuszii'), phylo:='Sousa chinensis',]
  d[Species=='Berardius arnuxii', phylo:='Berardius bairdii',]
  d[Species=='Hyperoodon planifrons', phylo:='Hyperoodon ampullatus',]
  d[Species=='Phocoenoides dalli', phylo:='Phocoena phocoena',]
  # Note: we are only changing phylo, still allowing species-specific random effects

  # log-transform body length
  d[,log_length_F:=log(length.mean_F),]
  d[,log_length_M:=log(length.mean_M),]
  
  # calculate sexual dimorphism in body length IF distinct F and M length measures exist
  d[!(Species %in% exclude_from_SSD),SSD:=length.mean_F/length.mean_M,]

  # Latitude and longitude
  d[,index:=.I,]
  d[,Latitude:=as.numeric(strsplit(Location,', ')[[1]][1]),by=index]
  d[,Longitude:=as.numeric(strsplit(Location,', ')[[1]][2]),by=index]
  
  # Study years
  d[!is.na(minYear),meanYear:=mean(c(minYear, maxYear)),by=index]
  
  # Calculate study duration
  d[,studyDuration:=(maxYear-minYear)+1,by=index]
  
  # Format network size as number
  d[,networkSize:=as.numeric(`Network size`),by=index]
  d[,`Network size`:=NULL,]
  
  # Format R and R_se as number
  d[,R:=as.numeric(R),by=index]
  d[,R_se:=as.numeric(R_se),by=index]
  
  return(d)
  
}
# d <- tar_read(data_c)
# exclude_from_SSD <- tar_read(exclude_from_SSD)



# Create IVSO dataframe ---------------------------------------------------
add_IVSO <- function(d) {
  
  # subset for Q
  q_data <- d[!(is.na(Title)) & !(is.na(Q)) & AssociationIndex %in% c('HWI','SRI','SAI','GAI (with HWI)') & !(Exclusions %in% c('Both','Q')) & SexFocus=='Assume mixed',,]
  q_data[, numObsQ:=.N, by=Species]
  q_data[, Q_IVSO:=calculate_cv(Q), by=Species]
  q_data <- unique(q_data[,c('Species','numObsQ','Q_IVSO')])
  
  # subset for S
  s_data <- d[!(is.na(Title)) & !(is.na(S)) & AssociationIndex %in% c('HWI','SRI','SAI','GAI (with HWI)') & !(Exclusions %in% c('Both','S')) & SexFocus=='Assume mixed',,]
  s_data <- s_data[S_method=='Likelihood',,]
  s_data[, numObsS:=.N, by=Species]
  s_data[, S_IVSO:=calculate_cv(S), by=Species]
  s_data <- unique(s_data[,c('Species','numObsS','S_IVSO')])
  s_data[, S_IVSO:=S_IVSO+0.001,] # avoids 0 for Globicephala macrorhynchus for modelling
   
  # merge back with life history data and remaining species
  u <- unique(d[,c('Species','phylo', 'lifespan_Post.Mean_F', 'log_length_F', 'age.mat_F', 'SSD')])
  u <- merge(u, q_data, by='Species', all.x=TRUE)
  u <- merge(u, s_data, by='Species', all.x=TRUE)
  
  return(u)
  
}
# d <- tar_read(data)



# Map_figure --------------------------------------------------------------
create_map <- function(d) {
  
  # remove rows without study
  d <- d[!is.na(Latitude),,]
  # create spatial object and load in world map
  spatial_data <- st_as_sf(d, coords = c("Longitude", "Latitude"), crs = 4326)
  world_map <- ne_countries(scale = "medium", returnclass = "sf")
  
  ggplot() +
    geom_sf(data = world_map, fill = "grey90", color = 'grey80') +
    geom_sf(data = spatial_data, aes(color = common, shape = common), alpha = 0.75, size = 3.5) +
    coord_sf() +
    scale_color_viridis(discrete = TRUE, option = proj_color, begin = 0.95, end = 0.0) +
    scale_shape_manual(values = rep(c(15,16,18), length.out = 23)) + # Recycle shapes to match 23 groups
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    ylim(-65, 90)
  
}
# d <- tar_read(data)



# Identify mismatches in species between tree and database ----------------
identify_mismatch <- function(species, tree, database) {
  
  tree_species <- unique(tree$tip.label)
  database_species <- unique(database$Species)
  
  test_that("All species in database are on official list", expect_true(all(database$Species %in% species$Species)))
  test_that("All species in tree are on species list", expect_true(all(tree$tip.label %in% species$Species)))
  
  # Messages to identify species mismatches between data and phylogeny
  Missing_from_tree <- species$Species[!(species$Species %in% tree_species)]
  
  return(c(Missing_from_tree))
  
}
# species <- tar_read(species)
# tree <- tar_read(tree)
# database <- tar_read(database)



# Coevolve_visual_summary -------------------------------------------------
GDPM_visual_summary <- function(fit, name1, name2) {
  
  g1 <- coev_plot_delta_theta(fit)
  # g2 <- coev_plot_flowfield(fit, var1=name1, var2=name2)
  g3 <- coev_plot_pred_series(fit)
  # g4 <- coev_plot_predictive_check(fit)
  g5 <- coev_plot_selection_gradient(fit, var1=name1, var2=name2)
  g6 <- coev_plot_trait_values(fit)
  
  (g1 | g3) / (g5 | g6)
  
}
# fit <- tar_read(coev_Q)
# name1 <- 'Q'
# name2 <- 'log_length_F'



# Identify species for valid SSD ------------------------------------------
SSD_exclusions <- function(species, ellisLifeHistory) {
  
  # Note: The summary of the marinelifehist length data includes what appear to be sex-specific measures
  # that may be have been assigned from mixed-sex groups or from the opposite sex when data for a given 
  # sex was unavailable. 
  # Here, we return to the raw data to identify a list of species for which both M and F length data
  # were actually available. It is for these species that calculating sexual size dimorphism will be appropriate.
  
  # We carried out two exclusions prior to calculating SSD:
  # 1. Any species lacking data from both Females and Males were excluded.
  # 2. Of species with both Female and Male data, if the sexes were ultimately assigned identical body length in the database,
  # we checked additional sources to ensure that it is plausible that they should be the same length. Any that were not plausibly
  # the same length across sexes but had identical listed mean lengths were excluded (N=3).

  # pull out raw size data, focusing on body length
  raw <- data.table(marine.lifehist.speciesdata$species_raw.sizes)
  raw <- raw[measurement=='length']
  
  # pull out species-sex-specific summaries of length and convert to wide: are any identical across sex? 
  length <- data.table(marine.lifehist.speciesdata$species_length)
  wide <- dcast(length, species ~ sex, value.var = c("mean", "sd"))
  
  # loop through species where derived mean length is equal for males and females: do any have F and M sex information?
  for (common in wide[mean_F==mean_M,species,]) {
    sex_info <- raw[species==common,unique(sex),]
    if ('M' %in% sex_info & 'F' %in% sex_info) print(common)
  }
  
  # Explore which species have data from both sexes
  
  raw[,M_data:=ifelse('M' %in% sex,1,0), by='species']
  raw[,F_data:=ifelse('F' %in% sex,1,0), by='species']
  
  raw[,MF_data:=sum(M_data, F_data), by=.I] # 2 if both
  raw[,MF_data:=sum(M_data, F_data), by=.I]
  
  # fix up species names
  species_key <- data.table(marine.lifehist.speciesdata$species_names.key)
  raw <- merge(raw, species_key[,c('species','latin.name')],by='species')
  raw[,Species:=latin.name,]

  # As of December 2024, the following species had M and F measures despite having the same ultimate mean value
  # Assess these species manually...
  # Can also check out bibliography: 
  bib <- marinelifehistdata::bibliography
  
  # CuviersBeakedWhale
  raw[species=='CuviersBeakedWhale',,] # max. values for M and F: duplicated, CHECK REFERENCE 131
  # From quick search, appears that Cuviers do tend to be similar in size
  # Accurate -- bizarrely, source does report separate findings of 670/700 for M/F AND vice-versa
  # Therefore Cuvier's can genuinely be assumed to have approx. the same length across sex -- no need to exclude
  
  # DuskyDolphin
  raw[species=='DuskyDolphin',,] # min. and max. values for both sexes, but appears to be using combined mean?
  # Trouble tracking down source, though it does appear that dusky dolphin sizes are pretty similar -- no need to exclude [but perhaps check with S.E.]
  
  # Franciscana
  raw[species=='Franciscana',,] # min. and max. values for both sexes
  # Same as above: using mean value from combined measure, despite some sex-specific data
  # However, both from database mean and max + online search, females should be larger -- EXCLUDE FOR NOW
  
  # GuianaDolphin
  raw[species=='GuianaDolphin',,] # max. values for both sexes
  # Using data from combined set, seems like they are similar in size -- no need to exclude, G.M.S. says similar [but perhaps check with S.E. re: database]
  
  # PacificWhiteSidedDolphin
  raw[species=='PacificWhiteSidedDolphin',,] # max. values for both sexes
  # Appears to be some sexual dimorphism based on an online search -- EXCLUDE FOR NOW
  
  # PygmyKillerWhale
  raw[species=='PygmyKillerWhale',,] # max. values for both sexes
  # Some evidence of sexual size dimorphism -- EXCLUDE FOR NOW
  
  # 1. Exclude species lacking data from both sexes
  exclusion_1 <- raw[MF_data<2, unique(Species),]
  
  # 2. Exclude species with sex-specific data, but equal mean value for each sex AND evidence of SSD in literature
  exclusion_2 <- c('Pontoporia blainvillei', 'Lagenorhynchus obliquidens', 'Feresa attenuata')

  # Combine exclusions
  exclusions <- c(exclusion_1, exclusion_2)

  # Species no longer excluded based on updated protocol [shouldn't make much of a difference, only measures of Q or S for first 2]
  # "Lagenorhynchus obscurus"
  # "Sotalia guianensis" 
  # "Ziphius cavirostris"  
  
  test_that("All excluded species match official species names", expect_true(all(exclusions %in% species$Species)))
  
  return(exclusions)

}
# species <- tar_read(species)
# ellisLifeHistory <- tar_read(ellisLifeHistory)

# Questions for Sam Ellis:

# Trying to figure out when F and M lengths can be reliably used to estimate SSD.
# Excluding species with ONLY combined-sex length data -- check -- is this sensible based on his protocol?
# However, some species with sex-specific data that have same mean length across sexes. 
# Some of these make sense (e.g., Cuvier's, Dusky), though others do seem to have SSD. I am excluding these -- is this sensible based on his protocol?
# Anything else I should be aware of or could be missing?

# Relatedly: why do killer whales have sex-specific mean lengths but not Guiana dolphins, for example?




# Ellis_Data --------------------------------------------------------------
extract_Ellis_data <- function(species) {
  
  # extract compiled data
  length <- data.table(marine.lifehist.speciesdata$species_length)
  maturity <- data.table(marine.lifehist.speciesdata$species_age.maturity)
  
  # merge into single dataframe
  combined <- merge(length, maturity, by=c('species', 'sex')) 
  combined[,length.mean:=mean,]
  combined[,length.sd:=sd,]
  combined[,c('mean','sd'):=NULL] # replace with more informative column names

  # fix up species names
  species_key <- data.table(marine.lifehist.speciesdata$species_names.key)
  combined <- merge(combined, species_key[,c('species','latin.name')],by='species')
  combined[,Species:=latin.name,]
  combined[,species:=NULL,] # remove common name
  
  # make wide for M and F trait values
  wide <- dcast(combined, latin.name + Species ~ sex, value.var = c("age.mat", "age.mat.range.min", "age.mat.range.max", "length.mean", "length.sd"))
  
  # confirm all species in official list
  test_that("All species in database are on official list", expect_true(all(wide$Species %in% species$Species)))
  
  return(wide)  
  
}
# species <- tar_read(species)



# Prepare McGowen trees ---------------------------------------------------
prepare_mcgowen_trees <- function(species, key, path) {
  
  # load data
  tree <- read.nexus(path)
  tree <- drop.tip(tree, tip='Dcap108471') # Genetic record for D. capensis excised, as lumped into D. delphis in official species list we are using
  
  # check conflicts between official species list and McGowen format
  setdiff(species$Species, key$Species)
  setdiff(key$Species, species$Species) # not toothed whales

  # merge species codes and species names
  tips <- data.frame(Code = tree$tip.label)
  tips <- merge(tips, key, by='Code', all.x=TRUE, sort = FALSE) # sort=FALSE retains order, but ensure to check resulting trees carefully 

  ###### gene keys all look sensible expect for sotalia, which may be S. fluviatilis. Would not change the analysis, but might be nice to confirm. 
  
  # assign species based on key    
  tree$tip.label <- tips$Species
  
  # exclude terrestrial mammals
  tree <- drop.tip(tree, tip = c('Equus caballus',
                                  'Vicugna pacos',
                                  'Camelus bactrianus',
                                  'Bos taurus',
                                  'Oryx leucoryx',
                                  'Tragelaphus eurycerus',
                                  'Choeropsis liberiensis',
                                  'Ovis aries',
                                  'Hippopotamus amphibius',
                                  'Panthalops hodgsonii',
                                  'Gazella arabica',
                                  'Sus scrofa'))
  
  # fix spelling mistake in gray whale species name
  tree$tip.label[tree$tip.label=='Eschrictius robustus'] <- 'Eschrichtius robustus'
  
  # exclude baleen whales
  tree <- drop.tip(tree, tip = c('Balaena mysticetus',
                                 'Eubalaena australis',
                                 'Eubalaena glacialis',
                                 'Eubalaena japonica',
                                 'Caperea marginata',
                                 'Balaenoptera acutorostrata',
                                 'Balaenoptera borealis',
                                 'Balaenoptera edeni',
                                 'Balaenoptera bonaerensis',
                                 'Balaenoptera musculus',
                                 'Balaenoptera physalus',
                                 'Eschrichtius robustus',
                                 'Megaptera novaeangliae'))
  
  # add tips for missing species
  tree <- add.cherry(tree, 'Berardius bairdii', new.tips = 'Berardius minimus')
  tree <- add.cherry(tree, 'Cephalorhynchus commersonii', new.tips = c('Cephalorhynchus hectori', 'Cephalorhynchus eutropia'))
  tree <- add.cherry(tree, 'Lagenorhynchus australis', new.tips = 'Lagenorhynchus cruciger')
  tree <- add.cherry(tree, 'Tasmacetus shepherdi', new.tips = 'Indopacetus pacificus')
  tree <- add.cherry(tree, 'Mesoplodon mirus', new.tips = 'Mesoplodon eueu')
  tree <- add.cherry(tree, 'Mesoplodon ginkgodens', new.tips = 'Mesoplodon hotaula')
  tree <- add.cherry(tree, 'Mesoplodon hectori', new.tips = 'Mesoplodon traversii')
  tree <- add.cherry(tree, 'Neophocaena phocaenoides', new.tips = 'Neophocaena asiaeorientalis')
  tree <- add.cherry(tree, 'Phocoena spinipinnis', new.tips = 'Phocoena sinus')
  tree <- add.cherry(tree, 'Sotalia guianensis', new.tips = 'Sotalia fluviatilis')
  tree <- add.cherry(tree, 'Sousa chinensis', new.tips = c('Sousa plumbea', 'Sousa sahulensis', 'Sousa teuszii'))
  
  # Berardius arnuxii - Michael says with other Berardius
  tree <- add.cherry(tree, 'Berardius bairdii', new.tips = 'Berardius arnuxii')
  # Hyperoodon planifrons - Michael says lump with other Hyperoodon
  tree <- add.cherry(tree, 'Hyperoodon ampullatus', new.tips = 'Hyperoodon planifrons')
  
  setdiff(tree$tip.label, species$Species) # Good - all species in tree are in official species list
  setdiff(species$Species, tree$tip.label) # Three species in official species list that are not in tree
  
  # Phocoenoides dalli 
  # Platanista gangetica
  # Platanista minor
  
  # ggtree(tree) + geom_tiplab(size=2) + xlim(0,0.4)
  return(tree)
  
}
# species <- tar_read(species)
# path <- 'Input/McGowen_Trees/FigTree_parts_6_mcmctree_AR.tre'
# key <- tar_read(gene_key)
# consensus <- 'FALSE'



# Prepare McGowen trees ---------------------------------------------------
prepare_mcgowen_trees_noSister <- function(species, key, path) {
  
  # load data
  tree <- read.nexus(path)
  tree <- drop.tip(tree, tip='Dcap108471') # slice out capensis -- maybe just justify this and go with it for final version
  
  # check conflicts between official species list and McGowen format
  setdiff(species$Species, key$Species)
  setdiff(key$Species, species$Species)
  
  
  tips <- data.frame(Code = tree$tip.label)
  tips <- merge(tips, key, by='Code', all.x=TRUE, sort = FALSE) # sort=FALSE retains order, but ensure to check resulting trees carefully ######
  
  # assign species based on key    
  tree$tip.label <- tips$Species
  
  # exclude terrestrial mammals
  tree <- drop.tip(tree, tip = c('Equus caballus',
                                 'Vicugna pacos',
                                 'Camelus bactrianus',
                                 'Bos taurus',
                                 'Oryx leucoryx',
                                 'Tragelaphus eurycerus',
                                 'Choeropsis liberiensis',
                                 'Ovis aries',
                                 'Hippopotamus amphibius',
                                 'Panthalops hodgsonii',
                                 'Gazella arabica',
                                 'Sus scrofa'))
  
  # fix spelling mistake in gray whale species name
  tree$tip.label[tree$tip.label=='Eschrictius robustus'] <- 'Eschrichtius robustus'
  
  # exclude baleen whales
  tree <- drop.tip(tree, tip = c('Balaena mysticetus',
                                 'Eubalaena australis',
                                 'Eubalaena glacialis',
                                 'Eubalaena japonica',
                                 'Caperea marginata',
                                 'Balaenoptera acutorostrata',
                                 'Balaenoptera borealis',
                                 'Balaenoptera edeni',
                                 'Balaenoptera bonaerensis',
                                 'Balaenoptera musculus',
                                 'Balaenoptera physalus',
                                 'Eschrichtius robustus',
                                 'Megaptera novaeangliae'))

  setdiff(tree$tip.label, species$Species) # Good - all species in tree are in official species list

  # ggtree(tree) + geom_tiplab(size=2) + xlim(0,0.4)
  return(tree)
  
}
# species <- tar_read(species)
# path <- 'Input/McGowen_Trees/FigTree_parts_6_mcmctree_AR.tre'
# key <- tar_read(gene_key)
# consensus <- 'FALSE'



# Function for adding cherry to tree --------------------------------------
add.cherry <- function(tree, tip, new.tips) {
  
  ## Find the edge leading to the tip
  tip_id <- match(tip, tree$tip.label)
  
  ## Create the new cherry
  tree_to_add <- ape::stree(length(c(tip, new.tips)))
  
  ## Naming the tips
  tree_to_add$tip.label <- c(tip, new.tips)
  
  ## Add 0 branch length
  tree_to_add$edge.length <- rep(0, Nedge(tree_to_add))
  
  ## Binding both trees
  return(bind.tree(tree, tree_to_add, where = tip_id))
}



# Plot univariate effect from static phylogenetic model -------------------
plot_univariate <- function(fit, traitName, ytraitName, x_lab, y_lab, ymax, color) {
  
  # Set up color
  cols <- viridis(6, option=proj_color)
  col=cols[color]
  
  # Extract data
  data <- fit$data
  
  # New data for prediction
  newdata <- expand.grid(trait = seq( min(data[,paste(traitName)]), max(data[,paste(traitName)]), length.out=30),
                         phylo = NA,
                         Species = NA)
  names(newdata)[names(newdata) == "trait"] <- paste(traitName)
  
  # Compute predictions
  pred <- data.table(epred_draws(fit, newdata, re_formula=NA)) # not incorporating varying effects in prediction here
  pred_means <- pred[, .(mean_epred=mean(.epred)),by=.(get(traitName))] 
  
  # Plotting
  p <- ggplot(pred, aes(x=get(traitName), y=.epred))+
    geom_point(inherit.aes = FALSE, data=data, aes(x=get(traitName), y=get(ytraitName)), alpha=0.75, size=1.5, color='grey40')+
    stat_lineribbon(fill=col, alpha=0.1, .width=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9))+
    labs(x='', y=y_lab)+
    # ylim(0,ymax)+
    ggtitle(x_lab)+
    theme_classic() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          legend.position = "none")
  p
  
  return(p)
  
}
# fit <- tar_read(m4q)
# traitName <- 'SSD' ######
# ytraitName <- 'Q'
# ymax <- 0.75
# x_lab <- 'SSD'
# y_lab <- 'Social modularity'
# color <- 1



# Predictions for ancestral state -----------------------------------------
ancestral_state <- function(fit, ancestral_length_raw, colnum, var, label_position, xlim) {
  
  # Using 250 cm as ancestral body length for prediction
  ancestral_length <- log(ancestral_length_raw)

  pred <- data.frame(posterior_epred(fit, 
                          newdata=data.frame(log_length_F = ancestral_length),
                          re_formula = NA))
  colnames(pred) <- 'anc_prediction'
  
  # colors
  cols <- viridis(10, option=proj_color)
  
  # calculate mean value
  med_val <- round(median(pred$anc_prediction),2)
  
  # plot labels
  xlab <- ifelse(var=='Q', 'Ancestral modularity', 'Ancestral social differentiation')
  
  # threshold to show
  threshold_int <- ifelse(var=='Q', 0.3, 0.5)
  
  ggplot(pred, aes(x=anc_prediction, y=0, fill=after_stat(x))) + 
    geom_density_ridges_gradient(alpha=0.5, color=NA) +
    scale_fill_viridis(option = proj_color, begin=0, end=0.9, direction=-1) +
    geom_vline(xintercept = threshold_int, color="grey", linetype='dashed', lwd=0.75) +
    labs(y='Density', x=xlab) +
    # annotate("text",label=bquote(bar(.(var))[Ancestral] == .(med_val)), x=label_position[1], y=label_position[2], color='black', size=5) +
    annotate("text", label = as.expression(bquote(bar(.(var))[Ancestral] == .(med_val))), x=label_position[1], y=label_position[2], color='black', size=5) +
    theme_classic() + 
    xlim(xlim[1], xlim[2])+
    theme(legend.position='none')
  
  
}
# fit <- tar_read(m3s)
# colnum <- 6
# var <- 'S'
# ancestral_length_raw <- 250
# label_position <- c(0.5, 2.5)
# label_position <- c(2, 1)
# xlim <- c(0,3)

# fit <- tar_read(m3q)
# colnum <- 6
# var <- 'Q'
# ancestral_length_raw <- 250
# label_position <- c(0.5, 2.5)
# label_position <- c(2, 1)
# xlim <- c(0,3)




# Calculate coefficient of variation --------------------------------------
calculate_cv <- function(measures) {
  
  # remove NAs and calculate coefficient of variation
  return(sd(measures, na.rm=TRUE) / mean(measures, na.rm=TRUE))
  
}
# measures <- c(1,2,3,4,5)
# Note this employs sample SD



# Plot Q vs. S ------------------------------------------------------------
plot_qs <- function(d) {
  
  # Note we are not incorporating any exclusions here
  
  # Calculate empirical 50% IQR of Q and S for each species
  d[,Q_lower:=quantile(Q,na.rm=TRUE,probs = 0.25), by='Species']
  d[,Q_upper:=quantile(Q,na.rm=TRUE,probs = 0.75), by='Species']
  d[,S_lower:=quantile(S,na.rm=TRUE,probs = 0.25), by='Species']
  d[,S_upper:=quantile(S,na.rm=TRUE,probs = 0.75), by='Species']
  
  # Calculate mean of Q and S by species
  d[,Q_mean:=mean(Q, na.rm=TRUE), by='Species']
  d[,S_mean:=mean(S, na.rm=TRUE), by='Species']

  # Subset to unique values
  d <- unique(d[,c('Species', 'common','Q_mean','S_mean','Q_lower','Q_upper','S_lower','S_upper')])
  
  ggplot(d, aes(x=S_mean, y=Q_mean, color=Species, label=common)) + 
    
    geom_hline(yintercept = 0.3, linetype='dashed', color='grey') +
    geom_vline(xintercept = 0.5, linetype='dashed', color='grey') +
    
    geom_errorbar(aes(ymin=Q_lower, ymax=Q_upper), width=0, linewidth=0.5) + 
    geom_errorbarh(aes(xmin=S_lower, xmax=S_upper), height=0, linewidth=0.5) +
    geom_point(pch=21, fill='white', size=2, stroke=1.6)  +
    geom_text_repel(size=2, segment.colour=NA)+
    
    scale_color_viridis(discrete=TRUE,option=proj_color, begin = 0, end=0.85)+
    labs(x='Social differentiation (S)', y='Modularity (Q)') +
    xlim(0,2)+ylim(0,1)+
    theme_classic()+
    theme(legend.position='none') 
  
}
# d <- tar_read(data)



# Plot tree with two traits on tips ---------------------------------------
plot_tree_two_adjacent_traits <- function(t, data, traitName1, traitName2, flip, just_labels) {
  
  # Restrict to species with a value for phylo
  d <- data[phylo %in% t$tip.label, , ]
  
  # Extract median values for both traits
  d[, medianTrait1 := median(get(traitName1), na.rm = TRUE), by = 'phylo']
  d[, medianTrait2 := median(get(traitName2), na.rm = TRUE), by = 'phylo']
  
  # Keep unique values for both traits
  d <- unique(d[, .(phylo, medianTrait1, medianTrait2)])
  
  # Create the initial tree plot
  p <- ggtree(t, colour = 'black') + theme(legend.position = "right")
  if (flip) {
    p <- p + scale_x_reverse() + theme(legend.position = c(0.8, 0.8))
  }
  
  # Merge the data into the ggtree data structure
  p <- p %<+% d
  
  # add background tip point
  # p <- p + geom_tippoint(size = 5.5, color='black') 
    # scale_color_viridis(name = paste(traitName1, "Median", sep = " "), option = proj_color, direction = -1, na.value = 'white')
  
  # Add tip points for the first trait (plot all points)
  p <- p + geom_tippoint(aes(color = medianTrait1), size = 5) +
    scale_color_viridis(name = paste(traitName1, "Median", sep = " "), option = proj_color, direction = -1, na.value = 'white')
  
  # Add adjacent tip points for the second trait, skipping NA values but keeping node info
  p <- p + geom_tippoint(data = p$data[!is.na(p$data$medianTrait2), ], 
                         aes(fill = medianTrait2), size = 5, shape = 21, color = "white", position = position_nudge(x = 0.02)) +
    scale_fill_viridis(name = paste(traitName2, "Median", sep = " "), option = proj_color, direction = -1)
  
  if (just_labels) {
    p <- ggplot(data = p$data, aes(y = node, x = 1, label = label)) +
      geom_text(size = 3) +
      scale_y_reverse(limits = c(nrow(p$data[!is.na(p$data$label), ]), 0)) +
      theme_void()
  }
  
  # format legends
  p <- p + theme(legend.position=c(0.25, 0.8), legend.direction = 'horizontal')
  
  # add clade labels
  clades <- data.frame(node = c(116, 100, 64, 98, 97, 94),
                       label=c('Sperm whales (Physeteroidea)',
                               'Beaked whales (Ziphiidae)',
                               'Oceanic dolphins (Delphinidae)',
                               'River dolphins (various)',
                               'Arctic whales (Monodontidae)','Porpoises (Phocoenidae)'))
  
  # Size scale for phylopic silhouettes
  size <- 0.009
  # Add phylopic silhouettes
  p <- p + geom_cladelab(node=clades$node, label=clades$label, offset=0.035) +
    xlim(0,0.6) +
    add_phylopic(name='Physeter macrocephalus', x=0.45, y=4.5, width = size*10, alpha=1, fill = 'grey30', horizontal=TRUE)+
    add_phylopic(name='Hyperoodon ampullatus', x=0.45, y=14, width = size*7, alpha=1, fill = 'grey30')+
    add_phylopic(name='Inia geoffrensis', x=0.45, y=22.75, width = size*2.5, alpha=1, fill = 'grey30')+
    add_phylopic(name='Phocoenoides dalli', x=0.45, y=28.75, width = size*2.1, alpha=1, fill = 'grey30')+
    add_phylopic(name='Orcinus orca', x=0.45, y=47, width = size*6, alpha=1, fill = 'grey30') +
    add_phylopic(name='Delphinapterus leucas', x=0.45, y=25.75, width = size*3.6, alpha=1, fill = 'grey30')
  
  # Annotations
  p <- p + annotate('text', label=expression(Q), x=0.3425, y=-1.75, size=6) +
    annotate('text', label=expression(S), x=0.361, y=-1.75, size=6)
  
  p
  
  # SW: 11m
  # HA: 7m
  # Inia: 2.5m
  # Dall's: 2.1
  # Orca: 6m
  # Narwhal: 4

  # ggtree(t) + geom_text(aes(label=node), hjust=-.3)
  # ggtree(t) + geom_tiplab()

  return(p)
}
# t <- tar_read(tree)
# data <- tar_read(data)
# traitName1 <- 'Q'
# traitName2 <- 'S'
# flip <- FALSE
# just_labels <- FALSE



# Trait change plot -------------------------------------------------------
trait_change_plot_exclude <- function(fit, exclude) {
  
  var1 <- names(fit$variables[1])
  var2 <- names(fit$variables[2])
  
  lab1 <- var1
  lab2 <- var2
  
  delta_theta_1 <- data.table(data.frame(coev_calculate_delta_theta(fit, response = var1, predictor=var2)))
  delta_theta_1[,response:=var1]
  delta_theta_1[,predictor:=var2]
  
  delta_theta_2 <- data.table(data.frame(coev_calculate_delta_theta(fit, response = var2, predictor=var1)))
  delta_theta_2[,response:=var2]
  delta_theta_2[,predictor:=var1]
  
  df <- rbindlist(list(delta_theta_1, delta_theta_2))
  # order factor for consistent plotting
  df$response <- factor(df$response, levels=c(paste(var2), paste(var1)))

  # exclusions if desired for making figures for presentations
  df <- df[predictor==exclude,X1.delta_theta:=NA,]
  
  # rename with nice labels for plotting
  if (var1=='log_length_F') (lab1<-'Length (log)')
  if (var2=='log_length_F') (lab2<-'Length (log)')
  
  if (var1=='length.mean_F') (lab1<-'Length')
  if (var2=='length.mean_F') (lab2<-'Length')
  
  if (var1=='age.mat_F') (lab1<-'Age maturity')
  if (var2=='age.mat_F') (lab2<-'Age maturity')
  
  if (var1=='lifespan_Post.Mean_F') (lab1<-'Lifespan')
  if (var2=='lifespan_Post.Mean_F') (lab2<-'Lifespan')
  
  # Set up plot colors
  plot_cols <- viridis(2, option=proj_color, begin=0.75, end=0.05)
  if (exclude=='Q') (plot_cols <- viridis(2, option=proj_color, begin=0.05, end=0.75)) # Just for Princeton plotting
  
  # Plotting
  p <- ggplot(df, aes(x=X1.delta_theta)) +
    
    geom_density(aes(fill=response, color=response), alpha=0.6) +
    geom_vline(xintercept = 0, color="black", linetype='dashed', lwd=0.8) +
    annotate("text",
             label = as.expression(bquote(bold(.(lab2) %->% .(lab1)) ~ ", " ~ italic(p)[dir] ~ "=" ~ .(round(as.numeric(p_direction(df[response == var1, X1.delta_theta])), 2)))),
             x = -Inf, y = Inf, hjust = -0.1, vjust = 2, size = 4, color = plot_cols[1]) +
    annotate("text",
             label = as.expression(bquote(bold(.(lab1) %->% .(lab2)) ~ ", " ~ italic(p)[dir] ~ "=" ~ .(round(as.numeric(p_direction(df[response == var2, X1.delta_theta])), 2)))),
             x = -Inf, y = Inf, hjust = -0.1, vjust = 4, size = 4, color = plot_cols[2]) +
    theme_minimal() +
    scale_x_continuous(limits=c(-20,20)) +
    ylim(0,0.38)+
    scale_color_manual(values=plot_cols)+
    scale_fill_manual(values=plot_cols)+
        scale_size_manual(values=c(1,1,1.5)) +
    theme_bw(base_size=14) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), legend.position = "none") +
    xlab(expression(paste(Delta, theta)["z"])) +
    ylab("")
  
  p
  
}
# fit <- tar_read(coev_Q_noTransform)
# exclude <- 'Q'




# Trait change plot -------------------------------------------------------
trait_change_plot <- function(fit) {
  
  var1 <- names(fit$variables[1])
  var2 <- names(fit$variables[2])
  
  lab1 <- var1
  lab2 <- var2
  
  delta_theta_1 <- data.table(data.frame(coev_calculate_delta_theta(fit, response = var1, predictor=var2)))
  delta_theta_1[,response:=var1]
  delta_theta_1[,predictor:=var2]
  
  delta_theta_2 <- data.table(data.frame(coev_calculate_delta_theta(fit, response = var2, predictor=var1)))
  delta_theta_2[,response:=var2]
  delta_theta_2[,predictor:=var1]
  
  df <- rbindlist(list(delta_theta_1, delta_theta_2))
  # order factor for consistent plotting
  df$response <- factor(df$response, levels=c(paste(var2), paste(var1)))
  
  # rename with nice labels for plotting
  if (var1=='log_length_F') (lab1<-'Length (log)')
  if (var2=='log_length_F') (lab2<-'Length (log)')
  
  if (var1=='length.mean_F') (lab1<-'Length')
  if (var2=='length.mean_F') (lab2<-'Length')
  
  if (var1=='age.mat_F') (lab1<-'Age maturity')
  if (var2=='age.mat_F') (lab2<-'Age maturity')
  
  if (var1=='lifespan_Post.Mean_F') (lab1<-'Lifespan')
  if (var2=='lifespan_Post.Mean_F') (lab2<-'Lifespan')
  
  # Set up plot colors
  plot_cols <- viridis(2, option=proj_color, begin=0.75, end=0.05)

  # Plotting
  p <- ggplot(df, aes(x=X1.delta_theta)) +
    
    geom_density(aes(fill=response, color=response), alpha=0.6) +
    geom_vline(xintercept = 0, color="black", linetype='dashed', lwd=0.8) +
    annotate("text",
             label = as.expression(bquote(bold(.(lab2) %->% .(lab1)) ~ ", " ~ italic(p)[dir] ~ "=" ~ .(round(as.numeric(p_direction(df[response == var1, X1.delta_theta])), 2)))),
             x = -Inf, y = Inf, hjust = -0.1, vjust = 2, size = 4, color = plot_cols[2]) +
    annotate("text",
             label = as.expression(bquote(bold(.(lab1) %->% .(lab2)) ~ ", " ~ italic(p)[dir] ~ "=" ~ .(round(as.numeric(p_direction(df[response == var2, X1.delta_theta])), 2)))),
             x = -Inf, y = Inf, hjust = -0.1, vjust = 4, size = 4, color = plot_cols[1]) +
    theme_minimal() +
    scale_x_continuous(limits=c(-20,20)) +
    # ylim(0,0.38)+
    scale_color_manual(values=plot_cols)+
    scale_fill_manual(values=plot_cols)+
    scale_size_manual(values=c(1,1,1.5)) +
    theme_bw(base_size=14) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), legend.position = "none") +
    xlab(expression(paste(Delta, theta)["z"])) +
    ylab("Density")
  
  p
  
}
# fit <- tar_read(coev_Q_noTransform)
# fit <- tar_read(coev_Q_lifespan)



# Describe sample size of brms model --------------------------------------
describe_sample_size <- function(model_fit) {
  
  # Extract the number of observations
  num_obs <- nobs(model_fit)
  
  # Extract the number of groups for each random effect
  random_effects <- ngrps(model_fit)
  
  # Handle cases where random effects are missing or empty
  if (length(random_effects) == 0) {
    effect_statements <- "no random effects specified"
  } else {
    effect_statements <- paste(sapply(names(random_effects), function(effect) {
      paste0(effect, " with ", random_effects[[effect]], " groups")
    }), collapse = ", ")
  }
  
  # Construct the final statement
  statement <- paste0(
    "The model was fitted using ", num_obs, 
    " observations. It includes the following random effects: ", 
    effect_statements, "."
  )
  
  return(statement)
}
# model_fit <- tar_read(m1q)



# Describe sample size of brms model, with rhat and div. trans. -----------
describe_sample_size_diagnostic <- function(model_fit) {
  # Extract the number of observations
  num_obs <- nobs(model_fit)
  
  # Extract the number of groups for each random effect
  random_effects <- ngrps(model_fit)
  
  # Handle cases where random effects are missing or empty
  if (length(random_effects) == 0) {
    effect_statements <- "no random effects specified"
  } else {
    effect_statements <- paste(sapply(names(random_effects), function(effect) {
      paste0(effect, " with ", random_effects[[effect]], " groups")
    }), collapse = ", ")
  }
  
  # Extract the number of divergent transitions
  num_divergent <- sum(subset(nuts_params(model_fit), Parameter == "divergent__")$Value)
  
  # Extract the maximum Rhat
  max_rhat <- max(summary(model_fit)$fixed[, "Rhat"], na.rm = TRUE) # Adjust for the right element if necessary
  
  # Construct the final statement
  statement <- paste0(
    "The model was fitted using ", num_obs, 
    " observations. It includes the following random effects: ", 
    effect_statements, ". The number of divergent transitions is ", 
    num_divergent, ". The maximum Rhat value is ", round(max_rhat, 3), "."
  )
  
  return(statement)
}
# model_fit <- tar_read(m1q)



# Table for GPDM (coevolutionary model) -----------------------------------
coevolve_table <- function(object, prob=0.9, robust=FALSE) {
  
  # adapted from "summary" function of coevolve package 
  
  # percentiles to use
  probs <- c(((1 - prob) / 2), 1 - ((1 - prob) / 2))
  
  # get summary of selection and drift parameters from cmdstanr
  s <-
    as.data.frame(
      object$fit$summary(
        NULL,
        Estimate = ifelse(robust, "median", "mean"),
        `Est.Error` = ifelse(robust, "mad", "sd"),
        ~quantile(.x, probs = probs),
        Rhat = "rhat",
        Bulk_ESS = "ess_bulk",
        Tail_ESS = "ess_tail"
      )
    )
  
  # summarise autoregressive selection effects
  A <- s[stringr::str_starts(s$variable, pattern = "A\\["),]
  equal <-
    readr::parse_number(
      stringr::str_extract(A$variable, pattern = "A\\[(\\d+)\\,")
    ) ==
    readr::parse_number(
      stringr::str_extract(A$variable, pattern = "\\,\\d+\\]")
    )
  auto <- A[equal,]
  rownames(auto) <- names(object$variables)[
    readr::parse_number(
      stringr::str_extract(auto$variable, pattern = "A\\[(\\d+)\\,")
    )
  ]
  auto <- auto[, 2:ncol(auto)]
  
  
  # summarise cross selection effects
  cross <- A[!equal,]
  rownames(cross) <-
    paste0(
      names(object$variables)[
        readr::parse_number(
          stringr::str_extract(cross$variable, pattern = "\\,\\d+\\]")
        )
      ],
      " \U27F6 ",
      names(object$variables)[
        readr::parse_number(
          stringr::str_extract(cross$variable, pattern = "A\\[(\\d+)\\,")
        )
      ]
    )
  cross <- cross[, 2:ncol(cross)]
  # only include cross selection effects that have been estimated in summary
  cross <- cross[!is.na(cross$Rhat),]
  
  # summarise drift sd parameters
  sd_drift <- s[stringr::str_starts(s$variable, pattern = "Q_sigma\\["),]
  rownames(sd_drift) <- paste0(
    "sd(",
    names(object$variables)[
      readr::parse_number(
        stringr::str_extract(sd_drift$variable, pattern = "Q_sigma\\[(\\d+)\\]")
      )
    ],
    ")"
  )
  sd_drift <- sd_drift[, 2:ncol(sd_drift)]
  
  # summarise drift cor parameters
  cor_drift <- NULL
  if (object$estimate_Q_offdiag) {
    cor_drift <- s[stringr::str_starts(s$variable, "cor_R"),]
    for (i in 1:length(object$variables)) {
      for (j in 1:length(object$variables)) {
        if (i >= j) {
          var <- paste0("cor_R[", i, ",", j, "]")
          cor_drift <- cor_drift[cor_drift$variable != var,]
        }
      }
    }
    rownames(cor_drift) <- paste0(
      "cor(",
      names(object$variables)[
        readr::parse_number(
          stringr::str_extract(
            cor_drift$variable,
            pattern = "cor\\_R\\[(\\d+)\\,"
          )
        )
      ],
      ",",
      names(object$variables)[
        readr::parse_number(
          stringr::str_extract(
            cor_drift$variable,
            pattern = "\\,(\\d+)\\]"
          )
        )
      ],
      ")"
    )
    cor_drift <- cor_drift[, 2:ncol(cor_drift)]
  }
  
  # summarise SDE intercepts
  sde_intercepts <- s[stringr::str_starts(s$variable, "b"),]
  rownames(sde_intercepts) <- names(object$variables)[
    readr::parse_number(sde_intercepts$variable)
  ]
  sde_intercepts <- sde_intercepts[, 2:ncol(sde_intercepts)]
  
  # summarise residual sds and correlations
  sd_residual <- NULL
  cor_residual <- NULL
  if (any(duplicated(object$data[,get(object$id)]))) { # Note: modification of original code, was not working without 'get'
    # sd parameters
    sd_residual <- s[stringr::str_starts(s$variable, "sigma_residual"),]
    rownames(sd_residual) <- paste0(
      "sd(",
      names(object$variables)[readr::parse_number(sd_residual$variable)],
      ")"
    )
    sd_residual <- sd_residual[, 2:ncol(sd_residual)]
    # correlation parameters
    cor_residual <- s[stringr::str_starts(s$variable, "cor_residual"),]
    for (i in 1:length(object$variables)) {
      for (j in 1:length(object$variables)) {
        if (i >= j) {
          var <- paste0("cor_residual[", i, ",", j, "]")
          cor_residual <- cor_residual[cor_residual$variable != var,]
        }
      }
    }
    rownames(cor_residual) <- paste0(
      "cor(",
      names(object$variables)[
        readr::parse_number(
          stringr::str_extract(
            cor_residual$variable,
            pattern = "cor\\_residual\\[(\\d+)\\,"
          )
        )
      ],
      ",",
      names(object$variables)[
        readr::parse_number(
          stringr::str_extract(
            cor_residual$variable,
            pattern = "\\,(\\d+)\\]"
          )
        )
      ],
      ")"
    )
    cor_residual <- cor_residual[, 2:ncol(cor_residual)]
  }
  
  # Function to add section headers as rows
  add_section_header <- function(df, section_name) {
    new <- rbind(NA, rownames_to_column(df, var='Parameter'))
    new$Type<- NA  # Create Section column
    new$Type[1] <- section_name  # Assign section name only to the first row
    new <- new[, c("Type", setdiff(names(new), "Type"))] # reorder columns
    return(new)
  }
  
  # combine drift parameters as in Coevolve summary format
  drift <- rbind(sd_drift, cor_drift)
  
  # Apply section headers to each dataset
  auto <- add_section_header(auto, "Autoregressive selection effects:")
  cross <- add_section_header(cross, "Cross selection effects:")
  drift <- add_section_header(drift, "Drift parameters:")
  sde_intercepts <- add_section_header(sde_intercepts, "Continuous time intercept parameters:")
  cor_residual <- add_section_header(cor_residual, "Residual parameters:")
  
  # Combine all tables
  all_tables <- bind_rows(auto, cross, drift, sde_intercepts, cor_residual)
  
  all_tables <- rbind(auto, cross, drift, sde_intercepts, cor_residual)
  
  # Generate a nicely formatted table
  nice_table(all_tables)
  
}
# object <- tar_read(coev_Q)
# prob <- 0.9
# robust <- FALSE




# Coevolutionary flowfield plot (fixed version) -------------------------------------------
custom_coev_plot_flowfield <- function(object, var1, var2, nullclines,
                                       limits, var1_lab, var2_lab) {

  # Check for required package
  if (!requireNamespace("ggplotify", quietly = TRUE)) {
    stop("Package 'ggplotify' is needed. Please install it with install.packages('ggplotify')")
  }

  # Use ggplotify to convert the base R plot to a ggplot object
  p <- ggplotify::as.ggplot(function() {

    # get IDs for variables
    id_var1 <- which(names(object$variables) == var1)
    id_var2 <- which(names(object$variables) == var2)
    # get posterior draws
    draws <- posterior::as_draws_rvars(object$fit)
    # medians and median absolute deviations for all variables
    eta  <- apply(
      draws$eta[,1:object$stan_data$N_tips,], 3, posterior::rvar_median
    )
    meds <- unlist(lapply(eta, stats::median))
    mads <- unlist(lapply(eta, stats::mad))
    lowers <- meds + limits[1]*mads
    uppers <- meds + limits[2]*mads
    # get median parameter values for A and b
    A <- stats::median(draws$A)
    b <- stats::median(draws$b)
    # function for flow field diagram
    OU <- function(t, y, parameters) {
      dy <- numeric(2)
      # variable 1
      dy[1] <- b[id_var1]
      for (j in 1:length(names(object$variables))) {
        if (j == id_var1) {
          # autoregressive effect
          dy[1] <- dy[1] + A[id_var1,j]*y[1]
        } else if (j == id_var2) {
          # cross-lagged effect of predictor
          dy[1] <- dy[1] + A[id_var1,j]*y[2]
        } else {
          # cross-lagged effects of other variables held at their median values
          dy[1] <- dy[1] + A[id_var1,j]*meds[j]
        }
      }
      # variable 2
      dy[2] <- b[id_var2]
      for (j in 1:length(names(object$variables))) {
        if (j == id_var2) {
          # autoregressive effect
          dy[2] <- dy[2] + A[id_var2,j]*y[2]
        } else if (j == id_var1) {
          # cross-lagged effect of predictor
          dy[2] <- dy[2] + A[id_var2,j]*y[1]
        } else {
          # cross-lagged effects of other variables held at their median values
          dy[2] <- dy[2] + A[id_var2,j]*meds[j]
        }
      }
      return(list(dy))
    }
    # create flow field diagram
    suppressWarnings({
      OU.flowField <-
        phaseR::flowField(
          OU,
          xlim = c(lowers[id_var1], uppers[id_var1]),
          ylim = c(lowers[id_var2], uppers[id_var2]),
          parameters = NA,
          add = FALSE,
          xlab = var1_lab,  # Set x-axis label here
          ylab = var2_lab,  # Set y-axis label here
          points = 12,
          col = "black",
          xaxt = 'n',
          yaxt = 'n',
          arrow.type = "proportional",
          frac = 1.5,
          xaxs = "i",
          yaxs = "i",
          axes = FALSE,
          lwd = 1.5
        )
    })

    # # Add axis labels directly using mtext
    # graphics::mtext(
    #   side = 1,
    #   text = var1_lab,
    #   line = 2.5,
    #   cex = 1.3
    # )
    # 
    # graphics::mtext(
    #   side = 2,
    #   text = var2_lab,
    #   line = 2.5,
    #   cex = 1.3
    # )

    # add nullclines to phase plane
    suppressWarnings({
      if (nullclines) {
        nc <-
          phaseR::nullclines(
            OU,
            xlim = c(lowers[id_var1], uppers[id_var1]),
            ylim = c(lowers[id_var2], uppers[id_var2]),
            parameters = NA,
            points = 20,
            axes = FALSE,
            col = c("#c55852","#5387b6"),
            add.legend = FALSE,
            lwd = 2
          )
      }
    })
    # add axes
    graphics::axis(
      side = 1,
      at = c(lowers[id_var1], meds[id_var1], uppers[id_var1]),
      labels =
        (c(lowers[id_var1], meds[id_var1], uppers[id_var1]) - meds[id_var1]) /
        mads[id_var1]
    )
    graphics::axis(
      side = 2,
      at = c(lowers[id_var2], meds[id_var2], uppers[id_var2]),
      labels =
        (c(lowers[id_var2], meds[id_var2], uppers[id_var2]) - meds[id_var2]) /
        mads[id_var2]
    )

    # Add bounding box to the plot
    # graphics::box()
  })

  # No need for these ggplot annotations as we're using mtext now
  # Remove these lines:
  # p <- p + annotate("text", x = 0.525, y = -Inf, label = paste(var1_lab), vjust = -1, size = 4) +
  #   annotate("text", x = -Inf, y = 0.525, label = paste(var2_lab), hjust = -1, size = 4)

  return(p)
}




# Extract posterior summary for a single fixed effect ----------------------
extract_coef_summary <- function(fit, coef) {
  
  draws <- as_draws_df(fit)[[coef]]
  if (is.null(draws)) stop(paste0("Coefficient '", coef, "' not found in model (check b_scale naming)."))
  
  data.table(
    median = median(draws),
    ll50   = as.numeric(quantile(draws, 0.25)),
    ul50   = as.numeric(quantile(draws, 0.75)),
    ll90   = as.numeric(quantile(draws, 0.05)),
    ul90   = as.numeric(quantile(draws, 0.95))
  )
  
}
# fit <- tar_read(m1q)
# coef <- 'b_scalelifespan_Post.Mean_F'



# Build table of coefficients across traits x robustness variants ----------
robustness_coef_table <- function(trait_vars, trait_labels, baseline_models, ns_models, dur_models) {
  
  variant_labels <- c('Baseline', '+ Network size', '+ Duration')
  
  out <- list()
  for (i in seq_along(trait_vars)) {
    
    coef_name <- paste0('b_scale', trait_vars[i])
    fits <- list(baseline_models[[i]], ns_models[[i]], dur_models[[i]])
    
    for (j in seq_along(fits)) {
      s <- extract_coef_summary(fits[[j]], coef_name)
      s[, trait := trait_labels[i]]
      s[, variant := variant_labels[j]]
      out[[length(out) + 1]] <- s
    }
  }
  
  d <- rbindlist(out)
  
  # Trait order top-to-bottom in facets = order given in trait_labels
  d[, trait := factor(trait, levels = trait_labels)]
  # Variant order bottom-to-top within a facet: Duration, Network size, Baseline (so Baseline ends up on top)
  d[, variant := factor(variant, levels = rev(variant_labels))]
  
  return(d)
  
}
# trait_vars <- c('lifespan_Post.Mean_F', 'age.mat_F')
# trait_labels <- c('Lifespan', 'Age at maturity')
# baseline_models <- list(tar_read(m1q), tar_read(m2q))
# ns_models <- list(tar_read(ns_m1q), tar_read(ns_m2q))
# dur_models <- list(tar_read(dur_m1q), tar_read(dur_m2q))



# Helper: symmetric x-axis limits centered on zero --------------------------
symmetric_x_limits <- function(d, pad = 1.05) {
  m <- max(abs(c(d$ll90, d$ul90)), na.rm = TRUE) * pad
  c(-m, m)
}



# Forest plot of coefficients across robustness checks -----------------------
plot_robustness_forest <- function(trait_vars, trait_labels, baseline_models, ns_models, dur_models,
                                   x_lab = 'Standardized effect', title = NULL) {
  
  d <- robustness_coef_table(trait_vars, trait_labels, baseline_models, ns_models, dur_models)
  
  cols <- viridis(3, option = proj_color, begin = 0.15, end = 0.85)
  names(cols) <- c('Baseline', '+ Network size', '+ Duration')
  
  ggplot(d, aes(y = variant, color = variant)) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey50') +
    geom_errorbarh(aes(xmin = ll90, xmax = ul90), height = 0, linewidth = 0.6) +
    geom_errorbarh(aes(xmin = ll50, xmax = ul50), height = 0, linewidth = 2.2) +
    geom_point(aes(x = median), size = 2.3, color = 'black') +
    scale_color_manual(values = cols, breaks = c('Baseline', '+ Network size', '+ Duration')) +
    scale_x_continuous(limits = symmetric_x_limits(d)) +
    facet_grid(trait ~ ., scales = 'free_y', space = 'free', switch = 'y') +
    labs(x = x_lab, y = NULL, title = title)+
    theme_classic() +
    theme(strip.placement = 'outside',
          strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0, face = 'bold', size = 11),
          axis.text.y = element_text(size = 10),
          panel.spacing = unit(1.2, 'lines'),
          legend.position = 'none',
          plot.caption = element_text(size = 8, color = 'grey40'))
  
}
# trait_vars <- c('lifespan_Post.Mean_F', 'age.mat_F', 'log_length_F', 'SSD')
# trait_labels <- c('Lifespan', 'Age at maturity', 'Body length (log)', 'SSD')
# baseline_models <- list(tar_read(m1q), tar_read(m2q), tar_read(m3q), tar_read(m4q))
# ns_models <- list(tar_read(ns_m1q), tar_read(ns_m2q), tar_read(ns_m3q), tar_read(ns_m4q))
# dur_models <- list(tar_read(dur_m1q), tar_read(dur_m2q), tar_read(dur_m3q), tar_read(dur_m4q))
# plot_robustness_forest(trait_vars, trait_labels, baseline_models, ns_models, dur_models, title = 'Q')




# Build table of control-variable coefficients across robustness models -----
control_coef_table <- function(control_labels, control_coefs, model_lists, row_labels) {
  
  out <- list()
  for (i in seq_along(control_labels)) {
    
    fits <- model_lists[[i]]
    for (j in seq_along(fits)) {
      s <- extract_coef_summary(fits[[j]], control_coefs[i])
      s[, control := control_labels[i]]
      s[, trait := row_labels[j]]
      out[[length(out) + 1]] <- s
    }
  }
  
  d <- rbindlist(out)
  
  # Control-variable blocks ordered top-to-bottom as given
  d[, control := factor(control, levels = control_labels)]
  # Within a block, first row_label ends up on top
  d[, trait := factor(trait, levels = rev(row_labels))]
  
  return(d)
  
}
# control_labels <- c('Network size', 'Duration')
# control_coefs  <- c('b_scalenetworkSize', 'b_scalestudyDuration')
# row_labels     <- c('Lifespan', 'Age at maturity', 'Body length', 'SSD')
# model_lists <- list(list(tar_read(ns_m1q), tar_read(ns_m2q), tar_read(ns_m3q), tar_read(ns_m4q)),
#                      list(tar_read(dur_m1q), tar_read(dur_m2q), tar_read(dur_m3q), tar_read(dur_m4q)))



# Forest plot of control-variable effects across robustness models ----------
plot_control_forest <- function(control_labels, control_coefs, model_lists, row_labels,
                                x_lab = 'Standardized effect', title = NULL) {
  
  d <- control_coef_table(control_labels, control_coefs, model_lists, row_labels)
  
  cols <- viridis(length(row_labels), option = proj_color, begin = 0.1, end = 0.9)
  names(cols) <- row_labels
  
  ggplot(d, aes(y = trait, color = trait)) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey50') +
    geom_errorbarh(aes(xmin = ll90, xmax = ul90), height = 0, linewidth = 0.6) +
    geom_errorbarh(aes(xmin = ll50, xmax = ul50), height = 0, linewidth = 2.2) +
    geom_point(aes(x = median), size = 2.3, color = 'black') +
    scale_color_manual(values = cols, breaks = row_labels) +
    scale_x_continuous(limits = symmetric_x_limits(d)) +
    facet_grid(control ~ ., scales = 'free_y', space = 'free', switch = 'y') +
    labs(x = x_lab, y = NULL, title = title)+
    theme_classic() +
    theme(strip.placement = 'outside',
          strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0, face = 'bold', size = 11),
          axis.text.y = element_text(size = 10),
          panel.spacing = unit(1.2, 'lines'),
          legend.position = 'none',
          plot.caption = element_text(size = 8, color = 'grey40'))
  
}


# Combined Q + S panel -------------------------------------------------------
plot_control_lifehistory_QS <- function(control_labels, control_coefs, row_labels,
                                        model_lists_Q, model_lists_S) {
  
  pQ <- plot_control_forest(control_labels, control_coefs, model_lists_Q, row_labels, title = 'Q (Modularity)')
  pS <- plot_control_forest(control_labels, control_coefs, model_lists_S, row_labels, title = 'S (Social differentiation)')
  
  pQ | pS
  
}




print('Cleared functions')
