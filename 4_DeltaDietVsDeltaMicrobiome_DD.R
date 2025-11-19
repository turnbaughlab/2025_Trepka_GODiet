# Load functions and directories
prefix <- "/mnt/tank/labmainshare/qb3share/" # "/Volumes/turnbaughlab/qb3share/"
dir <- paste0(prefix, "ktrepka/GODiet_v2/")
source(paste0(dir, "0_Functions.R"))

# Load diet data and normalize
fn <- paste0(dir, "Tables/TotalsMergedSequencedFilteredAverage_WithMetadata.csv")
t <- read.csv(fn)
t_norm <- normalize_all_to_kcal(t)
t_norm$Sample_ID = gsub("Final", "EOT", t_norm$Sample_ID)

t_filt = t_norm 
pt_keep = t_filt %>% group_by(PatientID) %>% dplyr::summarise(n = n()) %>% filter(n >= 2) %>% dplyr::select(PatientID) %>% pull()
t_filt = t_filt %>% filter(PatientID %in% pt_keep) #%>% filter(Cycle != "EOT") # exclude variable end-of-treatment timepoint

# Set variables to test
all_vars = c(colnames(t_filt)[diet_index_start:diet_index_end], colnames(t_filt)[hei_index_start:hei_index_end])
all_vars = setdiff(all_vars, feature_exclude_dd)
#vars = c("HEI", "AHEI", "G_REFINED", "PF_NUTSDS", "VK", "THEO", "ATOC", "CARB", "TFAT", "D_YOGURT")
vars = all_vars
min_clr = 0
rerun_diversity = FALSE
rerun_taxa = TRUE
rerun_fc = FALSE

if (rerun_fc){
  print("alpha_fc_shannon")
  alpha_dd(taxLevel = "ASVs", div_index = "shannon", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "Species", div_index = "shannon", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "Genus", div_index = "shannon", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "KO", div_index = "shannon", nsamp = 3, cutoff = 0.00001, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "Pathway", div_index = "shannon", nsamp = 3, cutoff = 0.00001, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "ASVs", div_index = "shannon", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "Species", div_index = "shannon", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "Genus", div_index = "shannon", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "KO", div_index = "shannon", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "Pathway", div_index = "shannon", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  
  print("alpha_fc_observed")
  alpha_dd(taxLevel = "ASVs", div_index = "observed", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "Species", div_index = "observed", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "Genus", div_index = "observed", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "KO", div_index = "observed", nsamp = 3, cutoff = 0.00001, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "Pathway", div_index = "observed", nsamp = 3, cutoff = 0.00001, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "ASVs", div_index = "observed", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "Species", div_index = "observed", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "Genus", div_index = "observed", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "KO", div_index = "observed", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd(taxLevel = "Pathway", div_index = "observed", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
}

if (rerun_diversity){
  print("alpha shannon")
  # Alpha diversity calculations
  alpha_dd_rank(taxLevel = "ASVs", div_index = "shannon", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "Species", div_index = "shannon", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "Genus", div_index = "shannon", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "KO", div_index = "shannon", nsamp = 3, cutoff = 0.00001, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "Pathway", div_index = "shannon", nsamp = 3, cutoff = 0.00001, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "ASVs", div_index = "shannon", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "Species", div_index = "shannon", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "Genus", div_index = "shannon", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "KO", div_index = "shannon", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "Pathway", div_index = "shannon", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  
  print("alpha observed")
  alpha_dd_rank(taxLevel = "ASVs", div_index = "observed", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "Species", div_index = "observed", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "Genus", div_index = "observed", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "KO", div_index = "observed", nsamp = 3, cutoff = 0.00001, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "Pathway", div_index = "observed", nsamp = 3, cutoff = 0.00001, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "ASVs", div_index = "observed", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "Species", div_index = "observed", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "Genus", div_index = "observed", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "KO", div_index = "observed", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  alpha_dd_rank(taxLevel = "Pathway", div_index = "observed", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  
  print("beta_fc_and_rank")
  #### Beta diversity calculations
  beta_dd_rank(taxLevel = "ASVs", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  beta_dd_rank(taxLevel = "Species", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  beta_dd_rank(taxLevel = "Genus", nsamp = 3, cutoff = 0.0005, vars = vars, t_filt = t_filt)
  beta_dd_rank(taxLevel = "KO", nsamp = 3, cutoff = 0.00001, vars = vars, t_filt = t_filt)
  beta_dd_rank(taxLevel = "Pathway", nsamp = 3, cutoff = 0.00001, vars = vars, t_filt = t_filt)
  beta_dd_rank(taxLevel = "ASVs", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  beta_dd_rank(taxLevel = "Species", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  beta_dd_rank(taxLevel = "Genus", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  beta_dd_rank(taxLevel = "KO", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
  beta_dd_rank(taxLevel = "Pathway", nsamp = 3, cutoff = 0, vars = vars, t_filt = t_filt)
}

if (rerun_taxa){
  ## every KO/pathway/taxa/etc
  print(1)
  feature_dd_rank(taxLevel = "Genus", nsamp = 25, cutoff = 0.0005, vars = vars, t_filt = t_filt, min_clr = min_clr)
  feature_dd_rank(taxLevel = "Species", nsamp = 25, cutoff = 0.0005, vars = vars, t_filt = t_filt, min_clr = min_clr)
  feature_dd_rank(taxLevel = "ASVs", nsamp = 25, cutoff = 0.0005, vars = vars, t_filt = t_filt, min_clr = min_clr)
  print(2)
  feature_dd_rank(taxLevel = "Genus", nsamp = 25, cutoff = 0, vars = vars, t_filt = t_filt, min_clr = min_clr)
  print(3)
  feature_dd_rank(taxLevel = "Species", nsamp = 25, cutoff = 0, vars = vars, t_filt = t_filt, min_clr = min_clr)
  print(4)
  feature_dd_rank(taxLevel = "ASVs", nsamp = 25, cutoff = 0, vars = vars, t_filt = t_filt, min_clr = min_clr)
  print(5)
  feature_dd_rank(taxLevel = "Pathway", nsamp = 25, cutoff = 0.00001, vars = vars, t_filt = t_filt, min_clr = min_clr)
  print(6)
  feature_dd_rank(taxLevel = "KO", nsamp = 25, cutoff = 0.00001, vars = vars, t_filt = t_filt, min_clr = min_clr)
  print(7)
  feature_dd_rank(taxLevel = "Pathway", nsamp = 25, cutoff = 0, vars = vars, t_filt = t_filt, min_clr = min_clr)
  print(8)
  feature_dd_rank(taxLevel = "KO", nsamp = 25, cutoff = 0, vars = vars, t_filt = t_filt, min_clr = min_clr)
}

