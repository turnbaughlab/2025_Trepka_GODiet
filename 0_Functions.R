# Libraries and dependencies
library(compositions)
library(ppcor) # partial correlations
library(reshape2)
library(tidyr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(nlme)
library(DescTools) # Cramer's V
library(polycor) 
library(vegan) # alpha
library(emmeans)
library(dplyr)


# Directories
prefix <- "/mnt/tank/labmainshare/qb3share/" # "/Volumes/turnbaughlab/qb3share/"
dir <- paste0(prefix, "ktrepka/GODiet_v2/")
rawdir <- paste0(dir, "RawData/")
figdir <- paste0(dir, "Figures/")
tabledir <- paste0(dir, "Tables/")
taxadir <- paste0(rawdir, "TaxaSummary/")
genedir <- paste0(rawdir, "GO_MGS/")
dir16s <- paste0(rawdir, "GO_16S/")
fig_meta <- paste0(figdir, "MetadataQC/")
fig_time <- paste0(figdir, "DietVsTime/")
table_time <- paste0(tabledir, "DietVsTime/")
fig_baseline <- paste0(figdir, "DietVsMicrobiomeBaseline/")
table_baseline <- paste0(tabledir, "DietVsMicrobiomeBaseline/")
fig_mt <- paste0(figdir, "DietVsMicrobiomeTime/")
table_mt <- paste0(tabledir, "DietVsMicrobiomeTime/")
table_mt_alpha <- paste0(tabledir, "DietVsMicrobiomeTime/Alpha/")
table_mt_beta <- paste0(tabledir, "DietVsMicrobiomeTime/Beta/")
table_mt_granular <- paste0(tabledir, "DietVsMicrobiomeTime/IndividualTaxaGenes/")
fig_mt_alpha <- paste0(figdir, "DietVsMicrobiomeTime/Alpha/")
fig_mt_beta <- paste0(figdir, "DietVsMicrobiomeTime/Beta/")
fig_mt_granular <- paste0(figdir, "DietVsMicrobiomeTime/IndividualTaxaGenes/")

# Create directories if they don't exist
for (d in c(figdir, tabledir, fig_meta, fig_time, table_time, fig_baseline, table_baseline,
            fig_mt, table_mt, table_mt_alpha, table_mt_beta, table_mt_granular, fig_mt_alpha,
            fig_mt_beta, fig_mt_granular)) {
  suppressWarnings(dir.create(d))
}

# Metadata location
fn_meta <- paste0(rawdir, "Metadata/Metadata250319.csv")

# Other global variables
diet_index_start = 4
diet_index_end = 105 
hei_index_start = 166
hei_index_end = 180 
metadata_index_start= 136
metadata_index_end = 162
duplicate_features <- c("V_LEGUMES", "LZ", "P182", "M181", "A_DRINKS", "V_REDOR_TOTAL", "FDFE", "V_STARCHY_POTATO", 
                        "SOLID_FATS", "S080", "S100", "S120", "S140", "CHOLN", "P204", "P205", "P225", "P184",
                        "VITE_ADD", "V_STARCHY_OTHER", "PF_ORGAN", "PF_SEAFD_HI", "PF_SEAFD_LOW", "PF_SOY") # duplicates + median of food group = 0; kept ALC
feature_exclude_dd = unique(c(duplicate_features))
metadata_exclude = c("Metadata_Presence_of_Primary_Tumor_.Yes.No.", "Metadata_Presence_of_Primary_Tumor", "Metadata_Procedure")

# Global palettes
enriched_palette <- c(
  "Enriched" = "blue", 
  "Depleted" = "orange", 
  "Not Significant (FDR > 0.2)" = "gray",
  "FDR > 0.2" = "gray",
  "Not Significant" = "gray"
)

# Functions
clr_transform <- function(data){
  pseudo_count = min(data[data > 0], na.rm = TRUE)/2 # Define pseudocount
  data_clr <- t(compositions::clr(data + pseudo_count)) %>% as.data.frame() # CLR transform
  data_clr <- data_clr[, sapply(data_clr, function(x) var(x, na.rm = TRUE) > 0)] # Remove zero-variance columns
  data_clr[data_clr == 0] <- min(data_clr, na.rm = TRUE) # Set 0s to minimum value
  data_clr[is.na(data_clr)] <- min(data_clr, na.rm = TRUE) # Set NA to minimum value
  return(data_clr)
}

load_seq_metadata <- function() {
  seq <- read.csv(fn_meta) %>% 
    dplyr::rename(PatientID = Patient_ID, Cycle = Treatment_Cycle)
  seq$Cycle[seq$Cycle == "Final"] <- "EOT"
  seq$PatientID[seq$PatientID == "9b"] <- "9"
  return(seq)
}

# method = "per1000" or "residuals"
normalize_column_to_kcal <- function(data, nutrient_col, method = "residuals") {
  if (nutrient_col %in% c("PROT", "CARB")){
    data[[nutrient_col]] <- 4*data[[nutrient_col]]/data[["KCAL"]]*100 # percentage
  } else if (nutrient_col == "TFAT"){
    data[[nutrient_col]] <- 9*data[[nutrient_col]]/data[["KCAL"]]*100 # percentage
  }
  else if (nutrient_col != "KCAL") {
    if (method == "residuals") {
      
      # Identify zero vs non-zero entries
      nonzero_idx <- data[[nutrient_col]] != 0
      
      # Subset for non-zero values
      y_nonzero <- data[[nutrient_col]][nonzero_idx]
      x_nonzero <- data[["KCAL"]][nonzero_idx]
      
      # First residuals approach
      mod <- lm(y_nonzero ~ x_nonzero)
      res <- as.numeric(mod$residuals) + mean(y_nonzero, na.rm = TRUE)
      
      # Check if any residuals are negative
      if (all(res >= 0)) {
        data[[nutrient_col]][nonzero_idx] <- res
      } else {
        # Log-transform fallback for non-zero entries
        y_log <- log(y_nonzero + min(y_nonzero[y_nonzero > 0]))
        x_log <- log(x_nonzero)
        mod <- lm(y_log ~ x_log)
        res_log <- exp(as.numeric(mod$residuals) + mean(y_log, na.rm = TRUE))
        data[[nutrient_col]][nonzero_idx] <- res_log
      }
      
      # Leave zeros untouched
      data[[nutrient_col]][!nonzero_idx] <- 0
    } else {
      data[[nutrient_col]] <- data[[nutrient_col]]/data[["KCAL"]]*1000 # per 1000 KCAL
    }
  }
  return(data)
}

read_file <- function(folder, dir, type = "Totals"){
  fold <- paste0(dir, "ASA24DataExports/", folder, "/")
  files <- list.files(fold)
  fn_totals <- files[grepl(type, files)]
  totals <- read.csv(paste0(fold, fn_totals))
  return(totals)
}

add_hei <- function(fn){
  t <- read.csv(file = fn)
  hei <- t
  
  hei$FWHOLEFRT <- hei$F_CITMLB + hei$F_OTHER 
  hei$MONOPOLY = hei$MFAT + hei$PFAT
  hei$VTOTALLEG = hei$V_TOTAL + hei$V_LEGUMES
  hei$VDRKGRLEG = hei$V_DRKGR + hei$V_LEGUMES
  hei$PFALLPROTLEG = hei$PF_MPS_TOTAL + hei$PF_EGGS + hei$PF_NUTSDS + hei$PF_SOY + hei$PF_LEGUMES
  hei$PFSEAPLANTLEG = hei$PF_SEAFD_HI + hei$PF_SEAFD_LOW + hei$PF_NUTSDS + hei$PF_SOY + hei$PF_LEGUMES
  
  totalFruits = hei$F_TOTAL/(hei$KCAL/1000)
  totalFruits_score = case_when(totalFruits < 0.8 ~ 5*totalFruits/0.8,
                                TRUE ~ 5)
  
  wholeFruits = hei$FWHOLEFRT/(hei$KCAL/1000)
  wholeFruits_score = case_when(wholeFruits < 0.4 ~ 5*wholeFruits/0.4,
                                TRUE ~ 5)
  
  totalVeg = hei$VTOTALLEG/(hei$KCAL/1000)
  totalVeg_score = case_when(totalVeg < 1.1 ~ 5*totalVeg/1.1,
                             TRUE ~ 5)
  
  greensBeans = hei$VDRKGRLEG/(hei$KCAL/1000)
  greensBeans_score = case_when(greensBeans < 0.2 ~ 5*greensBeans/0.2,
                                TRUE ~ 5)
  
  wholeGrains = hei$G_WHOLE/(hei$KCAL/1000)
  wholeGrains_score = case_when(wholeGrains < 1.5 ~ 10*wholeGrains/1.5,
                                TRUE ~ 10)
  
  dairy = hei$D_TOTAL/(hei$KCAL/1000)
  dairy_score = case_when(dairy < 1.3 ~ 10*dairy/1.3,
                          TRUE ~ 10)
  
  totalProtein = hei$PFALLPROTLEG/(hei$KCAL/1000)
  totalProtein_score = case_when(totalProtein < 2.5 ~ 5*totalProtein/2.5,
                                 TRUE ~ 5)
  
  seaPlantProtein = hei$PFSEAPLANTLEG/(hei$KCAL/1000)
  seaPlantProtein_score = case_when(seaPlantProtein < .8 ~ 5*seaPlantProtein/.8,
                                    TRUE ~ 5)
  
  refinedGrains = hei$G_REFINED/(hei$KCAL/1000) 
  refinedGrains_score = case_when(refinedGrains < 1.8 ~ 10,
                                  refinedGrains > 4.3 ~ 0,
                                  TRUE ~ 10*(4.3 - refinedGrains)/(4.3 - 1.8))
  
  addedSugar = hei$ADD_SUGARS*16/(hei$KCAL)*100 # percent of energy from sugar
  addedSugar_score = case_when(addedSugar < 6.5 ~ 10,
                               addedSugar > 26 ~ 0,
                               TRUE ~ 10*(26 - addedSugar)/(26 - 6.5)) 
  
  sodium = (hei$SODI/1000)/(hei$KCAL/1000) # convert to grams, normalize to diet
  sodium_score = case_when(sodium < 1.1 ~ 10,
                           sodium > 2 ~ 0,
                           TRUE ~ 10*(2 - sodium)/(2 - 1.1)) 
  
  satFats = (hei$SFAT*9)/(hei$KCAL)*100 # percent of energy from saturated fat
  satFats_score = case_when(satFats < 8 ~ 10,
                            satFats > 16 ~ 0,
                            TRUE ~ 10*(16 - satFats)/(16 - 8)) 
  
  fattyAcids = hei$MONOPOLY/hei$SFAT # ratio
  fattyAcids_score = case_when(fattyAcids > 2.5 ~ 10,
                               fattyAcids < 1.2 ~ 0,
                               TRUE ~ 10*(fattyAcids - 1.2)/(2.5 - 1.2))  
  
  hei_score = totalFruits_score + wholeFruits_score + totalVeg_score + greensBeans_score + wholeGrains_score + dairy_score + totalProtein_score + seaPlantProtein_score + refinedGrains_score + addedSugar_score + sodium_score + satFats_score + fattyAcids_score
  
  t$HEI = hei_score
  t$HEI_totalFruits_score = totalFruits_score
  t$HEI_wholeFruits_score = wholeFruits_score
  t$HEI_totalVeg_score = totalVeg_score
  t$HEI_greensBeans_score = greensBeans_score
  t$HEI_wholeGrains_score = wholeGrains_score
  t$HEI_dairy_score = dairy_score
  t$HEI_totalProtein_score = totalProtein_score
  t$HEI_seaPlantProtein_score = seaPlantProtein_score
  t$HEI_refinedGrains_score = refinedGrains_score
  t$HEI_addedSugar_score = addedSugar_score
  t$HEI_sodium_score = sodium_score
  t$HEI_satFats_score = satFats_score
  t$HEI_fattyAcids_score = fattyAcids_score
  write.csv(t, file = fn, row.names = FALSE)
}

add_ahei <- function(fn, return_df = FALSE){
  t <- read.csv(fn)
  hei <- t
  
  totalVeg = (hei$V_TOTAL + hei$V_LEGUMES - hei$V_STARCHY_POTATO) / (1/2) # convert to servings. 1/2 cup per serving
  totalVeg_score = case_when(totalVeg < 5 ~ 10*totalVeg/5,
                             TRUE ~ 10)
  
  totalFruit = (hei$F_CITMLB + hei$F_OTHER) / (1/2) # convert to servings. 1/2 cup per serving
  totalFruit_score = case_when(totalFruit < 4 ~ 10*totalFruit/5,
                               TRUE ~ 10)
  
  wholeGrains = hei$G_WHOLE * 16 # convert to grams. 16g/oz
  wholeGrains_score = case_when(hei$Metadata_Sex == "Male" & wholeGrains < 90 ~ 10*wholeGrains/90,
                                hei$Metadata_Sex == "Male" & wholeGrains >= 90 ~ 10,
                                hei$Metadata_Sex == "Female" & wholeGrains < 75 ~ 10*wholeGrains/75,
                                hei$Metadata_Sex == "Female" & wholeGrains >= 75 ~ 10,
                                TRUE ~ 0)
  
  sweetening = hei$F_JUICE + hei$ADD_SUGARS/6.6 # Assuming added sugars are in soda, and using Coke (10 tsp/1.5 cup, or 6.6 tsp/cup) to calculate cup equivalents
  sweetening_score = case_when(sweetening < 1 ~ 10*(1-sweetening),
                               TRUE ~ 0)
  
  nutsLeg = hei$PF_NUTSDS + hei$PF_LEGUMES # 1 serving = 1 oz
  nutsLeg_score = case_when(nutsLeg < 1 ~ 10*nutsLeg,
                            TRUE ~ 10)
  
  meat = hei$PF_MEAT/4 + hei$PF_ORGAN/4 + hei$PF_CUREDMEAT/1.5 # servings of red meat + cured meats
  meat_score = case_when(meat < 1.5 ~ 10*(1 - meat/1.5),
                         TRUE ~ 0)
  
  longChain = (hei$P205 + hei$P225 + hei$P226)*1000 # mg
  longChain_score = case_when(longChain < 250 ~ 10*longChain/250,
                              TRUE ~ 10)
  
  pufa = hei$PFAT*9/hei$KCAL * 100 # percent of energy
  pufa_score = case_when(pufa < 2 ~ 0,
                         pufa > 10 ~ 10,
                         TRUE ~ 10*(pufa - 2)/(10 - 2)) 
  
  sodium = hei$SODI # sodium, mg
  sodi_min = quantile(sodium, 0.1); sodi_max = quantile(sodium, 0.9)
  sodium_score = case_when(sodium < sodi_min ~ 10,
                           sodium > sodi_max ~ 0,
                           TRUE ~ 10*(sodi_max - sodium)/(sodi_max - sodi_min)) 
  
  alcohol = hei$A_DRINKS # number of drinks
  alcohol_score = case_when(hei$Metadata_Sex == "Male" & alcohol < 1 ~ 10,
                            hei$Metadata_Sex == "Male" & alcohol > 3.5 ~ 0,
                            hei$Metadata_Sex == "Male" ~ 10*(3.5 - alcohol)/(3.5 - 1),
                            hei$Metadata_Sex == "Female" & alcohol < 1 ~ 10,
                            hei$Metadata_Sex == "Female" & alcohol > 2.5 ~ 0,
                            hei$Metadata_Sex == "Female" ~ 10*(2.5 - alcohol)/(2.5 - 1),
                            TRUE ~ 0)
  
  ahei_score = totalVeg_score + totalFruit_score + wholeGrains_score + sweetening_score + nutsLeg_score + meat_score + longChain_score + pufa_score + sodium_score + alcohol_score
  
  if (return_df){
    t$AHEI = ahei_score
    t$AHEI_wholeGrains_score = wholeGrains_score
    t$AHEI_longChain_score = longChain_score
    t$AHEI_pufa_score = pufa_score 
    t$AHEI_sweetening_score = sweetening_score 
    t$AHEI_nutsLeg_score = nutsLeg_score 
    return(t)
  } else {
    write.csv(t, file = fn, row.names = FALSE)
  }
}

normalize_all_to_kcal <- function(t){
  t_norm = t
  fn = paste0(tabledir, "DietCategory.csv")
  cat_df <- read.csv(file = fn) %>% as.data.frame()
  norm_col <- cat_df %>% dplyr::filter(Type %in% c("Nutrient", "FoodGroup")) %>% 
    dplyr::select(Diet) %>% pull() %>% unique() 
  for (col in norm_col){
    t_norm <- normalize_column_to_kcal(t_norm, nutrient_col = col)
  }
  return(t_norm)
}

keep_k_pts <- function(t_norm, k = 2){
  pt_keep <- t_norm %>% dplyr::group_by(PatientID) %>% dplyr::summarise(n = n()) %>% 
    dplyr::filter(n >= k) %>% 
    dplyr::select(PatientID) %>% pull() 
  return(t_norm %>% dplyr::filter(PatientID %in% pt_keep))
}

fit_time_model <- function(dat, timevar = "Linear", logT = FALSE){
  fn = paste0(tabledir, "DietCategory.csv")
  cat_df <- read.csv(file = fn) %>% as.data.frame()
  
  if (logT & timevar != "Cat"){
    minval = min(dat[[timevar]][dat[[timevar]] > 0], na.rm = TRUE)/2
    dat$TimeVar = log(dat[[timevar]] + minval)
    } else {
    dat$TimeVar = dat[[timevar]]
  }
  
  res_all = NULL
  for (i in 1:length(cat_df$Diet)){
    var_name = cat_df$Diet[i]
    dat$Var <- dat[[var_name]]
    mod <- nlme::lme(Var ~ TimeVar, random = ~1|PatientID, data = dat)
    test <- summary(mod)
    res <- test$tTable %>% as.data.frame()
    res$Term = rownames(res)
    res$Variable = cat_df$Diet[i]
    res$Baseline = dat %>% dplyr::filter(Cycle == "C1D1") %>% dplyr::select(Var) %>% pull() %>% mean(na.rm = TRUE)
    res$Post = dat %>% dplyr::filter(Cycle != "C1D1") %>% dplyr::select(Var) %>% pull() %>% mean(na.rm = TRUE)
    res <- res %>% dplyr::filter(Term != "(Intercept)")
    if(length(res_all) == 0){res_all <- res}else{res_all <- rbind(res_all, res)}
  }
  res_all_corrected <- fdr_time(res_all)
  return(res_all_corrected)
}

# FDR correct, by category
fdr_time <- function(r){
  fn = paste0(tabledir, "DietCategory.csv")
  cat_df <- read.csv(file = fn) %>% as.data.frame()
  cat_df$Variable <- cat_df$Diet
  
  r <- r %>% dplyr::filter(!(Variable %in% feature_exclude_dd))
  
  r <- left_join(r, cat_df, by = "Variable")
  df_list <- split(r, r$Type)
  df_list <- lapply(df_list, function(sub_df) {
    sub_df$FDR <- p.adjust(sub_df$`p-value`, method = "BH")
    return(sub_df)
  })
  result_df <- bind_rows(df_list)
  return(result_df)
}

plot_hits <- function(box, hits, savedir = NULL){
  if (is.null(savedir)) savedir <- fig_time
  for (i in 1:length(hits$Diet)){
    hit = hits$Diet[i]
    p_value = hits$`p-value`[i]
    box$Var = box[[hit]]
    p <- ggplot(box, aes(x = Cycle, y = Var)) + 
      geom_point(aes(color = Cycle)) + 
      geom_boxplot(alpha = 0.2, color = "black", aes(fill = Cycle), outlier.size = -1) + 
      theme_pubr() + 
      ylab(hit) + 
      xlab("") + 
      annotate("text", x = 2, y = 1.1*max(box$Var), 
               label = paste("italic(p) == ", formatC(signif(p_value, digits=2), 
                                                      digits = 2, format="fg", flag="#")), parse = TRUE)
    if (hits$Value[i] < 0){
      p <- p + scale_color_manual(values = c("grey", "orange", "orange")) + 
        scale_fill_manual(values = c("grey", "orange", "orange")) + 
        theme(legend.position = "none") 
    } else {
      p <- p + scale_color_manual(values = c("grey", "blue", "blue")) + 
        scale_fill_manual(values = c("grey", "blue", "blue")) + 
        theme(legend.position = "none") 
    }
    fn <- paste0(savedir, "Boxplot", hit, ".pdf")
    ggsave(filename = fn, plot = p, width = 1.75, height = 3)  
  }
}

plot_hits_regimen <- function(box, hits){
  for (i in 1:length(hits$Diet)){
    hit = hits$Diet[i]
    box$Var = box[[hit]]
    p <- ggplot(box %>% dplyr::filter(!is.na(RegimenIntensityYesNo)), aes(x = Cycle, y = Var)) + 
      geom_point(aes(color = Cycle)) + 
      geom_boxplot(alpha = 0.2, color = "black", aes(fill = Cycle), outlier.size = -1) + 
      theme_pubr() + 
      ylab(hit) + 
      xlab("") + 
      facet_wrap(~RegimenIntensityYesNo)
    if (hits$Value[i] < 0){
      p <- p + scale_color_manual(values = c("grey", "orange", "orange")) + 
        scale_fill_manual(values = c("grey", "orange", "orange")) + 
        theme(legend.position = "none") 
    } else {
      p <- p + scale_color_manual(values = c("grey", "blue", "blue")) + 
        scale_fill_manual(values = c("grey", "blue", "blue")) + 
        theme(legend.position = "none") 
    }
    fn <- paste0(fig_time, "Boxplot", hit, "_RegimenFacet.pdf")
    ggsave(filename = fn, plot = p, width = 2, height = 3)  
  }
}

plot_medi_box_go <- function(medi, hits = c("PROT", "TFAT", "CARB"), norm = FALSE, 
                          cycles_stats = c("BL", "C1D1", "C1D3", "C1D7"), logT = FALSE,
                          xpos = 2.5){
  medi_box <- medi
  for (i in 1:length(hits)){
    hit = hits[i]
    medi_box$Cycle[medi_box$Cycle == "Baseline"] <- "BL" # easier plotting
    medi_box$Var = medi_box[[hit]]
    medi_box$Var[is.na(medi_box$Var)] <- 0 # fill in empty values with 0s for abundance
    medi_box <- medi_box %>% dplyr::filter(!is.na(Var))
    if (length(cycles_stats) == 4){
      medi_box$Linear <- case_when(medi_box$Cycle %in% c("BL", "C1D1") ~ 0,
                                   medi_box$Cycle %in% c("C1D3") ~ 3,
                                   TRUE ~ 7)
    } else {
      medi_box$Linear <- case_when(medi_box$Cycle %in% c("BL", "C1D1") ~ 0,
                                   TRUE ~ 1)
    }
    if (logT){
      medi_box$Var = log10(medi_box$Var + min(medi_box$Var[medi_box$Var > 0]))
    }
    test <- nlme::lme(Var ~ Linear, random = ~1|Patient_ID, data = medi_box %>% dplyr::filter(Cycle %in% cycles_stats))
    p_value <- summary(test)$tTable[2,5]
    p <- ggplot(medi_box, aes(x = Cycle, y = Var)) + 
      geom_point(aes(color = Cycle)) + 
      geom_boxplot(alpha = 0.2, color = "black", aes(fill = Cycle), outlier.size = -1) + 
      theme_pubr() + 
      ylab(hit) + 
      xlab("") + 
      scale_color_manual(values = c("grey", "grey", "orange", "orange", "orange", "orange", "orange")) + 
      scale_fill_manual(values = c("grey", "grey", "orange", "orange", "orange", "orange", "orange")) + 
      theme(legend.position = "none") +
      annotate("text", x = xpos, y = 1.1*max(medi_box$Var), 
               label = paste("italic(p) == ", formatC(signif(p_value, digits=2), 
                                                      digits = 2, format="fg", flag="#")), parse = TRUE) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

    fn <- paste0(fig_time, "MEDIGO", hit, ".pdf")
    ggsave(filename = fn, width = 3, height = 3)
  }
}

plot_medi_box_ne <- function(medi, hits = c("PROT", "TFAT", "CARB"), norm = FALSE, 
                             cycles_stats = c("T1", "T2", "T3"), logT = FALSE,
                             xpos = 2){
  medi_box <- medi
  for (i in 1:length(hits)){
    hit = hits[i]
    medi_box$Var = medi_box[[hit]]
    medi_box$Var[is.na(medi_box$Var)] <- 0 # fill in empty values with 0s for abundance
    medi_box <- medi_box %>% dplyr::filter(!is.na(Var))
    medi_box$Linear <- case_when(medi_box$Cycle %in% c("T1") ~ 0,
                                   medi_box$Cycle %in% c("T2") ~ 1,
                                   medi_box$Cycle %in% c("T3") ~ 1,
                                   TRUE ~ NA)
    if (logT){
      medi_box$Var = log10(medi_box$Var + min(medi_box$Var[medi_box$Var > 0]))
    }
    test <- nlme::lme(Var ~ Linear, random = ~1|Patient_ID, data = medi_box %>% dplyr::filter(Cycle %in% cycles_stats))
    p_value <- summary(test)$tTable[2,5]
    p <- ggplot(medi_box, aes(x = Cycle, y = Var)) + 
      geom_point(aes(color = Cycle)) + 
      geom_boxplot(alpha = 0.2, color = "black", aes(fill = Cycle), outlier.size = -1) + 
      theme_pubr() + 
      ylab(hit) + 
      xlab("") + 
      scale_color_manual(values = c("grey", "orange", "orange")) + 
      scale_fill_manual(values = c("grey", "orange", "orange")) + 
      theme(legend.position = "none") +
      annotate("text", x = xpos, y = 1.1*max(medi_box$Var), 
               label = paste("italic(p) == ", formatC(signif(p_value, digits=2), 
                                                      digits = 2, format="fg", flag="#")), parse = TRUE)
    
    fn <- paste0(fig_time, "MEDINE", hit, ".pdf")
    ggsave(filename = fn, width = 2, height = 3)
  }
}

mixed_corr_matrix <- function(data) {
  # Convert character columns to factors
  data[sapply(data, is.character)] <- lapply(data[sapply(data, is.character)], as.factor)
  
  n <- ncol(data)
  corr_matrix <- matrix(NA, n, n, dimnames = list(names(data), names(data)))
  
  for (i in 1:n) {
    for (j in 1:n) {
      var1 <- data[[i]]
      var2 <- data[[j]]
      
      # Skip if a variable has only one unique value (except NA)
      if (length(unique(na.omit(var1))) <= 1 || length(unique(na.omit(var2))) <= 1) {
        next  # Skip computation for constant columns
      }
      
      # Numeric-Numeric: Pearson correlation
      if (is.numeric(var1) && is.numeric(var2)) {
        corr_matrix[i, j] <- cor(var1, var2, method = "pearson", use = "complete.obs")
        
        # Categorical-Categorical: Cramér's V
      } else if (is.factor(var1) && is.factor(var2)) {
        corr_matrix[i, j] <- CramerV(table(var1, var2), bias.correct = TRUE)
        
        # Numeric-Categorical: Point-Biserial (if binary) or Polyserial (if multi-level)
      } else if (is.numeric(var1) && is.factor(var2)) {
        if (length(unique(na.omit(var2))) == 2) {  
          # Binary categorical -> Point-Biserial correlation
          corr_matrix[i, j] <- cor(var1, as.numeric(var2), method = "pearson", use = "complete.obs")
        } else {  
          # Multi-level categorical -> Polyserial correlation
          corr_matrix[i, j] <- polyserial(var1, as.numeric(var2))
        }
        
        # Categorical-Numeric: Same as above, but reversed
      } else if (is.factor(var1) && is.numeric(var2)) {
        if (length(unique(na.omit(var1))) == 2) {  
          # Binary categorical -> Point-Biserial correlation
          corr_matrix[i, j] <- cor(var2, as.numeric(var1), method = "pearson", use = "complete.obs")
        } else {  
          # Multi-level categorical -> Polyserial correlation
          corr_matrix[i, j] <- polyserial(var2, as.numeric(var1))
        }
      }
    }
  }
  
  return(as.data.frame(corr_matrix))
}


load_16s <- function(taxLevel = "Genus", cutoff = 0.0001, nsamp = 3, clr = TRUE, raw = FALSE,
                     min_clr = 0){
  fn <- paste0(rawdir, "GO_16S/TaxaSummary", taxLevel, ".csv")
  a <- read.csv(fn) %>% as.data.frame() %>% dplyr::rename(Taxa = X)
  colnames(a) <- gsub("\\.", "-", colnames(a)) # match formatting with metadata
  colnames(a) = gsub("Pt9b", "Pt9", colnames(a)) # matching with our data
  colnames(a) = gsub("Final", "EOT", colnames(a)) # matching with our data
  if(raw){
    return(a[-1])
  }
  a[-1] <- sweep(a[-1], 2, colSums(a[-1]), FUN = "/") # normalize to 1
  a_filt <- a[ rowSums(a[-1] >= cutoff) >= nsamp, ] # dplyr::filter: k relative abundance in n samples
  taxa_keep = a_filt$Taxa; rownames(a_filt) = a_filt$Taxa; a_filt <- a_filt %>% dplyr::select(-Taxa)
  rownames(a) <- a$Taxa; a <- a %>% dplyr::select(-Taxa)
  if(clr){
    a <- clr_transform(t(a)) %>% as.data.frame()
    a[a < min_clr] <- min_clr
    a_filt = a[taxa_keep,]
  }
  return(a_filt) 
}

load_mgs_taxa <- function(taxLevel = "Genus", cutoff = 0.001, nsamp = 3, clr = TRUE){
  fn <- paste0(rawdir, "GO_MGS/TaxaSummary", taxLevel, ".csv")
  m <- read.csv(fn) %>% as.data.frame() %>% dplyr::rename(Taxa = X)
  colnames(m) <- gsub("\\_", "-", colnames(m)) # match formatting with metadata
  colnames(m) = gsub("Pt9b", "Pt9", colnames(m)) # matching with our metadata
  colnames(m) = gsub("Final", "EOT", colnames(m)) # matching with our metadata
  m[-1] <- sweep(m[-1], 2, colSums(m[-1]), FUN = "/") # normalize to 1
  m <- m[ rowSums(m[-1] >= cutoff) >= nsamp, ] # filter: k relative abundance in n samples
  rownames(m) <- m$Taxa; m <- m %>% dplyr::select(-Taxa)
  if (clr){
    m <- clr_transform(t(m)) %>% as.data.frame()
  }
  return(m)
}

plot_16s_vs_mgs <- function(taxLevel = "Class"){
  taxLetter = tolower(substr(taxLevel, 1, 1))
  a <- load_16s(taxLevel = taxLevel, clr = TRUE)
  a <- a[!grepl("uncultured", rownames(a)),]
  a <- a[grepl(paste0(taxLetter, "__"), rownames(a)),]
  a$Taxa <- gsub(paste0(".*", taxLetter, "__"), "", rownames(a))
  a <- reshape2::melt(a, by = "Taxa")
  colnames(a) <- c("Taxa", "Patient_ID", "Abundance_16S")
  
  m <- load_mgs_taxa(taxLevel = taxLevel, clr = TRUE)
  m$Taxa <- gsub(paste0(".*", taxLetter, "__"), "", rownames(m))
  m <- reshape2::melt(m, by = "Taxa")
  colnames(m) <- c("Taxa", "Patient_ID", "Abundance_MGS")
  
  b <- inner_join(a, m, by = c("Taxa", "Patient_ID"))
  test <- cor.test(b$Abundance_16S, b$Abundance_MGS, method = "spearman")
  rho = test$estimate
  p_value = test$p.value
  p <- ggplot(b, aes(x = Abundance_16S, y = Abundance_MGS)) + 
    geom_point(alpha = 0.2, shape = 16) + 
    theme_pubr() + 
    xlab("16S abundance (CLR)") + 
    ylab("MGS abundance (CLR)") + 
    annotate("text", x = min(b$Abundance_16S), y = max(b$Abundance_MGS), label = paste("italic(r) == ", formatC(signif(rho, digits=2), digits = 2, format="fg", flag="#")), parse = TRUE, hjust = 0) + 
    annotate("text", x = min(b$Abundance_16S), y = 0.9*max(b$Abundance_MGS), label = paste("italic(p) == ", formatC(signif(p_value, digits=2), digits = 2, format="fg", flag="#")), parse = TRUE, hjust = 0)
  
  fn <- paste0(fig_meta, "MGSvs16S", taxLevel, ".pdf")
  ggsave(plot = p, filename = fn, height = 3, width = 3)
}

plot_pairs <- function(dat, diet1 = "V_DRKGR", diet2 = "VK", 
                       diet1label = "Dark leafy greens (cups)", diet2label = "Vitamin K1 (mcg)"){
  p <- ggplot(dat, aes(x = .data[[diet1]], y = .data[[diet2]])) + 
    geom_point() + 
    xlab(diet1label) + 
    ylab(diet2label) + 
    theme_pubr() + 
    geom_smooth(method = "lm") 
  fn <- paste0(fig_meta, "Corr_", diet1, "_", diet2, ".pdf")
  ggsave(plot = p, filename = fn, width = 3, height = 3)
}

load_mgs_genes <- function(taxLevel = "Pathway", cutoff = 0, min_clr = 0, nsamp = 20, clr = TRUE, raw = FALSE){
  if (grepl(taxLevel, "Pathways")){
    fn <- paste0(rawdir, "GO_MGS/PathAbundanceUnstratified.csv")
    m <- read.csv(fn) %>% as.data.frame() %>% dplyr::rename(Taxa = Pathway)
  } else if (grepl(taxLevel, "KOs")){
    fn <- paste0(rawdir, "GO_MGS/genefamilies_unstratified_nonzero0.csv")
    m <- read.csv(fn) %>% as.data.frame() %>% dplyr::rename(Taxa = KO) %>% dplyr::select(-X)
  } else {
    print("Wrong taxLevel! Select Pathway or KO")
    break()
  }
  colnames(m) <- gsub("\\_", "-", colnames(m)) # match formatting with metadata
  colnames(m) <- gsub("\\.x", "", colnames(m)) # remove duplicates
  if (grepl("KO", taxLevel)){
    colnames(m) = gsub("Pt9-", "Pt9c-", colnames(m)) # matching with our data
    colnames(m) = gsub("Pt9b", "Pt9", colnames(m)) # matching with our data
  }
  if (grepl("Pathway", taxLevel)){
    colnames(m) = gsub("Pt9-Baseline", "Pt9c-Baseline", colnames(m)) # matching with our data
    colnames(m) = gsub("Pt9b-Baseline", "Pt9-Baseline", colnames(m)) # matching with our data
  }
  colnames(m) = gsub("Final", "EOT", colnames(m)) # matching with our data
  if(raw){
    rownames(m) <- m$Taxa
    return(m[-1])
  }
  m[-1] <- sweep(m[-1], 2, colSums(m[-1]), FUN = "/") # normalize to 1
  m_filt <- m[ rowSums(m[-1] >= cutoff, na.rm = TRUE) >= nsamp, ] 
  taxa_keep = m_filt$Taxa; rownames(m_filt) = m_filt$Taxa; m_filt <- m_filt %>% dplyr::select(-Taxa)
  rownames(m) <- m$Taxa; m <- m %>% dplyr::select(-Taxa)
  m <- m[, !(apply(m, 2, function(x) all(is.na(x) | x == 0)))] # remove NA columns
  m_filt <- m_filt[, !(apply(m_filt, 2, function(x) all(is.na(x) | x == 0)))] # remove NA columns
  if(clr){
    m <- clr_transform(t(m)) %>% as.data.frame()
    m[m < min_clr] <- min_clr
    m_filt = m[taxa_keep,]
  }
  return(m_filt)
}

# alpha delta delta function
alpha_dd <- function(taxLevel = "ASVs", div_index = "shannon", nsamp = 20, cutoff = 0.0005, 
                     vars = c("HEI", "AHEI"), t_filt, return_alpha = FALSE){
  # Load data, calculate shannon index
  if (taxLevel %in% c("KO", "Pathway")){
    df <- load_mgs_genes(taxLevel = taxLevel, clr = FALSE, nsamp = nsamp, cutoff = cutoff) %>% t()
  } else {
    df <- load_16s(taxLevel = taxLevel, clr = FALSE, nsamp = nsamp, cutoff = cutoff) %>% t() 
  }
  
  if (div_index == "observed"){
    alpha_calculated = rowSums(df > 0)
  } else {
    alpha_calculated <- vegan::diversity(df, index = div_index)
  }
  
  alpha_df <- data.frame(Sample_ID = names(alpha_calculated), Alpha = as.numeric(alpha_calculated))
  alpha_df <- inner_join(t_filt, alpha_df, by = "Sample_ID") %>% dplyr::select(all_of(c(vars, "PatientID", "Alpha", "Cycle")))
  
  for (var in vars){
    threshold = (max(alpha_df[[var]]) - min(alpha_df[[var]]))/100
    alpha_df[[var]][alpha_df[[var]] < threshold] = threshold
  }
  
  if (return_alpha){
    return(alpha_df)
  }
  
  alpha_df <- alpha_df %>% mutate(Cycle = factor(Cycle, levels = sort(unique(Cycle))))
  alpha_df <- alpha_df %>% mutate(PatientID = factor(PatientID, levels = sort(unique(PatientID))))
  
  delta_abs <- alpha_df %>%
    arrange(PatientID, Cycle) %>%  # Ensure correct ordering
    group_by(PatientID) %>%
    mutate(across(where(is.numeric), ~ . - lag(.), .names = "{.col}")) %>%
    dplyr::filter(!is.na(HEI))  # Remove first row per patient (no prior cycle to compare)
  
  delta_rank <- alpha_df %>%
    mutate(across(.cols = where(is.numeric) & !any_of("Alpha"), .fns = ~ rank(.), .names = "{.col}")) %>% 
    arrange(PatientID, Cycle) %>%  # Ensure correct ordering
    group_by(PatientID) %>%
    mutate(across(where(is.numeric), ~ . - lag(.), .names = "{.col}")) %>%
    dplyr::filter(!is.na(HEI))  # Remove first row per patient (no prior cycle to compare)
  
  delta_fc <- alpha_df %>%
    arrange(PatientID, Cycle) %>%  # Ensure correct ordering
    group_by(PatientID) %>%
    mutate(
      across(
        .cols = where(is.numeric) & !any_of("Alpha"),
        .fns = ~ {
          # compute pseudocount: min positive value in column
          pos_vals <- .[. > 0 & !is.na(.)]
          pseudocount <- if (length(pos_vals) == 0) 1 else min(pos_vals)
          log2(. + pseudocount) - log2(lag(. + pseudocount))
        },
        .names = "{.col}"
      ),
      across(
        .cols = contains("Alpha"),
        .fns = ~ . - lag(.),
        .names = "{.col}"
      )
    ) %>%
    dplyr::filter(!is.na(HEI))  # Remove first row per patient
  
  res = data.frame(Var = character(0), Norm = character(0), p.value = numeric(0), Estimate = numeric(0))
  for (var in vars){
   
    est <- cor.test(delta_abs$Alpha, delta_abs[[var]], method = "spearman")$estimate[[1]]
    p <- cor.test(delta_abs$Alpha, delta_abs[[var]], method = "spearman")$p.value[[1]]
    row = data.frame(Var = var, Norm = "AbsDiff", p.value = p, Estimate = est)
    res <- rbind(res, row)
    
    est <- cor.test(delta_rank$Alpha, delta_rank[[var]], method = "spearman")$estimate[[1]]
    p <- cor.test(delta_rank$Alpha, delta_rank[[var]], method = "spearman")$p.value[[1]]
    row = data.frame(Var = var, Norm = "RankDiff", p.value = p, Estimate = est)
    res <- rbind(res, row)
    
    est <- cor.test(delta_fc$Alpha, delta_fc[[var]], method = "spearman")$estimate[[1]]
    p <- cor.test(delta_fc$Alpha, delta_fc[[var]], method = "spearman")$p.value[[1]]
    row = data.frame(Var = var, Norm = "logFC", p.value = p, Estimate = est)
    res <- rbind(res, row)
  }
  
  res <- res %>% arrange(p.value)
  
  fn = paste0(table_mt_alpha, "Alpha", div_index, taxLevel, "vsDiet_DeltaDelta_cutoff", cutoff, ".csv")
  write.csv(res, file = fn)
  
  fn = paste0(table_mt_alpha, "Alpha", div_index, taxLevel, "vsDiet_DeltaFC_cutoff", cutoff, ".csv")
  write.csv(delta_fc, file = fn)
  
  fn = paste0(table_mt_alpha, "Alpha", div_index, taxLevel, "vsDiet_DeltaRank_cutoff", cutoff, ".csv")
  write.csv(delta_rank, file = fn)
  
  fn = paste0(table_mt_alpha, "Alpha", div_index, taxLevel, "vsDiet_DeltaAbs_cutoff", cutoff, ".csv")
  write.csv(delta_abs, file = fn)
}

# Extract distances
convert_delta = function(delta, dist_matrix){
  delta$Beta = NA
  for (i in 1:length(delta$Name1)){
    try(delta$Beta[i] <- dist_matrix[delta$Name1[i], delta$Name2[i]], silent = TRUE) # fill in beta, ignoring errors
    try(delta$Beta[i] <- dist_matrix[delta$Name2[i], delta$Name1[i]], silent = TRUE) # fill in beta, ignoring errors
  }
  delta <- delta %>% dplyr::filter(!is.na(Beta)) 
  return(delta)
}

# Beta delta delta function
beta_dd <- function(taxLevel = "ASVs", nsamp = 20, cutoff = 0.0005, vars = c("HEI", "AHEI"), t_filt, min_clr = 0){
  if (taxLevel %in% c("KO", "Pathway")){
    df <- load_mgs_genes(taxLevel = taxLevel, clr = TRUE, nsamp = nsamp, cutoff = cutoff, min_clr = min_clr) %>% t()
  } else {
    df <- load_16s(taxLevel = taxLevel, clr = TRUE, nsamp = nsamp, cutoff = cutoff, min_clr = min_clr) %>% t() 
  }
  dist_matrix <- as.matrix(dist(df, method = "euclidean"))
  
  
  beta_df <- t_filt %>% dplyr::select(all_of(c(vars, "PatientID", "Cycle")))
  beta_df <- beta_df %>% mutate(Cycle = factor(Cycle, levels = sort(unique(Cycle))))
  beta_df <- beta_df %>% mutate(PatientID = factor(PatientID, levels = sort(unique(PatientID))))
  
  delta_abs <- beta_df %>%
    arrange(PatientID, Cycle) %>%  # Ensure correct ordering
    group_by(PatientID) %>%
    mutate(across(where(is.numeric), ~ . - lag(.), .names = "{.col}")) %>%
    dplyr::filter(!is.na(HEI))  # Remove first row per patient (no prior cycle to compare)
  
  delta_rank <- beta_df %>%
    mutate(across(.cols = where(is.numeric) & !any_of("Alpha"), .fns = ~ rank(.), .names = "{.col}")) %>% 
    arrange(PatientID, Cycle) %>%  # Ensure correct ordering
    group_by(PatientID) %>%
    mutate(across(where(is.numeric), ~ . - lag(.), .names = "{.col}")) %>%
    dplyr::filter(!is.na(HEI))  # Remove first row per patient (no prior cycle to compare)
  
  delta_fc <- beta_df %>%
    arrange(PatientID, Cycle) %>%  # Ensure correct ordering
    group_by(PatientID) %>%
    mutate(
      across(
        .cols = where(is.numeric) & !any_of("Alpha"),
        .fns = ~ {
          # compute pseudocount: min positive value in column
          pos_vals <- .[. > 0 & !is.na(.)]
          pseudocount <- if (length(pos_vals) == 0) 1 else min(pos_vals)
          log2(. + pseudocount) - log2(lag(. + pseudocount))
        },
        .names = "{.col}"
      ),
      across(
        .cols = contains("Alpha"),
        .fns = ~ . - lag(.),
        .names = "{.col}"
      )
    ) %>%
    dplyr::filter(!is.na(HEI))  # Remove first row per patient
  
  delta_fc$Cycle2 <- case_when(delta_fc$Cycle == "C2D1" ~ "C1D1", delta_fc$Cycle == "C3D1" ~ "C2D1", TRUE ~ "C3D1")
  delta_rank$Cycle2 <- case_when(delta_rank$Cycle == "C2D1" ~ "C1D1", delta_rank$Cycle == "C3D1" ~ "C2D1", TRUE ~ "C3D1")
  delta_abs$Cycle2 <- case_when(delta_abs$Cycle == "C2D1" ~ "C1D1", delta_abs$Cycle == "C3D1" ~ "C2D1", TRUE ~ "C3D1")
  
  delta_fc$Name1 <- paste0("GO-Pt", delta_fc$PatientID, "-", delta_fc$Cycle)
  delta_fc$Name2 <- paste0("GO-Pt", delta_fc$PatientID, "-", delta_fc$Cycle2)
  
  delta_rank$Name1 <- paste0("GO-Pt", delta_rank$PatientID, "-", delta_rank$Cycle)
  delta_rank$Name2 <- paste0("GO-Pt", delta_rank$PatientID, "-", delta_rank$Cycle2)
  
  delta_abs$Name1 <- paste0("GO-Pt", delta_abs$PatientID, "-", delta_abs$Cycle)
  delta_abs$Name2 <- paste0("GO-Pt", delta_abs$PatientID, "-", delta_abs$Cycle2)
  
  delta_abs = convert_delta(delta_abs, dist_matrix)
  delta_rank = convert_delta(delta_rank, dist_matrix)
  delta_fc = convert_delta(delta_fc, dist_matrix)
  
  res = data.frame(Var = character(0), Norm = character(0), p.value = numeric(0), Estimate = numeric(0))
  for (var in vars){
    est <- cor.test(delta_abs$Beta, abs(delta_abs[[var]]), method = "spearman")$estimate[[1]]
    p <- cor.test(delta_abs$Beta, abs(delta_abs[[var]]), method = "spearman")$p.value[[1]]
    row = data.frame(Var = var, Norm = "AbsDiff", p.value = p, Estimate = est)
    res <- rbind(res, row)
    
    est <- cor.test(delta_rank$Beta, abs(delta_rank[[var]]), method = "spearman")$estimate[[1]]
    p <- cor.test(delta_rank$Beta, abs(delta_rank[[var]]), method = "spearman")$p.value[[1]]
    row = data.frame(Var = var, Norm = "RankDiff", p.value = p, Estimate = est)
    res <- rbind(res, row)
    
    est <- cor.test(delta_fc$Beta, abs(delta_fc[[var]]), method = "spearman")$estimate[[1]]
    p <- cor.test(delta_fc$Beta, abs(delta_fc[[var]]), method = "spearman")$p.value[[1]]
    row = data.frame(Var = var, Norm = "logFC", p.value = p, Estimate = est)
    res <- rbind(res, row)
  }
  
  res <- res %>% arrange(p.value)
  fn = paste0(table_mt_beta, "Beta", taxLevel, "vsDiet_DeltaDelta_cutoff", cutoff, ".csv")
  write.csv(res, file = fn)
  
  fn = paste0(table_mt_beta, "Beta", taxLevel, "vsDiet_DeltaFC_cutoff", cutoff, ".csv")
  write.csv(delta_fc, file = fn)
  
  fn = paste0(table_mt_beta, "Beta", taxLevel, "vsDiet_DeltaRank_cutoff", cutoff, ".csv")
  write.csv(delta_rank, file = fn)
  
  fn = paste0(table_mt_beta, "Beta", taxLevel, "vsDiet_DeltaAbs_cutoff", cutoff, ".csv")
  write.csv(delta_abs, file = fn)
}

# Load fuso data
load_fuso = function(){
  fn = paste0(rawdir, "GO_qPCR/fuso_abundance.csv")
  f = read.csv(fn) %>% as.data.frame()
  
  # Set all samples below limit of detection to the limit of detection/2
  f$RelAb = 100*2^(f$rRNA - f$Fn1 - 1)
  f$RelAb[f$Fn1 >= 40] = min(f$RelAb[f$Fn1 < 40])

  mat = data.frame("ID" = f$Sample_ID, "Fn_RelAb" = f$RelAb, "Fn_ExpCq" = 2^(-f$Fn1), 
                   "Fn_logRelAb" = log(f$RelAb + min(f$RelAb[f$RelAb > 0])), # take natural log
                   "Fn_Cq" = f$Fn1)
  mat$Fn_logRelAb = mat$Fn_logRelAb - min(mat$Fn_logRelAb) # set min value -> 0 for modeling
  mat$ID = gsub("Final", "EOT", mat$ID)
  rownames(mat) = mat$ID; mat = mat %>% dplyr::select(-ID)
  return(mat)
}

# Features modeling
feature_dd = function(taxLevel = "ASVs", nsamp = 20, cutoff = 0.0005, vars = c("HEI", "AHEI"), t_filt, min_clr = 0){
  if (taxLevel %in% c("KO", "Pathway")){
    df <- load_mgs_genes(taxLevel = taxLevel, clr = TRUE, nsamp = nsamp, cutoff = cutoff, min_clr = min_clr) %>% t()
  } else if (grepl("Fuso|Fn", taxLevel)){
    df = load_fuso()
  } else {
    df <- load_16s(taxLevel = taxLevel, clr = TRUE, nsamp = nsamp, cutoff = cutoff, min_clr = min_clr) %>% t() 
  }

  
  samp_both = intersect(t_filt$Sample_ID, rownames(df))
  df_pt <- t_filt %>% dplyr::filter(Sample_ID %in% samp_both)
  df_tax = df[samp_both, ]
  
  output = data.frame(Var = character(0), Taxa = character(0), Norm = character(0), 
                   p.value = numeric(0), Estimate = numeric(0), stringsAsFactors = FALSE)
  for (taxa in colnames(df_tax)){
    alpha_df <- data.frame(Sample_ID = rownames(df_tax), Alpha = as.numeric(df_tax[,taxa]))
    alpha_df <- inner_join(t_filt, alpha_df, by = "Sample_ID") %>% dplyr::select(all_of(c(vars, "PatientID", "Alpha", "Cycle")))
    pt_keep = alpha_df %>% group_by(PatientID) %>% dplyr::summarise(s = sum(Alpha > min_clr/2, na.rm = TRUE)) %>% dplyr::filter(s > 0) %>% dplyr::select(PatientID) %>% pull()
    alpha_df <- alpha_df %>% dplyr::filter(PatientID %in% pt_keep) # dplyr::filter for only patients that have the taxa in at least 1 sample
    alpha_df <- alpha_df %>% mutate(Cycle = factor(Cycle, levels = sort(unique(Cycle))))
    alpha_df <- alpha_df %>% mutate(PatientID = factor(PatientID, levels = sort(unique(PatientID))))
    
    delta_abs <- alpha_df %>%
      arrange(PatientID, Cycle) %>%  # Ensure correct ordering
      group_by(PatientID) %>%
      mutate(across(where(is.numeric), ~ . - lag(.), .names = "{.col}")) %>%
      dplyr::filter(!is.na(HEI))  # Remove first row per patient (no prior cycle to compare)
    
    delta_rank <- alpha_df %>%
      mutate(across(.cols = where(is.numeric) & !any_of("Alpha"), .fns = ~ rank(.), .names = "{.col}")) %>% 
      arrange(PatientID, Cycle) %>%  # Ensure correct ordering
      group_by(PatientID) %>%
      mutate(across(where(is.numeric), ~ . - lag(.), .names = "{.col}")) %>%
      dplyr::filter(!is.na(HEI))  # Remove first row per patient (no prior cycle to compare)
    
    delta_fc <- alpha_df %>%
      arrange(PatientID, Cycle) %>%  # Ensure correct ordering
      group_by(PatientID) %>%
      mutate(across(.cols = where(is.numeric) & !any_of("Alpha"), .fns = ~ log2(.) - log2(lag(.)), .names = "{.col}"),
             across(.cols = contains("Alpha"), .fns = ~ . - lag(.), .names = "{.col}")) %>%
      dplyr::filter(!is.na(HEI))  # Remove first row per patient (no prior cycle to compare)
    

    for (var in vars){
      if(length(delta_abs$Alpha) >= 10){ # at least 10 observations
        maxval = max(alpha_df[[var]], na.rm = TRUE)
        pt_keep = alpha_df %>% group_by(PatientID) %>% dplyr::summarise(s = sum(.data[[var]] > maxval/100)) %>% dplyr::filter(s > 0) %>%dplyr::select(PatientID) %>% pull()
        
        delta_abs_test <- delta_abs %>% dplyr::filter(PatientID %in% pt_keep)
        res <- cor.test(delta_abs_test$Alpha, delta_abs_test[[var]], method = "spearman")
        est = res$estimate[[1]]
        p = res$p.value[[1]]
        row = data.frame(Var = var, Taxa = taxa, Norm = "AbsDiff", p.value = p, Estimate = est)
        output <- rbind(output, row)
        
        delta_rank_test <- delta_rank %>% dplyr::filter(PatientID %in% pt_keep)
        res <- cor.test(delta_rank_test$Alpha, delta_rank_test[[var]], method = "spearman")
        est = res$estimate[[1]]
        p = res$p.value[[1]]
        row = data.frame(Var = var, Taxa = taxa, Norm = "RankDiff", p.value = p, Estimate = est)
        output <- rbind(output, row)
        
        delta_fc_test <- delta_fc %>% dplyr::filter(PatientID %in% pt_keep)
        res <- cor.test(delta_fc_test$Alpha, delta_fc_test[[var]], method = "spearman")
        est = res$estimate[[1]]
        p = res$p.value[[1]]
        row = data.frame(Var = var, Taxa = taxa, Norm = "logFC", p.value = p, Estimate = est)
        output <- rbind(output, row)
      }
    }
  }
  
  output <- output %>% arrange(Norm, p.value)
  
  fn = paste0(table_mt_granular, taxLevel, "vsDiet_DeltaDelta_cutoff", cutoff, ".csv")
  write.csv(output, file = fn)
}


plot_all_volcanos = function(taxLevel, cutoff = 0.0005, fdr_cut = 0.2){
  fn <- paste0(table_mt_granular, taxLevel, "vsDiet_DeltaDelta_cutoff", cutoff, ".csv")
  v <- read.csv(fn) %>% as.data.frame()
  v_new = data.frame(matrix(ncol = ncol(v), nrow = 0))
  colnames(v_new) <- colnames(v)
  for (var in unique(v$Var)){
    for (norm in unique(v$Norm)){
      v_filt <- v %>% dplyr::filter(Var == var) %>% dplyr::filter(Norm == norm)
      v_filt$FDR <- p.adjust(v_filt$p.value, method = "BH")
      v_new <- rbind(v_new, v_filt)
    }
  }
  
  v_new$Significance <- case_when(v_new$FDR < fdr_cut & v_new$Estimate > 0 ~ "Enriched",
                                  v_new$FDR < fdr_cut & v_new$Estimate < 0 ~ "Depleted",
                                  TRUE ~ "FDR > 0.2")
  p <- ggplot(v_new, aes(x = Estimate, y = -log10(p.value), color = Significance)) + 
    geom_point() + 
    scale_color_manual(values = c("Enriched" = "blue", "Depleted" = "orange", "FDR > 0.2" = "grey"), name = NULL) + 
    facet_grid(Norm ~ Var) + 
    theme_pubr() +
    theme(legend.position = "none")
  
  fn <- paste0(fig_mt_granular, "Volcano", taxLevel, "cutoff_", cutoff, ".pdf")
  ggsave(filename = fn, plot = p, width = 12, height = 6)
}

plot_all_alphas <- function(taxLevel = "ASVs", cutoff = 0.0005){
  fn <- paste0(table_mt_alpha, "Alpha", taxLevel, "vsDiet_DeltaAbs_cutoff", cutoff, ".csv")
  t_abs <- read.csv(fn) %>%dplyr::select(-X)
  t_abs <- t_abs %>% tidyr::pivot_longer(cols = all_of(setdiff(colnames(t_abs), c("PatientID", "Alpha", "Cycle"))), 
                                         names_to = "Var", 
                                         values_to = "Diet")
  
  fn <- paste0(table_mt_alpha, "Alpha", taxLevel, "vsDiet_DeltaFC_cutoff", cutoff, ".csv")
  t_fc <- read.csv(fn) %>% dplyr::select(-X) 
  t_fc <- t_fc %>% tidyr::pivot_longer(cols = all_of(setdiff(colnames(t_fc), c("PatientID", "Alpha", "Cycle"))), 
                                       names_to = "Var", 
                                       values_to = "Diet")
  
  fn <- paste0(table_mt_alpha, "shannon", taxLevel, "_Raw_cutoff", cutoff, ".csv")
  t_rank <- read.csv(fn) %>% dplyr::select(-X) 
  t_rank <- t_rank %>% tidyr::pivot_longer(cols = all_of(setdiff(colnames(t_rank), c("PatientID", "Alpha", "Cycle"))), 
                                           names_to = "Var", 
                                           values_to = "Diet")
  
  p1 <- ggplot(t_abs, aes(x = Diet, y = Alpha)) + 
    geom_point() + 
    geom_smooth(method = "lm") + 
    ylab("Delta alpha") + 
    xlab("Delta diet (absolute units)") + 
    facet_wrap( ~ Var, scales = "free_x", nrow = 1) + 
    theme_pubr() 
  
  p2 <- ggplot(t_fc, aes(x = Diet, y = Alpha)) + 
    geom_point() + 
    geom_smooth(method = "lm") + 
    ylab("Delta alpha") + 
    xlab("Delta diet (log2 fold-change)") + 
    facet_wrap( ~ Var, scales = "free_x", nrow = 1) + 
    theme_pubr() 
  
  p3 <- ggplot(t_rank, aes(x = Diet, y = Alpha)) + 
    geom_point() + 
    geom_smooth(method = "lm") + 
    ylab("Delta alpha") + 
    xlab("Delta diet (rank difference)") + 
    facet_wrap( ~ Var, scales = "free_x", nrow = 1) + 
    theme_pubr() 
  
  p <- ggarrange(p1, p2, p3, nrow = 3)
  fn <- paste0(fig_mt_alpha, "Alpha", taxLevel, "cutoff_", cutoff, ".pdf")
  ggsave(filename = fn, plot = p, width = 12, height = 6)
}

plot_all_betas <- function(taxLevel = "ASVs", cutoff = 0.0005){
  fn <- paste0(table_mt_beta, "Beta", taxLevel, "vsDiet_DeltaAbs_cutoff", cutoff, ".csv")
  t_abs <- read.csv(fn) %>% dplyr::select(-all_of(c("X", "Cycle2", "Name1", "Name2")))
  t_abs <- t_abs %>% tidyr::pivot_longer(cols = all_of(setdiff(colnames(t_abs), c("PatientID", "Beta", "Cycle"))), 
                                         names_to = "Var", 
                                         values_to = "Diet")
  
  fn <- paste0(table_mt_beta, "Beta", taxLevel, "vsDiet_DeltaFC_cutoff", cutoff, ".csv")
  t_fc <- read.csv(fn) %>% dplyr::select(-all_of(c("X", "Cycle2", "Name1", "Name2")))
  t_fc <- t_fc %>% tidyr::pivot_longer(cols = all_of(setdiff(colnames(t_fc), c("PatientID", "Beta", "Cycle"))), 
                                       names_to = "Var", 
                                       values_to = "Diet")
  
  fn <- paste0(table_mt_beta, taxLevel, "_Raw_cutoff", cutoff, ".csv")
  t_rank <- read.csv(fn) %>% dplyr::select(-all_of(c("X", "Cycle2", "Name1", "Name2")))
  t_rank <- t_rank %>% tidyr::pivot_longer(cols = all_of(setdiff(colnames(t_rank), c("PatientID", "Beta", "Cycle"))), 
                                           names_to = "Var", 
                                           values_to = "Diet")
  
  p1 <- ggplot(t_abs, aes(x = Diet, y = Beta)) + 
    geom_point() + 
    geom_smooth(method = "lm") + 
    ylab("Delta beta") + 
    xlab("Delta diet (absolute units)") + 
    facet_wrap( ~ Var, scales = "free_x", nrow = 1) + 
    theme_pubr() 
  
  p2 <- ggplot(t_fc, aes(x = Diet, y = Beta)) + 
    geom_point() + 
    geom_smooth(method = "lm") + 
    ylab("Delta beta") + 
    xlab("Delta diet (log2 fold-change)") + 
    facet_wrap( ~ Var, scales = "free_x", nrow = 1) + 
    theme_pubr() 
  
  p3 <- ggplot(t_rank, aes(x = Diet, y = Beta)) + 
    geom_point() + 
    geom_smooth(method = "lm") + 
    ylab("Delta beta") + 
    xlab("Delta diet (rank difference)") + 
    facet_wrap( ~ Var, scales = "free_x", nrow = 1) + 
    theme_pubr() 
  
  p <- ggarrange(p1, p2, p3, nrow = 3)
  fn <- paste0(fig_mt_beta, "Beta", taxLevel, "cutoff_", cutoff, ".pdf")
  ggsave(filename = fn, plot = p, width = 12, height = 6)
}

## Manuscript plotting functions
plot_hit_manuscript <- function(box, var = "KCAL", ylab = "XXX", savedir, y2 = 1.15, beeswarm = FALSE, w = 1.75, h = 3){
  box$Var = box[[var]]
  res = nlme::lme(Var ~ as.numeric(CycleInt), random = ~1|PatientID, data = box %>% filter(CycleInt < 3))
  p_value = summary(res)$tTable[2,5]
  value = res$coefficients$fixed[2]
  print(p_value) # check other model p value
  
  # Different model
  box2 <- box
  box2$CycleInt <- as.factor(box2$CycleInt)
  model <- nlme::lme(Var ~ CycleInt, random = ~1 | PatientID, data = box2)
  emm <- emmeans::emmeans(model, "CycleInt")
  contrast_res <- contrast(emm, method = "pairwise", adjust = "none")
  contrast_df <- as.data.frame(contrast_res)
  
  
  stat_df <- contrast_df %>%
    tidyr::separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
    mutate(
      label = formatC(p.value, format = "fg", digits = 2, flag = "#")
    ) %>%
    dplyr::filter(!(group1 == "2" & group2 == "3")) #dplyr::filter(p.value < 0.05)
  y.position = c(1.05, y2)*max(box$Var, na.rm = TRUE)
 
  p <- ggplot(box %>% arrange(PatientID), aes(x = as.character(CycleInt), y = Var)) + 
    geom_point(aes(color = as.character(CycleInt))) +
    geom_boxplot(alpha = 0.2, color = "black", aes(fill = as.character(CycleInt)), outlier.size = -1) 
  if (beeswarm){
    p <- ggplot(box %>% arrange(PatientID), aes(x = as.character(CycleInt), y = Var)) + 
      geom_beeswarm(aes(color = as.character(CycleInt))) + 
      geom_violin(alpha = 0.2, color = "black", aes(fill = as.character(CycleInt)))
  }
  
  p = p + 
    theme_pubr() + 
    ylab(ylab) + 
    xlab("Cycle") 
  
  num_sig = length(stat_df$label)
  if (num_sig > 0){
    y.position = y.position[1:num_sig]
    stat_df$y.position = y.position
    stat_df$group1 <- gsub("CycleInt", "", stat_df$group1)
    stat_df$group2 <- gsub("CycleInt", "", stat_df$group2)
    p <- p + stat_pvalue_manual(stat_df, label = "label", tip.length = 0.01) +
      coord_cartesian(ylim = c(min(box$Var), max(box$Var) * (y2 + 0.05)))
  }
    
  
  if (value < 0){
    p <- p + scale_color_manual(values = c("grey", "orange", "orange")) + 
      scale_fill_manual(values = c("grey", "orange", "orange")) + 
      theme(legend.position = "none") 
  } else {
    p <- p + scale_color_manual(values = c("grey", "blue", "blue")) + 
      scale_fill_manual(values = c("grey", "blue", "blue")) + 
      theme(legend.position = "none") 
  }
  fn <- paste0(savedir, "Boxplot", var, ".pdf")
  ggsave(filename = fn, plot = p, width = w, height = h)  
}

# Function to load medi macro + micronutrients
load_medi_nutrients = function(names = c("PROT", "TFAT", "CARB", "ATOC", "LZ", "M181", "THEO", "VK", "BCAR"), 
                               kcal_norm = TRUE){
  # Load MEDI data
  fn <- paste0(rawdir, "GO_MEDI/food_content.csv")
  f <- read.csv(fn) %>% as.data.frame()
  names = unique(c(c("PROT", "TFAT", "CARB"), names))
  
  fn <- paste0(rawdir, "GO_MEDI/MEDI_ASA24_NamePairs.xlsx")
  map <- read_excel(fn) %>% as.data.frame() %>% dplyr::select(Name_ASA24, name = Name_MEDI_Putative) %>%
    tidyr::separate_rows(name, sep = ";")
  
  # Select variables of interest
  names_medi <- map %>% dplyr::filter(Name_ASA24 %in% names) %>% dplyr::select(name) %>% pull()
  
  # Format 
  medi <- f %>% dplyr::filter(name %in% names_medi) %>% dplyr::filter(unit == "mg/100g")
  medi <- left_join(medi, map, by = "name")
  medi <- medi %>% dplyr::select(sample_id, Name_ASA24, abundance) # consider using abundance for micro, abundance-max for macro
  medi <- medi %>% group_by(sample_id, Name_ASA24) %>% dplyr::summarise(abundance = sum(abundance)) %>% as.data.frame()
  medi <- medi %>% tidyr::pivot_wider(names_from = Name_ASA24, values_from = abundance) %>% as.data.frame()
  medi$KCAL = 4*medi$CARB + 9*medi$TFAT + 4*medi$PROT
  if (kcal_norm){
    for (var in union(c("CARB", "TFAT", "PROT"), names)){
      medi <- normalize_column_to_kcal(medi, var, method = "per1000")
    }
  }
  medi$sample_id <- gsub("_S.*", "", medi$sample_id)
  medi$sample_id <- gsub("_L.*", "", medi$sample_id)
  medi$Sample_ID <- gsub("_", "-", medi$sample_id)
  medi <- medi %>%
    tidyr::separate(Sample_ID, into = c(NA, "Patient_ID", "Cycle"), sep = "-", remove = FALSE)
  medi$Cycle[medi$Cycle == "Final"] <- "EOT"
  medi <- medi %>% dplyr::filter(!is.na(Cycle))
  return(medi)
}


# Load medi food group info. Can use either food_group or food_subgroup
# Units: reads per million (RPM), like MEDI paper
load_medi_food = function(method = "food_group", names = c("Nuts")){
  fn <- paste0(rawdir, "GO_MEDI/food_abundance.csv")
  f <- read.csv(file = fn) %>% as.data.frame()
  medi <- f %>% dplyr::filter(!!sym(method) %in% names)
  medi$sample_id <- gsub("_S.*", "", medi$sample_id)
  medi$sample_id <- gsub("_L.*", "", medi$sample_id)
  medi$Sample_ID <- gsub("_", "-", medi$sample_id)
  medi <- medi %>%
    tidyr::separate(Sample_ID, into = c(NA, "Patient_ID", "Cycle"), sep = "-", remove = FALSE)
  medi$Cycle[medi$Cycle == "Final"] <- "EOT"
  medi <- medi %>% dplyr::filter(!is.na(Cycle))
  medi[is.na(medi)] <- 0
  medi <- medi %>% group_by(Sample_ID, Patient_ID, Cycle) %>% dplyr::summarise(Abundance = 1000000*sum(relative_abundance, na.rm = TRUE))
  return(medi)
}

# Plot
plot_medi_nutrient = function(nutrient = "CARB", logT = FALSE, axis_name = "Nutrient",
                              cycle_keep = c("C1D1", "C1D3", "C1D7", "C2D1"),
                              asa24_label = "Carbohydrate (% kcal, ASA24)", 
                              w = 2.5, yval = 100, xval = 10, 
                              axis_factor = 2.75, return_df = TRUE){
  medi_all = load_medi_nutrients(names = nutrient) 
  medi = medi_all %>% filter(Cycle %in% cycle_keep)
  medi_filt = medi %>% dplyr::select(Sample_ID, MEDI = !!sym(nutrient), Cycle)
  medi_all_filt = medi_all %>% dplyr::select(Sample_ID, MEDI = !!sym(nutrient), Cycle)
  t_filt <- t_norm %>% dplyr::select(Sample_ID, ASA = !!sym(nutrient))
  merged_filt <- inner_join(medi_all_filt, t_filt, by = "Sample_ID")
  medi_filt[is.na(medi_filt)] <- 0
  merged_filt[is.na(merged_filt)] = 0
  
  medi_filt$PatientID <- gsub("-.*", "", gsub("GO-", "", medi_filt$Sample_ID))
  medi_filt <- medi_filt %>% 
    group_by(PatientID) %>%
    dplyr::filter(sum(MEDI) > 0) %>%
    ungroup()
  #medi_filt <- medi_filt %>% dplyr::select(-PatientID)
  medi_filt$Day = case_when(medi_filt$Cycle %in% c("Baseline", "C1D1") ~ 0,
                    medi_filt$Cycle == "C1D3" ~ 2,
                    medi_filt$Cycle == "C1D7" ~ 6,
                    medi_filt$Cycle == "C2D1" ~ 20,
                    medi_filt$Cycle == "C3D1" ~ 41,
                    medi_filt$Cycle == "EOT" ~ 48)
  if (logT){
    p_cor <- ggplot(merged_filt, aes(x = log10(ASA + min(ASA[ASA > 0])), y = log10(MEDI + min(MEDI[MEDI > 0])))) + 
      xlab(bquote(.(asa24_label))) +
      ylab(bquote(.(axis_name) ~ "(log" [10] ~ "mg/100g, MEDI)")) + 
      coord_cartesian(ylim = c(NA, log10(axis_factor*max(merged_filt$MEDI))))
    p_medi <- ggplot(medi_filt, aes(x = Cycle, y = log10(MEDI + min(MEDI[MEDI > 0])))) + 
      ylab(bquote(.(axis_name) ~ "(log" [10] ~ "mg/100g)"))
    p_medi_day <- ggplot(medi_filt, aes(x = Day, y = log10(MEDI + min(MEDI[MEDI > 0])))) + 
      ylab(bquote(.(axis_name) ~ "(log" [10] ~ "mg/100g)"))
  } else if (nutrient %in% c("CARB", "TFAT", "PROT")) {
    p_cor <- ggplot(merged_filt, aes(x = ASA, y = MEDI)) +
      xlab(bquote(.(asa24_label))) +
      ylab(paste0(axis_name, " (% kcal, MEDI)")) +
      coord_cartesian(ylim = c(NA, 1.1*max(merged_filt$MEDI))) 
    p_medi <- ggplot(medi_filt, aes(x = Cycle, y = MEDI)) +
      ylab(paste0(axis_name, " (% kcal)"))
    p_medi_day <- ggplot(medi_filt, aes(x = Day, y = MEDI)) +
      ylab(paste0(axis_name, " (% kcal)"))
  } else {
    p_cor <- ggplot(merged_filt, aes(x = ASA, y = MEDI)) +
      xlab(bquote(.(asa24_label))) + 
      ylab(paste0(axis_name, " (mg/100g, MEDI)")) +
      coord_cartesian(ylim = c(NA, 1.1*max(merged_filt$MEDI)))
    p_medi <- ggplot(medi_filt, aes(x = Cycle, y = MEDI)) +
      ylab(paste0(axis_name, " (mg/100g)"))
    p_medi_day <- ggplot(medi_filt, aes(x = Day, y = MEDI)) +
      ylab(paste0(axis_name, " (mg/100g)"))
  }
  
  p_cor = p_cor + 
    geom_point(color = "black") + 
    geom_smooth(method = "lm", color = "black", se = TRUE) + 
    theme_pubr() + 
    stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 1, vjust = -1)
  
  p_medi_box = p_medi + 
    geom_point(aes(color = Cycle)) +  #, shape = 16, stroke = 0, size = 2
    geom_boxplot(aes(fill = Cycle), alpha = 0.2, outlier.size = -1) + 
    xlab("Time") + 
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values = c("Baseline" = "grey", 
                                 "C1D1" = "grey",
                                 "C1D3" = "orange", 
                                 "C1D7" = "orange", 
                                 "C2D1" = "orange", 
                                 "C3D1" = "orange", 
                                 "EOT" = "orange")) + 
    scale_color_manual(values = c("Baseline" = "grey", 
                                  "C1D1" = "grey",
                                  "C1D3" = "orange", 
                                  "C1D7" = "orange", 
                                  "C2D1" = "orange", 
                                  "C3D1" = "orange", 
                                  "EOT" = "orange")) +
    theme(legend.position = "none")
  

  # p value
  medi_filt$Day_grouped <- factor(ifelse(medi_filt$Day == 0, "Day0", "OtherDays"))
  model <- lme(
    MEDI ~ Day_grouped,
    random = ~ 1 | PatientID,
    data = medi_filt,
    method = "REML"
  )
  coefs <- summary(model)$tTable
  pval <- coefs["Day_groupedOtherDays", "p-value"]
  print(pval)
  label_df <- data.frame(
    x = xval,  # adjust to desired position
    y = yval,  # slightly above top
    label = paste0("italic(p) == ", signif(pval, 2))
  )
  
  p_medi_line = p_medi_day +
    stat_summary(fun = mean, geom = "line", aes(group = 1)) +
    stat_summary(fun = mean,  geom = "point") +
    stat_summary(fun.data = mean_se,  geom = "errorbar", width = 1.5) + 
    theme_pubr() + 
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = label),
      parse = TRUE
    )
  
  fn <- paste0(figdir, "MEDI", nutrient, "_Time_Line.pdf")
  ggsave(filename = fn, plot = p_medi_line, height = 3, width = w)
  
  fn <- paste0(figdir, "MEDI", nutrient, "_corr.pdf")
  ggsave(filename = fn, plot = p_cor, height = 3, width = w)
  
  fn <- paste0(figdir, "MEDI", nutrient, "_Time_Boxplot.pdf")
  ggsave(filename = fn, plot = p_medi_box, height = 3, width = 2)
  
  if (return_df){return(medi_filt)}
}

## Plot MEDI food groups
plot_medi_food <- function(method = "food_group", name_asa24 = "PF_NUTSDS",
                           names = c("Nuts"), axis_name = "Nuts", 
                           cycle_keep = c("C1D1", "C1D3", "C1D7", "C2D1"),
                           asa24_label = "Nuts (oz, ASA24)",
                           read_cut = 2000, w = 2.5, yval = 1.6, xval = 10){
  medi_all <- load_medi_food(method = method, names = names) 
  medi_all$Abundance[medi_all$Abundance > read_cut] = read_cut 
  medi = medi_all %>% filter(Cycle %in% cycle_keep)

  medi$Day = case_when(medi$Cycle %in% c("Baseline", "C1D1") ~ 0,
                       medi$Cycle == "C1D3" ~ 2,
                       medi$Cycle == "C1D7" ~ 6,
                       medi$Cycle == "C2D1" ~ 20,
                       medi$Cycle == "C3D1" ~ 41,
                       medi$Cycle == "EOT" ~ 48)
  minval = min(medi_all$Abundance[medi_all$Abundance > 0])
  medi$Abundance = log10(medi$Abundance + minval)
  medi_all$Abundance = log10(medi_all$Abundance + minval)
  
  p <- ggplot(medi, aes(x = Cycle, y = Abundance)) + 
    geom_point(aes(color = Cycle)) +  #, shape = 16, stroke = 0, size = 2
    geom_boxplot(aes(fill = Cycle), alpha = 0.2, outlier.size = -1) + 
    xlab("Time") + 
    ylab(paste0(axis_name, " (log10 rpm)")) + 
    theme_pubr() + 
    scale_fill_manual(values = c("Baseline" = "grey", 
                                  "C1D1" = "grey",
                                  "C1D3" = "orange", 
                                  "C1D7" = "orange", 
                                  "C2D1" = "orange", 
                                  "C3D1" = "orange", 
                                  "EOT" = "orange")) + 
    scale_color_manual(values = c("Baseline" = "grey", 
                                 "C1D1" = "grey",
                                 "C1D3" = "orange", 
                                 "C1D7" = "orange", 
                                 "C2D1" = "orange", 
                                 "C3D1" = "orange", 
                                 "EOT" = "orange")) +
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  fn <- paste0(figdir, "MEDI_", axis_name, "_Time_Boxplot.pdf")
  ggsave(filename = fn, plot = p, width = 2, height = w)
  
  # add p values
  medi$Day_grouped <- factor(ifelse(medi$Day == 0, "Day0", "OtherDays"))
  model <- lme(
    Abundance ~ Day_grouped,
    random = ~ 1 | Patient_ID,
    data = medi,
    method = "REML"
  )
  coefs <- summary(model)$tTable
  pval <- coefs["Day_groupedOtherDays", "p-value"]
  print(pval)
  label_df <- data.frame(
    x = xval,  # adjust to desired position
    y = yval,  # slightly above top
    label = paste0("italic(p) == ", signif(pval, 2))
  )
  
  p_medi_line = ggplot(medi, aes(x = Day, y = Abundance)) +
    stat_summary(fun = mean, geom = "line", aes(group = 1)) +
    stat_summary(fun = mean,  geom = "point") +
    xlab("Time (days)") + 
    ylab(bquote(.(axis_name) ~ "(log" [10] ~ "rpm)")) + 
    stat_summary(fun.data = mean_se,  geom = "errorbar", width = 1.5) + 
    theme_pubr() + 
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = label),
      parse = TRUE
    )
  
  fn <- paste0(figdir, "MEDI_", axis_name, "_Time.pdf")
  ggsave(filename = fn, plot = p_medi_line, height = 3, width = w)
  
  medi_bl <- medi %>%
    group_by(Patient_ID) %>%
    mutate(
      baseline_mean = mean(Abundance[Day == 0], na.rm = TRUE),
      Abundance_norm = Abundance - baseline_mean
    ) %>%
    ungroup()
  
  # Correlation plot
  t_filt <- t_norm %>% dplyr::select(Sample_ID, ASA = !!sym(name_asa24))
  t_filt$ASA = log10(t_filt$ASA + min(t_filt$ASA[t_filt$ASA > 0]))
  merged = full_join(t_filt, medi_all, by = "Sample_ID")
  merged = merged %>% filter(!is.na(ASA))
  merged[is.na(merged)] = 0
  
  p_cor = ggplot(merged, aes(x = ASA, y = Abundance)) + 
    xlab(bquote(.(asa24_label))) + 
    ylab(bquote(.(axis_name) ~ "(log" [10] ~ "rpm, MEDI)")) + 
    geom_point(color = "black") + 
    geom_smooth(method = "lm", color = "black", se = TRUE) + 
    theme_pubr() +
    coord_cartesian(ylim = c(NA, 1.1*max(merged$Abundance))) + 
    stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 1, vjust = -1) 
  
  fn <- paste0(figdir, "MEDI_", axis_name, "_corr.pdf")
  ggsave(filename = fn, plot = p_cor, width = w, height = 3)
}

# Function to load medi macro + micronutrients for NE study
load_medi_nutrients_ne = function(names = c("PROT", "TFAT", "CARB"), kcal_norm = TRUE){
  # Load MEDI data
  fn <- paste0(rawdir, "NE_MEDI/food_content.csv")
  f <- read.csv(fn) %>% as.data.frame()
  names = unique(c(c("PROT", "TFAT", "CARB"), names))
  
  fn <- paste0(rawdir, "GO_MEDI/MEDI_ASA24_NamePairs.xlsx")
  map <- read_excel(fn) %>% as.data.frame() %>% dplyr::select(Name_ASA24, name = Name_MEDI_Putative) %>%
    tidyr::separate_rows(name, sep = ";")
  
  # Select variables of interest
  names_medi <- map %>% dplyr::filter(Name_ASA24 %in% names) %>% dplyr::select(name) %>% pull()
  
  # Format 
  medi <- f %>% dplyr::filter(name %in% names_medi) %>% dplyr::filter(unit == "mg/100g")
  medi <- left_join(medi, map, by = "name")
  medi <- medi %>% dplyr::select(sample_id, Name_ASA24, abundance_max) # consider using abundance for micro, abundance-max for macro. This goes in this line and the below line
  medi <- medi %>% group_by(sample_id, Name_ASA24) %>% dplyr::summarise(abundance = sum(abundance_max)) %>% as.data.frame()
  medi <- medi %>% tidyr::pivot_wider(names_from = Name_ASA24, values_from = abundance) %>% as.data.frame()
  medi$KCAL = 4*medi$CARB + 9*medi$TFAT + 4*medi$PROT
  if (kcal_norm){
    for (var in union(c("CARB", "TFAT", "PROT"), names)){
      medi <- normalize_column_to_kcal(medi, var, method = "per1000")
    }
  }
  medi$Sample <- gsub("_S.*", "", medi$sample_id)
  medi <- medi %>%
    tidyr::separate(Sample, into = c(NA, "Patient_ID", "Cycle"), sep = "_", remove = FALSE)
  medi <- medi %>% dplyr::filter(!is.na(Cycle))
  return(medi)
}


# Load medi food group info. Can use either food_group or food_subgroup. This is for NE study
# Units: reads per million (RPM), like MEDI paper
load_medi_food_ne = function(method = "food_group", names = c("Nuts")){
  fn <- paste0(rawdir, "NE_MEDI/food_abundance.csv")
  f <- read.csv(file = fn) %>% as.data.frame()
  medi <- f %>% dplyr::filter(!!sym(method) %in% names)
  medi$Sample <- gsub("_S.*", "", medi$sample_id)
  medi <- medi %>%
    tidyr::separate(Sample, into = c(NA, "Patient_ID", "Cycle"), sep = "_", remove = FALSE)
  medi <- medi %>% dplyr::filter(!is.na(Cycle))
  medi[is.na(medi)] <- 0
  medi <- medi %>% group_by(Sample, Patient_ID, Cycle) %>% dplyr::summarise(Abundance = 1000000*sum(relative_abundance, na.rm = TRUE))
  return(medi)
}



load_go_tox <- function(t_in){
  # From REDCap
  hfs <- c("25", "1", "2", "11", "14", "20", "21") # "hfs", "hand", "neuropathy", "erythrodys*" - ensure it's new during treatment
  stomatitis <- c("1", "25") # "stomatitis" 
  reduc <- c("7", "11", "25", "20", "38", "1") # "reduc", "decreas"
  gi <- c("1", "7", "20", "21", "54", "51", "22", "33", "31") #  "diarrhea" or "constipation" - ensure it's new during treatment
  nausea <- c("19", "33", "4", "54", "1") # "nausea" 
  fatigue <- c("23", "4", "11", "35", "45", "54") # "fatigue"
  neutropenia <- c("4") # "neutropenia"
  
  # From validated clinical trial spreadsheet, "PCB". Note that the Patient 6 toxicities are all several months after the initial cycles of treatment, so we exclude these
  hfs_validated <- c("1", "2", "13", "14", "18") # "Definite"
  
  # Final tox of interest
  tox_all <- unique(c(stomatitis, hfs, reduc, gi, fatigue, nausea, neutropenia, hfs_validated)) 
  tox_hfs = unique(c(hfs, hfs_validated))
  
  # Patients with progressive disease (PD)
  pd <- c("20", "30", "47", "14", "6", "17", # "progress*" in REDCAP
          "2", "4", "7", "8", "13", "15", "44") # GO_Samples_TLab_DS spreadsheet for discontinued (clinical decline) or deceased
  
  # Add to dataframe
  t_pt <- t_in
  t_pt$AnyTox <- case_when(t_pt$PatientID %in% tox_all ~ 1, TRUE ~ 0)
  t_pt$HFS <- case_when(t_pt$PatientID %in% tox_hfs ~ 1, TRUE ~ 0)
  t_pt$Dose <- case_when(t_pt$PatientID %in% reduc ~ 1, TRUE ~ 0)
  t_pt$PD <- case_when(t_pt$PatientID %in% pd ~ 1, TRUE ~ 0)
  t_pt$GI <- case_when(t_pt$PatientID %in% c(nausea, gi) ~ 1, TRUE ~ 0)
  return(t_pt)
}

## Plot GO Tox vs Diet. outcome_var options: PD, HFS, Dose, AnyTox
plot_go_tox <- function(t_pt, diet_var = "F_JUICE", outcome_var = "PD", xlabel = "Outcome", 
                        ylabel = "Diet", method = "single", ylabel_fc = NULL, w = 2.25,
                        stat_test = "t.test"){
  if (is.null(ylabel_fc)){
    ylabel_fc <- bquote("Change in" ~ .(ylabel))
  }
  t_pt_filt = t_pt
  t_pt_filt$Diet = t_pt_filt[[diet_var]]
  t_pt_filt$Outcome = t_pt_filt[[outcome_var]]
  t_pt_filt = t_pt_filt %>% dplyr::select(Outcome, Cycle, Diet, PatientID)
  if (method == "single"){t_pt_single = t_pt_filt %>% filter(Cycle == "C1D1")}
  if (method == "average"){t_pt_single = t_pt_filt  %>% filter(Cycle %in% c("C1D1", "C2D1")) %>% dplyr::group_by(PatientID) %>% 
    dplyr::summarise(Diet = mean(Diet), Outcome = max(Outcome, na.rm = TRUE))}
  t_pt_single$Outcome <- case_when(t_pt_single$Outcome == 0 ~ "No",
                                 TRUE ~ "Yes")
  t_pt_filt$Outcome <- case_when(t_pt_filt$Outcome == 0 ~ "No",
                                 TRUE ~ "Yes")
  
  p1 <- ggplot(t_pt_single, aes(x = Outcome, y = Diet)) + 
    geom_boxplot(alpha = 0.2, color = "black", outlier.size = -1, aes(fill = Outcome)) + 
    geom_point(aes(color = Outcome)) + 
    scale_color_manual(values = c("No" = "grey", "Yes" = "red")) + 
    scale_fill_manual(values = c("No" = "grey", "Yes" = "red")) +     theme_pubr() + 
    xlab(xlabel) + 
    ylab(ylabel) + 
    stat_compare_means(label.x = 1.25, method = stat_test, 
                       aes(label = paste0("italic(p) ==", after_stat(p.format))), parse = TRUE,
                       label.y.npc = 1) + 
    coord_cartesian(clip = "off")+
    theme(legend.position = "none")
  
  pt_pair <- t_pt_filt %>% dplyr::filter(Cycle %in% c("C1D1", "C2D1")) %>% dplyr::group_by(PatientID) %>% dplyr::summarise(n = n()) %>% dplyr::filter(n == 2) %>% dplyr::select(PatientID) %>% pull()
  t_pt_delta <- t_pt_filt %>% dplyr::filter(PatientID %in% pt_pair)
  
  df_delta <- t_pt_delta %>% 
    group_by(PatientID) %>% 
    mutate(
      DeltaDiet = Diet[Cycle == "C2D1"] - Diet[Cycle == "C1D1"]
    ) %>% 
    ungroup() %>% 
    filter(Cycle == "C1D1")
  
  p2 <- ggplot(df_delta, aes(x = Outcome, y = DeltaDiet)) + 
    geom_point(aes(color = Outcome)) + 
    geom_boxplot(alpha = 0.2, color = "black", outlier.size = -1, aes(fill = Outcome)) + 
    scale_color_manual(values = c("No" = "grey", "Yes" = "red")) + 
    scale_fill_manual(values = c("No" = "grey", "Yes" = "red")) + 
    xlab(xlabel) + 
    ylab(ylabel_fc) +
    theme_pubr() + 
    stat_compare_means(label.x = 1.25, method = "t.test", 
                       aes(label = paste0("italic(p) ==", after_stat(p.format))), parse = TRUE,
                       label.y.npc = 1) +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none")
  
  fn <- paste0(figdir, "Baseline_", diet_var, "_", outcome_var, ".pdf")
  ggsave(filename = fn, plot = p1, width = w, height = 3)
  
  fn <- paste0(figdir, "Delta_", diet_var, "_", outcome_var, ".pdf")
  ggsave(filename = fn, plot = p2, width = w, height = 3)
}

# test = "wilcox" or "t.test"
run_baseline_vs_tox <- function(t_pt, test = "wilcox"){
  t_pt_filt = t_pt %>% filter(Cycle == "C1D1")
  diet_vars = colnames(t_pt_filt[c(diet_index_start:diet_index_end,
                                   hei_index_start:hei_index_end)])
  outcome_vars = c("PD", "HFS", "Dose", "AnyTox", "GI")
  
  df = NULL
  for (diet_var in diet_vars){
    t_pt_filt$DietTest = t_pt_filt[[diet_var]]
    for (outcome_var in outcome_vars){ 
      t_pt_filt$OutcomeTest = t_pt_filt[[outcome_var]]
      x1 = t_pt_filt %>% filter(OutcomeTest == 0) %>% dplyr::select(DietTest) %>% pull()
      x2 = t_pt_filt %>% filter(OutcomeTest == 1) %>% dplyr::select(DietTest) %>% pull()
      if (test == "wilcox"){res = wilcox.test(x1, x2)}
      if (test == "t.test"){res = t.test(x1, x2)}
      pv = res$p.value
      new_row = data.frame("Diet"= diet_var, "Outcome" = outcome_var, "p.value" = pv, 
                           "Mean_0" = mean(x1, na.rm = TRUE), "Mean_1" = mean(x2, na.rm = TRUE),
                           "Statistic" = res$statistic[[1]])
      if (length(df) == 0){
        df = new_row
      } else {
        df = rbind(df, new_row)
      }
    }
  }

  # None survive FDR correction... but many are individually significant! 
  df$FDR <- p.adjust(df$p.value, method = "BH")
  
  # Save
  fn <- paste0(tab_ms, "OutcomesVsDiet_GO_ASA24.csv")
  write.csv(df %>% arrange(p.value), file = fn, quote = FALSE, row.names = FALSE)
}


# test = "wilcox" or "t.test"
run_delta_vs_tox <- function(t_pt, test = "wilcox"){
  t_pt_filt = t_pt %>% filter(Cycle %in% c("C1D1", "C2D1"))
  pt_keep = t_pt_filt %>% dplyr::group_by(PatientID) %>% dplyr::summarise(n = n()) %>% dplyr::filter(n == 2) %>% dplyr::select(PatientID) %>% pull()
  t_pt_filt <- t_pt_filt %>% filter(PatientID %in% pt_keep)
  diet_vars = colnames(t_pt_filt[c(diet_index_start:diet_index_end,
                                   hei_index_start:hei_index_end)])
  outcome_vars = c("PD", "HFS", "Dose", "AnyTox", "GI")
  
  df = NULL
  for (diet_var in diet_vars){
    t_pt_filt$DietTest = t_pt_filt[[diet_var]]
    for (outcome_var in outcome_vars){ 
      t_pt_filt$OutcomeTest = t_pt_filt[[outcome_var]]
      d1_baseline = t_pt_filt %>% arrange(PatientID) %>% filter(OutcomeTest == 0) %>% dplyr::filter(Cycle == "C1D1") %>% dplyr::select(DietTest) %>% pull()
      d1_post = t_pt_filt %>% arrange(PatientID) %>% filter(OutcomeTest == 0) %>% dplyr::filter(Cycle == "C2D1") %>% dplyr::select(DietTest) %>% pull()
      d2_baseline = t_pt_filt %>% arrange(PatientID) %>% filter(OutcomeTest == 1) %>% dplyr::filter(Cycle == "C1D1") %>% dplyr::select(DietTest) %>% pull()
      d2_post = t_pt_filt %>% arrange(PatientID) %>% filter(OutcomeTest == 1) %>% dplyr::filter(Cycle == "C2D1") %>% dplyr::select(DietTest) %>% pull()
      delta1 = d1_post - d1_baseline
      delta2 = d2_post - d2_baseline
      if (test == "wilcox"){res = wilcox.test(delta1, delta2)}
      if (test == "t.test"){res = t.test(delta1, delta2)}
      pv = res$p.value
      new_row = data.frame("Diet"= diet_var, "Outcome" = outcome_var, "p.value" = pv,
                           "Mean_0" = mean(delta1, na.rm = TRUE), "Mean_1" = mean(delta2, na.rm = TRUE),
                           "Statistic" = res$statistic[[1]])
      if (length(df) == 0){
        df = new_row
      } else {
        df = rbind(df, new_row)
      }
    }
  }
  
  # None survive FDR correction... but many are individually significant! 
  df$FDR <- p.adjust(df$p.value, method = "BH")
  
  # Save
  fn <- paste0(tab_ms, "OutcomesVsDiet_GO_ASA24Delta.csv")
  write.csv(df %>% arrange(p.value), file = fn, quote = FALSE, row.names = FALSE)
}

load_tox_medi_go <- function(names = c("PROT", "TFAT", "CARB", "ATOC", "COPP"), kcal_norm = TRUE){
  pd = c("20", "30", "47", "14", "6", "17", "2", "4", "7", "8", "13", "15", "44")
  tox_all = c("1",  "25", "2", "11", "14", "51", "54", "7", "20", "38", "22", "21", "33", "31", "19", "4", "23", "35", "45", "13", "18")
  tox_hfs = c("25", "1", "2", "11", "14", "20", "21")
  reduc = c("7", "11", "25", "20", "38", "1")
  gi = c("19", "33", "4", "54", "1", "7", "20", "21", "51", "22", "31")
  
  dat <- load_medi_nutrients(names = names, kcal_norm = kcal_norm)
  dat$PatientID = gsub("Pt", "", dat$Patient_ID)
  dat$AnyTox <- case_when(dat$PatientID %in% tox_all ~ 1, TRUE ~ 0)
  dat$HFS <- case_when(dat$PatientID %in% tox_hfs ~ 1, TRUE ~ 0)
  dat$Dose <- case_when(dat$PatientID %in% reduc ~ 1, TRUE ~ 0)
  dat$PD <- case_when(dat$PatientID %in% pd ~ 1, TRUE ~ 0)
  dat$GI <- case_when(dat$PatientID %in% gi ~ 1, TRUE ~ 0)
  return(dat)
}


plot_tox_medi_go <- function(dat, diet_var = "ATOC", outcome_var = "AnyTox", 
                             xlabel = "Outcome", ylabel = "Diet", ylabel_fc = NULL){
  if (is.null(ylabel_fc)){
    ylabel_fc = paste0("Change in ", ylabel)
  }
  t_pt_filt = dat
  t_pt_filt$Diet = t_pt_filt[[diet_var]]
  t_pt_filt$Outcome = t_pt_filt[[outcome_var]]
  t_pt_filt$Outcome <- case_when(t_pt_filt$Outcome == 0 ~ "No",
                                 TRUE ~ "Yes")
  
  p1 <- ggplot(t_pt_filt %>% filter(Cycle == "C1D1"), aes(x = Outcome, y = Diet)) + 
    geom_boxplot(alpha = 0.2, color = "black", outlier.size = -1, aes(fill = Outcome)) + 
    geom_point(aes(color = Outcome)) + 
    scale_color_manual(values = c("No" = "grey", "Yes" = "red")) + 
    scale_fill_manual(values = c("No" = "grey", "Yes" = "red")) +     theme_pubr() + 
    xlab(xlabel) + 
    ylab(ylabel) + 
    stat_compare_means(label.x = 1.25, method = "t.test", 
                       aes(label = paste0("italic(p) ==", after_stat(p.format))), parse = TRUE,
                       label.y.npc = 1) + 
    coord_cartesian(clip = "off")+
    theme(legend.position = "none")
  
  pt_pair <- t_pt_filt %>% dplyr::filter(Cycle %in% c("C1D1", "C2D1")) %>% dplyr::group_by(PatientID) %>% dplyr::summarise(n = n()) %>% dplyr::filter(n == 2) %>% dplyr::select(PatientID) %>% pull()
  t_pt_delta <- t_pt_filt %>% dplyr::filter(PatientID %in% pt_pair)
  
  df_delta <- t_pt_delta %>% 
    group_by(PatientID) %>% 
    mutate(
      DeltaDiet = Diet[Cycle == "C2D1"] - Diet[Cycle == "C1D1"]
    ) %>% 
    ungroup() %>% 
    filter(Cycle == "C1D1")
  
  p2 <- ggplot(df_delta, aes(x = Outcome, y = DeltaDiet)) + 
    geom_boxplot(alpha = 0.2, color = "black", outlier.size = -1, aes(fill = Outcome)) + 
    geom_point(aes(color = Outcome)) + 
    scale_color_manual(values = c("No" = "grey", "Yes" = "red")) + 
    scale_fill_manual(values = c("No" = "grey", "Yes" = "red")) + 
    xlab(xlabel) + 
    ylab(ylabel_fc) +
    theme_pubr() + 
    stat_compare_means(label.x = 1.25, method = "t.test", 
                       aes(label = paste0("italic(p) ==", after_stat(p.format))), parse = TRUE,
                       label.y.npc = 1) +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none")
  
  # Save
  fn <- paste0(figdir, "Baseline_", diet_var, "_", outcome_var, "_MEDI_GO.pdf")
  ggsave(filename = fn, plot = p1, width = 2.25, height = 3)
  
  fn <- paste0(figdir, "Delta_", diet_var, "_", outcome_var, "_MEDI_GO.pdf")
  ggsave(filename = fn, plot = p2, width = 2.25, height = 3)
}


plot_tox_medi_ne <- function(dat, diet_var = "ATOC", outcome_var = "AnyTox_T2", 
                             xlabel = "Outcome", ylabel = "Diet", ylabel_fc = NULL){
  if (is.null(ylabel_fc)){ylabel_fc = paste0("Change in ", ylabel)}
  t_pt_filt = dat
  t_pt_filt$Diet = t_pt_filt[[diet_var]]
  t_pt_filt$Outcome = t_pt_filt[[outcome_var]]
  t_pt_filt <- t_pt_filt %>% 
    filter(!is.na(Outcome)) %>%
    filter(Outcome < 100) 
  t_pt_filt$Outcome <- case_when(t_pt_filt$Outcome == 0 ~ "No",
                                 TRUE ~ "Yes")
  
  p1 <- ggplot(t_pt_filt %>% filter(Time == "Baseline"), aes(x = Outcome, y = Diet)) + 
    geom_boxplot(alpha = 0.2, color = "black", outlier.size = -1, aes(fill = Outcome)) + 
    geom_point(aes(color = Outcome)) + 
    scale_color_manual(values = c("No" = "grey", "Yes" = "red")) + 
    scale_fill_manual(values = c("No" = "grey", "Yes" = "red")) +     
    theme_pubr() + 
    xlab(xlabel) + 
    ylab(ylabel) + 
    stat_compare_means(label.x = 1.25, method = "t.test", 
                       aes(label = paste0("italic(p) ==", after_stat(p.format))), parse = TRUE,
                       label.y.npc = 1) + 
    coord_cartesian(clip = "off")+
    theme(legend.position = "none")
  
  pt_pair <- t_pt_filt %>% dplyr::filter(Time %in% c("Baseline", "Cycle 3")) %>% 
    dplyr::group_by(Patient_ID) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::filter(n == 2) %>% 
    dplyr::select(Patient_ID) %>% pull()
  t_pt_delta <- t_pt_filt %>% dplyr::filter(Patient_ID %in% pt_pair)
  t_pt_delta$Diet[is.na(t_pt_delta$Diet)] <- 0
  
  df_delta <- t_pt_delta %>% 
    group_by(Patient_ID) %>% 
    mutate(
      DeltaDiet = Diet[Time == "Cycle 3"] - Diet[Time == "Baseline"]
    ) %>% 
    ungroup() %>% 
    filter(Time == "Baseline")
  
  
  p2 <- ggplot(df_delta, aes(x = Outcome, y = DeltaDiet)) + 
    geom_boxplot(alpha = 0.2, color = "black", outlier.size = -1, aes(fill = Outcome)) + 
    geom_point(aes(color = Outcome)) + 
    scale_color_manual(values = c("No" = "grey", "Yes" = "red")) + 
    scale_fill_manual(values = c("No" = "grey", "Yes" = "red")) + 
    xlab(xlabel) + 
    ylab(ylabel_fc) +
    theme_pubr() + 
    stat_compare_means(label.x = 1.25, method = "t.test", 
                       aes(label = paste0("italic(p) ==", after_stat(p.format))), parse = TRUE,
                       label.y.npc = 1) +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none")
  
  # Save
  fn <- paste0(figdir, "Baseline_", diet_var, "_", outcome_var, "_MEDI_NE.pdf")
  ggsave(filename = fn, plot = p1, width = 2.25, height = 3)
  
  fn <- paste0(figdir, "Delta_", diet_var, "_", outcome_var, "_MEDI_NE.pdf")
  ggsave(filename = fn, plot = p2, width = 2.25, height = 3)
}

load_tox_medi_ne <- function(names = c("PROT", "TFAT", "CARB", "ATOC", "COPP"), kcal_norm = TRUE){
  dat <- load_medi_nutrients_ne(names = names, kcal_norm = kcal_norm)
  fn <- paste0(rawdir, "NE_MEDI/NEMetadata.csv")
  tox <- read.csv(fn)
  merged_all <- inner_join(tox %>% dplyr::select(-Patient_ID), dat, by = "Sample")
  return(merged_all)
}

# Transform: "abs" or "fc"
plot_diet_beta_dd <- function(nutrient = "HEI", nutrient_name = "HEI", taxLevel = "Genus", figdir,
                              transform = "abs", w = 3, h = 3, return_plot = FALSE){
  if (transform == "fc"){
    fn <- paste0(datadir, "Tables/DietVsMicrobiomeTime/Beta/", taxLevel, "_RawFC_cutoff0.csv")
  } else if (transform == "abs"){
    fn <- paste0(datadir, "Tables/DietVsMicrobiomeTime/Beta/", taxLevel, "_RawAbs_cutoff0.csv")
  }
  t <- read.csv(fn) 
  t$Nutrient = t[[nutrient]]
  
  p1_1 <- ggplot(t %>% filter(abs(Nutrient) > 0), aes(x = abs(Nutrient), y = Beta)) + 
    geom_point(color = "black") +
    theme_pubr() + 
    xlab(bquote("| " * Delta ~ .(nutrient_name) * " |")) + 
    geom_smooth(method = "lm", color = "black") + 
    ylab(expression(Delta ~ beta ~ "diversity (CLR-Euclidean)")) +
    stat_cor(method = "pearson", label.x = min(abs(t$Nutrient), na.rm = TRUE), label.y = 1.05*max(t$Beta, na.rm = TRUE)) + 
    ylim(c(min(t$Beta, na.rm = TRUE), 1.05*max(t$Beta, na.rm = TRUE)))
  p1_2 <- ggplot(t, aes(x = Nutrient, y = Beta)) + 
    geom_point(color = "black") +
    theme_pubr() + 
    geom_smooth(method = "lm", color = "black") + 
    stat_cor(method = "pearson", label.x = min(abs(t$Nutrient), na.rm = TRUE), label.y = 1.05*max(t$Beta, na.rm = TRUE))
  
  if (transform == "fc"){p1_1 = p1_1 + xlab(bquote("| " * Delta ~ .(nutrient_name) * " (log"[2]*"FC) |"))}
  
  p1 <- ggarrange(p1_2, p1_1)
  fn <- paste0(figdir, "BetaCorrAll", nutrient, ".pdf")
  ggsave(filename = fn, plot = p1, width = 7, height = 3)

  fn <- paste0(figdir, "BetaCorr", nutrient, ".pdf")
  ggsave(filename = fn, plot = p1_1, width = w, height = h)
  
  t$Nutrient = abs(t$Nutrient)
  t_bin = t
  med = median(t_bin$Nutrient)
  t_bin$Nutrient <- case_when(t_bin$Nutrient <= med ~ "Low", 
                              TRUE ~ "High")
  t_bin$Nutrient <- factor(t_bin$Nutrient, levels = c("Low", "High"))
  
  p2 <- ggplot(t_bin, aes(x = Nutrient, y = Beta)) + 
    geom_boxplot(aes(fill = Nutrient), alpha = 0.2, outlier.size = -1, width = 0.6) + 
    geom_point(aes(color = Nutrient)) + 
    theme_pubr() + 
    xlab(bquote("| " * Delta ~ .(nutrient_name) * " |")) + 
    #xlab(bquote(Delta ~ .(nutrient_name))) + 
    ylab(expression(Delta ~ beta ~ "diversity (CLR-Euclidean)")) + 
    stat_compare_means(label.x = 1, label.y = 1.05*max(t_bin$Beta), 
                       method = "wilcox.test", 
                       aes(label = paste0("italic(p) ==", after_stat(p.format))),
                       parse = TRUE) + 
    ylim(c(min(t_bin$Beta), 1.15*max(t_bin$Beta))) +
    theme(legend.position = "none") + 
    scale_color_manual(values = c("High" = "blue", "Low" = "orange")) + 
    scale_fill_manual(values = c("High" = "blue", "Low" = "orange")) 
  
  if (return_plot){return(p2)}
  
  fn <- paste0(figdir, "BetaBoxplot", nutrient, ".pdf")
  ggsave(filename = fn, plot = p2, width = 2, height = 3)
}

plot_diet_alpha_dd = function(nutrient = "THEO", nutrient_name = "Theobromine", taxLevel = "Genus", figdir, 
                              div_index = "shannon", return_plot = FALSE){
  tax_name = case_when(taxLevel == "Genus" ~ "genera",
                       taxLevel == "Species" ~ "species",
                       taxLevel == "ASVs" ~ "ASVs",
                       TRUE ~ "taxa")
  fn <- paste0(datadir, "Tables/DietVsMicrobiomeTime/Alpha/Alpha", div_index, taxLevel, "vsDiet_DeltaFC_cutoff0.csv")
  t <- read.csv(fn)
  t$Nutrient = t[[nutrient]]
  #t <- t %>% filter(Cycle == "C2D1") # for exploration only
  
  t_bin = t
  t_bin$Nutrient <- case_when(t_bin$Nutrient < 0 ~ "Down", 
                              t_bin$Nutrient > 0 ~ "Up")
  t_bin = t_bin %>% filter(!is.na(Nutrient))
  p2 <- ggplot(t_bin, aes(x = Nutrient, y = Alpha)) + 
    geom_boxplot(aes(fill = Nutrient), alpha = 0.2, outlier.size = -1) + 
    geom_point(aes(color = Nutrient)) + 
    theme_pubr() + 
    xlab(bquote(Delta ~ .(nutrient_name))) + 
    stat_compare_means(label.x = 1.1, label.y = 1.05*max(t_bin$Alpha), 
                       method = "wilcox", 
                       aes(label = paste0("italic(p) ==", after_stat(p.format))),
                       parse = TRUE) + 
    ylim(c(min(t_bin$Alpha), 1.15*max(t_bin$Alpha))) +
    theme(legend.position = "none") + 
    scale_color_manual(values = c("Up" = "blue", "Down" = "orange")) + 
    scale_fill_manual(values = c("Up" = "blue", "Down" = "orange")) 
  
  if (div_index == "shannon"){p2 = p2+ylab(expression(Delta ~ alpha ~ "diversity (Shannon)"))}
  if (div_index == "observed"){p2 = p2+ylab(bquote(Delta ~ alpha ~ "(" * observed ~ .(tax_name) * ")"))}
  
  fn <- paste0(figdir, "AlphaBoxplot", nutrient, "_", div_index, ".pdf")
  ggsave(filename = fn, plot = p2, width = 2, height = 3)
  
  t_plot = t %>% filter(abs(Nutrient) > max(abs(Nutrient))/100)
  p3 = ggplot(t_plot, aes(x = Nutrient, y = Alpha)) + 
    geom_point(color = "black") + 
    theme_pubr() + 
    geom_smooth(method = "lm", color = "black") + 
    stat_cor(method = "pearson", label.x = min(t_plot$Nutrient), label.y = max(t_plot$Alpha)) + 
    xlab(bquote(Delta ~ .(nutrient_name) ~ "(log"[2]*"FC)")) + 
    ylim(min(t_plot$Alpha), 1.05*max(t_plot$Alpha))
  
  if (div_index == "shannon"){p3 = p3+ylab(expression(Delta ~ alpha ~ "diversity (Shannon)"))}
  if (div_index == "observed"){p3 = p3+ylab(bquote(Delta ~ alpha ~ "(" * observed ~ .(tax_name) * ")"))}
  
  if (return_plot){return(p3)}
  
  fn <- paste0(figdir, "AlphaCorrplot", nutrient, "_", div_index, ".pdf")
  ggsave(filename = fn, plot = p3, width = 3, height = 3)
}

plot_alpha_bargraph <- function(taxLevel = "Genus", n_hits = 10, figdir = figdir, fdr_cut = 0.1, 
                                div_index = "shannon", h = 3, w = 4, ylabel = ""){
  fn <- paste0(datadir, "Tables/DietVsMicrobiomeTime/Alpha/", div_index, taxLevel, "_Res_cutoff0.csv")
  t <- read.csv(fn)
  
  t <- t %>% filter(Norm == "RankDiff") %>% arrange(p.value)
  t <- t %>% filter(!grepl("HEI_|AHEI", Var))
  t = t %>% filter(!(Var %in% feature_exclude_dd))
  
  t$FDR <- p.adjust(t$p.value, method = "BH")
  
  # Plot bargraph
  df_sorted <- t %>%
    arrange(FDR) %>%
    head(n_hits)
  
  df_sorted$Significance <- case_when(df_sorted$FDR < fdr_cut ~ "FDR < 0.2", 
                                      TRUE ~ "n.s.")
  
  # Need to find a better way to code this section
  df_sorted = rename_var(df_sorted)
  
  ylabs <- setNames(df_sorted$Var, df_sorted$Var)
  ylabs["Vitamin B12"] <- "Vitamin~B[12]"
  
  p <- ggplot(df_sorted, aes(x = Estimate, y = reorder(Var, Estimate), fill = Significance)) +
    geom_col() +
    scale_fill_manual(values = c("FDR < 0.2" = "purple", "n.s." = "grey"), name = NULL) + 
    #xlab(expression("-log10(" * italic("p") * ")")) + 
    xlab(expression(rho * " (" * Delta ~ "diet vs " * Delta ~ alpha ~ "diversity)")) + 
    theme_pubr() + 
    scale_y_discrete(labels = function(x) {
      # Parse only the one that needs it
      ifelse(x == "Vitamin B12", 
             parse(text = "Vitamin~B[12]"), 
             x)
    }) + 
    ylab(ylabel)
  
  fn <- paste0(figdir, "AlphaBargraph_", taxLevel, "_", div_index, ".pdf")
  ggsave(filename = fn, plot = p, width = w, height = h)
}

plot_beta_bargraph = function(taxLevel = "Genus", n_hits = 10, figdir, fdr_cut = 0.1, ylabel = "", h = 3, w = 4){
  fn <- paste0(datadir, "Tables/DietVsMicrobiomeTime/Beta/", taxLevel, "_Res_cutoff0.csv")
  t <- read.csv(fn)
  
  t <- t %>% filter(Norm == "RankDiff")
  t <- t %>% filter(!(grepl("AHEI|HEI_", Var))) %>% arrange(p.value)
  t = t %>% filter(!(Var %in% feature_exclude_dd))
  t$FDR <- p.adjust(t$p.value, method = "BH")
  
  # Plot bargraph
  df_sorted <- t %>%
    arrange(FDR) %>%
    head(n_hits)
  
  df_sorted$Significance <- case_when(df_sorted$FDR < fdr_cut ~ "FDR < 0.2", 
                                      TRUE ~ "n.s.")
  
  
  df_sorted <- rename_var(df_sorted) 
  
  p <- ggplot(df_sorted, aes(x = Estimate, y = reorder(Var, Estimate), fill = Significance)) +
    geom_col() +
    scale_fill_manual(values = c("FDR < 0.2" = "purple", "n.s." = "grey"), name = NULL) + 
    #xlab(expression("-log10(" * italic("p") * ")")) + 
    xlab(expression(rho * " (|" * Delta ~ "diet| vs " * Delta ~ beta ~ "diversity)")) + 
    ylab(ylabel) + 
    theme_pubr()
  
  fn <- paste0(figdir, "BetaBargraph_", taxLevel, ".pdf")
  ggsave(filename = fn, plot = p, width = w, height = h)
}



## Plot change in diet vs change in taxa/gene 
# useMedian: whether to use the median or 0 as the cutoff for up/down
# raw_genes: whether to load raw genes (RPKG) or normalize to 1
# clr_genes: whether to clr-transform the genes 
# set t_filt = t_norm, figdir = figdir
plot_dd_taxa_diet = function(t_filt, figdir, var = "D_YOGURT", var_name = "Yogurt", 
                             taxa = "fd496fd32dc8c08ade2e8b6c9d8ee13d", taxa_name = "test",
                             taxLevel = "ASVs", logT_diet = TRUE, nsamp = 20, cutoff = 0.0005, 
                             useMedian = FALSE, raw_genes = FALSE, clr_genes = TRUE, samp_exclude = "",
                             w = 2.75, return_df = FALSE, return_raw_df = FALSE, return_plot = FALSE, 
                             method = "pearson", wbox = 2){
  if (grepl("Pathway|KO", taxLevel)){cutoff = 0}
  
  # Load data
  if (grepl("Pathway|KO", taxLevel)){
    df <- load_mgs_genes(taxLevel = taxLevel, raw = raw_genes, clr = clr_genes, nsamp = nsamp, cutoff = cutoff) %>% t()
    if((!clr_genes) | raw_genes){df <- log10(df + 0.0001)}
  } else if (grepl("Fuso|Fn", taxLevel)){
    df = load_fuso()
  } else if (grepl("MGSGenus", taxLevel)){
    df = load_mgs_taxa(taxLevel = "Genus", cutoff = 0, nsamp = 0, clr = TRUE) %>% t()
  } else {
    df <- load_16s(taxLevel = taxLevel, clr = TRUE, nsamp = nsamp, cutoff = cutoff) %>% t()
  }
  samp_both = intersect(t_filt$Sample_ID, rownames(df))
  
  df_pt <- t_filt %>% filter(Sample_ID %in% samp_both)
  df_tax = df[samp_both, ]
  df_pt$Taxa = df_tax[,taxa]
  df_pt$Var <- df_pt[[var]] 
  df_pt <- df_pt %>% filter(!(Sample_ID %in% samp_exclude))
  
  if (return_raw_df){return(df_pt)}
  
  ### Plots
  
  # Plot 1
  df_pt_plot = df_pt
  if (logT_diet){
    minval = min(df_pt_plot[[var]][df_pt_plot[[var]] > 0], na.rm = TRUE)/2
    df_pt_plot$Var <- log2(df_pt_plot[[var]] + minval)
  } 
  p4 <- ggplot(df_pt_plot, aes(x = Var, y = Taxa)) + 
    geom_point(color = "black") + 
    #geom_path(aes(group = PatientID), arrow = arrow(type = "closed", length = unit(0.1, "inches"))) + 
    xlab(var_name) + 
    ylab("Abundance (CLR)") + 
    theme_pubr() +
    geom_smooth(method = "lm", level = 0.95, color = "black") + 
    stat_cor(method = "spearman", label.x = min(df_pt_plot$Var), label.y = max(df_pt_plot$Taxa))
  
  # Delta delta plot
  pt_keep = df_pt %>% group_by(PatientID) %>% dplyr::summarise(n = n()) %>% filter(n > 1) %>% dplyr::select(PatientID) %>% pull()
  
  taxa_threshold = (max(df_pt$Taxa) - min(df_pt$Taxa))/100
  df_pt$Taxa[df_pt$Taxa < taxa_threshold] = taxa_threshold
  
  diet_threshold = (max(df_pt$Var) - min(df_pt$Var))/100
  df_pt$Var[df_pt$Var < diet_threshold] = diet_threshold
  
  if (logT_diet){
    minval = min(df_pt[[var]][df_pt[[var]] > 0], na.rm = TRUE)/2
    df_pt$Var <- log2(df_pt[[var]] + minval)
  }
  
  df_pair = df_pt %>% 
    filter(PatientID %in% pt_keep) %>% 
    dplyr::group_by(PatientID) %>%
    arrange(Cycle) %>% 
    dplyr::reframe(
      Cycle1 = combn(Cycle, 2, FUN = function(x) x[1]),
      Cycle2 = combn(Cycle, 2, FUN = function(x) x[2]),
      DeltaTaxa = combn(Taxa, 2, FUN = function(x) x[2] - x[1]),
      DeltaVar = combn(Var, 2, FUN = function(x) x[2] - x[1])
    )
  df_pair = df_pair %>% filter(abs(DeltaTaxa) > 0) %>% filter(abs(DeltaVar) > 0)
  df_pair$DeltaTaxa = log2(exp(df_pair$DeltaTaxa)) # convert to log2FC
  
  df_pair <- df_pair %>% filter((Cycle1 == "C1D1" & Cycle2 == "C2D1") | (Cycle1 == "C2D1" & Cycle2 == "C3D1") | (Cycle1 == "C3D1" & Cycle2 == "EOT")) # filter (optional, this is to match the model)
  
  if (return_df){return(df_pair)}
  
  p5 <- ggplot(df_pair, aes(x = DeltaVar, y = DeltaTaxa)) + 
    geom_point(color = "black") + 
    ylab(bquote(Delta ~ italic(.(taxa_name)) ~ "(" * log[2] * "FC)")) +
    xlab(bquote(Delta ~ .(var_name))) + 
    theme_pubr() + 
    geom_smooth(method = "lm", level = 0.95, color = "black") + 
    stat_cor(method = method, label.x = min(df_pair$DeltaVar), label.y = 1.25*max(df_pair$DeltaTaxa))
  
  if (grepl("Pathway|KO", taxLevel)){
    p5 = p5 + ylab(bquote(Delta ~ .(taxa_name) ~ "(" * log[2] * "FC)")) 
  }
  
  # Boxplot
  df_pair$DeltaVarCat = case_when(df_pair$DeltaVar > 0 ~ "Up",
                                  TRUE ~ "Down")
  p6 <- ggplot(df_pair, aes(x = DeltaVarCat, y = DeltaTaxa)) + 
    geom_point(aes(color = DeltaVarCat)) + 
    geom_boxplot(aes(fill = DeltaVarCat), alpha = 0.2) + 
    scale_color_manual(values = c("orange", "blue")) + 
    scale_fill_manual(values = c("orange", "blue")) + 
    xlab(bquote(Delta~.(var_name))) + 
    ylab(bquote(Delta ~ italic(.(taxa_name)) ~ "(" * log[2] * "FC)")) +
    theme_pubr() + 
    stat_compare_means(method = "wilcox.test", label.x = 1.2, label.y = 1.15*max(df_pair$DeltaTaxa),
                       aes(label = paste0("italic(p) ==", after_stat(p.format))), parse = TRUE) + 
    ylim(min(df_pair$DeltaTaxa), 1.25*max(df_pair$DeltaTaxa)) + 
    theme(legend.position = "none")
  
  if (logT_diet){
    p4 <- p4 + xlab(bquote(.(var_name) ~ "(" * log[2] * ")"))
    p5 <- p5 + xlab(bquote(Delta ~ .(var_name) ~ "(" * log[2] * "FC)"))
  }
  
  if(return_plot){return(p5)}
  
  # Save
  fn = paste0(figdir, "DD_Raw_", taxLevel, "_", var, "_", taxa, ".pdf")
  ggsave(filename = fn, plot = p4, width = 3, height = 3)
  
  fn = paste0(figdir, "DD_Delta_", taxLevel, "_", var, "_", taxa, ".pdf")
  ggsave(filename = fn, plot = p5, width = w, height = 3)
  
  fn = paste0(figdir, "DD_Box_", taxLevel, "_", var, "_", taxa, ".pdf")
  ggsave(filename = fn, plot = p6, width = wbox, height = 3)
}

# Plot timecourse, stratified by change C1D1
plot_time_taxa_diet = function(t_filt, figdir, var = "D_YOGURT", var_name = "Yogurt", 
                               taxLevel = "ASVs", taxa = "fd496fd32dc8c08ade2e8b6c9d8ee13d", 
                               taxa_name = "S. thermophilus", cutoff = 0.0005, 
                               useMedian = FALSE, 
                               cycle_keep = c("C1D1", "C1D3", "C1D7", "C2D1", "C3D1"),
                               plot_lines = FALSE, 
                               c1 = "C1D1", c2 = "C2D1", 
                               exclude_zeros = TRUE, return_df = FALSE, 
                               ylabel = NULL, 
                               return_plot = FALSE){
  # Read in data
  if (grepl("Pathway|KO", taxLevel)){cutoff = 0}
  
  if (grepl("KO|Pathway", taxLevel)){
    t <- load_mgs_genes(taxLevel = taxLevel, clr = TRUE, raw = FALSE)
    if (sum(t) > 1){t <- log2(t+0.0001)} else {t <- log2(exp(t))}
  } else if (grepl("Fuso|Fn", taxLevel)){
    t = load_fuso() %>% t()
  } else {
    t <- load_16s(taxLevel = taxLevel, clr = TRUE)
  }
  t_all <- t
  
  t <- t_all[taxa,]
  tax <- data.frame("Sample_ID" = names(t), "Abundance" = as.numeric(t))
  tax <- tax %>% filter(!grepl("9c", tax$Sample_ID)) # Patient 9 = Patient 9b
  tax <- tax %>% filter(!(grepl("Batch|Zymo|duplicate|Tox|blank", Sample_ID))) # remove controls before regex
  tax$PatientID <- as.numeric(unlist(regmatches(tax$Sample_ID, gregexpr("(?<=Pt)\\d+", tax$Sample_ID, perl = TRUE))))
  tax$Cycle <- gsub(".*-", "", tax$Sample_ID)
  tax <- tax %>% filter(Cycle %in% cycle_keep)
  
  # Filter only for patients that have above 0 CLR in at least one samples and at least 2 total samples
  pt_keep = tax %>% group_by(PatientID) %>% 
    dplyr::summarise(s = sum(Abundance > 0)) %>% 
    filter(s > 0) %>% dplyr::select(PatientID) %>% pull()
  if (exclude_zeros){
    tax <- tax %>% filter(PatientID %in% pt_keep)
  }
  pt_keep_2 = tax %>% group_by(PatientID) %>% 
    dplyr::summarise(n = n()) %>% filter(n >= 2) %>% dplyr::select(PatientID) %>% pull()
  tax <- tax %>% filter(PatientID %in% pt_keep_2)
  
  # Segregate by T2-T1 diet, plot 
  t_c12 <- t_filt %>% filter(Cycle %in% c(c1, c2))
  t_c12$Var = t_c12[[var]]
  pt_keep <- t_c12 %>% group_by(PatientID) %>% dplyr::summarise(n = n()) %>% filter(n == 2) %>% dplyr::select(PatientID) %>% pull()
  t_c12 <- t_c12 %>% filter(PatientID %in% pt_keep)
  t_delta <- t_c12 %>%
    dplyr::group_by(PatientID) %>%
    dplyr::summarise(Difference = Var[Cycle == c2] - Var[Cycle == c1], .groups = "drop") %>% 
    as.data.frame()
  if (useMedian){
    t_delta$IntakeShift = case_when(t_delta$Difference >= median(t_delta$Difference) ~ "> median", #median(t_delta$Difference)
                                    TRUE ~ "< median")
  } else {
    t_delta$IntakeShift = case_when(t_delta$Difference >= 0 ~ "Increase", 
                                    TRUE ~ "Decrease")
  }
  
  
  tax_merged <- left_join(tax, t_delta, by = "PatientID")
  
  tax_merged <- tax_merged %>% filter(!is.na(IntakeShift))
  tax_merged <- tax_merged %>% # baseline-normalize (C1D1)
    group_by(PatientID) %>%
    mutate(Abundance_Orig = Abundance,
           Abundance = Abundance - Abundance[Cycle == "C1D1"]) %>% 
    ungroup()
  
  if (return_df){return(tax_merged)}
  
  if (grepl("KO|Pathway", taxLevel)){
    tax_merged <- tax_merged %>% filter(Cycle != "Baseline") # filter(Cycle %in% c("C1D1", "C1D7", "C2D1"))
  } else {
    tax_merged <- tax_merged %>% filter(Cycle != "Baseline")
  }
  
  legend_title = var_name # paste0(var_name, "\n(C2-C1)")
  tax_merged <- tax_merged %>% filter(abs(Abundance) < 10)
  p <- ggplot(tax_merged, aes(x = Cycle, y = Abundance, color = IntakeShift, fill = IntakeShift)) + 
    theme_pubr() +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) + 
    stat_summary(fun = mean, geom = "line", aes(group = IntakeShift), linewidth = 1) +  # Connect means
    scale_color_manual(values = c("< median" = "orange", "> median" = "blue", "Decrease" = "orange", "Increase" = "blue")) + 
    scale_fill_manual(values = c("< median" = "orange", "> median" = "blue", "Decrease" = "orange", "Increase" = "blue")) + 
    labs(color = legend_title, fill = legend_title)
  if (is.null(ylabel)){
    p = p + ylab(bquote(italic(.(taxa_name)) ~ "vs C1D1 (CLR)")) 
  } else {
    p = p + ylab(ylabel)
  }
  
  if(plot_lines){p = p + geom_line(aes(group = PatientID), alpha = 0.3)}
  if(return_plot){return(p)}
  
  # Save
  fn = paste0(figdir, "DD_Time_", taxLevel, "_", var_name, "_", taxa_name, ".pdf")
  ggsave(filename = fn, plot = p, width = 3.5, height = 3)
}

# Plot NE. method = "nutrient" or "food_group"
plot_ne_medi_box = function(method = "nutrient", axis_name = "Carbohydrate", names = c("CARB"), kcal_norm = TRUE, logT = FALSE, figdir){
  if (method == "food_group"){
    medi = load_medi_food_ne(method = method, names = names)
  } else {
    medi = load_medi_nutrients_ne(names = names, kcal_norm = kcal_norm)
    medi$Abundance = medi[[names]]
  }
  
  medi$Abundance[is.na(medi$Abundance)] = 0 
  if (logT){
    minval = min(medi$Abundance[medi$Abundance > 0], na.rm = TRUE)/2
    medi$Abundance = log10(medi$Abundance + minval)
  }
  
  medi$CycleInt = case_when(medi$Cycle == "T1" ~ 1,
                            medi$Cycle == "T2" ~ 2, 
                            TRUE ~ 3)
  medi$CycleInt = as.factor(medi$CycleInt)

  # Keep pt with all 3
  pt_keep = names(table(medi$Patient_ID)[table(medi$Patient_ID) >= 3])
  medi = medi %>% filter(Patient_ID %in% pt_keep)
  
  # Stats
  model <- nlme::lme(Abundance ~ CycleInt, random = ~1 | Patient_ID, data = medi)
  emm <- emmeans::emmeans(model, "CycleInt")
  contrast_res <- contrast(emm, method = "pairwise", adjust = "none")
  contrast_df <- as.data.frame(contrast_res)
  
  
  stat_df <- contrast_df %>%
    tidyr::separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
    mutate(
      label = formatC(p.value/2, format = "fg", digits = 2, flag = "#") # one-sided testing
    ) %>% 
    dplyr::filter(!(group1 == "T2" & group2 == "T3"))
  y.position = c(1.05, 1.15)*max(medi$Abundance, na.rm = TRUE)
  
  # Plot
  p = ggplot(medi, aes(x = CycleInt, y = Abundance, group = Cycle)) +
    geom_point(aes(color = Cycle)) + 
    geom_boxplot(outlier.size = -1, aes(fill = Cycle), alpha = 0.2) + 
    xlab("Timepoint") + 
    theme_pubr()
  
  # Label
  if (method == "food_group"){
    if (logT){p = p + ylab(paste0(axis_name, " abundance (log10 rpm)"))}else{p = p + ylab(paste0(axis_name, " abundance (rpm)"))}
  } else if (names %in% c("CARB", "TFAT", "PROT")){
    p = p + ylab(paste0(axis_name, " (% kcal)"))
  } else {
    if(logT){p = p + ylab(paste0(axis_name, " (log10 mg/100g)"))}else{p = p + ylab(paste0(axis_name, " (mg/100g)"))}
  }
  
  # Add stats
  num_sig = length(stat_df$label)
  if (num_sig > 0){
    y.position = y.position[1:num_sig]
    stat_df$y.position = y.position
    stat_df$group1 <- gsub("CycleInt", "", stat_df$group1)
    stat_df$group2 <- gsub("CycleInt", "", stat_df$group2)
    p <- p + stat_pvalue_manual(stat_df, label = "label", tip.length = 0.01) +
      coord_cartesian(ylim = c(min(medi$Abundance), max(medi$Abundance) * (1 + 0.1*num_sig)))
  }
  
  # Add color
  if (contrast_df[1,2] > 0){
    p <- p + scale_color_manual(values = c("grey", "orange", "orange")) + 
      scale_fill_manual(values = c("grey", "orange", "orange")) + 
      theme(legend.position = "none") 
  } else {
    p <- p + scale_color_manual(values = c("grey", "blue", "blue")) + 
      scale_fill_manual(values = c("grey", "blue", "blue")) + 
      theme(legend.position = "none") 
  }
  
  # Save
  fn = paste0(figdir, "NE_", axis_name, ".pdf")
  ggsave(filename = fn, plot = p, width = 2, height = 3)
}


### Rerun alpha, rank only
alpha_dd_rank <- function(taxLevel = "ASVs", div_index = "shannon", nsamp = 20, cutoff = 0.0005, 
                     vars = c("HEI", "AHEI"), t_filt, return_alpha = FALSE){
  # Load data, calculate shannon index
  if (taxLevel %in% c("KO", "Pathway")){
    df <- load_mgs_genes(taxLevel = taxLevel, clr = FALSE, nsamp = nsamp, cutoff = cutoff) %>% t()
  } else {
    df <- load_16s(taxLevel = taxLevel, clr = FALSE, nsamp = nsamp, cutoff = cutoff) %>% t() 
  }
  
  if (div_index == "observed"){
    alpha_calculated = rowSums(df > 0)
  } else {
    alpha_calculated <- vegan::diversity(df, index = div_index)
  }
  
  alpha_df <- data.frame(Sample_ID = names(alpha_calculated), Alpha = as.numeric(alpha_calculated))
  alpha_df <- inner_join(t_filt, alpha_df, by = "Sample_ID") %>% dplyr::select(all_of(c(vars, "PatientID", "Alpha", "Cycle")))
  
  if (return_alpha){
    return(alpha_df)
  }
  
  alpha_df <- alpha_df %>% mutate(Cycle = factor(Cycle, levels = sort(unique(Cycle))))
  alpha_df <- alpha_df %>% mutate(PatientID = factor(PatientID, levels = sort(unique(PatientID))))
  
  # Set threshold that allows 0 change in diet
  for (var in vars){
    threshold = (max(alpha_df[[var]]) - min(alpha_df[[var]]))/100
    alpha_df[[var]][alpha_df[[var]] < threshold] = threshold
  }
  
  delta_rank_in <- alpha_df %>%
    mutate(across(.cols = where(is.numeric), .fns = ~ rank(.), .names = "{.col}")) %>% 
    arrange(PatientID, Cycle) %>%  # Ensure correct ordering
    group_by(PatientID) %>%
    mutate(across(where(is.numeric), ~ . - lag(.), .names = "{.col}")) %>%
    dplyr::filter(!is.na(HEI))  # Remove first row per patient (no prior cycle to compare)
  
  res = data.frame(Var = character(0), Norm = character(0), p.value = numeric(0), Estimate = numeric(0))
  for (var in vars){
    # Remove 0 comparisons
    delta_rank = delta_rank_in
    delta_rank$Var = delta_rank[[var]]
    delta_rank = delta_rank %>% filter(abs(Var) > 0)

    est <- cor.test(delta_rank$Alpha, delta_rank[[var]], method = "spearman")$estimate[[1]]
    p <- cor.test(delta_rank$Alpha, delta_rank[[var]], method = "spearman")$p.value[[1]]
    row = data.frame(Var = var, Norm = "RankDiff", p.value = p, Estimate = est)
    res <- rbind(res, row)
  }
  
  res <- res %>% arrange(p.value)
  
  fn = paste0(table_mt_alpha, div_index, taxLevel, "_Res_cutoff", cutoff, ".csv")
  write.csv(res, file = fn)
  
  fn = paste0(table_mt_alpha, div_index, taxLevel, "_Raw_cutoff", cutoff, ".csv")
  write.csv(delta_rank_in, file = fn)
}

# Beta delta delta function - rank only
beta_dd_rank <- function(taxLevel = "ASVs", nsamp = 20, cutoff = 0.0005, vars = c("HEI", "AHEI"), t_filt, min_clr = 0){
  if (taxLevel %in% c("KO", "Pathway")){
    df <- load_mgs_genes(taxLevel = taxLevel, clr = TRUE, nsamp = nsamp, cutoff = cutoff, min_clr = min_clr) %>% t()
  } else {
    df <- load_16s(taxLevel = taxLevel, clr = TRUE, nsamp = nsamp, cutoff = cutoff, min_clr = min_clr) %>% t() 
  }
  dist_matrix <- as.matrix(dist(df, method = "euclidean"))
  
  beta_df <- t_filt %>% dplyr::select(all_of(c(vars, "PatientID", "Cycle")))
  beta_df <- beta_df %>% mutate(Cycle = factor(Cycle, levels = sort(unique(Cycle))))
  beta_df <- beta_df %>% mutate(PatientID = factor(PatientID, levels = sort(unique(PatientID))))
  
  # Set threshold that allows 0 change in diet
  for (var in vars){
    threshold = (max(beta_df[[var]]) - min(beta_df[[var]]))/100
    beta_df[[var]][beta_df[[var]] < threshold] = threshold
  }
  
  delta_rank <- beta_df %>%
    mutate(across(.cols = where(is.numeric), .fns = ~ rank(.), .names = "{.col}")) %>% # add  & !any_of("Alpha") for old behavior
    arrange(PatientID, Cycle) %>%  # Ensure correct ordering
    group_by(PatientID) %>%
    mutate(across(where(is.numeric), ~ . - lag(.), .names = "{.col}")) %>%
    dplyr::filter(!is.na(HEI))  # Remove first row per patient (no prior cycle to compare)
  
  delta_fc <- beta_df %>%
    arrange(PatientID, Cycle) %>%  # Ensure correct ordering
    group_by(PatientID) %>%
    mutate(
      across(
        .cols = where(is.numeric) & !any_of("Alpha"),
        .fns = ~ {
          # compute pseudocount: min positive value in column
          pos_vals <- .[. > 0 & !is.na(.)]
          pseudocount <- if (length(pos_vals) == 0) 1 else min(pos_vals)
          log2(. + pseudocount) - log2(lag(. + pseudocount))
        },
        .names = "{.col}"
      ),
      across(
        .cols = contains("Alpha"),
        .fns = ~ . - lag(.),
        .names = "{.col}"
      )
    ) %>%
    dplyr::filter(!is.na(HEI))  # Remove first row per patient
  
  delta_abs <- beta_df %>%
    arrange(PatientID, Cycle) %>%  # Ensure correct ordering
    group_by(PatientID) %>%
    mutate(across(where(is.numeric), ~ . - lag(.), .names = "{.col}")) %>%
    dplyr::filter(!is.na(HEI))  # Remove first row per patient (no prior cycle to compare)
  
  delta_rank$Cycle2 <- case_when(delta_rank$Cycle == "C2D1" ~ "C1D1", delta_rank$Cycle == "C3D1" ~ "C2D1", TRUE ~ "C3D1")
  delta_fc$Cycle2 <- case_when(delta_fc$Cycle == "C2D1" ~ "C1D1", delta_fc$Cycle == "C3D1" ~ "C2D1", TRUE ~ "C3D1")
  delta_abs$Cycle2 <- case_when(delta_abs$Cycle == "C2D1" ~ "C1D1", delta_abs$Cycle == "C3D1" ~ "C2D1", TRUE ~ "C3D1")
  
  delta_rank$Name1 <- paste0("GO-Pt", delta_rank$PatientID, "-", delta_rank$Cycle)
  delta_rank$Name2 <- paste0("GO-Pt", delta_rank$PatientID, "-", delta_rank$Cycle2)
  delta_fc$Name1 <- paste0("GO-Pt", delta_fc$PatientID, "-", delta_fc$Cycle)
  delta_fc$Name2 <- paste0("GO-Pt", delta_fc$PatientID, "-", delta_fc$Cycle2)
  delta_abs$Name1 <- paste0("GO-Pt", delta_abs$PatientID, "-", delta_abs$Cycle)
  delta_abs$Name2 <- paste0("GO-Pt", delta_abs$PatientID, "-", delta_abs$Cycle2)
  
  delta_rank_in = convert_delta(delta_rank, dist_matrix)
  delta_fc_in = convert_delta(delta_fc, dist_matrix)
  delta_abs_in = convert_delta(delta_abs, dist_matrix)
  
  res = data.frame(Var = character(0), Norm = character(0), p.value = numeric(0), Estimate = numeric(0))
  for (var in vars){
    # Remove 0 comparisons
    delta_rank = delta_rank_in
    delta_rank$Var = delta_rank[[var]]
    delta_rank = delta_rank %>% filter(abs(Var) > 0)
    
    est <- cor.test(delta_rank$Beta, abs(delta_rank[[var]]), method = "spearman")$estimate[[1]]
    p <- cor.test(delta_rank$Beta, abs(delta_rank[[var]]), method = "spearman")$p.value[[1]]
    row = data.frame(Var = var, Norm = "RankDiff", p.value = p, Estimate = est)
    res <- rbind(res, row)
  }
  
  res <- res %>% arrange(p.value)
  fn = paste0(table_mt_beta, taxLevel, "_Res_cutoff", cutoff, ".csv")
  write.csv(res, file = fn)
  
  fn = paste0(table_mt_beta, taxLevel, "_Raw_cutoff", cutoff, ".csv")
  write.csv(delta_rank_in, file = fn)
  
  fn = paste0(table_mt_beta, taxLevel, "_RawFC_cutoff", cutoff, ".csv")
  write.csv(delta_fc_in, file = fn)
  
  fn = paste0(table_mt_beta, taxLevel, "_RawAbs_cutoff", cutoff, ".csv")
  write.csv(delta_abs_in, file = fn)
}


# Features modeling - rank only
feature_dd_rank = function(taxLevel = "ASVs", nsamp = 20, cutoff = 0.0005, vars = c("HEI", "AHEI"), t_filt, min_clr = 0){
  if (taxLevel %in% c("KO", "Pathway")){
    df <- load_mgs_genes(taxLevel = taxLevel, clr = TRUE, nsamp = nsamp, cutoff = cutoff, min_clr = min_clr) %>% t()
  } else if (grepl("Fuso|Fn", taxLevel)){
    df = load_fuso()
  } else {
    df <- load_16s(taxLevel = taxLevel, clr = TRUE, nsamp = nsamp, cutoff = cutoff, min_clr = min_clr) %>% t() 
  }
  
  samp_both = intersect(t_filt$Sample_ID, rownames(df))
  df_pt <- t_filt %>% dplyr::filter(Sample_ID %in% samp_both)
  df_tax = df[samp_both, ]
  
  output = data.frame(Var = character(0), Taxa = character(0), Norm = character(0), 
                      p.value = numeric(0), Estimate = numeric(0), stringsAsFactors = FALSE)
  for (taxa in colnames(df_tax)){
    alpha_df <- data.frame(Sample_ID = rownames(df_tax), Alpha = as.numeric(df_tax[,taxa]))
    alpha_df <- inner_join(t_filt, alpha_df, by = "Sample_ID") %>% dplyr::select(all_of(c(vars, "PatientID", "Alpha", "Cycle")))
    alpha_df <- alpha_df %>% mutate(Cycle = factor(Cycle, levels = sort(unique(Cycle))))
    alpha_df <- alpha_df %>% mutate(PatientID = factor(PatientID, levels = sort(unique(PatientID))))

    # Set threshold that allows 0 change in diet, ASV
    for (var in vars){
      threshold = (max(alpha_df[[var]]) - min(alpha_df[[var]]))/100
      alpha_df[[var]][alpha_df[[var]] < threshold] = threshold
      threshold = (max(alpha_df[["Alpha"]]) - min(alpha_df[["Alpha"]]))/100
      alpha_df[["Alpha"]][alpha_df[["Alpha"]] < threshold] = threshold
    }
    
    delta_rank_in <- alpha_df %>%
      mutate(across(.cols = where(is.numeric), .fns = ~ rank(.), .names = "{.col}")) %>% # add  & !any_of("Alpha") to first argument for old behavior
      arrange(PatientID, Cycle) %>%  # Ensure correct ordering
      group_by(PatientID) %>%
      mutate(across(where(is.numeric), ~ . - lag(.), .names = "{.col}")) %>%
      dplyr::filter(!is.na(HEI))  # Remove first row per patient (no prior cycle to compare)
    
    for (var in vars){
      # Remove 0 comparisons
      delta_rank = delta_rank_in
      delta_rank$Var = delta_rank[[var]]
      delta_rank = delta_rank %>% filter(abs(Var) > 0) %>% filter(abs(Alpha) > 0)
      
      if(length(delta_rank$Alpha) >= 10){ # at least 10 observations
        res <- cor.test(delta_rank$Alpha, delta_rank[[var]], method = "spearman")
        est = res$estimate[[1]]
        p = res$p.value[[1]]
        row = data.frame(Var = var, Taxa = taxa, Norm = "RankDiff", p.value = p, Estimate = est)
        output <- rbind(output, row)
      }
    }
  }
  
  output <- output %>% arrange(Norm, p.value)
  
  fn = paste0(table_mt_granular, taxLevel, "vsDiet_DeltaDelta_cutoff", cutoff, ".csv")
  write.csv(output, file = fn)
}

## Rename variables
rename_var = function(df_sorted){
  df_sorted$Var <- case_when(df_sorted$Var == "THEO" ~ "Theobromine", 
                             df_sorted$Var == "Metadata_Cancer_stage" ~ "Cancer stage", 
                             df_sorted$Var == "CALC" ~ "Calcium", 
                             df_sorted$Var == "S180" ~ "SFA 18:0", 
                             df_sorted$Var == "V_TOTAL" ~ "Total vegetables", 
                             df_sorted$Var == "Metadata_Colostomy" ~ "Colostomy", 
                             df_sorted$Var == "Metadata_Ethnicity" ~ "Ethnicity", 
                             df_sorted$Var == "Metadata_Sex" ~ "Sex", 
                             df_sorted$Var == "P226" ~ "PUFA 22:6", 
                             df_sorted$Var == "CHOLE" ~ "Cholesterol", 
                             df_sorted$Var == "Metadata_Procedure_Cat" ~ "Prior surgery", 
                             df_sorted$Var == "Metadata_Concurrent_Tx" ~ "Concurrent treatment", 
                             df_sorted$Var == "CARB" ~ "Carbohydrates", 
                             df_sorted$Var == "G_REFINED" ~ "Refined grains", 
                             df_sorted$Var == "VITD" ~ "Vitamin D",
                             df_sorted$Var == "V_DRKGR" ~ "Dark, leafy greens", 
                             df_sorted$Var == "PF_TOTAL" ~ "Protein foods", 
                             df_sorted$Var == "PF_MPS_TOTAL" ~ "Meat, poultry, seafood", 
                             df_sorted$Var == "F_JUICE" ~ "Fruit juice", 
                             df_sorted$Var == "OILS" ~ "Oils", 
                             df_sorted$Var == "F_TOTAL" ~ "Total fruits", 
                             df_sorted$Var == "ADD_SUGARS" ~ "Added sugars", 
                             df_sorted$Var == "CAFF" ~ "Caffeine", 
                             df_sorted$Var == "V_OTHER" ~ "Other vegetables", 
                             df_sorted$Var == "M201" ~ "MUFA 20:1", 
                             df_sorted$Var == "PFAT" ~ "Polyunsaturated fat", 
                             df_sorted$Var == "VB6" ~ "Vitamin B6", 
                             df_sorted$Var == "ZINC" ~ "Zinc", 
                             df_sorted$Var == "FA" ~ "Folic acid", 
                             df_sorted$Var == "VB12" ~ "Vitamin B12", 
                             df_sorted$Var == "PF_EGGS" ~ "Eggs", 
                             df_sorted$Var == "D_CHEESE" ~ "Cheese", 
                             df_sorted$Var == "PF_POULT" ~ "Poultry", 
                             df_sorted$Var == "PHOS" ~ "Phosphorus", 
                             df_sorted$Var == "V_STARCHY_TOTAL" ~ "Starchy vegetables", 
                             df_sorted$Var == "VC" ~ "Vitamin C",
                             df_sorted$Var == "HEI" ~ "HEI", 
                             df_sorted$Var == "KCAL" ~ "Energy (kcal)", 
                             df_sorted$Var == "AHEI" ~ "AHEI", 
                             df_sorted$Var == "FIBE" ~ "Fiber", 
                             df_sorted$Var == "SFAT" ~ "Saturated fat", 
                            df_sorted$Var == "ALC" ~ "Alcohol", 
                            df_sorted$Var == "S160" ~ "SFA 16:0", 
                            df_sorted$Var == "LYCO" ~ "Lycopene", 
                             df_sorted$Var == "IRON" ~ "Iron",
                            df_sorted$Var == "F_CITMLB" ~ "CMB fruits",
                            df_sorted$Var == "D_MILK" ~ "Milk",
                            df_sorted$Var == "SODI" ~ "Sodium",
                            df_sorted$Var == "SELE" ~ "Selenium",
                            df_sorted$Var == "PF_LEGUMES" ~ "Legumes",
                            df_sorted$Var == "V_REDOR_OTHER" ~ "Red/orange vegetables",
                            df_sorted$Var == "VARA" ~ "Vitamin A",
                            df_sorted$Var == "PF_NUTSDS" ~ "Nuts",
                            df_sorted$Var == "MAGN" ~ "Magnesium",
                            df_sorted$Var == "D_YOGURT" ~ "Yogurt",
                            df_sorted$Var == "CRYP" ~ "Cryptoxanthine",
                            df_sorted$Var == "ATOC" ~ "Vitamin E",
                            df_sorted$Var == "VK" ~ "Vitamin K1",
                            df_sorted$Var == "FF" ~ "Free folate",
                            df_sorted$Var == "V_REDOR_TOMATO" ~ "Tomatoes",
                            df_sorted$Var == "TFAT" ~ "Fat",
                            df_sorted$Var == "COPP" ~ "Copper",
                            df_sorted$Var == "G_WHOLE" ~ "Whole grains",
                            df_sorted$Var == "SUGR" ~ "Sugar",
                            df_sorted$Var == "POTA" ~ "Potassium",
                            df_sorted$Var == "S040" ~ "Butanoic acid",
                            df_sorted$Var == "P183" ~ "PUFA 18:3",
                            df_sorted$Var == "B12_ADD" ~ "Added B12",
                            df_sorted$Var == "M161" ~ "MUFA 16:1",
                             TRUE ~ df_sorted$Var)
  return(df_sorted)
}


plot_pcoa = function(dist, t_dist, var = "CHOLE", title = "Cholesterol intake",
                     ellipse = TRUE, max2 = 10000, quartiles = FALSE, norm = TRUE, 
                     axis1 = 1, axis2 = 2){
  pcoa_res <- ape::pcoa(dist)
  
  # Extract first 4 coordinates
  pcoa_df <- data.frame(
    Sample = rownames(pcoa_res$vectors),
    Axis1 = pcoa_res$vectors[,1],
    Axis2 = pcoa_res$vectors[,2],
    Axis3 = pcoa_res$vectors[,3],
    Axis4 = pcoa_res$vectors[,4]
  )
  
  # Merge with metadata
  pcoa_df <- pcoa_df %>%
    left_join(t_dist, by = "Sample")
  
  # % variance explained
  var_exp <- round(100 * pcoa_res$values$Relative_eig[axis1:axis2], 1)
  
  # Plot with ggplot2
  if (norm){
    pcoa_df$Var_Cat = case_when(pcoa_df[[var]] <= median(pcoa_df[[var]]) ~ "Low", TRUE ~ "Hi")
    if (quartiles){
      pcoa_df <- pcoa_df %>%
        mutate(
          Var_Cat = ntile(.data[[var]], 4),       # 4 = quartiles; 3 = tertiles
          Var_Cat = factor(
            Var_Cat,
            levels = 1:4,
            labels = c("Q1", "Q2", "Q3", "Q4")
          )
        )    
    }
  } else {pcoa_df$Var_Cat = pcoa_df[[var]]}

  if (max2 < 10000){print(pcoa_df %>% filter(Axis2 > max2) %>% dplyr::select(Sample) %>% pull())}
  
  p <- ggplot(pcoa_df %>% filter(Axis2 < max2), aes(x = .data[[paste0("Axis", axis1)]], y = .data[[paste0("Axis", axis2)]], color = Var_Cat)) +
    geom_point(size = 2, alpha = 0.8) +
    labs(
      x = paste0("PCoA", axis1, " (", var_exp[1], "%)"),
      y = paste0("PCoA", axis2, " (", var_exp[2], "%)"), 
      color = title) + 
    scale_color_manual(values = c("Hi" = "blue", "Low" = "darkorange", 
                                  "Yes" = "blue", "No" = "darkorange", 
                                  "Yes surgery" = "blue", "No surgery" = "darkorange", 
                                  "Q1" = "lightblue", "Q2" = "blue", "Q3" = "darkblue", "Q4" = "black")) + 
    theme_pubr() 
  if (ellipse){p = p + stat_ellipse(level = 0.95, linetype = 2)}
  
  return(p)
}

# Plot volcano by diet
plot_volcano_by_diet = function(t_in, figdir, diet_var = "COPP", diet_name = "Copper", h = 3, w = 4, return_df = FALSE){
  var = diet_var
  if (!(var %in% c("HEI", "AHEI"))){
    diet_name <- tools::toTitleCase(tolower(diet_name))
  }
  taxa_name = taxaLevel
  if(taxa_name == "ASVs"){taxa_name = "ASV"}
  
  # Filter, FDR correct
  v = t_in %>% filter(Var == var)
  v$FDR <- p.adjust(v$p.value, method = "BH")
  v$Association <- with(v, ifelse(FDR < fdr_cut & Estimate < 0, "(-)",
                                  ifelse(FDR < fdr_cut & Estimate > 0, "(+)",
                                         "n.s.")))
  v$Association = factor(v$Association, levels = c("n.s.", "(-)", "(+)"))
  
  # Plot
  p <- ggplot(v, aes(x = Estimate, y = -log10(p.value))) + 
    geom_point(aes(color = Association)) + 
    theme_pubr() + 
    xlab(bquote(rho ~ "(" * Delta ~ .(diet_name) ~ "vs " * Delta ~ "Gene)")) +
    ylab(expression("Significance (-log"[10]*italic(p)*")")) + 
    scale_color_manual(values = c("(+)" = "blue", "(-)" = "orange", "n.s." = "grey")) + 
    theme(legend.position = "top") + 
    labs(color = NULL)
  
  # Save
  fn = paste0(figdir, "Volcano", var, ".pdf")
  ggsave(plot = p, filename = fn, height = h, width = w)
  
  if (return_df){return(v)}
}

## Read in fuso qpcr data
read_cq <- function(date = "050223", projdir){
  # Set parameters
  loaddir <- paste0(projdir, "GOqPCR", date, "/")
  
  # Read in Cq file
  fn <- paste0(loaddir, "GO_qPCR_Fn_", date, "_QuantificationCqResults1000.csv") # Default single threshold quantification
  dat <- read.csv(fn)
  
  # Read in and format metadata file
  fn <- paste0(loaddir, "Layout_", date, ".csv")
  meta <- read.csv(fn) %>% as.data.frame
  meta <- reshape2::melt(meta, id = "X")
  meta$variable <- gsub("X", "", meta$variable)
  meta$Well <- paste0(meta$X, meta$variable)
  meta$Sample <- meta$value
  meta <- meta %>% dplyr::select(c("Well", "Sample"))
  meta[meta == ""] <- NA
  meta <- na.omit(meta)
  
  # Merge
  dat <- dat %>% dplyr::select(c("Well", "Cq"))
  df <- left_join(meta, dat, by = "Well") %>% as.data.frame()
  
  # Add additional parameters of interest
  df$Patient_ID <- gsub("-.*", "", df$Sample)
  df$Treatment_Cycle <- gsub(".*-", "", df$Sample)
  df$Primers <- case_when(grepl("01|02|03|04|05|06|07|08|09|10|11|12", df$Well) ~ "Fn1",
                          TRUE ~ "rRNA") 
  df <- df %>% dplyr::select(-"Well")
  df$Control <- case_when(grepl("Water|Fuso", df$Patient_ID) ~ "Control",
                          TRUE ~ "Sample")
  df$Concentration <- case_when(grepl("undil", df$Patient_ID) ~ 1,
                                grepl("e1", df$Patient_ID) ~ 10^(-1),
                                grepl("e2", df$Patient_ID) ~ 10^(-2),
                                grepl("e3", df$Patient_ID) ~ 10^(-3),
                                grepl("e4", df$Patient_ID) ~ 10^(-4),
                                grepl("e5", df$Patient_ID) ~ 10^(-5),
                                grepl("e6", df$Patient_ID) ~ 10^(-6),
                                grepl("e7", df$Patient_ID) ~ 10^(-7),
                                grepl("e8", df$Patient_ID) ~ 10^(-8),
                                grepl("Water", df$Patient_ID) ~ 0)
  df$Dilution <- case_when(grepl("undil", df$Patient_ID) ~ 0,
                           grepl("e1", df$Patient_ID) ~ -1,
                           grepl("e2", df$Patient_ID) ~ -2,
                           grepl("e3", df$Patient_ID) ~ -3,
                           grepl("e4", df$Patient_ID) ~ -4,
                           grepl("e5", df$Patient_ID) ~ -5,
                           grepl("e6", df$Patient_ID) ~ -6,
                           grepl("e7", df$Patient_ID) ~ -7,
                           grepl("e8", df$Patient_ID) ~ -8,
                           grepl("Water", df$Patient_ID) ~ -10)
  
  # Summarize
  diff <- df %>% 
    distinct() %>% 
    tidyr::pivot_wider(names_from = Primers, values_from = Cq, values_fn = mean)
  diff$CqNorm <- diff$rRNA - diff$Fn1
  diff$Date = date
  return(diff)
}

## Plot outcome vs diet-associated taxa
plot_outcome = function(taxa_var = "581a55014641e1cd55cdc272c0365a28", outcome_var = "Dose", min_clr = 0, 
                        cycles = c("C1D1", "Baseline", "C1D3", "C1D7"), pt_keep =NULL, method = "t.test",
                        clr= TRUE, summary = "mean"){
  t = load_16s(taxLevel = "ASVs", min_clr = min_clr, clr = clr)
  tf = t[taxa_var,]
  tf = data.frame(SampleID = names(tf), Abundance = as.numeric(tf))
  tf <- tf %>%
    mutate(SampleID = trimws(SampleID)) %>%  # remove leading spaces
    separate(SampleID, into = c("drop", "PatientID", "Cycle"), sep = "-") %>%
    dplyr::select(-drop)   # remove the "GO" part
  tf$PatientID = gsub("Pt", "", tf$PatientID)
  
  tf = load_go_tox(tf)
  
  tf = tf %>% filter(Cycle %in% cycles)
  tf$Outcome = tf[[outcome_var]]
  if (summary == "min"){
    tf = tf %>% group_by(PatientID) %>% 
      dplyr::summarise(Outcome = mean(Outcome), Abundance = min(Abundance))
  }
  if (summary == "mean"){
    tf = tf %>% group_by(PatientID) %>% 
      dplyr::summarise(Outcome = mean(Outcome), Abundance = mean(Abundance))
  }
  if (summary == "delta"){
    tf <- tf %>%
      filter(Cycle %in% cycles) %>%                       # keep only relevant cycles
      group_by(PatientID) %>%
      filter(n_distinct(Cycle) == 2) %>%                  # remove patients missing one of the cycles
      dplyr::summarise(Outcome = mean(Outcome),
                Abundance = Abundance[Cycle == cycles[2]] - Abundance[Cycle == cycles[1]])
  }
  
  if (!(is.null(pt_keep))){
    tf = tf %>% filter(PatientID %in% pt_keep)
  }
  
  tf$Outcome = case_when(tf$Outcome == 0 ~ "No", TRUE ~ "Yes")
  
  p = ggplot(tf, aes(x = Outcome, y = Abundance)) +
    geom_beeswarm(aes(color = Outcome)) + 
    geom_boxplot(alpha = 0.2, aes(fill = Outcome), outlier.size = -1) + 
    theme_pubr() +
    stat_compare_means(label.x = 1.25, method = method, 
                       aes(label = paste0("italic(p) ==", after_stat(p.format))), parse = TRUE,
                       label.y.npc = 0.97) + 
    scale_color_manual(values = c("grey", "darkred")) + 
    scale_fill_manual(values = c("grey", "darkred")) 
    
  return(p)
}


