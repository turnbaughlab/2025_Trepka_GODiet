# Load functions and directories
prefix <- "/mnt/tank/labmainshare/qb3share/" # "/Volumes/turnbaughlab/qb3share/"
dir <- paste0(prefix, "ktrepka/GODiet_v2/")
source(paste0(dir, "0_Functions.R"))

# Load diet data and save different categories
fn <- paste0(dir, "Tables/TotalsMergedSequencedFilteredAverage_WithMetadata.csv")
t <- read.csv(fn)

all_names <- colnames(t)[diet_index_start:diet_index_end]
foodgroup = c("F_TOTAL", "F_CITMLB", "F_OTHER", "F_JUICE", "V_TOTAL", "V_DRKGR", "V_REDOR_TOTAL", "V_REDOR_TOMATO", "V_REDOR_OTHER", "V_STARCHY_TOTAL", "V_STARCHY_POTATO", "V_STARCHY_OTHER", "V_OTHER", "V_LEGUMES", "G_TOTAL", "G_WHOLE", "G_REFINED", "PF_TOTAL", "PF_MPS_TOTAL", "PF_MEAT", "PF_CUREDMEAT", "PF_ORGAN", "PF_POULT", "PF_SEAFD_HI", "PF_SEAFD_LOW", "PF_EGGS", "PF_SOY", "PF_NUTSDS", "PF_LEGUMES", "D_TOTAL", "D_MILK", "D_YOGURT", "D_CHEESE", "OILS", "SOLID_FATS", "ADD_SUGARS", "A_DRINKS")
nutrient = setdiff(all_names, foodgroup)
hei = colnames(t)[hei_index_start:hei_index_end]

cat_df = data.frame("Diet" = c(foodgroup, nutrient, hei),
                    "Type" = c(rep("FoodGroup", length(foodgroup)), 
                               rep("Nutrient", length(nutrient)), 
                               rep("HEI", length(hei))))
fn = paste0(tabledir, "DietCategory.csv")
write.csv(cat_df, file = fn, row.names = FALSE)

# Normalize by KCAL (if relevant)
t_norm = normalize_all_to_kcal(t)
4*t_norm$CARB + 4*t_norm$PROT + 9*t_norm$TFAT # double check that normalization worked - this should equal 1

# Filter to patients with at least k samples
t_norm <- keep_k_pts(t_norm, k = 2)

# Modeling - for nutrients, normalize to KCAL. For foods, adjust for KCAL. For HEI, nothing
# 3 models:
# 1. Mixed effects model, pre vs post
# 2. Mixed effects model, time
# 3. Mixed effects model, every comparison vs baseline
dat <- t_norm %>% filter(Cycle != "EOT")
dat$PrePost = case_when(dat$Cycle == "C1D1" ~ 0,
                        dat$Cycle %in% c("C2D1", "C3D1", "EOT") ~ 1,
                        TRUE ~ NA)
dat$Linear = case_when(dat$Cycle == "C1D1" ~ 0,
                       dat$Cycle == "C2D1" ~ 1,
                       dat$Cycle == "C3D1" ~ 2,
                       dat$Cycle == "EOT" ~ 3,
                       TRUE ~ NA)
dat$Cat <- factor(dat$Cycle, levels = c("C1D1", "C2D1", "C3D1", "EOT"))

res_prepost <- fit_time_model(dat, timevar = "PrePost")
fn <- paste0(table_time, "PrePost.csv")
write.csv(res_prepost %>% arrange(FDR), fn)

res_linear <- fit_time_model(dat, timevar = "Linear")
fn <- paste0(table_time, "Linear.csv")
write.csv(res_linear %>% arrange(FDR), fn)

res_cat <- fit_time_model(dat, timevar = "Cat")
fn <- paste0(table_time, "Cat.csv")
write.csv(res_cat %>% arrange(FDR), fn)


