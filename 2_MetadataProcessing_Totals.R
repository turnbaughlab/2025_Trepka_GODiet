############### Libraries and directories ############### 
# Load functions and directories
prefix <- "/mnt/tank/labmainshare/qb3share/" # "/Volumes/turnbaughlab/qb3share/"
dir <- paste0(prefix, "ktrepka/GODiet_v2/")
source(paste0(dir, "0_Functions.R"))

# Global constants
enriched_palette <- c(
  "Enriched" = "blue", 
  "Depleted" = "orange", 
  "Not Significant (FDR > 0.2)" = "gray",
  "FDR > 0.2" = "gray",
  "Not Significant" = "gray"
)

# Merge and update metadata
fn <- paste0(rawdir, "Metadata/RegimenIntensity.csv")
reg <- read.csv(fn) %>% dplyr::rename(Patient_ID = PatientID)
reg$Patient_ID = as.character(reg$Patient_ID)
fn <- paste0(rawdir, "Metadata/Metadata170523.csv")
meta_old = read.csv(fn) %>% dplyr::select(-X.1)
meta_new = left_join(meta_old, reg, by = "Patient_ID")
fn <- paste0(rawdir, "Metadata/Metadata250319.csv")
write.csv(meta_new, file = fn)

############### Map GO to ASA24 ############### 
# Exploration of available data / column names differences 
d_2016 = read_file(folder = "ASA24_2016", dir = rawdir)
d_2018 = read_file(folder = "ASA24_2018/GODIET_2021-03-04_52418", dir = rawdir)
d_2020 = read_file(folder = "ASA24_2020", dir = rawdir)
diff1 = setdiff(colnames(d_2018), colnames(d_2016)) # just need to fix capitalization for this one
diff2 = setdiff(colnames(d_2020), colnames(d_2018))

# Read in mapping
fn <- paste0(rawdir, "REDCapData/ASA_GO_Map.xlsx")
map <- read_excel(fn, sheet = 1) %>% as.data.frame()
colnames(map) <- c("UserName", "PatientID")

# Make sure all diet records have a mapping
compiled_data <- NULL
folders = c("ASA24_2016", "ASA24_2018/GODIET_2021-03-04_52418", "ASA24_2018/GODIET_2021-08-05_57721", "ASA24_2020")
for (folder in folders){
  fold <- paste0(rawdir, "ASA24DataExports/", folder, "/")
  files <- list.files(fold)
  fn_totals <- files[grepl("Totals", files)]
  totals <- read.csv(paste0(fold, fn_totals))
  if ("AMTUsual" %in% colnames(totals)){
    totals <- totals %>% dplyr::rename(AmtUsual = AMTUsual)
  }
  totals$Version = folder
  if (length(compiled_data) == 0){
    compiled_data <- totals
  } else {
    colnames = intersect(colnames(totals), colnames(compiled_data))
    compiled_data <- compiled_data[,colnames]
    totals <- totals[,colnames]
    compiled_data <- rbind(compiled_data, totals)
  }
}
compiled_data$AmtUsual <- case_when(compiled_data$AmtUsual == 3 ~ 1,
                                    compiled_data$AmtUsual == 2 ~ 0,
                                    compiled_data$AmtUsual == 0 ~ -1,
                                    TRUE ~ NA)
compiled_data <- left_join(compiled_data, map, by = "UserName")
compiled_data <- compiled_data %>% filter(!is.na(PatientID))
fn <- paste0(tabledir, "TotalsMerged.csv")
write.csv(file = fn, compiled_data, row.names = FALSE)

# Add cycle data manually (based on DietRecordCompletion)
re_annotate = FALSE
if (re_annotate){
  to_annotate <- compiled_data %>% dplyr::select(UserName, PatientID, RecordNo, RecordDayNo, ReportingDate)
  write.csv(file = paste0(rawdir, "REDCapData/TimeMap_ToAnnotate.csv"), to_annotate, row.names = FALSE)
}

# After annotating manually, load the annotated data and merge on it
m <- load_seq_metadata()

fn <- paste0(rawdir, "REDCapData/TimeMap_Annotated.csv")
time <- read.csv(fn) %>% as.data.frame() %>% dplyr::select(-CycleConfidence) %>% unique()
time$PtCycle = paste0(time$PatientID, time$Cycle); m$PtCycle <- paste0(m$PatientID, m$Cycle)
pt_cycle_keep = m %>% dplyr::select(PtCycle) %>% pull() %>% unique()
time <- time %>% filter(PtCycle %in% pt_cycle_keep)
time <- time %>% dplyr::select(PatientID, RecordNo, RecordDayNo, Cycle, ReportingDate)
time$ReportingDate = format(as.Date(time$ReportingDate, "%m/%d/%y"), "%m/%d/%Y")

fn <- paste0(tabledir, "TotalsMerged.csv")
t <- read.csv(fn) %>% as.data.frame()
t <- left_join(t, time, by = c("PatientID", "RecordNo", "RecordDayNo", "ReportingDate"), relationship = "many-to-many") # Add cycle info
t <- t %>% filter(!is.na(Cycle)) # Remove duplicates
t <- t %>% distinct(PatientID, Cycle, RecordDayNo, .keep_all = TRUE)
t$KCAL = 4*t$PROT + 9*t$TFAT + 4*t$CARB # Fix KCAL - from macronutrients (back-calculated)
t <- t %>% filter(Cycle != "Baseline") # remove spurious baseline sample
fn <- paste0(tabledir, "TotalsMergedSequenced.csv")
write.csv(t, file = fn, row.names = FALSE)

############### Filter based on calorie percentiles ############### 
fn <- paste0(tabledir, "TotalsMergedSequenced.csv")
id_keep <- c("2_C2D1_RecordDay1", "6_C3D1_RecordDay3", "15_C2D1_RecordDay1", "54_C1D1_RecordDay3") # list of patients with KCAL <= 458 to manually keep
t <- read.csv(fn) %>% as.data.frame()
t$ID <- paste0("Pt" = t$PatientID, 
               "_", t$Cycle, 
               "_RecordDay", t$RecordDayNo)
t_filt <- t %>%
  filter(KCAL > 458 | ID %in% id_keep) %>% 
  filter(KCAL < 5105)
t <- t %>% dplyr::select(-ID)
fn <- paste0(tabledir, "TotalsMergedSequencedFiltered.csv")
write.csv(t_filt, file = fn, row.names = FALSE)

############### Add metadata ############### 
# Every day is separate
fn <- paste0(tabledir, "TotalsMergedSequencedFiltered.csv")
d <- read.csv(fn) %>% as.data.frame()

fn <- fn_meta
seq <- read.csv(fn) %>% dplyr::rename(PatientID = Patient_ID, Cycle = Treatment_Cycle)
seq$Cycle[seq$Cycle == "Final"] <- "EOT"; seq$PatientID[seq$PatientID == "9b"] <- "9"; seq$PatientID <- as.numeric(seq$PatientID)
d <- inner_join(d, seq, by = c("PatientID", "Cycle"))
fn <- paste0(dir, "Tables/TotalsMergedSequencedFiltered_WithMetadata.csv")
write.csv(d, file = fn, row.names = FALSE)

# Averaged over 3 day period
fn <- paste0(tabledir, "TotalsMergedSequencedFiltered.csv")
d <- read.csv(fn) 
d <- d[,12:length(colnames(d))]
d <- d %>% dplyr::select(-DataComp) %>% dplyr::select(-Version)
d_surv <- d %>%
  group_by(PatientID, Cycle) %>% dplyr::summarise(NumberValidSurveys = n())
d_sum <- d %>%
  group_by(PatientID, Cycle) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = 'drop') %>%
  as.data.frame()

fn <- fn_meta
seq <- read.csv(fn) %>% dplyr::rename(PatientID = Patient_ID, Cycle = Treatment_Cycle)
seq$Cycle[seq$Cycle == "Final"] <- "EOT"; seq$PatientID[seq$PatientID == "9b"] <- "9"; seq$PatientID <- as.numeric(seq$PatientID)
d_sum <- inner_join(d_sum, seq, by = c("PatientID", "Cycle"))
d_sum <- left_join(d_sum, d_surv, by = c("PatientID", "Cycle"))

fn <- paste0(tabledir, "TotalsMergedSequencedFilteredAverage_WithMetadata.csv")
write.csv(d_sum, file = fn, row.names = FALSE)

fn <- paste0(dir, "Tables/TotalsMergedSequencedFilteredAverage_WithMetadata.csv")
t <- read.csv(fn)

############### Add healthy eating index (HEI) and alternative healthy eating index (AHEI) ############### 
# Healthy eating index (HEI)
fn <- paste0(tabledir, "TotalsMergedSequencedFilteredAverage_WithMetadata.csv")
add_hei(fn)

fn <- paste0(tabledir, "TotalsMergedSequencedFiltered_WithMetadata.csv")
add_hei(fn)

# Alternative healthy eating index (AHEI)
fn <- paste0(tabledir, "TotalsMergedSequencedFilteredAverage_WithMetadata.csv")
add_ahei(fn)

fn <- paste0(tabledir, "TotalsMergedSequencedFiltered_WithMetadata.csv")
add_ahei(fn)

############### Metadata summaries for Erin ############### 
fn <- paste0(tabledir, "TotalsMergedSequencedFilteredAverage_WithMetadata.csv")
metadata <- read.csv(fn) %>% as.data.frame() %>% filter(Cycle == "C1D1")
diet_cats = colnames(metadata)[diet_index_start:diet_index_end]

summary = NULL
for (diet_cat in diet_cats){
  diet <- metadata[[diet_cat]]
  mu = mean(diet)
  stdev = sd(diet)
  qs <- as.numeric(quantile(diet))
  row = data.frame("Variable" = diet_cat, "Mean" = mu, "StDev" = stdev, "Percentile0" = qs[1],
                   "Percentile25" = qs[2], "Percentile50" = qs[3], "Percentile75" = qs[4], "Percentile100" = qs[5])
  if(length(summary) == 0){
    summary = row
  } else {
    summary = rbind(summary, row)
  }
}

fn <- paste0(tabledir, "DietVariableStats.csv")
write.csv(summary, file = fn, row.names = FALSE, quote = FALSE)