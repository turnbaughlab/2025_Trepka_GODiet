############### Libraries and directories ############### 
# Load functions and directories
prefix <- "/mnt/tank/labmainshare/qb3share/" # "/Volumes/turnbaughlab/qb3share/"
dir <- paste0(prefix, "ktrepka/GODiet_v2/")
source(paste0(dir, "0_Functions.R"))

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
d_2016 = read_file(folder = "ASA24_2016", dir = rawdir, type = "Items")
d_2018 = read_file(folder = "ASA24_2018/GODIET_2021-03-04_52418", dir = rawdir, type = "Items")
d_2020 = read_file(folder = "ASA24_2020", dir = rawdir, type = "Items")
diff1 = setdiff(colnames(d_2018), colnames(d_2016)) 
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
  fn_totals <- files[grepl("Items", files)]
  totals <- read.csv(paste0(fold, fn_totals))

  totals$Version = folder
  if (is.null(compiled_data)){
    compiled_data <- totals
  } else {
    colnames = intersect(colnames(totals), colnames(compiled_data))
    compiled_data <- compiled_data[,colnames]
    totals <- totals[,colnames]
    compiled_data <- rbind(compiled_data, totals)
  }
}

compiled_data <- left_join(compiled_data, map, by = "UserName")
compiled_data <- compiled_data %>% filter(!is.na(PatientID))
fn <- paste0(tabledir, "ItemsMerged.csv")
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

fn <- paste0(tabledir, "ItemsMerged.csv")
t <- read.csv(fn) %>% as.data.frame()
t <- left_join(t, time, by = c("PatientID", "RecordNo", "RecordDayNo", "ReportingDate"), relationship = "many-to-many") # Add cycle info
t <- t %>% filter(!is.na(Cycle)) # Remove duplicates
#t <- t %>% distinct(PatientID, Cycle, RecordDayNo, .keep_all = TRUE) # doesn't work for items
t <- t %>%
  group_by(PatientID, Cycle, RecordDayNo, Food_Description, Version) %>%
  mutate(row_in_group = row_number()) %>%      # number rows within the group
  filter(row_in_group %% 2 == 1) %>%          # keep 1st, 3rd, 5th, etc.
  ungroup() %>%
  dplyr::select(-row_in_group)
t$KCAL = 4*t$PROT + 9*t$TFAT + 4*t$CARB # Fix KCAL - from macronutrients (back-calculated)
t <- t %>% filter(Cycle != "Baseline") # remove spurious baseline sample
fn <- paste0(tabledir, "ItemsMergedSequenced.csv")
write.csv(t, file = fn, row.names = FALSE)

################ Filter by KCAL cutoff ################ 
fn <- paste0(tabledir, "TotalsMergedSequenced.csv")
id_keep <- c("2_C2D1_RecordDay1", "6_C3D1_RecordDay3", "15_C2D1_RecordDay1", "54_C1D1_RecordDay3") # list of patients with KCAL <= 458 to manually keep
totals <- read.csv(fn) %>% as.data.frame()
totals$ID <- paste0("Pt" = totals$PatientID, 
               "_", totals$Cycle, 
               "_RecordDay", totals$RecordDayNo)
totals_filt <- totals %>%
  filter(KCAL > 458 | ID %in% id_keep) %>% 
  filter(KCAL < 5105)
totals_filt = totals_filt %>% dplyr::select(PatientID, Cycle, RecordNo, RecordDayNo, IntakeStartDateTime, Version)
tf = inner_join(t, totals_filt, by = c("PatientID", "Cycle", "RecordDayNo", "RecordNo", "IntakeStartDateTime", "Version"))
fn <- paste0(tabledir, "ItemsMergedSequencedFiltered.csv")
write.csv(tf, file = fn, row.names = FALSE)

############### Add metadata ############### 
# Every day is separate
fn <- paste0(tabledir, "ItemsMergedSequencedFiltered.csv")
d <- read.csv(fn) %>% as.data.frame()

fn <- fn_meta
seq <- read.csv(fn) %>% dplyr::rename(PatientID = Patient_ID, Cycle = Treatment_Cycle)
seq$Cycle[seq$Cycle == "Final"] <- "EOT"; seq$PatientID[seq$PatientID == "9b"] <- "9"; seq$PatientID <- as.numeric(seq$PatientID)
d <- inner_join(d, seq, by = c("PatientID", "Cycle"))
fn <- paste0(dir, "Tables/ItemsMergedSequencedFiltered_WithMetadata.csv")
write.csv(d, file = fn, row.names = FALSE)


## Sanity check - there shouldn't be duplicate entries
d %>% filter(PatientID == "31") %>% 
  filter(Cycle == "C1D1") %>% 
  filter(THEO > 0) %>% 
  arrange(-THEO) %>% 
  dplyr::select(RecordDayNo, Version, Food_Description, THEO)

