### ---------------------------------------------- RT predictor benchmark for IP  ----------------------------------------------
# description:  benchmark retention time (RT) predictors for small scale immunopeptidome data
#               
# input:        1. MHC-I peptide dataset
# output:       
#               - Performance plots
#               
# author:       YH, JL

library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(vroom)

wd = "/home/yhorokh/data/Results_reports/wd/RT_benchmarks/"

# Benchmark parameters
train_split = 500
n_folds = 10
peptide_data <- read_xlsx("data/Table_S1_SPI-ART_rev-v31.xlsx", sheet = 2)

# Create folders
{
  suppressWarnings(dir.create("results/train"))
  suppressWarnings(dir.create("results/train/achrom"))
  suppressWarnings(dir.create("results/train/AutoRT"))
  suppressWarnings(dir.create("results/train/DeepLC"))
  
  suppressWarnings(dir.create("results/test"))
  suppressWarnings(dir.create("results/test/achrom"))
  suppressWarnings(dir.create("results/test/AutoRT"))
  suppressWarnings(dir.create("results/test/DeepLC"))
  
  suppressWarnings(dir.create("results/plots"))
}
# ---------------------------------------- (1) Pre-processing ------------------------------------------------
train <- peptide_data %>%
  separate(orig_raw_file, sep = "_", into = c("run", "person", "date", "SPL", "replicate")) %>%
  select(orig_raw_scan, run, date, SPL, Mascot_seq, RT) %>%
  mutate(SPL = str_split_fixed(SPL, "-", 2)[,1]) %>%
  group_by(date, SPL) %>%
  unique() %>%
  split(f = ~ date + SPL, drop = T)

for (i in seq_along(train)) {
  
  # Split the data into train and test sets
  # suppressWarnings(dir.create(paste0("results/predict/", names(train)[[i]])))
  # suppressWarnings(dir.create(paste0("results/train/", names(train)[[i]])))
  
  for (j in 1:n_folds) {
    train_i = train[[i]] %>%
      ungroup() %>%
      select(Mascot_seq, RT) %>%
      unique() %>%
      slice_sample(n = train_split)
    
    test_i = train[[i]] %>%
      filter(!Mascot_seq %in% train_i) %>%
      ungroup() %>%
      select(Mascot_seq, RT) %>%
      unique()
    
    ### Achrom
    vroom_write(train_i, file = paste0("results/train/achrom/", names(train)[[i]], ".sample", j, ".csv"), num_threads = 8, append = F)
    vroom_write(test_i, file = paste0("results/test/achrom/", names(train)[[i]], ".sample", j, ".csv"), num_threads = 8, append = F)
    
    ### DeepLC
    train_i %>%
      mutate(modifications = "") %>%
      rename(seq=Mascot_seq, tr=RT) %>%
      select(seq, modifications, tr) %>%
      vroom_write(file = paste0("results/train/DeepLC/", names(train)[[i]], ".sample", j, ".csv"), num_threads = 8, append = F)
    
    test_i %>%
      mutate(modifications = "") %>%
      rename(seq=Mascot_seq, tr=RT) %>%
      select(seq, modifications, tr) %>%
      vroom_write(file = paste0("results/test/DeepLC/", names(train)[[i]], ".sample", j, ".csv"), num_threads = 8, append = F)
    
    ### AutoRT
    train_i %>%
      rename(x=Mascot_seq, y=RT) %>%
      select(x, y) %>%
      vroom_write(file = paste0("results/train/AutoRT/", names(train)[[i]], ".sample", j, ".csv"), num_threads = 8, append = F)
    
    test_i %>%
      rename(x=Mascot_seq, y=RT) %>%
      select(x, y) %>%
      vroom_write(file = paste0("results/test/AutoRT/", names(train)[[i]], ".sample", j, ".csv"), num_threads = 8, append = F)
  }
}

# ---------------------------------------- (2) Calibrate and predict ------------------------------------------------
### pytheomics Achrom
library(reticulate)
# ### Create a new environment 
# conda_create("r-reticulate")
# use_condaenv("r-reticulate")
# conda_install(envname = "r-reticulate", 
#               packages = c("r-reticulate", "pyteomics", "numpy", "matplotlib", "pandas"), 
#               forge = T, 
#               channel = c("conda-forge","bioconda"))

use_condaenv("r-reticulate")
pyteomics <- import("pyteomics")




# DEV: ON --------------------------------

pep <- r_to_py(train_i$Mascot_seq)
RTs <- r_to_py(train_i$RT)


out <- list()
for (i in seq_along(pep)) {
  out <- pyteomics$achrom$get_RCs(sequences = pep, RTs = RTs)
}





# >>> from pyteomics import achrom
# >>> RCs = achrom.get_RCs(sequences, RTs)
# >>> achrom.calculate_RT('PEPTIDE', RCs)




### DeepLC 
# deeplc --file_pred  <path/to/peptide_file.csv> --file_cal <path/to/peptide_file_with_tr.csv>
# seq,modifications,tr



# ---------------------------------------- (3) Performance metrics ------------------------------------------------

