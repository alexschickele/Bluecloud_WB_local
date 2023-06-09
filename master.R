#' ============================= MASTER SCRIPT =================================
#' For the Blue-Cloud 2026 project workbench
#' To be translated to D4-Science and notebook at the end
#' A. Schickele 2023
#' =============================================================================
#' ============================== TO DO LIST ===================================
#' - check data access service as input
#' - check blue cloud data miner compatibility of the functions
#' - include MBTR in the prototypes for proportion data
#' =============================================================================

# --- 0. Start up and load functions
# All will be called in the config file later
rm(list=ls())
closeAllConnections()
setwd("/net/meso/work/aschickele/Bluecloud_WB_local")
source(file = "./code/00_config.R")
run_name <- "mbtr_test"

# --- 1a. List the available species
# Within the user defined selection criteria
list_bio <- list_bio_wrapper(FOLDER_NAME = run_name,
                             DATA_SOURCE = "cont",
                             SAMPLE_SELECT = list(MIN_SAMPLE = 50, MIN_DEPTH = 0, MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016))

# Define the list of species to consider
# sp_list <- c("5820", "9760") # random OTU short selection
sp_list <- list_bio %>% 
  # dplyr::filter(grepl("Calanus | Calanoides ", scientificname)) %>%
  dplyr::filter(grepl("Thalassiosira ", scientificname)) %>%
  # dplyr::filter(grepl("Chaetoceros ", scientificname)) %>%
  dplyr::select(worms_id) %>% 
  unique() %>% pull()

# --- 1b. Create the output folder, initialize parallelisation and parameters
# (1) Create an output folder containing all species-level runs, (2) Stores the 
# global parameters in an object, (3) Checks for environmental correlated variables
subfolder_list <- run_init(FOLDER_NAME = run_name,
                           SP_SELECT = sp_list,
                           LOAD_FROM = NULL,
                           DATA_TYPE = "cont",
                           ENV_VAR = NULL,
                           ENV_PATH = c("/net/meso/work/aschickele/Bluecloud_WB_local/data/bio_oracle", 
                                        "/net/meso/work/aschickele/Bluecloud_WB_local/data/features_mean_from_monthly"),
                           ENV_COR = 0.8,
                           NFOLD = 3,
                           FOLD_METHOD = "lon")

# --- 2a. Query biological data
# Get the biological data of the species we wish to model
mcmapply(FUN = query_bio_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 2b. Query environmental data
# This functions returns an updated subfolder_list object to avoid computing
# species with less than the user defined minimum occurrence number
subfolder_list <- mcmapply(FUN = query_env,
                  FOLDER_NAME = run_name,
                  SUBFOLDER_NAME = subfolder_list,
                  mc.cores = min(length(subfolder_list), MAX_CLUSTERS)) %>% 
  na.omit(subfolder_list) %>% 
  as.vector()

# --- 3. Generate pseudo-absences if necessary
mcmapply(FUN = pseudo_abs,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         METHOD_PA = "density",
         PER_RANDOM = 0.25,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 4. Outliers, Environmental predictor and MESS check 
mcmapply(FUN = query_check,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         OUTLIER = TRUE,
         UNIVARIATE = TRUE,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 5. Generate split and re sampling folds
mcmapply(FUN = folds,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 6. Hyper parameters to train
hyperparameter(FOLDER_NAME = run_name,
               MODEL_LIST = c("GLM","GAM","RF","MLP","SVM","BRT"),
               LEVELS = 3)

# --- 7. Model fit -- FIX : RF is very long for big data
mcmapply(FUN = model_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 8. Model evaluation
# Performance metric and variable importance
mcmapply(FUN = eval_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         ENSEMBLE = TRUE,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 9. Model projections
mcmapply(FUN = proj_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         N_BOOTSTRAP = 10,
         CUT = 0.1,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 10. Output plots
# --- 10.1. Standard maps per algorithms
mcmapply(FUN = standard_maps,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         ENSEMBLE = TRUE,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 10.2. Partial dependency plots - TAKES AGES FOR LARGE OCCURRENCE NUMBER
mcmapply(FUN = pdp,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         N_BOOTSTRAP = 10,
         ENSEMBLE = TRUE,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 10.3 Diversity
diversity_maps(FOLDER_NAME = run_name,
               SUBFOLDER_NAME = subfolder_list,
               BUFFER = 1,
               N_BOOTSTRAP = 10)
# --- END --- 