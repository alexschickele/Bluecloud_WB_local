#' ============================= MASTER SCRIPT =================================
#' For the Blue-Cloud 2026 project workbench
#' To be translated to D4-Science and notebook at the end
#' A. Schickele 2023
#' =============================================================================
#' ============================== TO DO LIST ===================================
#' - data access service
#' - make sure the models are sending dataminer jobs
#' - have an easy and common architecture of files
#' - check the compatibility between dataminer and Biomod2...
#' - have a prototype wrapper for the 3 sub-pipelines
#' =============================================================================

# --- START UP
# All will be called in the config file later
rm(list=ls())
setwd("/net/meso/work/aschickele/Diversity")
source(file = "./code/00_config.R")

# --- 1. Query biological data
# First check which species are available -- TOO LONG : add key on DB
list_bio <- list_bio(DATA_TYPE = "pres",
                     SAMPLE_SELECT = list(MIN_SAMPLE = 50, MIN_DEPTH = 0, MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016))

# Then query the species of interest -- TOO LONG : add key on DB
query <- query_bio(DATA_TYPE = "pres",
                   SP_SELECT = 102,
                   SAMPLE_SELECT = list(MIN_SAMPLE = 50, MIN_DEPTH = 0, MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016))

# --- 2. Query environmental data
query <- query_env(QUERY_BIO = query,
                   ENV_VAR = NULL,
                   ENV_PATH = "/net/meso/work/aschickele/Diversity/data/features_monthly")

# --- 3. Generate pseudo-absences if necessary
if(query$CALL$DATA_TYPE == "pres"){
  query <- pseudo_abs(QUERY = query,
                      METHOD_PA = "disk")
}

# --- 4. Further data processing
query <- query_check(QUERY = query,
                     OUTLIER = TRUE,
                     ENV_COR = 0.8,
                     MESS = TRUE)

# --- 5. Generate split and re sampling folds
query <- folds(QUERY = query,
               NFOLD = 5,
               FOLD_METHOD = "lon")

# --- 6. Hyper parameters definition











# --- END --- 