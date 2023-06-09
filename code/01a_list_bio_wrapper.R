#' =============================================================================
#' @name list_bio_wrapper
#' @description wrapper around functions to list the available species within the
#' different data types and access service sources.
#' @param FOLDER_NAME name of the folder to create, corresponding to the run
#' @param DATA_SOURCE type of data to query from : string among "cont", "pres" or
#' "prop" corresponding to continuous (e.g. Ecotaxa), presence-only (e.g. OBIS) 
#' or proportions (e.g. omics) data
#' @param SAMPLE_SELECT list of sample selection criteria, including :
#' MIN_SAMPLE : minimum number of geographical points to consider the selection (i.e. non-zero records)
#' MIN_DEPTH and MAX_DEPTH: minimum and maximum depth levels in m
#' START_YEAR and STOP_YEAR: start and stop years for the considered period
#' @return for "cont" and "pres" data, returns a data frame with the number
#' of occurrences per worms_ID available
#' @return for "prop" data, returns a complete list of metadata, samples, taxonomic
#' annotations available.
#' @return the returned object is saved in the run file to avoid re-running the query

list_bio_wrapper <- function(FOLDER_NAME = "test_run",
                             DATA_SOURCE = "cont",
                             SAMPLE_SELECT = list(MIN_SAMPLE = 50, MIN_DEPTH = 0, MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016)){
  
  # --- 1. Initialize
  # --- 1.1. Parameter checking
  if(DATA_SOURCE != "cont" & DATA_SOURCE != "pres" & DATA_SOURCE != "omic"){
    stop("The specified data source should be 'cont', 'pres' or 'omic'")
  }
  
  # --- 1.2. Folder creation
  # If the directory exist, stops the creation and inform the user to manually delete the
  # previous runs, or name the directory differently
  folderpath <- paste0(project_wd,"/output/",FOLDER_NAME)
  if(file.exists(folderpath)==TRUE){
    stop("--- This foldername is already used")
  } else {
    dir.create(folderpath)
  }
  
  # --- 2. Redirection to ATLANTECO data access
  # For continuous and presence source data
  if(DATA_SOURCE == "cont" | DATA_SOURCE == "pres"){
    # --- 2.1. Load function
    source(file = paste0(project_wd, "/code/01c_list_atlanteco.R"))
    
    # --- 2.2. Run function
    LIST_BIO <- list_atlanteco(DATA_SOURCE = DATA_SOURCE,
                               SAMPLE_SELECT = SAMPLE_SELECT)
  } # End ATLANTECO redirection
  
  # --- 3. Redirection to MGNIFY data access
  # For continuous and presence source data
  if(DATA_SOURCE == "omic"){
    # --- 3.1. Load function
    source(file = paste0(project_wd, "/code/01b_list_mgnify.R"))
    
    # --- 3.2. Run function
    LIST_BIO <- list_mgnify(SAMPLE_SELECT = SAMPLE_SELECT)
  } # End ATLANTECO redirection
  
  # --- 4. Wrap up and save
  CALL <- list(DATA_SOURCE = DATA_SOURCE, 
               SAMPLE_SELECT = SAMPLE_SELECT, 
               LIST_BIO = LIST_BIO)
  save(CALL, file = paste0(project_wd,"/output/", FOLDER_NAME, "/CALL.RData"))
  return(CALL$LIST_BIO)
  
} # END FUNCTION
