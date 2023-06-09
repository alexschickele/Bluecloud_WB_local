#' =============================================================================
#' @name model_wrapper
#' @description wrapper function redirecting towards the sub-pipeline
#' corresponding to the type of data. 
#' @description In case of proportion data, an input converting section is run, 
#' to properly pass the inputs to python library MBTR
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @param MODEL_LIST vector of string model names. It can be different from the
#' list passed to the hyperparameter function in the previous step
#' @return a model list object containing the different model objects
#' @return in case of proportion data, the model list object contains the path
#' to the model files as it cannot be passed as an object in memory
#' @return outputs are saved in a MODEL.RData object

# TO DO : implement the input converter for MBTR

model_wrapper <- function(FOLDER_NAME = NULL,
                          SUBFOLDER_NAME = NULL,
                          MODEL_LIST = NULL){
  
  # --- 1. Initialize function
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : model_wrapper ********************"))
  
  # --- 1.2. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  HP <- CALL$HP
  if(is.null(MODEL_LIST)){
    MODEL_LIST <- HP$CALL$MODEL_LIST
  }
  
  # --- 2. Redirection to PRESENCE model
  if(CALL$DATA_TYPE == "pres"){
    # --- 2.1. Load function
    source(file = paste0(project_wd, "/code/07b_model_pres.R"))
    
    # --- 2.2. Run function
    MODEL <- model_pres(CALL,
                        QUERY = QUERY,
                        HP = HP,
                        MODEL_LIST = MODEL_LIST)
  } # END if pres
  
  # --- 3. Redirection to CONTINUOUS model
  if(CALL$DATA_TYPE == "cont"){
    # --- 3.1. Load function
    source(file = paste0(project_wd, "/code/07c_model_cont.R"))
    
    # --- 3.2. Run function
    MODEL <- model_cont(CALL,
                        QUERY = QUERY,
                        HP = HP,
                        MODEL_LIST = MODEL_LIST)
  } # END if pres
  
  
  # --- 4. Redirection to PROPORTION model
  # TO BE IMPLEMENTED
  
  # --- 5. Wrap up and save
  # --- 5.1. Save file(s)
  save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  # --- 5.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)

} # END FUNCTION
