#' =============================================================================
#' @name query_env
#' @description appends the query_bio output with a list of environmental
#' values at the sampling stations and a path to the environmental raster that
#' will be used for projections
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @return X: a data frame of environmental values at the sampling stations and 
#' @return ENV_PATH: the path to the environmental .nc or raster
#' @return Updates the output in a QUERY.RData and CALL.Rdata files

query_env <- function(FOLDER_NAME = NULL,
                      SUBFOLDER_NAME = NULL){
  
  # --- 1. Initialize function
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : query_env ********************"))
  # --- 1.2. Load the run metadata and query
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  
  # --- 2. Open features gridded data and names
  # /!\ To change later for a direct query in a .nc file : avoid 2 different raster files
  features <- stack(CALL$ENV_PATH) %>% readAll()
  features_name <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
    names()

  # --- 3. Re-grid sample on the raster resolution and filter
  # (1) The cell centers are at .5, thus it is re-gridded to the nearest .5 value
  # (2) /!\ Depth is not taken into account for now, neither year (1990-2016 = WOA)
  res <- res(features)[[1]]
  digit <- nchar(sub('^0+','',sub('\\.','',res)))-1
  sample <- QUERY$S %>% 
    cbind(QUERY$Y) %>% 
    mutate(decimallatitude = round(decimallatitude+0.5*res, digits = digit)-0.5*res) %>% 
    mutate(decimallongitude = round(decimallongitude+0.5*res, digits = digit)-0.5*res) 
  
  # --- 4. Select one sample per group of identical coordinates x month
  # Among each group of identical lat, long and month, one random point is selected
  S <- sample %>% 
    dplyr::select(-names(QUERY$Y)) %>% 
    group_by(decimallongitude, decimallatitude, month) %>% 
    slice_sample(n = 1) %>% 
    dplyr::ungroup() 
  
  # --- 5. Average measurement value per group of identical coordinates x month
  Y <- sample %>% 
    dplyr::select(decimallongitude, decimallatitude, month, names(QUERY$Y)) %>% 
    group_by(decimallongitude, decimallatitude, month) %>% 
    summarize_at(names(QUERY$Y), mean) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(names(QUERY$Y))
  
  # --- 6. Extract the environmental data in the data frame
  # If there is an NA, extract from nearest non-NA cells
  X <- NULL
  for(j in 1:nrow(S)){
    id <- grep(pattern = S$month[j], names(features))
    xy <- S[j,] %>% dplyr::select(x = decimallongitude, y = decimallatitude)
    
    tmp <- raster::extract(features[[id]], xy) %>% 
      as.data.frame()
    
    # if(is.na(sum(tmp))){
    #   tmp <- mclapply(features[[id]]@layers, function(a_layer) sample_raster_NA(a_layer, xy), mc.cores = 1) %>% 
    #     lapply(mean) %>% # average between two points if a coordinate is EXACTLY at an integer value (i.e., between two cells)
    #     as.data.frame()
    # } # If extract is NA
    if(is.na(sum(tmp))){
      tmp <- lapply(features[[id]]@layers, function(a_layer) sample_raster_NA(a_layer, xy)) %>% 
        lapply(mean) %>% 
        as.data.frame()
    }

    colnames(tmp) <- features_name
    X <- rbind(X, tmp)
  } # End for j
  
  # --- 7. Wrap up and save
  # --- 7.1. Append QUERY with the environmental values and save
  # And updated Y and S tables with dupplicate coordinate removed
  QUERY[["Y"]] <- Y
  QUERY[["S"]] <- S
  QUERY[["X"]] <- X
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  
  # --- 7.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  
} # END FUNCTION