#' =============================================================================
#' @name query_env
#' @description appends the query_bio output with a list of environmental
#' values at the sampling stations and a path to the environmental raster that
#' will be used for projections
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @return X: a data frame of environmental values at the sampling stations
#' @return Y and S updated with duplicate stations removed
#' @return Updates the output in a QUERY.RData file
#' @return an updated list of subfolders according to the minimum number of occurrence criteria

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
  features <- stack(CALL$ENV_PATH) %>% readAll()
  features_name <- features %>% names()

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
  # Among each group of identical lat and long, concatenates description
  # Updates the ID as the previous one is overwritten (we do not keep the raw data)
  S <- sample %>%
    dplyr::select(-names(QUERY$Y)) %>%
    group_by(decimallongitude, decimallatitude) %>%
    reframe(across(everything(), ~ str_flatten(unique(.x), collapse = ";"))) %>% 
    mutate(ID = row_number())
  
  # --- 5. Average measurement value per group of identical coordinates x month
  # The corresponding biological value is averaged across all samples of identical lat and long
  Y0 <- sample %>% 
    dplyr::select(decimallongitude, decimallatitude, names(QUERY$Y)) 
  
  Y <- NULL
  for(n in 1:nrow(S)){
    tmp <- Y0 %>% 
      inner_join(S[n,], by = c("decimallongitude", "decimallatitude")) %>% 
      dplyr::select(names(QUERY$Y))
    tmp <- apply(tmp, 2, mean)
    Y <- rbind(Y, tmp)
  }
  Y <- as.data.frame(Y)

  # --- 6. Extract the environmental data in the data frame
  # If there is an NA, extract from nearest non-NA cells
  X <- NULL
  for(j in 1:nrow(S)){
    xy <- S[j,] %>% dplyr::select(x = decimallongitude, y = decimallatitude)

    tmp <- raster::extract(features, xy) %>%
      as.data.frame()

    if(is.na(sum(tmp))){
      r_dist <- distanceFromPoints(features, xy) # Compute distance to NA point
      r_dist <- synchroniseNA(stack(r_dist, features[[1]]))[[1]] # Synchronize NA
      min_dist <- which.min(getValues(r_dist)) # Get closest non-NA point ID
      tmp <- features[min_dist] %>%
        as.data.frame()
    }

    colnames(tmp) <- features_name
    X <- rbind(X, tmp)
  } # End for j
  
  # --- 7. Wrap up and save
  # --- 7.1. Append QUERY with the environmental values and save
  # And updated Y and S tables with duplicate coordinate removed
  QUERY[["Y"]] <- Y
  QUERY[["S"]] <- S
  QUERY[["X"]] <- X
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  
  # --- 7.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  
  # --- 7.3. Update list of SUBFOLDER_NAME
  if(nrow(S) >= CALL$SAMPLE_SELECT$MIN_SAMPLE){
    return(SUBFOLDER_NAME)
  } else {
    return(NA)
  }
  
} # END FUNCTION
