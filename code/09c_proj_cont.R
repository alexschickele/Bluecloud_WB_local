#' =============================================================================
#' @name proj_cont
#' @description computes spatial projections for the continuous data sub-pipeline
#' @param QUERY the query object from the master pipeline
#' @param MODELS the models object from the master pipeline
#' @param N_BOOTSTRAP number of bootstrap to do for the projections
#' @param PROJ_PATH (optional) path to a environmental raster, potentially 
#' different than the one given in the QUERY object. This is the case for 
#' supplementary projections in other time and climate scenarios for example. 
#' To your own risk and only for expert users !
#' @return an updated model list object containing the projections objects
#' embedded in each model sub-list.

proj_cont <- function(QUERY,
                      MODELS,
                      N_BOOTSTRAP,
                      PROJ_PATH){
  
  # --- 1. Load environmental data - TO FIX DYNAMICALLY
  features <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
    readAll() %>% 
    raster::subset(QUERY$CALL$ENV_VAR) %>% 
    rasterToPoints() %>% 
    as.data.frame() %>% 
    dplyr::select(-c(x, y))
  
  # --- 2. Define bootstraps
  # --- 2.1. Re-assemble all query tables
  tmp <- cbind(QUERY$Y, QUERY$X, QUERY$S)
  
  # --- 2.2. Run the bootstrap generation from tidy models
  boot_split <- bootstraps(tmp, times = N_BOOTSTRAP)
  
  # --- 3. Fit model on bootstrap
  # fit_resamples() does not save models by default. Thus the control_resamples()
  boot_fit <- MODELS$GAM$final_wf %>% 
    fit_resamples(resamples = boot_split,
                  control = control_resamples(extract = function (x) extract_fit_parsnip(x))) %>% 
    unnest(.extracts)
  
  # --- 4. Compute one prediction per bootstrap
  # As we extracted the model information in a supplementary column, we can
  # directly compute the bootstrap within the synthetic resample object.
  boot_proj <- boot_fit %>% 
    mutate(proj = map(.extracts, function(x)(x = predict(x, features))))
  
  # --- 5. Compute average and CV across bootstraps
  # First transform the object into a cell x bootstrap matrix
  # /!\ Need to create a unique row identifier for pivot_wider to work...
  tmp <- boot_proj %>% 
    dplyr::select(id, proj) %>% 
    unnest(c(id, proj)) %>% 
    as.data.frame() %>% 
    group_by(id) %>%
    mutate(row = row_number()) %>%
    pivot_wider(names_from = id, values_from = .pred) %>%
    dplyr::select(-row)
  
  # Open a raster to have the list of cells
  r_val <- raster(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
    getValues()
  
  # Assign the desired values to the non-NA cells in the list
  y_hat <- apply(tmp, 2, function(x){
    r <- r_val
    r[!is.na(r)] <- x
    x <- r
  })
  
  # --- 6. Append the MODEL object
  MODELS[[i]][["proj"]][["y_hat"]] <- y_hat
  
  
  return()
  
} # END FUNCTION