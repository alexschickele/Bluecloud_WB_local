#' =============================================================================
#' @name proj_binary
#' @description computes spatial projections for the binary data sub-pipeline
#' @param CALL the call object from the master pipeline
#' @param QUERY the query object from the master pipeline
#' @param MODEL the models object from the master pipeline
#' @return an updated model list object containing the projections objects
#' embedded in each model sub-list.

proj_binary <- function(QUERY,
                        MODEL,
                        CALL){


  # --- 1. Initialize function
  # --- 1.1. Open base raster and values
  r0 <- CALL$ENV_DATA[[1]][[1]]
  r_val <- getValues(r0)

  # --- 1.2. Define the projections to compute
  # All algorithms if FAST == FALSE; only the ones that passed QC otherwise
  if(CALL$FAST == FALSE){
    loop_over <- CALL$HP$MODEL_LIST
  } else {
    loop_over <- MODEL$MODEL_LIST
  }

  # --- 2. Define bootstraps
  # --- 2.1. Re-assemble all query tables
  tmp <- cbind(QUERY$Y, QUERY$X)#, QUERY$S)

  # --- 2.2. Run the bootstrap generation from tidy models
  boot_split <- bootstraps(tmp, times = CALL$N_BOOTSTRAP, strata = "measurementvalue")

  # --- 3. Start the loop over algorithms
  for(i in loop_over){
    # --- 3.1. Fit model on bootstrap
    # fit_resamples() does not save models by default. Thus the control_resamples()
    boot_fit <- MODEL[[i]][["final_wf"]] %>%
      fit_resamples(resamples = boot_split,
                    control = control_resamples(extract = function (x) extract_fit_parsnip(x),
                    verbose = FALSE)) %>%
      unnest(.extracts)


    # --- 4. Loop over month for predictions
    y_hat <- NULL
    for(m in 1:length(CALL$ENV_DATA)){

      # --- 4.1. Load the right features
      features <- CALL$ENV_DATA[[m]] %>%
        raster::subset(QUERY$SUBFOLDER_INFO$ENV_VAR) %>%
        rasterToPoints() %>%
        as.data.frame() %>%
        dplyr::select(-c(x, y))

      # --- 4.2. Compute one prediction per bootstrap
      # As we extracted the model information in a supplementary column, we can
      # directly compute the bootstrap within the synthetic resample object.
      boot_proj <- boot_fit %>%
        mutate(proj = map(.extracts, function(x)(x = predict(x, features))))

      # --- 4.3. First transform the object into a cell x bootstrap matrix
      # /!\ Need to create a unique row identifier for pivot_wider to work...
      tmp <- boot_proj %>%
        dplyr::select(id, proj) %>%
        unnest(c(id, proj)) %>%
        as.data.frame() %>%
        group_by(id) %>%
        mutate(row = row_number()) %>%
        pivot_wider(names_from = id, values_from = .pred) %>%
        dplyr::select(-row)

      # --- 4.4. Assign the desired values to the non-NA cells in the base raster
      tmp <- apply(tmp, 2, function(x){
        r <- r_val
        r[!is.na(r)] <- x
        x <- r
      })

      # --- 4.5. Concatenate with previous month
      y_hat <- abind(y_hat, tmp, along = 3)
      #message(paste("--- PROJ : month", m, "done \t"))
    } # for m month

    # --- 5. Cut spatial discontinuities
    if(!is.null(CALL$CUT)){
      tmp <- apply(y_hat, c(2,3), function(x){
        # --- 5.1. Open presence data
        xy <- QUERY$S %>%
          dplyr::select(decimallongitude, decimallatitude)
        xy <- xy[which(QUERY$Y == 1),] # specific to presence data

        # --- 5.2. Cut y_hat
        x[x < CALL$CUT*max(x, na.rm = TRUE)] <- 0
        r <- setValues(r0, x)

        # --- 5.3. Define patches and overlap with presence points
        r_patch <- clump(r)
        id_patch <- r_patch %>%
          raster::extract(xy) %>% unique() %>%
          .[!is.na(.)]

        # --- 5.4. Subset values from a patch overlapping with presences
        r_patch[!(r_patch %in% id_patch)] <- 0
        r_patch <- getValues(r_patch)
        r_patch[r_patch > 0] <- 1
        x <- x*r_patch
        return(x)
      })
      y_hat <- tmp
    } # if CUT

    # --- 6. Compute the average CV across bootstrap runs as a QC
    if(dim(y_hat)[[2]] == CALL$N_BOOTSTRAP) {
      NSD <- apply(y_hat, c(1,3), function(x)(x = sd(x, na.rm = TRUE))) %>%
        mean(na.rm = TRUE)
      NSD <- NSD/mean(y_hat, na.rm = TRUE)
    } else {
      NSD <- NA
      message(paste("--- PROJ: Model", i, " discarded, bootstrap did not complete"))
    }
    # --- 7. Append the MODEL object
    # --- 7.1. Save the evaluation metric and projections
    MODEL[[i]][["proj"]][["y_hat"]] <- y_hat
    MODEL[[i]][["eval"]][["NSD"]] <- NSD

    # --- 7.2. Discard low quality models according to NSD
    # Fixed at 50% average across the projection
    if(MODEL[[i]][["eval"]][["NSD"]] > 0.5 | is.na(MODEL[[i]][["eval"]][["NSD"]]) | is.null(MODEL[[i]][["eval"]][["NSD"]]){
      MODEL$MODEL_LIST <- MODEL$MODEL_LIST[MODEL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to NSD =", MODEL[[i]][["eval"]][["NSD"]], "> 0.5 \n"))
    }

    # --- 7.3. Discard if the PRE_VIP QC in the query is too low
    # In case of FAST = TRUE; the species was kept until here but we do not want
    # to compute an ensemble for it.
    if(QUERY$eval$PRE_VIP < 0.05){
      MODEL$MODEL_LIST <- MODEL$MODEL_LIST[MODEL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to PRE_VIP =", QUERY$eval$PRE_VIP, "< 0.05 \n"))
    }
  } # for i model loop

  # --- 8. Build the ensemble NSD if there is an ensemble
  if(CALL$ENSEMBLE == TRUE & (length(MODEL$MODEL_LIST) > 1)){
    MODEL[["ENSEMBLE"]][["eval"]][["NSD"]] <- lapply(MODEL$MODEL_LIST,
                                                     FUN = function(x){
                                                       x <- MODEL[[x]]$eval$NSD
                                                     }) %>%
      unlist() %>% mean()
  } # End ENSEMBLE QC

  return(MODEL)

} # END FUNCTION


