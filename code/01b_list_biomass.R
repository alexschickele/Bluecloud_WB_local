#' =============================================================================
#' @name list_biomass
#' @description extracts available Aphia_ID corresponding to user defined criteria
#' among the available data within the Atlanteco database, build locally.
#' @param DATA_SOURCE parameter passed from the wrapper function
#' @param SAMPLE_SELECT parameter passed from the wrapper function
#' @return a list of available Worms ID or Aphia ID and number of occurrences
#' within the data type and sample criteria

list_biomass <- function(DATA_SOURCE,
                           SAMPLE_SELECT){

  # --- 1. Connect to database
  # This database is only available locally at ETH for now. It contains the data
  # built by Fabio during the Atlanteco project.
  db <- dbConnect(RSQLite::SQLite(), paste0(project_wd,"/data/DB_clean.sqlite"))

  # --- 2. Query the database
  # Returns the Worms ID, or Aphia ID and number of occurrences within the data
  # type and sample criteria
  message("--- LIST BIO : retrieving species available in ATLANTECO")
  data_list <- tbl(db, paste0(DATA_SOURCE,"_data")) %>%
    dplyr::filter(depth >= !!SAMPLE_SELECT$TARGET_MIN_DEPTH &
                  depth <= !!SAMPLE_SELECT$TARGET_MAX_DEPTH &
                  year >= !!SAMPLE_SELECT$START_YEAR &
                  year <= !!SAMPLE_SELECT$STOP_YEAR &
                  measurementvalue != "Absence" &
                  measurementvalue != 0,
                  !is.na(measurementvalue),
                  worms_id != "Not found") %>%
    group_by(worms_id) %>%
    mutate(nb_occ = n()) %>%
    dplyr::select(-c(decimallongitude, decimallatitude, month, depth, year, measurementvalue, measurementunit)) %>%
    distinct() %>%
    dplyr::filter(nb_occ >= !!SAMPLE_SELECT$MIN_SAMPLE) %>%
    collect()

  # --- 3. Disconnect from database
  dbDisconnect(db)

  # --- 4. Wrap up and save
  return(data_list)

} # END FUNCTION
