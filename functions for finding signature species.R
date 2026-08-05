species_by_hub_matrix <- function(survey_lookup, survey_species_matrix) {
  
  library(dplyr)
  library(tidyr)
  
  survey_species_matrix %>%
    pivot_longer(
      cols = -survey_id,
      names_to = "species",
      values_to = "present"
    ) %>%
    filter(present == 1) %>%
    left_join(
      survey_lookup %>%
        select(survey_id, hub),
      by = "survey_id"
    ) %>%
    count(species, hub, name = "records") %>%
    pivot_wider(
      names_from = hub,
      values_from = records,
      values_fill = 0
    ) %>%
    arrange(species)
}

library(dplyr)
library(tibble)
library(vegan)

calculate_hub_species_ses <- function(
    survey_lookup,
    survey_species_matrix,
    nsim = 999,
    ses_threshold = 2,
    null_method = "curveball",
    burnin = 50000,
    thin = 500,
    seed = NULL
) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Separate survey IDs from the presence/absence matrix
  survey_ids <- survey_species_matrix$survey_id
  
  species_matrix <- survey_species_matrix %>%
    select(-survey_id) %>%
    as.matrix()
  
  # Overall number of surveys in which each species occurs
  overall_observations <- colSums(species_matrix)
  
  storage.mode(species_matrix) <- "numeric"
  
  # Check that the species matrix is binary
  if (!all(species_matrix %in% c(0, 1))) {
    stop(
      "survey_species_matrix must contain only 0s and 1s ",
      "apart from the survey_id column."
    )
  }
  
  # Match each matrix row to its hub
  hub_lookup <- survey_lookup$hub[
    match(survey_ids, survey_lookup$survey_id)
  ]
  
  if (anyNA(hub_lookup)) {
    
    missing_ids <- survey_ids[is.na(hub_lookup)]
    
    stop(
      "The following survey IDs have no matching hub: ",
      paste(missing_ids, collapse = ", ")
    )
  }
  
  # Fix hub order
  hub_factor <- factor(
    hub_lookup,
    levels = unique(hub_lookup)
  )
  
  hub_names <- levels(hub_factor)
  
  # Design matrix: surveys x hubs
  hub_design <- model.matrix(
    ~ hub_factor - 1
  )
  
  colnames(hub_design) <- hub_names
  
  # Observed hub-by-species counts
  observed_hub_species <- crossprod(
    hub_design,
    species_matrix
  )
  
  # Convert to species x hubs
  observed_species_hub <- t(observed_hub_species)
  

  # Null model preserving row and column totals
  
  message("Creating null model...")
  
  null_model <- vegan::nullmodel(
    species_matrix,
    method = null_method
  )
  
  if (null_method == "curveball") {
    message("Starting burn-in...")
    
    null_model <- update(
      null_model,
      nsim = burnin
    )
    
    message("Burn-in complete.")
  }
  
  message("Null model created.")
  
  n_species <- ncol(species_matrix)
  n_hubs <- length(hub_names)
  
  # Store hub-level counts only:
  # species x hubs x simulations
  null_counts <- array(
    NA_real_,
    dim = c(n_species, n_hubs, nsim),
    dimnames = list(
      species = colnames(species_matrix),
      hub = hub_names,
      simulation = paste0("sim_", seq_len(nsim))
    )
  )
  
  for (i in seq_len(nsim)) {
    
    if (i <= 5) {
      message("Starting simulation ", i)
    }
    
    if (null_method == "curveball") {
      
      random_matrix <- simulate(
        null_model,
        nsim = 1,
        thin = thin
      )[, , 1]
      
    } else {
      
      random_matrix <- simulate(
        null_model,
        nsim = 1
      )[, , 1]
    }
    
    random_hub_species <- crossprod(
      hub_design,
      random_matrix
    )
    
    null_counts[, , i] <- t(random_hub_species)
    
    if (i %% 50 == 0) {
      message("Completed ", i, " of ", nsim, " simulations")
    }
   
  }
  
  message("Calculating expected values and SES...")
  
  expected_species_hub <- rowMeans(
    null_counts,
    dims = 2
  )
  
  message("Calculating null standard deviations...")
  
  null_sd_species_hub <- apply(
    null_counts,
    c(1, 2),
    sd
  )
  
  message("Calculating SES...")
  
  # Standardised effect size
  ses_species_hub <- (
    observed_species_hub - expected_species_hub
  ) / null_sd_species_hub
  
  ses_species_hub[null_sd_species_hub == 0] <- NA_real_  
  # SES output matrix
  ses_matrix <- ses_species_hub %>%
    as.data.frame(check.names = FALSE) %>%
    rownames_to_column("species") %>%
    as_tibble()
  
  all_results <- bind_rows(
    lapply(seq_along(hub_names), function(h) {
      
      tibble(
        species = rownames(observed_species_hub),
        hub = hub_names[h],
        observations = observed_species_hub[, h],
        overall_observations = overall_observations[rownames(observed_species_hub)],
        expected_observations = expected_species_hub[, h],
        null_sd = null_sd_species_hub[, h],
        ses = ses_species_hub[, h]
      )
    })
  ) %>%
    arrange(hub, desc(ses), species)
  
  
  hub_tables <- all_results %>%
    filter(
      !is.na(ses),
      ses > ses_threshold
    ) %>%
    select(
      hub,
      species,
      observations,
      overall_observations,
      expected_observations,
      ses
    ) %>%
    mutate(
      expected_observations = round(expected_observations, 2),
      ses = round(ses, 2)
    ) %>%
    group_split(hub) %>%
    setNames(
      all_results %>%
        filter(
          !is.na(ses),
          ses > ses_threshold
        ) %>%
        distinct(hub) %>%
        pull(hub)
    ) %>%
    lapply(function(x) {
      x %>%
        select(-hub) %>%
        arrange(desc(ses), desc(observations), species)
    })  
  return(
    list(
      ses_matrix = ses_matrix,
      hub_tables = hub_tables,
      all_results = all_results,
      
      # Diagnostic outputs
      observed_counts = observed_species_hub,
      expected_counts = expected_species_hub,
      null_sd = null_sd_species_hub,
      null_counts = null_counts
    )
  )
}
