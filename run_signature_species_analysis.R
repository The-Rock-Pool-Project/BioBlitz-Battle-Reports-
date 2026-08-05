library(dplyr)
library(tibble)
library(vegan)
library(purrr)
library(tidyr)

inputs <- readRDS(
  "Laura/Biodiversity report Q2 2026/analysis_data/signature_species_inputs.rds"
)

survey_lookup <- inputs$survey_lookup
survey_species_matrix <- inputs$survey_species_matrix
taxon_data <- inputs$taxon_data

taxon_lookup <- taxon_data %>%
  transmute(
    species = taxon.name,
    rank = taxon.rank,
    taxon_id = taxon.id,
    common_name = taxon.preferred_common_name
  ) %>%
  filter(
    !is.na(species),
    species != ""
  ) %>%
  distinct(
    species,
    .keep_all = TRUE
  )

allowed_photo_licences <- c(
  "cc0",
  "cc-by",
  "cc-by-sa",
  "cc-by-nc",
  "cc-by-nc-sa"
)


normalise_licence <- function(x) {
  x <- tolower(trimws(x))
  x[x == ""] <- NA_character_
  x
}

# Convert an iNaturalist square-image URL to a larger version
make_medium_url <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    sub(
      "/square\\.(jpg|jpeg|png)$",
      "/medium.\\1",
      x,
      ignore.case = TRUE
    )
  )
}
source("Laura/Biodiversity report Q2 2026/functions for finding signature species.R")

greedyqswap_results <- calculate_hub_species_ses(
  survey_lookup = survey_lookup,
  survey_species_matrix = survey_species_matrix,
  nsim = 999,
  null_method = "greedyqswap",
  seed = 123
)

signature_species <- greedyqswap_results$all_results %>%
  left_join(
    taxon_lookup,
    by = "species"
  ) %>%
  filter(
    rank == "species",
    overall_observations >= 10,
    !is.na(ses),
    ses > 2
  ) %>%
  group_by(hub) %>%
  slice_max(
    order_by = ses,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    hub,
    species,
    common_name,
    taxon_id,
    observations,
    overall_observations,
    expected_observations,
    ses
  ) %>%
  mutate(
    expected_observations = round(expected_observations, 2),
    ses = round(ses, 2)
  ) %>%
  arrange(hub)

# -------------------------------------------------------------------
# EXTRACT PHOTOS FROM BRPC OBSERVATIONS
# -------------------------------------------------------------------

extract_observation_photos <- function(
    photos,
    species,
    hub
) {
  
  if (
    is.null(photos) ||
    !is.data.frame(photos) ||
    nrow(photos) == 0
  ) {
    return(tibble())
  }
  
  photos %>%
    as_tibble() %>%
    transmute(
      species = species,
      hub = hub,
      photo_id = id,
      image_url = make_medium_url(url),
      image_attribution = attribution,
      image_licence = normalise_licence(license_code),
      image_width = original_dimensions.width,
      image_height = original_dimensions.height,
      hidden = hidden
    )
}

brpc_species_photos <- pmap_dfr(
  list(
    photos = taxon_data$photos,
    species = taxon_data$taxon.name,
    hub = taxon_data$hub
  ),
  extract_observation_photos
) %>%
  filter(
    !is.na(species),
    species != "",
    !is.na(image_url),
    image_url != "",
    !is.na(image_licence),
    image_licence %in% allowed_photo_licences,
    is.na(hidden) | hidden == FALSE
  ) %>%
  distinct(
    photo_id,
    .keep_all = TRUE
  )


# -------------------------------------------------------------------
# CHOOSE THE BEST BRPC PHOTO
# -------------------------------------------------------------------

# Prefer landscape or near-landscape images for the wide cards.
# Within that preference, select the largest available image.
brpc_species_photos <- brpc_species_photos %>%
  mutate(
    aspect_ratio = image_width / image_height,
    landscape_preference = case_when(
      is.na(aspect_ratio) ~ 3L,
      aspect_ratio >= 1.2 ~ 1L,
      aspect_ratio >= 0.8 ~ 2L,
      TRUE ~ 3L
    ),
    image_area = image_width * image_height
  )

# First choice: same species and same hub
hub_species_photos <- brpc_species_photos %>%
  arrange(
    species,
    hub,
    landscape_preference,
    desc(image_area),
    photo_id
  ) %>%
  group_by(species, hub) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    species,
    hub,
    hub_image_url = image_url,
    hub_image_attribution = image_attribution,
    hub_image_licence = image_licence,
    hub_image_photo_id = photo_id
  )

# Second choice: same species from any BRPC hub
general_species_photos <- brpc_species_photos %>%
  arrange(
    species,
    landscape_preference,
    desc(image_area),
    photo_id
  ) %>%
  group_by(species) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    species,
    general_image_url = image_url,
    general_image_attribution = image_attribution,
    general_image_licence = image_licence,
    general_image_photo_id = photo_id,
    general_image_hub = hub
  )


# -------------------------------------------------------------------
# SUITABLY LICENCED DEFAULT TAXON PHOTOS
# -------------------------------------------------------------------

default_taxon_photos <- taxon_data %>%
  transmute(
    species = taxon.name,
    default_image_url = taxon.default_photo.medium_url,
    default_image_attribution = taxon.default_photo.attribution,
    default_image_licence = normalise_licence(
      taxon.default_photo.license_code
    )
  ) %>%
  filter(
    !is.na(species),
    species != ""
  ) %>%
  distinct(
    species,
    .keep_all = TRUE
  ) %>%
  mutate(
    default_image_allowed =
      !is.na(default_image_url) &
      default_image_url != "" &
      !is.na(default_image_licence) &
      default_image_licence %in% allowed_photo_licences
  )


# -------------------------------------------------------------------
# BUILD FINAL CARD DATA
# -------------------------------------------------------------------

signature_species_cards <- signature_species %>%
  left_join(
    hub_species_photos,
    by = c("species", "hub")
  ) %>%
  left_join(
    general_species_photos,
    by = "species"
  ) %>%
  left_join(
    default_taxon_photos,
    by = "species"
  ) %>%
  mutate(
    image_url = case_when(
      default_image_allowed ~ default_image_url,
      !is.na(hub_image_url) ~ hub_image_url,
      !is.na(general_image_url) ~ general_image_url,
      TRUE ~ NA_character_
    ),
    
    image_attribution = case_when(
      default_image_allowed ~ default_image_attribution,
      !is.na(hub_image_url) ~ hub_image_attribution,
      !is.na(general_image_url) ~ general_image_attribution,
      TRUE ~ NA_character_
    ),
    
    image_licence = case_when(
      default_image_allowed ~ default_image_licence,
      !is.na(hub_image_url) ~ hub_image_licence,
      !is.na(general_image_url) ~ general_image_licence,
      TRUE ~ NA_character_
    ),
    
    image_photo_id = case_when(
      default_image_allowed ~ NA_integer_,
      !is.na(hub_image_url) ~ hub_image_photo_id,
      !is.na(general_image_url) ~ general_image_photo_id,
      TRUE ~ NA_integer_
    ),
    
    image_source = case_when(
      default_image_allowed ~
        "iNaturalist default taxon image",
      
      !is.na(hub_image_url) ~
        "BRPC record from the signature hub",
      
      !is.na(general_image_url) ~
        paste0(
          "BRPC record from ",
          general_image_hub
        ),
      
      TRUE ~
        "No suitably licensed image found"
    )
  ) %>%
  select(
    hub,
    species,
    common_name,
    taxon_id,
    observations,
    overall_observations,
    expected_observations,
    ses,
    image_url,
    image_attribution,
    image_licence,
    image_photo_id,
    image_source
  ) %>%
  arrange(hub)


# -------------------------------------------------------------------
# CHECK PHOTO SELECTIONS
# -------------------------------------------------------------------

print(
  signature_species_cards %>%
    select(
      hub,
      species,
      image_source,
      image_licence
    ),
  n = Inf
)

missing_images <- signature_species_cards %>%
  filter(
    is.na(image_url) |
      image_url == ""
  )

if (nrow(missing_images) > 0) {
  warning(
    "No suitably licensed image was found for: ",
    paste(
      paste0(
        missing_images$hub,
        " — ",
        missing_images$species
      ),
      collapse = "; "
    )
  )
}

signature_species_outputs <- list(
  signature_species = signature_species,
  signature_species_cards = signature_species_cards,
  all_results = greedyqswap_results$all_results,
  hub_tables = greedyqswap_results$hub_tables,
  
  metadata = list(
    null_method = "greedyqswap",
    nsim = 999,
    seed = 123,
    ses_threshold = 2,
    minimum_overall_observations = 10,
    allowed_photo_licences = allowed_photo_licences,
    generated_at = Sys.time()
  )
)


saveRDS(
  signature_species_outputs,
  file = "Laura/Biodiversity report Q2 2026/analysis_data/signature_species_results.rds"
)

message(
  "Signature-species results saved to ",
  "Laura/Biodiversity report Q2 2026/analysis_data/signature_species_results.rds"
)


