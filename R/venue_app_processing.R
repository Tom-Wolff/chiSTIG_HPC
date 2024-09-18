library(tidyverse)

# Set directory for files
this_dir <- "~/Desktop/may10_runs/"

# Get names of `.rds` files
this_rds <- list.files(this_dir)[stringr::str_detect(list.files(this_dir), ".rds$")]

for (i in 1:length(this_rds)) {

  apps_treat <- FALSE
  venues_treat <- FALSE

  # Which treatment is this?
  if (stringr::str_detect(this_rds[[i]], "^apps")) {
    apps_treat <- TRUE
  }
  if (stringr::str_detect(this_rds[[i]], "venues")) {
    venues_treat <- TRUE
  }
  if (stringr::str_detect(this_rds[[i]], "both")) {
    apps_treat <- TRUE
    venues_treat <- TRUE
  }


this_sim <- readRDS(paste(this_dir, this_rds[[i]], sep = ""))

sim_list <- as.list(this_sim)

# Combine pertinent elements of `sim_list` to a single data frame
full_data <- dplyr::bind_rows(lapply(sim_list$raw.records[[1]], function(x){x$object}))

# Create a separate object containing just the infection acts
inf_events <- full_data %>% dplyr::filter(type == 4) %>%
  ### Select only relevant variables
  select(head_uid, tail_uid,
         ### Time of infection might be useful
         inf_time = time) %>%
  ### Infection indicator to confirm merge
  mutate(infection = 1)

# Heads and tails of the act list are arranged by each node's relative position
# in a sexual act. To maximize merging, we need to create a second `inf_events`
# data frame where head and tail are reversed, then combine this with the original
# `inf_events`
inf_events2 <- inf_events %>% dplyr::select(head_uid = tail_uid,
                                            tail_uid = head_uid,
                                            inf_time, infection)
inf_events <- dplyr::bind_rows(inf_events, inf_events2)

# Now let's store our actual edgelist in its own object
edges <- full_data %>% dplyr::filter(type != 4) %>%
  group_by(head_uid, tail_uid) %>%
  slice(1)

# Merge in `inf_events`
merged_data <- edges %>% dplyr::left_join(inf_events, by = c("head_uid", "tail_uid")) %>%
  # Keep just last ten years
  dplyr::filter(time > 3120)

# What percentage of cases in the infection act list successfully merged?
# Note that we use `nrow(inf_events2)` here for the accurate denominator

# merged_data %>% group_by(type) %>%
#   summarize(sum(infection, na.rm = TRUE))


# Let's read in the agent log
head_demo <- read_delim(paste(this_dir, "agent_log_apps_mar20test.txt", sep = ""),
           delim = ",") %>%
  select(time = tick, head_uid = agent_id, head_age = age, head_agegroup = age_group,
         head_race = race_ethnicity, head_hiv = hiv_status)
tail_demo <- head_demo %>%
  rename(tail_uid = head_uid, tail_age = head_age, tail_agegroup = head_agegroup,
         tail_race = head_race, tail_hiv = head_hiv)

merged_data <- merged_data %>%
  dplyr::left_join(head_demo, by = c("head_uid", "time")) %>%
  dplyr::left_join(tail_demo, by = c("tail_uid", "time"))


# Discordant edges at time of partnership formation
test <- merged_data %>%
  dplyr::mutate(hiv_disc = (head_hiv + tail_hiv == 1)) %>%
  dplyr::filter(hiv_disc == TRUE)

test %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(count = n())

# Now extract just the partnerships that led to infection
inf_edges <- merged_data %>%
  filter(infection == 1)

# How many by partnership type?
inf_edges %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(count = n())

# Loop over infection events
for (j in 1:nrow(inf_edges)) {

  # Get this infection event
  this_inf <- inf_edges[j,]

# APPS

  # if (apps_treat == TRUE) {
  # # Which apps did both use at tie formation?
  #     shared_apps <- unlist(stringr::str_split(this_inf$head_apps, "\\|"))[unlist(stringr::str_split(this_inf$head_apps, "\\|")) %in% unlist(stringr::str_split(this_inf$tail_apps, "\\|"))]
  #
      # if (length(shared_apps) == 0) {
      #   shared_apps <- "None"
      # }
  # } else {
  #   shared_apps <- "None"
  # }
  shared_apps <- unlist(stringr::str_split(this_inf$head_apps, "\\|"))[unlist(stringr::str_split(this_inf$head_apps, "\\|")) %in% unlist(stringr::str_split(this_inf$tail_apps, "\\|"))]

  if (length(shared_apps) == 0) {
    shared_apps <- "None"
  }

  # Store apps from this infection event
  this_apps <- data.frame(apps = shared_apps,
                          type = this_inf$type,
                          time = this_inf$time,
                          head_age = this_inf$head_age,
                          head_agegroup = this_inf$head_agegroup,
                          head_race = this_inf$head_race,
                          head_hiv = this_inf$head_hiv,
                          tail_age = this_inf$tail_age,
                          tail_agegroup = this_inf$tail_agegroup,
                          tail_race = this_inf$tail_race,
                          tail_hiv = this_inf$tail_hiv)

# VENUES
  # if (venues_treat == TRUE) {
  #   # Which apps did both use at tie formation?
  #   shared_venues <- unlist(stringr::str_split(this_inf$head_venues, "\\|"))[unlist(stringr::str_split(this_inf$head_venues, "\\|")) %in% unlist(stringr::str_split(this_inf$tail_venues, "\\|"))]
  #
  #   if (length(shared_venues) == 0) {
  #     shared_venues <- "None"
  #   }
  # } else {
  #   shared_venues <- "None"
  # }

  shared_venues <- unlist(stringr::str_split(this_inf$head_venues, "\\|"))[unlist(stringr::str_split(this_inf$head_venues, "\\|")) %in% unlist(stringr::str_split(this_inf$tail_venues, "\\|"))]

  if (length(shared_venues) == 0) {
        shared_venues <- "None"
       }

  # Store venues from this infection event
  this_venues <- data.frame(venues = shared_venues,
                          type = this_inf$type,
                          time = this_inf$time,
                          head_age = this_inf$head_age,
                          head_agegroup = this_inf$head_agegroup,
                          head_race = this_inf$head_race,
                          head_hiv = this_inf$head_hiv,
                          tail_age = this_inf$tail_age,
                          tail_agegroup = this_inf$tail_agegroup,
                          tail_race = this_inf$tail_race,
                          tail_hiv = this_inf$tail_hiv)

  if (j == 1) {
    inf_apps <- this_apps
    inf_venues <- this_venues
  } else {
    inf_apps <- dplyr::bind_rows(inf_apps, this_apps)
    inf_venues <- dplyr::bind_rows(inf_venues, this_venues)
  }
# End loop over infection events
}

# Filter out cases without matches
inf_apps <- inf_apps %>% dplyr::filter(apps != "None")
inf_venues <- inf_venues %>% dplyr::filter(venues != "None")

}

