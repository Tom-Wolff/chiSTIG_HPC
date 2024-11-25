library(tidyverse)

edgelist_overlap <- function(file) {

  # browser()

  # Read in file
  sim <- readRDS(file)

  full_main <- data.frame()
  full_casl <- data.frame()
  full_inst <- data.frame()

  # Extract new edges from `raw.records`. These edges have been stored to
  # keep track of the venues and apps associated with nodes week to week
  for (i in 1:length(sim$raw.records[[1]])) {

    if (sim$raw.records[[1]][[i]]$label == "infection_act") {
      next
    }

    this_obj <- sim$raw.records[[1]][[i]]$object
    this_time <- sim$raw.records[[1]][[i]]$time

    new_edges <- this_obj %>% filter(time == this_time)
    main <- new_edges %>% filter(type == 1)
    casl <- new_edges %>% filter(type == 2)
    inst <- new_edges %>% filter(type == 3)


      full_main <- dplyr::bind_rows(full_main, main)
      full_casl <- dplyr::bind_rows(full_casl, casl)
      full_inst <- dplyr::bind_rows(full_inst, inst)

  }

  # Create pasted "combined" variables of venues and apps. This
  # sets up for the sub-function I use to to detect overlap
  full_main <- full_main %>%
    mutate(venues_combined = paste(head_venues, tail_venues, sep = "|"),
           apps_combined = paste(head_apps, tail_apps, sep = "|"))

  full_casl <- full_casl %>%
    mutate(venues_combined = paste(head_venues, tail_venues, sep = "|"),
           apps_combined = paste(head_apps, tail_apps, sep = "|"))

  full_inst <- full_inst %>%
    mutate(venues_combined = paste(head_venues, tail_venues, sep = "|"),
           apps_combined = paste(head_apps, tail_apps, sep = "|"))

  full_main$venue_overlap <- detect_overlap(full_main$venues_combined)
  full_main$app_overlap <- detect_overlap(full_main$apps_combined)

  full_casl$venue_overlap <- detect_overlap(full_casl$venues_combined)
  full_casl$app_overlap <- detect_overlap(full_casl$apps_combined)

  full_inst$venue_overlap <- detect_overlap(full_inst$venues_combined)
  full_inst$app_overlap <- detect_overlap(full_inst$apps_combined)

  summary <- data.frame(file = file,
                        main_venue_overlap = mean(full_main$venue_overlap),
                        casl_venue_overlap = mean(full_casl$venue_overlap),
                        inst_venue_overlap = mean(full_inst$venue_overlap),
                        main_app_overlap = mean(full_main$app_overlap),
                        casl_app_overlap = mean(full_casl$app_overlap),
                        inst_app_overlap = mean(full_inst$app_overlap))

}

detect_overlap <- function(var) {
  split_vals <- stringr::str_split(var, "\\|")
  # This line sees whether the number of venues in the combined variable equals the number of unique venues. If it doesn't, we have overlap.
  has_overlap <- unlist(lapply(split_vals, length)) != unlist(lapply(split_vals, function(x){length(unique(x))}))
  return(has_overlap)
}

# Directory containing runs
this_dir <- "~/Desktop/sept12_runs"

files <- list.files(this_dir)
files <- files[stringr::str_detect(files, "rds$")]

for (i in 1:length(files)) {
  over_sums <- edgelist_overlap(paste(this_dir, files[[i]], sep = "/"))

  if (i == 1) {
    overlap_summaries <- over_sums
  } else {
    overlap_summaries <- dplyr::bind_rows(overlap_summaries, over_sums)
  }

}


overlap_summaries <- overlap_summaries %>%
  # mutate(treat = case_when(str_detect(file, "control") ~ "Control Model",
  #                          str_detect(file, "apps") ~ "Apps Only",
  #                          str_detect(file, "venues") ~ "Venues Only",
  #                          str_detect(file, "both") ~ "Venues and Apps",
  #                          TRUE ~ NA),
  #        target_shift = case_when(str_detect(file, "plus") ~ "Target Stats Increased",
  #                                 str_detect(file, "minus") ~ "Target Stats Decreased",
  #                                 TRUE ~ NA))
mutate(treat = case_when(str_detect(file, "apps") ~ "Apps",
                  str_detect(file, "both") ~ "Venues and Apps",
                  str_detect(file, "control") ~ "Control",
                  str_detect(file, "venues") ~ "Venues",
                  TRUE ~ NA),
      trial = str_extract(file, "\\d*_venue"))

overlap_summaries$trial <- as.numeric(str_replace_all(overlap_summaries$trial, "_venue", ""))

target_info <- read_delim("~/Downloads/edge_target_calibration_vals_test26sep.csv", delim = "\t") %>%
  rename(trial = fit_no)

overlap_summaries <- overlap_summaries %>%
  dplyr::left_join(target_info, by = "trial")

overlap_summaries %>%
  ggplot(aes(x = main_venue_overlap, fill = treat)) +
  geom_density() +
  # facet_grid(cols = vars(trial)) +
  theme_classic() +
  labs(x = "Proportion of Partnerships with Venue Overlap",
       title = "Main Partnerships")

overlap_summaries %>%
  ggplot(aes(x = casl_venue_overlap, fill = treat)) +
  geom_density() +
  # facet_grid(cols = vars(trial)) +
  theme_classic() +
  labs(x = "Proportion of Partnerships with Venue Overlap",
       title = "Casual Partnerships")

overlap_summaries %>%
  ggplot(aes(x = inst_venue_overlap, fill = treat)) +
  geom_density() +
  # facet_grid(cols = vars(trial)) +
  theme_classic() +
  labs(x = "Proportion of Partnerships with Venue Overlap",
       title = "One-Time Partnerships")

# APP

overlap_summaries %>%
  ggplot(aes(x = apps_main, y = main_app_overlap, color = treat)) +
  geom_point() +
  # facet_grid(cols = vars(trial)) +
  theme_classic() +
  labs(x = "Target Statistic",
       title = "Main Partnerships")

overlap_summaries %>%
  ggplot(aes(x = apps_casual, y = casl_app_overlap, color = treat)) +
  geom_point() +
  # facet_grid(cols = vars(trial)) +
  theme_classic() +
  labs(x = "Target Statistic",
       title = "Casual Partnerships")

overlap_summaries %>%
  ggplot(aes(x = inst_app_overlap, fill = treat)) +
  geom_density() +
  # facet_grid(cols = vars(trial)) +
  theme_classic() +
  labs(x = "Proportion of Partnerships with App Overlap",
       title = "One-Time Partnerships")
