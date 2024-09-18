library(tidyverse)

this_dir <- "~/Desktop/may10_runs/"
this_logs <- list.files(this_dir)[stringr::str_detect(list.files(this_dir), "^agent_log")]

for (i in 1:length(this_logs)) {

  if (stringr::str_detect(this_logs[[i]], "control") == TRUE) {
    next
  }

  this_demo <- read_delim(paste(this_dir, this_logs[[i]], sep = ""),
                          delim = ",") %>%
    dplyr::mutate(num_venues = ifelse(stringr::str_detect(venues_attended, "^s"), NA, venues_attended),
                  num_venues = nchar(stringr::str_replace_all(num_venues, "[a-z]|\\d", ""))+1,
                  num_venues = ifelse(is.na(num_venues), 0, num_venues),
                  num_apps = ifelse(stringr::str_detect(apps_used, "^s"), NA, apps_used),
                  num_apps = nchar(stringr::str_replace_all(num_apps, "[a-z]|\\d", ""))+1,
                  num_apps = ifelse(is.na(num_apps), 0, num_apps),
                  age_group2 = case_when(age < 20 ~ "16-19",
                                       age >= 20 & age < 25 ~ "20-24",
                                       TRUE ~ "25+"))

  age_props <- this_demo %>%
    group_by(age_group2) %>%
    summarize(prop_venues = mean(num_venues > 0),
              mean_venues = mean(num_venues),
              prop_apps = mean(num_apps > 0),
              mean_apps = mean(num_apps)) %>%
    mutate(sim = this_logs[[i]])

  race_props <- this_demo %>%
    group_by(race_ethnicity) %>%
    summarize(prop_venues = mean(num_venues > 0),
              mean_venues = mean(num_venues),
              prop_apps = mean(num_apps > 0),
              mean_apps = mean(num_apps)) %>%
    mutate(sim = this_logs[[i]])

  if (i == 1) {
    age_summaries <- age_props
    race_summaries <- race_props
  } else {
    age_summaries <- dplyr::bind_rows(age_summaries, age_props)
    race_summaries <- dplyr::bind_rows(race_summaries, race_props)
  }
}

age_summaries %>%
  group_by(age_group2) %>%
  summarize(prop_venues = mean(prop_venues),
            mean_venues = mean(mean_venues),
            prop_apps = mean(prop_apps),
            mean_apps = mean(mean_apps))

race_summaries %>%
  group_by(race_ethnicity) %>%
  summarize(prop_venues = mean(prop_venues),
            mean_venues = mean(mean_venues),
            prop_apps = mean(prop_apps),
            mean_apps = mean(mean_apps))
