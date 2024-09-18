# Venue and App Assignment/Attendance

library(tidyverse)
library(foreach)

# Set directory for files
this_dir <- "~/Desktop/20mar_test/19mar_20runs/"

# Get names of `.rds` files
logs <- list.files(this_dir)[stringr::str_detect(list.files(this_dir), "^agent_log")]

for (i in 1:length(logs)) {

  print(i)

  if (stringr::str_detect(logs[[i]], "control") == TRUE) {
    next
  }

  this_log <- read_delim(paste(this_dir, logs[[i]], sep = ""),
                         delim = ",")


  # Get list of unique venues and apps
  unique_venues <- unique(unlist(stringr::str_split(this_log$venues_attended, pattern = "\\|")))
  ### Remove anything starting with `"s"`
  unique_venues <- sort(unique_venues[!stringr::str_detect(unique_venues, "^s")])

  unique_apps <- unique(unlist(stringr::str_split(this_log$apps_used, pattern = "\\|")))
  ### Remove anything starting with `"s"`
  unique_apps <- sort(unique_apps[!stringr::str_detect(unique_apps, "^s")])

  # To avoid exceedingly long data, let's subset by strata first
  # ### Race
  # black <- this_log %>% filter(race_ethnicity == "blackNH")
  # hisp <- this_log %>% filter(race_ethnicity == "hispanic")
  # other <- this_log %>% filter(race_ethnicity == "otherNH")
  # white <- this_log %>% filter(race_ethnicity == "whiteNH")
  # ### Age Group
  # age16_20 <- this_log %>% filter(age_group == "16to20")
  # age21_29 <- this_log %>% filter(age_group == "21to29")
  # ### HIV Status
  # hivpos <- this_log %>% filter(hiv_status == 1)
  # hivneg <- this_log %>% filter(hiv_status == 0)

  ### Race x Age
  black16 <- this_log %>% filter(race_ethnicity == "blackNH") %>% filter(age_group == "16to20")
  black21 <- this_log %>% filter(race_ethnicity == "blackNH") %>% filter(age_group == "21to29")
  hisp16 <- this_log %>% filter(race_ethnicity == "hispanic") %>% filter(age_group == "16to20")
  hisp21 <- this_log %>% filter(race_ethnicity == "hispanic") %>% filter(age_group == "21to29")
  other16 <- this_log %>% filter(race_ethnicity == "otherNH") %>% filter(age_group == "16to20")
  other21 <- this_log %>% filter(race_ethnicity == "otherNH") %>% filter(age_group == "21to29")
  white16 <- this_log %>% filter(race_ethnicity == "whiteNH")  %>% filter(age_group == "16to20")
  white21 <- this_log %>% filter(race_ethnicity == "whiteNH")  %>% filter(age_group == "21to29")
  # ### Age Group
  # age16_20 <- this_log %>% filter(age_group == "16to20")
  # age21_29 <- this_log %>% filter(age_group == "21to29")
  # ### HIV Status
  # hivpos <- this_log %>% filter(hiv_status == 1)
  # hivneg <- this_log %>% filter(hiv_status == 0)

  # Setup parallelization
  cl <- parallel::makePSOCKcluster(8)
  doParallel::registerDoParallel(cl)


  # RECORDING PREVALENCE OF APPS BY STRATUM (PARALLEL)
  app_prev <- foreach(j = 1:length(unique_apps),
                  .combine = "rbind") %dopar% {

                    data.frame(app = unique_apps[[j]],

                               # prev_black = sum(stringr::str_detect(black$apps_used, unique_apps[[j]]))/nrow(black),
                               # prev_hisp = sum(stringr::str_detect(hisp$apps_used, unique_apps[[j]]))/nrow(hisp),
                               # prev_other = sum(stringr::str_detect(other$apps_used, unique_apps[[j]]))/nrow(other),
                               # prev_white = sum(stringr::str_detect(white$apps_used, unique_apps[[j]]))/nrow(white),
                               #
                               # prev_16to20 = sum(stringr::str_detect(age16_20$apps_used, unique_apps[[j]]))/nrow(age16_20),
                               # prev_21to29 = sum(stringr::str_detect(age21_29$apps_used, unique_apps[[j]]))/nrow(age21_29),
                               #
                               # prev_hivpos = sum(stringr::str_detect(hivpos$apps_used, unique_apps[[j]]))/nrow(hivpos),
                               # prev_hivneg = sum(stringr::str_detect(hivneg$apps_used, unique_apps[[j]]))/nrow(hivneg))

                               prev_black16 = sum(stringr::str_detect(black16$apps_used, unique_apps[[j]]))/nrow(black16),
                               prev_black21 = sum(stringr::str_detect(black21$apps_used, unique_apps[[j]]))/nrow(black21),
                               prev_hisp16 = sum(stringr::str_detect(hisp16$apps_used, unique_apps[[j]]))/nrow(hisp16),
                               prev_hisp21 = sum(stringr::str_detect(hisp21$apps_used, unique_apps[[j]]))/nrow(hisp21),
                               prev_other16 = sum(stringr::str_detect(other16$apps_used, unique_apps[[j]]))/nrow(other16),
                               prev_other21 = sum(stringr::str_detect(other21$apps_used, unique_apps[[j]]))/nrow(other21),
                               prev_white16 = sum(stringr::str_detect(white16$apps_used, unique_apps[[j]]))/nrow(white16),
                               prev_white21 = sum(stringr::str_detect(white21$apps_used, unique_apps[[j]]))/nrow(white21))


                  }

  # overall_app <- data.frame(app = "Any apps",
  #
  #                           prev_black = sum(!str_detect(black$apps_used, "^s"))/nrow(black),
  #                           prev_hisp =  sum(!str_detect(hisp$apps_used,  "^s"))/nrow(hisp),
  #                           prev_other = sum(!str_detect(other$apps_used, "^s"))/nrow(other),
  #                           prev_white = sum(!str_detect(white$apps_used, "^s"))/nrow(white),
  #
  #                           prev_16to20 = sum(!str_detect(age16_20$apps_used, "^s"))/nrow(age16_20),
  #                           prev_21to29 = sum(!str_detect(age21_29$apps_used, "^s"))/nrow(age21_29),
  #
  #                           prev_hivpos = sum(!str_detect(hivpos$apps_used, "^s"))/nrow(hivpos),
  #                           prev_hivneg = sum(!str_detect(hivneg$apps_used, "^s"))/nrow(hivneg))
  #
  #
  #
  # app_prev <- dplyr::bind_rows(overall_app, app_prev)
  app_prev$sim <- logs[[i]]


  venue_prev <- foreach(j = 1:length(unique_venues),
                         .combine = "rbind") %dopar% {

                           data.frame(venue = unique_venues[[j]],

                                      # prev_black =  sum(stringr::str_detect(black$venues_attended, unique_venues[[j]]))/nrow(black),
                                      # prev_hisp =   sum(stringr::str_detect(hisp$venues_attended, unique_venues[[j]]))/nrow(hisp),
                                      # prev_other =  sum(stringr::str_detect(other$venues_attended, unique_venues[[j]]))/nrow(other),
                                      # prev_white =  sum(stringr::str_detect(white$venues_attended, unique_venues[[j]]))/nrow(white),
                                      #
                                      # prev_16to20 = sum(stringr::str_detect(age16_20$venues_attended, unique_venues[[j]]))/nrow(age16_20),
                                      # prev_21to29 = sum(stringr::str_detect(age21_29$venues_attended, unique_venues[[j]]))/nrow(age21_29),
                                      #
                                      # prev_hivpos = sum(stringr::str_detect(hivpos$venues_attended, unique_venues[[j]]))/nrow(hivpos),
                                      # prev_hivneg = sum(stringr::str_detect(hivneg$venues_attended, unique_venues[[j]]))/nrow(hivneg))

                           prev_black16 = sum(stringr::str_detect(black16$venues_attended, unique_venues[[j]]))/nrow(black16),
                           prev_black21 = sum(stringr::str_detect(black21$venues_attended, unique_venues[[j]]))/nrow(black21),
                           prev_hisp16 = sum(stringr::str_detect(hisp16$venues_attended, unique_venues[[j]]))/nrow(hisp16),
                           prev_hisp21 = sum(stringr::str_detect(hisp21$venues_attended, unique_venues[[j]]))/nrow(hisp21),
                           prev_other16 = sum(stringr::str_detect(other16$venues_attended, unique_venues[[j]]))/nrow(other16),
                           prev_other21 = sum(stringr::str_detect(other21$venues_attended, unique_venues[[j]]))/nrow(other21),
                           prev_white16 = sum(stringr::str_detect(white16$venues_attended, unique_venues[[j]]))/nrow(white16),
                           prev_white21 = sum(stringr::str_detect(white21$venues_attended, unique_venues[[j]]))/nrow(white21))
                         }


  # overall_venue <- data.frame(venue = "Any venues",
  #
  #                           prev_black = sum(!str_detect(black$venues_attended, "^s"))/nrow(black),
  #                           prev_hisp =  sum(!str_detect(hisp$venues_attended,  "^s"))/nrow(hisp),
  #                           prev_other = sum(!str_detect(other$venues_attended, "^s"))/nrow(other),
  #                           prev_white = sum(!str_detect(white$venues_attended, "^s"))/nrow(white),
  #
  #                           prev_16to20 = sum(!str_detect(age16_20$venues_attended, "^s"))/nrow(age16_20),
  #                           prev_21to29 = sum(!str_detect(age21_29$venues_attended, "^s"))/nrow(age21_29),
  #
  #                           prev_hivpos = sum(!str_detect(hivpos$venues_attended, "^s"))/nrow(hivpos),
  #                           prev_hivneg = sum(!str_detect(hivneg$venues_attended, "^s"))/nrow(hivneg))
  #
  # venue_prev <- dplyr::bind_rows(overall_venue, venue_prev)
  venue_prev$sim <- logs[[i]]

  parallel::stopCluster(cl)

  if (i == 1) {
    all_apps <- app_prev
    all_venues <- venue_prev
  } else {
    all_apps <- dplyr::bind_rows(all_apps, app_prev)
    all_venues <- dplyr::bind_rows(all_venues, venue_prev)
  }

}


# Aggregating together in a single table

final_apps <- data.frame(app = rep(unique(sort(all_apps$app)), each = length(unique(all_apps$sim))),
                         sim = rep(sort(unique(all_apps$sim)), length(unique(all_apps$app)))) %>%
  dplyr::left_join(all_apps, by = c("app", "sim"))
final_apps[is.na(final_apps)] <- 0


final_venues <- data.frame(venue = rep(unique(sort(all_venues$venue)), each = length(unique(all_venues$sim))),
                           sim = rep(sort(unique(all_venues$sim)), length(unique(all_venues$venue)))) %>%
  dplyr::left_join(all_venues, by = c("venue", "sim"))
final_venues[is.na(final_venues)] <- 0

final_means_apps <- final_apps %>%
  dplyr::select(-sim) %>%
  dplyr::group_by(app) %>%
  dplyr::summarize_all(mean)

final_means_venues <- final_venues %>%
  dplyr::select(-sim) %>%
  dplyr::group_by(venue) %>%
  dplyr::summarize_all(mean)

write.csv(final_means_apps, "~/Desktop/app_prev_19mar.csv")
write.csv(final_means_venues, "~/Desktop/venue_prev_19mar.csv")


library(tidyverse)

app_matrix <- read.csv("~/Downloads/empop_ego_app_matrix_binary.csv")
demos <- read.csv("~/Downloads/empop_egoid_demo_def.csv")

app_long <- app_matrix %>%
  pivot_longer(starts_with("a"),
               names_to = "app",
               values_to = "use") %>%
  filter(use == 1) %>%
  rename(egoid = X) %>%
  left_join(demos, by = "egoid") %>%
  select(app, demo8, use) %>%
  group_by(app, demo8) %>%
  summarize(num_using = sum(use))

demo_sizes <- demos %>%
  group_by(demo8) %>%
  summarize(demo_total = n())

demo_labels <- demos %>%
  group_by(demo8) %>%
  slice(1) %>%
  mutate(stratum = paste(race_ethnicity, age, sep = "_")) %>%
  select(demo8, stratum)

app_merge <- app_long %>%
  left_join(demo_sizes, by = "demo8") %>%
  mutate(prev = num_using/demo_total) %>%
  left_join(demo_labels, by = "demo8") %>%
  select(app, prev, stratum) %>%
  pivot_wider(id_cols = app,
              names_from = stratum,
              values_from = prev#,
              #values_fill = 0
              )

# Rename final_means_apps columns
final_means_apps2 <- final_means_apps
colnames(final_means_apps2) <- stringr::str_replace_all(colnames(final_means_apps2), "21", "21to29")
colnames(final_means_apps2) <- stringr::str_replace_all(colnames(final_means_apps2), "16", "under21")
colnames(final_means_apps2) <- stringr::str_replace_all(colnames(final_means_apps2), "prev_black", "blackNH_")
colnames(final_means_apps2) <- stringr::str_replace_all(colnames(final_means_apps2), "prev_white", "whiteNH_")
colnames(final_means_apps2) <- stringr::str_replace_all(colnames(final_means_apps2), "prev_hisp", "hispanic_")
colnames(final_means_apps2) <- stringr::str_replace_all(colnames(final_means_apps2), "prev_other", "otherNH_")
colnames(final_means_apps2) <- paste(colnames(final_means_apps2), "_actual", sep = "")
colnames(final_means_apps2)[[1]] <- "app"

app_merge2 <- app_merge %>%
  dplyr::left_join(final_means_apps2, by = "app")

app_merge3 <- app_merge2 %>%
  filter(app %in% final_apps$app) %>%
  mutate(blackNH_under21 = ifelse(blackNH_under21 == 0, NaN, blackNH_under21_actual - blackNH_under21),
         blackNH_21to29 =   ifelse(blackNH_21to29 == 0, NaN, blackNH_21to29_actual - blackNH_21to29),
         hispanic_under21 = ifelse(hispanic_under21 == 0, NaN, hispanic_under21_actual - hispanic_under21),
         hispanic_21to29 =  ifelse(hispanic_21to29  == 0, NaN, hispanic_21to29_actual - hispanic_21to29),
         otherNH_under21 =  ifelse(otherNH_under21  == 0, NaN, otherNH_under21_actual - otherNH_under21),
         otherNH_21to29 =   ifelse(otherNH_21to29 == 0, NaN, otherNH_21to29_actual - otherNH_21to29),
         whiteNH_under21 =  ifelse(whiteNH_under21  == 0, NaN, whiteNH_under21_actual - whiteNH_under21),
         whiteNH_21to29 =   ifelse(whiteNH_21to29 == 0, NaN, whiteNH_21to29_actual - whiteNH_21to29)
         ) %>%
  select(-ends_with("actual"))

app_merge4 <- app_merge2 %>%
  filter(app %in% final_apps$app) %>%
  mutate(blackNH_under21 =    (blackNH_under21_actual - blackNH_under21)/blackNH_under21,
         blackNH_21to29 =       (blackNH_21to29_actual - blackNH_21to29)/blackNH_21to29,
         hispanic_under21 = (hispanic_under21_actual - hispanic_under21)/hispanic_under21,
         hispanic_21to29 =    (hispanic_21to29_actual - hispanic_21to29)/hispanic_21to29,
         otherNH_under21 =    (otherNH_under21_actual - otherNH_under21)/otherNH_under21,
         otherNH_21to29 =       (otherNH_21to29_actual - otherNH_21to29)/otherNH_21to29,
         whiteNH_under21 =    (whiteNH_under21_actual - whiteNH_under21)/whiteNH_under21,
         whiteNH_21to29 =   (whiteNH_21to29_actual - whiteNH_21to29)/whiteNH_21to29) %>%
  select(-ends_with("actual"))


final_means_apps

means_heat <- final_means_apps2 %>%
  tidyr::pivot_longer(cols = ends_with("actual"),
                      names_to = "stratum",
                      values_to = "value") %>%
  dplyr::mutate(stratum = str_replace(stratum, "_actual", ""),
                stratum = str_replace(stratum, "under21", "16to20")) %>%
  ggplot2::ggplot(ggplot2::aes(x = as.factor(stratum),
                               y = forcats::fct_rev(as.factor(app)),
                               fill = value)) +
  ggplot2::geom_tile(color = "black") +
  ggplot2::geom_text(ggplot2::aes(label = round(value, digits = 2)), color = "white", size = 4) +
  ggplot2::theme_minimal() +
  # scale_x_continuous(breaks = 1:max(cluster_edgelist$cluster)) +
  #   scale_y_continuous(breaks = 1:max(cluster_edgelist$cluster)) +
  # scale_y_reverse() +
  ggplot2::labs(x = "\nStratum", y = "App\n", title = "Proportion of Person-Weeks Using Apps, by Race x Age Group\n") +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(hjust = 0.5))


diffs_heat <- app_merge3 %>%
  tidyr::pivot_longer(cols = contains("_"),
                      names_to = "stratum",
                      values_to = "value") %>%
  dplyr::mutate(stratum = str_replace(stratum, "_actual", ""),
                stratum = str_replace(stratum, "under21", "16to20")) %>%
  # dplyr::mutate(value = ifelse(is.nan(value), 99, value)) %>%
  # dplyr::filter(!is.na(value)) %>%
  # dplyr::mutate(value = ifelse(value == 99, NA, value)) %>%
  ggplot2::ggplot(ggplot2::aes(x = as.factor(stratum),
                               y = forcats::fct_rev(as.factor(app)),
                               fill = value)) +
  ggplot2::geom_tile(color = "black") +
  ggplot2::geom_text(ggplot2::aes(label = round(value, digits = 2)), color = "white", size = 4) +
  ggplot2::theme_minimal() +
  # scale_x_continuous(breaks = 1:max(cluster_edgelist$cluster)) +
  #   scale_y_continuous(breaks = 1:max(cluster_edgelist$cluster)) +
  # scale_y_reverse() +
  ggplot2::labs(x = "\nStratum", y = "App\n", title = "Difference in Proportion of Person-Weeks Using Apps by Race x Age Group\n(Simulation Average - Average in Synthetic Population)\n") +
ggplot2::scale_fill_gradient2(low="#ad0206", mid = "grey", high="blue", #colors in the scale
                              midpoint=0,    #same midpoint for plots (mean of the range)
                              breaks=seq(-.15,
                                         .15, (.3/5)), #breaks in the scale bar
                              limits=c(-.15,
                                       .15))

diffs_heat2 <- app_merge3 %>%
  tidyr::pivot_longer(cols = contains("_"),
                      names_to = "stratum",
                      values_to = "value") %>%
  dplyr::mutate(stratum = str_replace(stratum, "_actual", ""),
                stratum = str_replace(stratum, "under21", "16to20")) %>%
  filter(!is.na(value))

std_diffs_heat <- app_merge4 %>%
  tidyr::pivot_longer(cols = contains("_"),
                      names_to = "stratum",
                      values_to = "value") %>%
  dplyr::mutate(stratum = str_replace(stratum, "_actual", ""),
                stratum = str_replace(stratum, "under21", "16to20")) %>%
  # dplyr::mutate(value = ifelse(is.nan(value), 99, value)) %>%
  # dplyr::filter(!is.na(value)) %>%
  # dplyr::mutate(value = ifelse(value == 99, NA, value)) %>%
  ggplot2::ggplot(ggplot2::aes(x = as.factor(stratum),
                               y = forcats::fct_rev(as.factor(app)),
                               fill = value)) +
  ggplot2::geom_tile(color = "black") +
  ggplot2::geom_text(ggplot2::aes(label = round(value, digits = 2)), color = "white", size = 4) +
  ggplot2::theme_minimal() +
  # scale_x_continuous(breaks = 1:max(cluster_edgelist$cluster)) +
  #   scale_y_continuous(breaks = 1:max(cluster_edgelist$cluster)) +
  # scale_y_reverse() +
  ggplot2::labs(x = "\nStratum", y = "App\n", title = "Percent Difference in Proportion of Person-Weeks Using Apps by Race x Age Group\n(Simulation Average - Average in Synthetic Population)/Average in Synthetic Population\n") +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(hjust = 0.5)) +
ggplot2::scale_fill_gradient2(low="#ad0206", mid = "grey", high="blue", #colors in the scale
                              midpoint=0,    #same midpoint for plots (mean of the range)
                              breaks=seq(-1,
                                         1, .5), #breaks in the scale bar
                              limits=c(-1,
                                       1))


###################################
#    A N Y   A P P   U S A G E    #
###################################


for (i in 1:length(logs)) {

  print(i)

  if (stringr::str_detect(logs[[i]], "control") == TRUE) {
    next
  }

  this_log <- read_delim(paste(this_dir, logs[[i]], sep = ""),
                         delim = ",")


  # Get list of unique venues and apps
  unique_venues <- unique(unlist(stringr::str_split(this_log$venues_attended, pattern = "\\|")))
  ### Remove anything starting with `"s"`
  unique_venues <- sort(unique_venues[!stringr::str_detect(unique_venues, "^s")])

  unique_apps <- unique(unlist(stringr::str_split(this_log$apps_used, pattern = "\\|")))
  ### Remove anything starting with `"s"`
  unique_apps <- sort(unique_apps[!stringr::str_detect(unique_apps, "^s")])

  # To avoid exceedingly long data, let's subset by strata first
  # ### Race
  # black <- this_log %>% filter(race_ethnicity == "blackNH")
  # hisp <- this_log %>% filter(race_ethnicity == "hispanic")
  # other <- this_log %>% filter(race_ethnicity == "otherNH")
  # white <- this_log %>% filter(race_ethnicity == "whiteNH")
  # ### Age Group
  # age16_20 <- this_log %>% filter(age_group == "16to20")
  # age21_29 <- this_log %>% filter(age_group == "21to29")
  # ### HIV Status
  # hivpos <- this_log %>% filter(hiv_status == 1)
  # hivneg <- this_log %>% filter(hiv_status == 0)

  ### Race x Age
  black16 <- this_log %>% filter(race_ethnicity == "blackNH") %>% filter(age_group == "16to20")
  black21 <- this_log %>% filter(race_ethnicity == "blackNH") %>% filter(age_group == "21to29")
  hisp16 <- this_log %>% filter(race_ethnicity == "hispanic") %>% filter(age_group == "16to20")
  hisp21 <- this_log %>% filter(race_ethnicity == "hispanic") %>% filter(age_group == "21to29")
  other16 <- this_log %>% filter(race_ethnicity == "otherNH") %>% filter(age_group == "16to20")
  other21 <- this_log %>% filter(race_ethnicity == "otherNH") %>% filter(age_group == "21to29")
  white16 <- this_log %>% filter(race_ethnicity == "whiteNH")  %>% filter(age_group == "16to20")
  white21 <- this_log %>% filter(race_ethnicity == "whiteNH")  %>% filter(age_group == "21to29")
  # ### Age Group
  # age16_20 <- this_log %>% filter(age_group == "16to20")
  # age21_29 <- this_log %>% filter(age_group == "21to29")
  # ### HIV Status
  # hivpos <- this_log %>% filter(hiv_status == 1)
  # hivneg <- this_log %>% filter(hiv_status == 0)

  overall_app <- data.frame(app = "Any apps",

                            prev_black16 = sum(!str_detect(black16$apps_used, "^s"))/nrow(black16),
                            prev_hisp16 =  sum(!str_detect(hisp16$apps_used,  "^s"))/nrow(hisp16),
                            prev_other16 = sum(!str_detect(other16$apps_used, "^s"))/nrow(other16),
                            prev_white16 = sum(!str_detect(white16$apps_used, "^s"))/nrow(white16),

                            prev_black21 = sum(!str_detect(black21$apps_used, "^s"))/nrow(black21),
                            prev_hisp21 =  sum(!str_detect(hisp21$apps_used,  "^s"))/nrow(hisp21),
                            prev_other21 = sum(!str_detect(other21$apps_used, "^s"))/nrow(other21),
                            prev_white21 = sum(!str_detect(white21$apps_used, "^s"))/nrow(white21))



  # app_prev <- dplyr::bind_rows(overall_app, app_prev)
  overall_app$sim <- logs[[i]]





  overall_venue <- data.frame(app = "Any venues",

                            prev_black16 = sum(!str_detect(black16$venues_attended, "^s"))/nrow(black16),
                            prev_hisp16 =  sum(!str_detect(hisp16$venues_attended,  "^s"))/nrow(hisp16),
                            prev_other16 = sum(!str_detect(other16$venues_attended, "^s"))/nrow(other16),
                            prev_white16 = sum(!str_detect(white16$venues_attended, "^s"))/nrow(white16),

                            prev_black21 = sum(!str_detect(black21$venues_attended, "^s"))/nrow(black21),
                            prev_hisp21 =  sum(!str_detect(hisp21$venues_attended,  "^s"))/nrow(hisp21),
                            prev_other21 = sum(!str_detect(other21$venues_attended, "^s"))/nrow(other21),
                            prev_white21 = sum(!str_detect(white21$venues_attended, "^s"))/nrow(white21))
  overall_venue$sim <- logs[[i]]


  if (i == 1) {
    all_apps <- overall_app
    all_venues <- overall_venue
  } else {
    all_apps <- dplyr::bind_rows(all_apps, overall_app)
    all_venues <- dplyr::bind_rows(all_venues, overall_venue)
  }

}

final_means_apps <- all_apps %>%
  dplyr::select(-sim) %>%
  dplyr::group_by(app) %>%
  dplyr::summarize_all(mean) %>%
  select(-app)

colnames(final_means_apps) <- stringr::str_replace_all(colnames(final_means_apps), "21", "21to29")
colnames(final_means_apps) <- stringr::str_replace_all(colnames(final_means_apps), "16", "under21")
colnames(final_means_apps) <- stringr::str_replace_all(colnames(final_means_apps), "prev_black", "blackNH_")
colnames(final_means_apps) <- stringr::str_replace_all(colnames(final_means_apps), "prev_white", "whiteNH_")
colnames(final_means_apps) <- stringr::str_replace_all(colnames(final_means_apps), "prev_hisp", "hispanic_")
colnames(final_means_apps) <- stringr::str_replace_all(colnames(final_means_apps), "prev_other", "otherNH_")

final_means_apps$measure <- "prev_simulation"

synthpop_usage <- data.frame(egoid = app_matrix[,1],
                             use_apps = rowSums(app_matrix[,2:ncol(app_matrix)]) > 0)

synthpop_usage2 <- demos %>%
  left_join(synthpop_usage, by = "egoid") %>%
  mutate(use_apps = ifelse(is.na(use_apps), FALSE, TRUE)) %>%
  group_by(race_ethnicity, age) %>%
  summarize(use_apps = mean(use_apps)) %>%
  ungroup() %>%
  mutate(stratum = paste(race_ethnicity, age, sep = "_")) %>%
  select(stratum, use_apps) %>%
  pivot_wider(names_from = stratum, values_from = use_apps) %>%
  mutate(measure = "prev_synthpop") %>%
  select(measure,
         blackNH_under21, blackNH_21to29,
         hispanic_under21, hispanic_21to29,
         otherNH_under21, otherNH_21to29,
         whiteNH_under21, whiteNH_21to29)

app_prev_values <- bind_rows(synthpop_usage2, final_means_apps)

app_diffs <- app_prev_values[2, 2:9] - app_prev_values[1, 2:9]
app_diffs$measure <- "difference_in_proportion"

pct_diffs <- (app_prev_values[2, 2:9] - app_prev_values[1, 2:9])/app_prev_values[1, 2:9]
pct_diffs$measure <- "pct_difference_in_proportion"

app_prev_values2 <- bind_rows(app_prev_values, app_diffs, pct_diffs)[c(2, 1, 3, 4),]

