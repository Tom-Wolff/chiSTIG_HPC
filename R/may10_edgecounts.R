library(tidyverse)

# Set directory for files
this_dir <- "~/Desktop/21jun_runs/"

# Get names of `.rds` files
this_rds <- list.files(this_dir)[stringr::str_detect(list.files(this_dir), ".rds$")]
#### Make a dataframe for merging
rds_df <- data.frame(rds = this_rds) %>%
  dplyr::mutate(apps_treat = stringr::str_detect(rds, "^apps_"),
                venues_treat = stringr::str_detect(rds, "^venues_"),
                appsvenues_treat = stringr::str_detect(rds, "^both"),
                num = as.numeric(stringr::str_extract(rds, "\\d+")))


# Get names of agent log files
this_log <- list.files(this_dir)[stringr::str_detect(list.files(this_dir), "^agent_log")]
#### Make a dataframe for merging
log_df <- data.frame(log = this_log) %>%
  dplyr::mutate(apps_treat = stringr::str_detect(log, "^agent_log_a"),
                venues_treat = stringr::str_detect(log, "^agent_log_b"),
                appsvenues_treat = stringr::str_detect(log, "^agent_log_c"),
                num = as.numeric(stringr::str_extract(log, "\\d+")))

# Merge
file_index <- rds_df %>% dplyr::left_join(log_df)




for (i in 1:nrow(file_index)) {

  print(i)

  file_info <- file_index[i,]

  apps_treat <- file_info$apps_treat | file_info$appsvenues_treat
  venues_treat <- file_info$venues_treat | file_info$appsvenues_treat

  this_sim <- readRDS(paste(this_dir, file_info$rds, sep = ""))

  sim_list <- as.list(this_sim)

  epi_df <- as.data.frame(sim_list$epi)
  epi_df$time <- 1:nrow(epi_df)
  colnames(epi_df) <- names(sim_list$epi)
  epi_df$apps_treat <- apps_treat
  epi_df$venues_treat <- venues_treat
  epi_df$sim <- file_info$num

  if (i == 1) {
    full_epi <- epi_df
  } else {
    full_epi <- dplyr::bind_rows(full_epi, epi_df)
  }

}

full_epi <- full_epi %>%
  mutate(treatment = case_when((venues_treat == FALSE & apps_treat == FALSE) ~ "Control",
                               (venues_treat == TRUE & apps_treat == FALSE) ~ "Venues Only",
                               (venues_treat == FALSE & apps_treat == TRUE) ~ "Apps Only",
                               (venues_treat == TRUE & apps_treat == TRUE) ~ "Full Model",
                               TRUE ~ NA)#,
         # time = `...249`


) %>%
  group_by(treatment) %>%
  mutate(time = row_number()) %>%
  ungroup()



# GENERAL NUMBER OF EDGES
### Read in `netstats` to get some targets
netstats <- readRDS("data/intermediate/estimates/netstats-level-local.rds")

full_epi %>%
  ggplot(aes(x = time,
             y = n_edges_main,
             color = treatment)) +
  geom_smooth(se = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 1630.057344, linetype = "dotted") +
  geom_hline(yintercept = netstats$main$edges) +
  geom_hline(yintercept = 1172.186555, linetype = "dotted") +
  labs(title = "Number of Main Partnerships",
       x = "Time",
       y = NULL)


full_epi %>%
  ggplot(aes(x = time,
             y = n_edges_casual,
             color = treatment)) +
  geom_smooth(se = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 2386.663237, linetype = "dotted") +
  geom_hline(yintercept = netstats$casl$edges) +
  geom_hline(yintercept = 1591.149932, linetype = "dotted") +
  labs(title = "Number of Casual Partnerships",
       x = "Time",
       y = NULL)

full_epi %>%
  ggplot(aes(x = time,
             y = n_edges_onetime,
             color = treatment)) +
  geom_smooth(se = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 95.26738158, linetype = "dotted") +
  geom_hline(yintercept = netstats$inst$edges) +
  geom_hline(yintercept = 60.6132855832735, linetype = "dotted") +
  labs(title = "Number of One-Time Partnerships",
       x = "Time",
       y = NULL)

# POPULATION SIZE

full_epi %>%
  ggplot(aes(x = time,
             y = num,
             color = treatment)) +
  geom_smooth(se = TRUE) +
  theme_classic() +
  geom_hline(yintercept = length(netstats$attr$numeric.id)) +
  labs(title = "Population Size",
       x = "Time",
       y = NULL)

# PROPORTION BY RACE

full_epi %>%
  mutate(prop_B = n__B/num) %>%
ggplot(aes(x = time,
           y = prop_B,
           color = treatment)) +
geom_smooth(se = TRUE) +
  theme_classic() +
  geom_hline(yintercept = mean(netstats$attr$race == 1)) +
  labs(title = "Proportion of Population Black",
       x = "Time",
       y = "Proportion") +
  scale_y_continuous(breaks = seq(0, .4, .1),
                     limits = c(0, .4))

full_epi %>%
  mutate(prop_H = n__H/num) %>%
  ggplot(aes(x = time,
             y = prop_H,
             color = treatment)) +
  geom_smooth(se = TRUE) +
  theme_classic() +
  geom_hline(yintercept = mean(netstats$attr$race == 2)) +
  labs(title = "Proportion of Population Hispanic",
       x = "Time",
       y = "Proportion") +
  scale_y_continuous(breaks = seq(0, .4, .1),
                     limits = c(0, .4))

full_epi %>%
  mutate(prop_O = n__O/num) %>%
  ggplot(aes(x = time,
             y = prop_O,
             color = treatment)) +
  geom_smooth(se = TRUE) +
  theme_classic() +
  geom_hline(yintercept = mean(netstats$attr$race == 3)) +
  labs(title = "Proportion of Population Other",
       x = "Time",
       y = "Proportion") +
  scale_y_continuous(breaks = seq(0, .4, .1),
                     limits = c(0, .4))

full_epi %>%
  mutate(prop_W = n__W/num) %>%
  ggplot(aes(x = time,
             y = prop_W,
             color = treatment)) +
  geom_smooth(se = TRUE) +
  theme_classic() +
  geom_hline(yintercept = mean(netstats$attr$race == 4)) +
  labs(title = "Proportion of Population White",
       x = "Time",
       y = "Proportion") +
  scale_y_continuous(breaks = seq(0, .4, .1),
                     limits = c(0, .4))

# PROPORTION BY AGE GROUP

full_epi %>%
mutate(prop = num.16to20/num) %>%
ggplot(aes(x = time,
           y = prop,
           color = treatment)) +
geom_smooth(se = TRUE) +
  theme_classic() +
  geom_hline(yintercept = mean(netstats$attr$age.grp == 1)) +
  labs(title = "Proportion of Population Aged 16-20",
       x = "Time",
       y = "Proportion") +
  scale_y_continuous(breaks = seq(0, .4, .1),
                     limits = c(0, .4))


full_epi %>%
mutate(prop = num.21plus/num) %>%
ggplot(aes(x = time,
           y = prop,
           color = treatment)) +
geom_smooth(se = TRUE) +
  theme_classic() +
  geom_hline(yintercept = mean(netstats$attr$age.grp == 2)) +
  labs(title = "Proportion of Population Aged 21-29",
       x = "Time",
       y = "Proportion")

# CONCURRENT TIES

data.frame(y = sim_list[[1]]$stats$nwstats[[1]][, "concurrent"]) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  geom_hline(yintercept = 118.0182859, linetype = "dotted") +
  geom_hline(yintercept = 28.29061478) +
  geom_hline(yintercept = 6.741955174, linetype = "dotted") +
  labs(title = "Number of Concurrent Partnerships (Main)",
       x = "Time",
       y = NULL)

data.frame(y = sim_list[[1]]$stats$nwstats[[2]][, "concurrent"]) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  geom_hline(yintercept = 1055.447181, linetype = "dotted") +
  geom_hline(yintercept = 735.5733031) +
  geom_hline(yintercept = 508.0176154, linetype = "dotted") +
  labs(title = "Number of Concurrent Partnerships (Casual)",
       x = "Time",
       y = NULL)

# Coefficients
full_epi %>%
  ggplot(aes(x = time,
             y = coef_main_nodefactor.age.grp.1,
             color = treatment)) +
  # geom_smooth(se = TRUE) +
  geom_line()+
  theme_classic() +
  labs(title = "Coefficient",
       x = "Time",
       y = NULL)



# NUMBER OF EDGES (AGE GROUP)
### Main


  data.frame(time = full_epi$time,
             under21 = full_epi$n_edges_main.under21,
             over21 = full_epi$n_edges_main.over21,
             sameage = full_epi$n_edges_main.sameage,
             treatment = full_epi$treatment) %>%
  tidyr::pivot_longer(cols = under21:sameage,
                      names_to = "age_group",
                      values_to = "n_edges") %>%
  ggplot(aes(x = time,
             y = n_edges,
             color = age_group,
             linetype = treatment)) +
  geom_smooth() +
  theme_classic() +
  # Over 21 limits
  geom_hline(yintercept = 997.555715, linetype = "dotted", color = scales::hue_pal()(3)[[1]]) +
  geom_hline(yintercept = 763.2791946, color = scales::hue_pal()(3)[[1]]) +
  geom_hline(yintercept = 538.9211228, linetype = "dotted", color = scales::hue_pal()(3)[[1]]) +
  # Same age group limits
  geom_hline(yintercept = 1345.816152, linetype = "dotted", color = scales::hue_pal()(3)[[2]]) +
  geom_hline(yintercept = 979.2628343, color = scales::hue_pal()(3)[[2]]) +
  geom_hline(yintercept = 649.7071525, linetype = "dotted", color = scales::hue_pal()(3)[[2]]) +
  # Under 21 limits
  geom_hline(yintercept = 348.2604366, linetype = "dotted", color = scales::hue_pal()(3)[[3]]) +
  geom_hline(yintercept = 215.9836397, color = scales::hue_pal()(3)[[3]]) +
  geom_hline(yintercept = 110.7860297, linetype = "dotted", color = scales::hue_pal()(3)[[3]]) +
  labs(title = "Number of Main Partnerships by Age Group Composition",
       x = "Time",
       y = "Number of Edges")

### Casual
  data.frame(time = full_epi$time,
             under21 = full_epi$n_edges_casual.under21,
             over21 = full_epi$n_edges_casual.over21,
             sameage = full_epi$n_edges_casual.sameage,
             treatment = full_epi$treatment) %>%
    tidyr::pivot_longer(cols = under21:sameage,
                        names_to = "age_group",
                        values_to = "n_edges") %>%
    ggplot(aes(x = time,
               y = n_edges,
               color = age_group,
               linetype = treatment)) +
    geom_smooth() +
    theme_classic() +
  # Over 21 limits
  geom_hline(yintercept = 1391.18333, linetype = "dotted", color = scales::hue_pal()(3)[[1]]) +
  geom_hline(yintercept = 1036.860988, color = scales::hue_pal()(3)[[1]]) +
  geom_hline(yintercept = 694.8788677, linetype = "dotted", color = scales::hue_pal()(3)[[1]]) +
  # Same age group limits
  geom_hline(yintercept = 1959.644774, linetype = "dotted", color = scales::hue_pal()(3)[[2]]) +
  geom_hline(yintercept = 1383.753222, color = scales::hue_pal()(3)[[2]]) +
  geom_hline(yintercept = 878.2011497, linetype = "dotted", color = scales::hue_pal()(3)[[2]]) +
  # Under 21 limits
  geom_hline(yintercept = 568.4614443, linetype = "dotted", color = scales::hue_pal()(3)[[3]]) +
  geom_hline(yintercept = 346.8922334, color = scales::hue_pal()(3)[[3]]) +
  geom_hline(yintercept = 183.322282, linetype = "dotted", color = scales::hue_pal()(3)[[3]]) +
  labs(title = "Number of Casual Partnerships by Age Group Composition",
       x = "Time",
       y = "Number of Edges")


  data.frame(time = full_epi$time,
             under21 = full_epi$n_edges_onetime.under21,
             over21 = full_epi$n_edges_onetime.over21,
             sameage = full_epi$n_edges_onetime.sameage,
             treatment = full_epi$treatment) %>%
    tidyr::pivot_longer(cols = under21:sameage,
                        names_to = "age_group",
                        values_to = "n_edges") %>%
    ggplot(aes(x = time,
               y = n_edges,
               color = age_group,
               linetype = treatment)) +
    geom_smooth() +
    theme_classic() +
  # Over 21 limits
  geom_hline(yintercept = 54.46798642, linetype = "dotted", color = scales::hue_pal()(3)[[1]]) +
  geom_hline(yintercept = 38.61964212, color = scales::hue_pal()(3)[[1]]) +
  geom_hline(yintercept = 22.91260422, linetype = "dotted", color = scales::hue_pal()(3)[[1]]) +
  # Same age group limits
  geom_hline(yintercept = 82.68117582, linetype = "dotted", color = scales::hue_pal()(3)[[2]]) +
  geom_hline(yintercept = 56.19889793, color = scales::hue_pal()(3)[[2]]) +
  geom_hline(yintercept = 32.46456615, linetype = "dotted", color = scales::hue_pal()(3)[[2]]) +
  # Under 21 limits
  geom_hline(yintercept = 28.2131894, linetype = "dotted", color = scales::hue_pal()(3)[[3]]) +
  geom_hline(yintercept = 17.57925581, color = scales::hue_pal()(3)[[3]]) +
  geom_hline(yintercept = 9.551961931, linetype = "dotted", color = scales::hue_pal()(3)[[3]]) +
  labs(title = "Number of One-Time Partnerships by Age Group Composition",
       x = "Time",
       y = "Number of Edges")


# Prevalence
  data.frame(time = full_epi$time,
             prevB = full_epi$i.prev.B,
             prevH = full_epi$i.prev.H,
             prevO = full_epi$i.prev.O,
             prevW = full_epi$i.prev.W,
             treatment = full_epi$treatment) %>%
    tidyr::pivot_longer(cols = prevB:prevW,
                        names_to = "race",
                        values_to = "prev") %>%
    ggplot(aes(x = time,
               y = prev,
               color = race,
               linetype = treatment)) +
    geom_smooth() +
    theme_classic()

  # Disparities
  data.frame(time = full_epi$time,
             dispB = full_epi$i.prev.B - full_epi$i.prev.W,
             dispH = full_epi$i.prev.H - full_epi$i.prev.W,
             dispO = full_epi$i.prev.O - full_epi$i.prev.W,
             treatment = full_epi$treatment) %>%
    tidyr::pivot_longer(cols = dispB:dispO,
                        names_to = "race",
                        values_to = "disparity") %>%
    ggplot(aes(x = time,
               y = disparity,
               color = race,
               linetype = treatment)) +
    geom_smooth() +
    theme_classic()


  full_epi$time
  full_epi$i.prev.B


for (i in 1:nrow(file_index)) {

  file_info <- file_index[i,]

  apps_treat <- file_info$apps_treat | file_info$appsvenues_treat
  venues_treat <- file_info$venues_treat | file_info$appsvenues_treat

  this_sim <- readRDS(paste(this_dir, file_info$rds, sep = ""))

  sim_list <- as.list(this_sim)

  # Combine pertinent elements of `sim_list` to a single data frame
  full_data <- dplyr::bind_rows(lapply(sim_list$raw.records[[1]], function(x){x$object}))

  full_data2 <- full_data %>%
    dplyr::filter(type != 4) %>%
    dplyr::group_by(head_uid, tail_uid) %>% dplyr::mutate(dyad_id = dplyr::cur_group_id(),
                                                          end = max(current_time)) %>%
    dplyr::rename(start = time) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()


  # full_data2 <- full_data %>%
  #   dplyr::group_by(head_uid, tail_uid) %>% dplyr::mutate(dyad_id = dplyr::cur_group_id()) %>%
  #   dplyr::ungroup() %>%
  #   arrange(dyad_id, time)

  # full_data3 <- full_data2 %>%
  #   dplyr::group_by(dyad_id) %>%
  #   dplyr::summarize(count = n(),
  #                    min_time = min(time),
  #                    max_time = max(time))



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
  edges <- full_data2 %>%
    select(dyad_id, type, head_uid, tail_uid, head_pid, tail_pid, start, end, head_venues, tail_venues, head_apps, tail_apps, current_time)

    # dplyr::filter(type != 4) %>%
    # group_by(head_uid, tail_uid) %>%
    # slice(1)

  # Merge in `inf_events`
  merged_data <- edges %>% dplyr::left_join(inf_events, by = c("head_uid", "tail_uid")) #%>%
    # Keep just last ten years
    # dplyr::filter(time > 3120)

  # What percentage of cases in the infection act list successfully merged?
  # Note that we use `nrow(inf_events2)` here for the accurate denominator

  # merged_data %>% group_by(type) %>%
  #   summarize(sum(infection, na.rm = TRUE))


  # Let's read in the agent log
  head_demo <- read_delim(paste(this_dir, file_info$log, sep = ""),
                          delim = ",") %>%
    select(start = tick, head_uid = agent_id, head_age = age, head_agegroup = age_group,
           head_race = race_ethnicity, head_hiv = hiv_status)
  tail_demo <- head_demo %>%
    rename(tail_uid = head_uid, tail_age = head_age, tail_agegroup = head_agegroup,
           tail_race = head_race, tail_hiv = head_hiv)

  merged_data <- merged_data %>%
    dplyr::left_join(head_demo, by = c("head_uid", "start")) %>%
    dplyr::left_join(tail_demo, by = c("tail_uid", "start"))

  #
  # # Discordant edges at time of partnership formation
  # test <- merged_data %>%
  #   dplyr::mutate(hiv_disc = (head_hiv + tail_hiv == 1)) %>%
  #   dplyr::filter(hiv_disc == TRUE)
  #
  # test %>%
  #   dplyr::group_by(type) %>%
  #   dplyr::summarize(count = n())
  #
  # # Now extract just the partnerships that led to infection
  # inf_edges <- merged_data %>%
  #   filter(infection == 1)
  #
  # # How many by partnership type?
  # inf_edges %>%
  #   dplyr::group_by(type) %>%
  #   dplyr::summarize(count = n())

  # Loop over infection events
  for (j in 1:nrow(merged_data)) {

    # Get this infection event
    this_inf <- merged_data[j,]

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
    # shared_apps <- unlist(stringr::str_split(this_inf$head_apps, "\\|"))[unlist(stringr::str_split(this_inf$head_apps, "\\|")) %in% unlist(stringr::str_split(this_inf$tail_apps, "\\|"))]
    this_inf$shared_apps <- sum(unlist(stringr::str_split(this_inf$head_apps, "\\|")) %in% unlist(stringr::str_split(this_inf$tail_apps, "\\|")))
    this_inf$shared_venues <- sum(unlist(stringr::str_split(this_inf$head_venues, "\\|")) %in% unlist(stringr::str_split(this_inf$tail_venues, "\\|")))

    if (j == 1) {
      va_edges <- this_inf
    } else {
      va_edges <- dplyr::bind_rows(va_edges, this_inf)
    }

  }

  for (k in 3041:max(va_edges$current_time, na.rm = T)) {

    this_time <- va_edges %>% filter(start <= k & k <= end)
    this_main <- this_time %>% filter(type == 1) %>%
      summarize(time = k,
                type = "main",
                edge_count = n(),
                venue_overlap = sum(shared_venues > 0),
                app_overlap = sum(shared_apps > 0),
                prop_venue = venue_overlap/edge_count,
                prop_app = app_overlap/edge_count)
    this_casl <- this_time %>% filter(type == 2) %>%
      summarize(time = k,
                type = "casual",
                edge_count = n(),
                venue_overlap = sum(shared_venues > 0),
                app_overlap = sum(shared_apps > 0),
                prop_venue = venue_overlap/edge_count,
                prop_app = app_overlap/edge_count)
    this_onetime <- this_time %>% filter(type == 3) %>%
      summarize(time = k,
                type = "onetime",
                edge_count = n(),
                venue_overlap = sum(shared_venues > 0),
                app_overlap = sum(shared_apps > 0),
                prop_venue = venue_overlap/edge_count,
                prop_app = app_overlap/edge_count)

    if (k == 3041) {
      diag_main <- this_main
      diag_casl <- this_casl
      diag_onetime <- this_onetime
    } else {
      diag_main <- dplyr::bind_rows(diag_main, this_main)
      diag_casl <- dplyr::bind_rows(diag_casl, this_casl)
      diag_onetime <- dplyr::bind_rows(diag_onetime, this_onetime)
    }

  }

  diag_main$apps_treat <- apps_treat
  diag_casl$apps_treat <- apps_treat
  diag_onetime$apps_treat <- apps_treat

  diag_main$venues_treat <- venues_treat
  diag_casl$venues_treat <- venues_treat
  diag_onetime$venues_treat <- venues_treat

  diag_main$file <- file_info$rds
  diag_casl$file <- file_info$rds
  diag_onetime$file <- file_info$rds

  if (i == 1) {
    full_diag_main <- diag_main
    full_diag_casl <- diag_casl
    full_diag_onetime  <- diag_onetime
  } else {
    full_diag_main <-     dplyr::bind_rows(full_diag_main, diag_main)
    full_diag_casl <-     dplyr::bind_rows(full_diag_casl, diag_casl)
    full_diag_onetime  <- dplyr::bind_rows(full_diag_onetime, diag_onetime)
  }

}



### Visualizations
############# Total Edges
full_diag_main %>%
  ggplot(aes(x = time, y = edge_count, color = file)) +
  geom_line() +
  theme_classic() +
  geom_hline(aes(yintercept = 1390)) +
  labs(title = "Number of Main Partnerships in Network, by Treatment")

full_diag_casl %>%
  ggplot(aes(x = time, y = edge_count, color = file)) +
  geom_line() +
  theme_classic() +
  geom_hline(aes(yintercept = 1974)) +
  labs(title = "Number of Casual Partnerships in Network, by Treatment")

full_diag_onetime %>%
  ggplot(aes(x = time, y = edge_count, color = file)) +
  geom_line() +
  theme_classic() +
  geom_hline(aes(yintercept = 77)) +
  labs(title = "Number of One-time Partnerships in Network, by Treatment")

############# Venue Overlap
full_diag_main %>%
  ggplot(aes(x = time, y = venue_overlap, color = file)) +
  geom_line() +
  theme_classic() +
  geom_hline(aes(yintercept = 210)) +
  labs(title = "Number of Main Partnerships with any Venue Overlap, by Treatment")

full_diag_casl %>%
  ggplot(aes(x = time, y = venue_overlap, color = file)) +
  geom_line() +
  theme_classic() +
  geom_hline(aes(yintercept = 285)) +
  labs(title = "Number of Casual Partnerships with any Venue Overlap, by Treatment")

full_diag_onetime %>%
  ggplot(aes(x = time, y = venue_overlap, color = file)) +
  geom_line() +
  theme_classic() +
  geom_hline(aes(yintercept = 13)) +
  labs(title = "Number of One-time Partnerships with any Venue Overlap, by Treatment")

############# App Overlap
full_diag_main %>%
  ggplot(aes(x = time, y = app_overlap, color = file)) +
  geom_line() +
  theme_classic() +
  geom_hline(aes(yintercept = 779)) +
  labs(title = "Number of Main Partnerships with any App Overlap, by Treatment")

full_diag_casl %>%
  ggplot(aes(x = time, y = app_overlap, color = file)) +
  geom_line() +
  theme_classic() +
  geom_hline(aes(yintercept = 1074)) +
  labs(title = "Number of Casual Partnerships with any App Overlap, by Treatment")

full_diag_onetime %>%
  ggplot(aes(x = time, y = app_overlap, color = file)) +
  geom_line() +
  theme_classic() +
  geom_hline(aes(yintercept = 51)) +
  labs(title = "Number of One-time Partnerships with any App Overlap, by Treatment")


# Adjusted ERGM Coefficients

full_epi %>%
  ggplot(aes(x = time,
             y = coef_main_edges,
             color = treatment)) +
  geom_smooth(se = TRUE) +
  theme_classic() +
  # geom_hline(yintercept = 1630.057344, linetype = "dotted") +
  # geom_hline(yintercept = netstats$main$edges) +
  # geom_hline(yintercept = 1172.186555, linetype = "dotted") +
  labs(title = "Edges Coefficient (Main)",
       x = "Time",
       y = NULL)


full_epi %>%
  ggplot(aes(x = time,
             y = coef_casual_edges,
             color = treatment)) +
  geom_smooth(se = TRUE) +
  theme_classic() +
  # geom_hline(yintercept = 1630.057344, linetype = "dotted") +
  # geom_hline(yintercept = netstats$main$edges) +
  # geom_hline(yintercept = 1172.186555, linetype = "dotted") +
  labs(title = "Edges Coefficient (Casual)",
       x = "Time",
       y = NULL)

full_epi %>%
  ggplot(aes(x = time,
             y = coef_onetime_edges,
             color = treatment)) +
  geom_smooth(se = TRUE) +
  theme_classic() +
  # geom_hline(yintercept = 1630.057344, linetype = "dotted") +
  # geom_hline(yintercept = netstats$main$edges) +
  # geom_hline(yintercept = 1172.186555, linetype = "dotted") +
  labs(title = "Edges Coefficient (Onetime)",
       x = "Time",
       y = NULL)



full_epi %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time,
             y = start_edges_main,
             color = treatment)) +
  geom_smooth(se = TRUE) +
  theme_classic() +
  # geom_hline(yintercept = 1630.057344, linetype = "dotted") +
  # geom_hline(yintercept = netstats$main$edges) +
  # geom_hline(yintercept = 1172.186555, linetype = "dotted") +
  labs(title = "New Edges Created (Main)",
       x = "Time",
       y = NULL)
