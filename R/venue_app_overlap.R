agent_log <- read.delim("~/Desktop/sept3_runs/agent_log_a1_aug30plus.txt", sep = ",")

timesteps <- unique(agent_log$tick)

for (i in 3601:length(timesteps)) {

  print(i)

  # Get all nodes for this particular time
  this_time <- agent_log %>% filter(tick == timesteps[[i]])

  # What's the population size?
  num_nodes <- nrow(this_time)
  # How many possible dyads exist?
  num_potential <- num_nodes*(num_nodes-1)/2

  # Split venues attended
  venue_vec <- stringr::str_split(this_time$venues_attended, "\\|")

  # Pivot into a two-mode adjacency matrix
  venue_el <- data.frame(ego = rep(this_time$ego_id, unlist(lapply(venue_vec, length))),
                         venue = unlist(venue_vec)) %>%
    filter(!stringr::str_detect(venue, "^s")) %>%
    mutate(attend = 1) %>%
    tidyr::pivot_wider(names_from = "venue",
                       values_from = "attend",
                       values_fill = 0)

  # Convert to matrix object
  venue_mat <- as.matrix(venue_el[,2:ncol(venue_el)])

  # Multiple two-mode matrix by its transpose to get person-level projection matrix
  venue_personmode <- venue_mat %*% t(venue_mat)
  # Set diagnonal to zero
  diag(venue_personmode) <- 0
  # We're only interested in whether ANY overlap exists, so values above 1 should be set to 1
  venue_personmode <- venue_personmode > 0

  # How many dyads have any overlap?
  venue_overlap <- sum(venue_personmode)/2

  # What proportion of potential dyads have overlap in venues?
  prop_venue <- venue_overlap/num_potential

  ##### Same thing for apps
  # Split apps attended
  app_vec <- stringr::str_split(this_time$apps_used, "\\|")

  # Pivot into a two-mode adjacency matrix
  app_el <- data.frame(ego = rep(this_time$ego_id, unlist(lapply(app_vec, length))),
                         app = unlist(app_vec)) %>%
    filter(!stringr::str_detect(app, "^s")) %>%
    mutate(attend = 1) %>%
    tidyr::pivot_wider(names_from = "app",
                       values_from = "attend",
                       values_fill = 0)

  # Convert to matrix object
  app_mat <- as.matrix(app_el[,2:ncol(app_el)])

  # Multiple two-mode matrix by its transpose to get person-level projection matrix
  app_personmode <- app_mat %*% t(app_mat)
  # Set diagnonal to zero
  diag(app_personmode) <- 0
  # We're only interested in whether ANY overlap exists, so values above 1 should be set to 1
  app_personmode <- app_personmode > 0

  # How many dyads have any overlap?
  app_overlap <- sum(app_personmode)/2

  # What proportion of potential dyads have overlap in venues?
  prop_app <- app_overlap/num_potential

  if (i == 1) {
    va_summary <- data.frame(time = timesteps[[i]],
                             num_nodes = num_nodes,
                             num_potential = num_potential,
                             venue_overlap = venue_overlap,
                             prop_venue = prop_venue,
                             app_overlap = app_overlap,
                             prop_app = prop_app)
  } else {
    this_summary <- data.frame(time = timesteps[[i]],
                               num_nodes = num_nodes,
                               num_potential = num_potential,
                               venue_overlap = venue_overlap,
                               prop_venue = prop_venue,
                               app_overlap = app_overlap,
                               prop_app = prop_app)

    va_summary <- dplyr::bind_rows(va_summary, this_summary)
  }

}
