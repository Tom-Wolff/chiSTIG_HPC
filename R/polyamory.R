# Identifying "polycules," defined as components of concurrent partnerships whose ties endure
# over a set period of time

file = "~/Desktop/25sept_runs/apps_1_test13sep_1.rds"

library(EpiModelHIV)
library(tidyverse)

sim <- readRDS(file)


component_summaries <- function(file,
                           interval = 4,
                           start_time = 3000) {



  el_cuml <- sim[[1]]$el.cuml
  el_cuml1 <- el_cuml[[1]]
  el_cuml1$type <- 1
  el_cuml2 <- el_cuml[[2]]
  el_cuml2$type <- 2
  onetime_el <- el_cuml[[3]]
  onetime_el$type <- 3

  # el_cuml <- sim$el.cuml[[1]]
  # el_cuml1 <- el_cuml[[1]]
  # el_cuml1$type <- 1
  # el_cuml2 <- el_cuml[[2]]
  # el_cuml2$type <- 2
  # onetime_el <- el_cuml[[3]]
  # onetime_el$type <- 3

  full_el <- dplyr::bind_rows(el_cuml1, el_cuml2)

  # Get full list of timesteps
  timesteps <- sort(unique(c(full_el$start, full_el$stop)))

  for (i in 1:max(timesteps)) {

    # Isolate edgelists within time interval
    this_time_main <- el_cuml1 %>%
      dplyr::filter(start <= i & stop > (i + (interval-1)))

    this_time_casl <- el_cuml2 %>%
      dplyr::filter(start <= i & stop > (i + (interval-1)))

    this_time_onetime <- onetime_el %>%
      dplyr::filter(start >= i & start <= (i + (interval-1)))

    this_time_full <- full_el %>%
      dplyr::filter(start <= i & stop > (i + (interval-1)))

    this_igraph_main <- igraph::graph_from_data_frame(this_time_main, directed = FALSE)
    this_igraph_casl <- igraph::graph_from_data_frame(this_time_casl, directed = FALSE)
    this_igraph_onetime <- igraph::graph_from_data_frame(this_time_onetime, directed = FALSE)
    this_igraph_full <- igraph::graph_from_data_frame(this_time_full, directed = FALSE)

    these_components <- igraph::components(this_igraph_full, mode = "weak")
    these_components_main <- igraph::components(this_igraph_main, mode = "weak")

    component_counts <- as.data.frame(table(these_components$csize[these_components$csize > 2]))
    component_counts_main <- as.data.frame(table(these_components_main$csize[these_components_main$csize > 2]))

    # If no components, skip timestep
    if (nrow(component_counts) == 0) {
      next
    }
    component_counts$time <- i
    colnames(component_counts) <- c("size", "count", "time")
    component_counts <- component_counts %>%
      tidyr::pivot_wider(names_from = "size", names_prefix = "size", values_from = "count")

    # From Levine et al. (2018)

    # MONOGAMY
    main_degree <- data.frame(node_id = names(igraph::degree(this_igraph_main)),
                              main_degree = igraph::degree(this_igraph_main))


    # OPEN RELATIONSHIPS
    ### Couples typically retain emotional intimacy within a primary relationship
    ### and pursue additional casual and/or sexual partnerships
    ### Operationalized as a single, central "main/serious" partnership with pendants of casual/onetime partnerships
    ###### Simplest way to do this, probably, would be to take a look at individual nodes' degree values across the
    ###### three partnership types:
    ############ If ego has a single serious tie but no other ties, we can consider them MONOGAMOUS
    ############ If ego has a single serious tie but additional casual/onetime partnerships,
    ############ we can consider them to be in OPEN RELATIONSHIPS (or cheating)


    casl_degree <- data.frame(node_id = names(igraph::degree(this_igraph_casl)),
                              casl_degree = igraph::degree(this_igraph_casl))

    onetime_degree <- data.frame(node_id = names(igraph::degree(this_igraph_onetime)),
                                 onetime_degree = igraph::degree(this_igraph_onetime))

    total_degree_df <- data.frame(node_id = sort(unique(c(main_degree$node_id,
                                                          casl_degree$node_id,
                                                          onetime_degree$node_id)))) %>%
      dplyr::left_join(main_degree, by = "node_id") %>%
      dplyr::left_join(casl_degree, by = "node_id") %>%
      dplyr::left_join(onetime_degree, by = "node_id")

    total_degree_df[is.na(total_degree_df)] <- 0

    total_degree_df <- total_degree_df %>%
      dplyr::mutate(monogamous = main_degree == 1 & casl_degree == 0 & onetime_degree == 0,
                    open = main_degree == 1 & (casl_degree > 0 | onetime_degree > 0))



    # POLYAMORY
    ### Individuals are open to the possibility of forming loving relationships with multiple partners
    ### Operationalized as components of multiple "main/serious" partnerships

    ##### This tells us how many components in the network represent the presence of polyamory
    component_counts_main$time <- i
    colnames(component_counts_main) <- c("size", "count", "time")
    component_counts_main <- component_counts_main %>%
      tidyr::pivot_wider(names_from = "size", names_prefix = "size", values_from = "count")

    ##### Next, we'll want to count how many nodes are in polyamorous relationships.
    ##### If we have a component of three nodes where one node bridges two others,
    ##### we're observing one node engaged in polyamory while the other two aren't
    these_components_main$membership
    these_components_main$csize
    these_components_main$no

    ######## Which components contain more than two nodes?
    poly_components <- which(these_components_main$csize > 2)
    ######## Which node IDs belong to these components?
    poly_nodes <- these_components_main$membership[these_components_main$membership %in% poly_components]
    poly_nodes <- data.frame(node_id = names(poly_nodes),
                             component = poly_nodes) %>%
      dplyr::arrange(component)

    for (j in 1:length(poly_components)) {
      poly_ids <- poly_nodes %>%
        dplyr::filter(component == poly_components[[j]])
      this_subgraph <- igraph::subgraph(this_igraph_main, vids = poly_ids$node_id)
      this_degree <- data.frame(degree = igraph::degree(this_subgraph),
                              node_id = names(igraph::degree(this_subgraph)))
      if (j == 1) {
        degree_df <- this_degree
      } else {
        degree_df <- dplyr::bind_rows(degree_df, this_degree)
      }
    }

    ######## Merge degree counts back into `poly_nodes`
    poly_nodes <- poly_nodes %>% dplyr::left_join(degree_df, by = "node_id")
    ######## How many nodes are exhibiting polyamory?
    poly_nodes <- poly_nodes %>%
      dplyr::mutate(polyamory = degree > 1)
    actual_poly_nodes <- poly_nodes %>% filter(degree > 1)
    ######## What's the distribution of partner counts among polyamorous nodes?
    poly_dist <- as.data.frame(table(actual_poly_nodes$degree))
    colnames(poly_dist) <- c("partners", "count")
    poly_dist <- poly_dist %>%
      tidyr::pivot_wider(names_from = "partners", values_from = "count", names_prefix = "partners") %>%
      dplyr::mutate(time = i)

    # POLYFIDELITY
    ### Three or more individuals form a closed romantic partnership
    ### Operationalized as a fully dense component of "main/serious" partnerships
    ##### In other words, we're looking for k-cores within the main partnership network,
    ##### specifically 1-cores
    ##### After talking with Jim, identifying cliques is probably a better way to go
    # First thing to do is remove isolated dyads
    # Or maybe the thing to do is take `just_poly_components` and see what cores we get

    just_poly_components <- igraph::delete.vertices(this_igraph_main, v = !(igraph::V(this_igraph_main)$name %in% poly_nodes$node_id))



    just_poly_components <- igraph::graph_from_data_frame(data.frame(ego = c(1, 1, 2, 4, 4, 4, 5, 4, 5, 6, 8, 8, 9, 3, 3, 11),
                                                     alter = c(2, 3, 3, 1, 5, 6, 6, 7, 7, 7, 9, 10, 10, 11, 12, 12)), directed = FALSE)

    clique_list <- igraph::max_cliques(just_poly_components, min = 3)


    if (length(clique_list) > 0) {

        for (k in 1:length(clique_list)) {
          node_clique <- data.frame(node_id = names(clique_list[[k]]),
                                    clique = k,
                                    clique_size = length(names(clique_list[[k]])))

          if (k == 1) {
            clique_mems <- node_clique
          } else {
            clique_mems <- dplyr::bind_rows(clique_mems, node_clique)
          }
        }

      degree_df <- data.frame(node_id = names(igraph::degree(just_poly_components)),
                              degree = igraph::degree(just_poly_components))

      clique_mems <- clique_mems %>% dplyr::left_join(degree_df, by = "node_id") %>%
        mutate(polyfidelity = degree == (clique_size-1))

      polyfidelity_nodes <- clique_mems %>%
                                dplyr::group_by(node_id) %>%
                                dplyr::summarize(num_cliques = n(),
                                                 min_clique_size = min(clique_size),
                                                 max_clique_size = max(clique_size),
                                                 polyfidelity = min(polyfidelity)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(time = i)

      clique_summary <- clique_mems %>%
                            dplyr::summarize(num_cliques = length(unique(clique)),
                                             min_clique_size = min(clique_size),
                                             max_clique_size = max(clique_size),
                                             mean_clique_size = mean(clique_size),
                                             med_clique_size = median(clique_size),
                                             num_in_cliques = length(unique(node_id)),
                                             num_polyfidelity = sum(polyfidelity)) %>%
                            dplyr::mutate(time = i)


    } else {

      polyfidelity_nodes <- data.frame()

      clique_summary <- data.frame(time = i,
                                   num_cliques = 0,
                                   min_clique_size =  0,
                                   max_clique_size =  0,
                                   mean_clique_size = 0,
                                   med_clique_size =  0,
                                   num_in_cliques = 0,
                                   num_polyfidelity = 0)

    }



    # SWINGING
    ### Couples pursue extradyadic sex
    ### Operationalized as four cycles where two "main/serious" partnerships are connected by 2-3 casual/onetime partnerships


    if (i == 1) {
      component_df <- component_counts
    } else {
      component_df <- dplyr::bind_rows(component_df, component_counts)
    }

  }




  component_df[is.na(component_df)] <- 0

  # Read in file
  this_file <- readRDS(file)

  raw <- this_file$raw.records[[1]]

  full_el <- data.frame()

  for (i in 1:length(raw)) {

    if (raw[[i]]$label == "infection_act") {
      next
    } else {
      full_el <- dplyr::bind_rows(full_el, raw[[i]]$object)
    }
  }

  main_cas <- full_el %>%
    filter(type != 3) %>%
    group_by(head_uid, tail_uid) %>%
    summarize(count = n(),
              start = min(time),
              stop = max(time)) %>%
    ungroup() %>%
    mutate(duration = stop - start)


  full_edgelist <- this_file$el.cuml
  edgelist1 <- full_edgelist[[1]][[1]]

  hm <- edgelist1 %>%
    group_by(head, tail) %>%
    summarize(start = min(start, na.rm = T),
              stop = max(start, na.rm = T),
              duration = stop - start)

  edgelist2 <- full_edgelist[[1]][[2]]

  edgelist <- this_file$el[[1]][[1]]



}
