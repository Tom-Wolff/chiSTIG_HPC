##
## 11. Epidemic Model Parameter Calibration, Processing of the simulation files
##

# Read in data

# Venues Model
#
# sim_files <- list.files(paste(getwd(), "/data/intermediate/calibration", sep = ""))
# sim_dir <- paste(getwd(), "/data/intermediate/", sep = "")
# sim_files <- paste(sim_dir, sim_files, sep = "")
#
#
# for (i in 1:length(sim_files)) {
#
#   this_sim <- readRDS(sim_files[[i]])
#
#   this_sim <- this_sim %>% as_tibble() %>%
#     mutate_calibration_targets() %>%
#     mutate(dep.total = dep.gen + dep.AIDS.off.tx + dep.AIDS.on.tx) %>%
#     mutate(total.incid = incid + exo.incid) %>%
#     mutate(prop_unprotected_discordant = num_unprotected_disc_acts/num_disc_acts) %>%
#     mutate(cc.dx.B = ifelse(is.nan(cc.dx.B), 0, cc.dx.B),
#            cc.dx.H = ifelse(is.nan(cc.dx.H), 0, cc.dx.H),
#            cc.dx.O = ifelse(is.nan(cc.dx.O), 0, cc.dx.O),
#            cc.dx.W = ifelse(is.nan(cc.dx.W), 0, cc.dx.W),
#
#            num_diagnosed.B = cc.dx.B*i.num.B,
#            num_diagnosed.H = cc.dx.H*i.num.H,
#            num_diagnosed.O = cc.dx.O*i.num.O,
#            num_diagnosed.W = cc.dx.W*i.num.W,
#            num_diagnosed = num_diagnosed.B + num_diagnosed.H + num_diagnosed.O
#            + num_diagnosed.W,
#            cc.dx = num_diagnosed/i.num) %>%
#     mutate(cc.vsupp.B = ifelse(is.nan(cc.vsupp.B), 0, cc.vsupp.B),
#            cc.vsupp.H = ifelse(is.nan(cc.vsupp.H), 0, cc.vsupp.H),
#            cc.vsupp.O = ifelse(is.nan(cc.vsupp.O), 0, cc.vsupp.O),
#            cc.vsupp.W = ifelse(is.nan(cc.vsupp.W), 0, cc.vsupp.W),
#
#            num_supp.B = num_diagnosed.B*cc.vsupp.B,
#            num_supp.H = num_diagnosed.H*cc.vsupp.H,
#            num_supp.O = num_diagnosed.O*cc.vsupp.O,
#            num_supp.W = num_diagnosed.W*cc.vsupp.W,
#            num_supp = num_supp.B + num_supp.H + num_supp.O
#            + num_supp.W,
#            cc.vsupp = num_supp/num_diagnosed) %>%
#     mutate(n_edges_total = n_edges_main + n_edges_casual + n_edges_onetime) %>%
#     mutate(n_edges_main.RHom = n_edges_main.B + n_edges_main.H + n_edges_main.O + n_edges_main.W,
#            n_edges_casual.RHom = n_edges_casual.B + n_edges_casual.H + n_edges_casual.O + n_edges_casual.W,
#            n_edges_onetime.RHom = n_edges_onetime.B + n_edges_onetime.H + n_edges_onetime.O + n_edges_onetime.W,
#            n_edges_total.RHom = n_edges_total.B + n_edges_total.H + n_edges_total.O + n_edges_total.W,
#            n_edges_total = n_edges_main + n_edges_casual + n_edges_onetime,
#
#            prop_edges_main.RHom = n_edges_main.RHom/n_edges_main,
#            prop_edges_casual.RHom = n_edges_casual.RHom/n_edges_casual,
#            prop_edges_onetime.RHom = n_edges_onetime.RHom/n_edges_onetime,
#            prop_edges_total.RHom = n_edges_total.RHom/n_edges_total) %>%
#     mutate(gen_dep_rate = dep.gen / num,
#            AIDStx_dep_rate = dep.AIDS.on.tx/num,
#            AIDSnotx_dep_rate = dep.AIDS.off.tx/num) %>%
#     mutate(sim = sim_files[[i]], venues_treat = TRUE) %>%
#     select(sim, venues_treat, everything())
#
#   if (i == 1) {
#     sim_data <- this_sim
#   } else {
#     sim_data <- dplyr::bind_rows(sim_data, this_sim)
#   }
#
# }
#
# sim_data <- sim_data %>%
#   mutate(treat = case_when(stringr::str_detect(sim_data$sim, "/app_") == TRUE ~ "Apps Only",
#                            stringr::str_detect(sim_data$sim, "/venue_sim") == TRUE ~ "Venues Only",
#                            stringr::str_detect(sim_data$sim, "/venue_app_") == TRUE ~ "Venues and Apps",
#                            TRUE ~ NA)) %>%
#   dplyr::select(sim, treat, everything())
#
# sim_data1 <- sim_data
#
# # Adrien says to just look at the last 52 weeks
# # sim_data <- sim_data %>%
# #   filter(time > (max(time)-52))
#
# saveRDS(sim_data, paste(sim_dir, "venues_runs_oct24.rds", sep = ""))
#
#
#
# # Control Model
#
# sim_files <- list.files(paste(getwd(), "/calib1/control/", sep = ""))
# sim_dir <- paste(getwd(), "/calib1/control/", sep = "")
# sim_files <- paste(sim_dir, sim_files, sep = "")
#
# for (i in 1:length(sim_files)) {
#
#   this_sim <- readRDS(sim_files[[i]])
#
#   this_sim <- this_sim %>% as_tibble() %>%
#     mutate_calibration_targets() %>%
#     mutate(dep.total = dep.gen + dep.AIDS.off.tx + dep.AIDS.on.tx) %>%
#     mutate(total.incid = incid + exo.incid) %>%
#     mutate(prop_unprotected_discordant = num_unprotected_disc_acts/num_disc_acts) %>%
#     mutate(cc.dx.B = ifelse(is.nan(cc.dx.B), 0, cc.dx.B),
#            cc.dx.H = ifelse(is.nan(cc.dx.H), 0, cc.dx.H),
#            cc.dx.O = ifelse(is.nan(cc.dx.O), 0, cc.dx.O),
#            cc.dx.W = ifelse(is.nan(cc.dx.W), 0, cc.dx.W),
#
#            num_diagnosed.B = cc.dx.B*i.num.B,
#            num_diagnosed.H = cc.dx.H*i.num.H,
#            num_diagnosed.O = cc.dx.O*i.num.O,
#            num_diagnosed.W = cc.dx.W*i.num.W,
#            num_diagnosed = num_diagnosed.B + num_diagnosed.H + num_diagnosed.O
#            + num_diagnosed.W,
#            cc.dx = num_diagnosed/i.num) %>%
#     mutate(cc.vsupp.B = ifelse(is.nan(cc.vsupp.B), 0, cc.vsupp.B),
#            cc.vsupp.H = ifelse(is.nan(cc.vsupp.H), 0, cc.vsupp.H),
#            cc.vsupp.O = ifelse(is.nan(cc.vsupp.O), 0, cc.vsupp.O),
#            cc.vsupp.W = ifelse(is.nan(cc.vsupp.W), 0, cc.vsupp.W),
#
#            num_supp.B = num_diagnosed.B*cc.vsupp.B,
#            num_supp.H = num_diagnosed.H*cc.vsupp.H,
#            num_supp.O = num_diagnosed.O*cc.vsupp.O,
#            num_supp.W = num_diagnosed.W*cc.vsupp.W,
#            num_supp = num_supp.B + num_supp.H + num_supp.O
#            + num_supp.W,
#            cc.vsupp = num_supp/num_diagnosed) %>%
#     mutate(n_edges_total = n_edges_main + n_edges_casual + n_edges_onetime) %>%
#     mutate(n_edges_main.RHom = n_edges_main.B + n_edges_main.H + n_edges_main.O + n_edges_main.W,
#            n_edges_casual.RHom = n_edges_casual.B + n_edges_casual.H + n_edges_casual.O + n_edges_casual.W,
#            n_edges_onetime.RHom = n_edges_onetime.B + n_edges_onetime.H + n_edges_onetime.O + n_edges_onetime.W,
#            n_edges_total.RHom = n_edges_total.B + n_edges_total.H + n_edges_total.O + n_edges_total.W,
#            n_edges_total = n_edges_main + n_edges_casual + n_edges_onetime,
#
#            prop_edges_main.RHom = n_edges_main.RHom/n_edges_main,
#            prop_edges_casual.RHom = n_edges_casual.RHom/n_edges_casual,
#            prop_edges_onetime.RHom = n_edges_onetime.RHom/n_edges_onetime,
#            prop_edges_total.RHom = n_edges_total.RHom/n_edges_total) %>%
#     mutate(gen_dep_rate = dep.gen / num,
#            AIDStx_dep_rate = dep.AIDS.on.tx/num,
#            AIDSnotx_dep_rate = dep.AIDS.off.tx/num) %>%
#     mutate(sim = sim_files[[i]], treat = "Control Model", venues_treat = FALSE) %>%
#     select(sim, treat, venues_treat, everything())
#
#   if (i == 1) {
#     sim_data <- this_sim
#   } else {
#     sim_data <- dplyr::bind_rows(sim_data, this_sim)
#   }
#
# }
#
# sim_data2 <- sim_data


# sim_data <- sim_data %>%
#   filter(time > (max(time)-52))

saveRDS(sim_data, paste(sim_dir, "novenues_runs_oct23.rds", sep = ""))


sim_targets <- tibble::as_tibble(sim) %>%
  mutate_calibration_targets() %>%
  mutate(cc.dx.B = ifelse(is.nan(cc.dx.B), 0, cc.dx.B),
         cc.dx.H = ifelse(is.nan(cc.dx.H), 0, cc.dx.H),
         cc.dx.O = ifelse(is.nan(cc.dx.O), 0, cc.dx.O),
         cc.dx.W = ifelse(is.nan(cc.dx.W), 0, cc.dx.W),

         num_diagnosed.B = cc.dx.B*i.num.B,
         num_diagnosed.H = cc.dx.H*i.num.H,
         num_diagnosed.O = cc.dx.O*i.num.O,
         num_diagnosed.W = cc.dx.W*i.num.W,
         num_diagnosed = num_diagnosed.B + num_diagnosed.H + num_diagnosed.O
         + num_diagnosed.W,
         cc.dx = num_diagnosed/i.num)

# Libraries --------------------------------------------------------------------
library("tidyverse")
library("future.apply")
library("EpiModelHIV")

# data_oct24 <- readRDS("~/Desktop/chistig_debug/all_data_oct24.rds")


# Custom function for generating summary plots
target_plot <- function(data, group, var, benchmark = NULL, title = NULL) {

  # Create placeholder of variable we need
  data2 <- data[, c("time", group, var)]
  colnames(data2) <- c("time", "venues_treat", "this_var")

  this_plot <- data2 %>%
    group_by(venues_treat, time) %>%
    summarize(this_var = median(this_var)) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x = time, y = this_var, color = venues_treat)) +
    geom_line() +
    theme_minimal() +
    labs(y = var)

  if (!is.null(benchmark)) {
    this_plot <- this_plot +
      geom_hline(aes(yintercept = benchmark), color = "red")
  }

  if (!is.null(title)) {
    this_plot <- this_plot +
      ggtitle(title)
  }

  return(this_plot)
}

# Settings ---------------------------------------------------------------------
source("./R/utils-0_project_settings.R")
context <- if (!exists("context")) "local" else "hpc"

if (context == "local") {
  plan(sequential)
} else if (context == "hpc") {
  plan(multisession, workers = ncores)
} else  {
  stop("The `context` variable must be set to either 'local' or 'hpc'")
}

# ------------------------------------------------------------------------------
# Necessary files
source("R/utils-chistig_basic_inputs.R") # generate `path_to_est`, `param` and `init`
path_to_est <- "./data/intermediate/estimates/basic_netest-local.rds"
source("./R/utils-targets.R")
batches_infos <- EpiModelHPC::get_scenarios_batches_infos(calib_dir)

# How many simulations are we comparing
num_batches <- 2
# Later on we'll use `list.files` to automate this:

# file_vec <- c("./data/intermediate/calibration/sim__1__1.rds",
#               "./data/intermediate/calibration/sim__2__1.rds",
#              "./data/intermediate/calibration/sim__3__1.rds",
#                "./data/intermediate/calibration/sim__4__1.rds",
#                "./data/intermediate/calibration/sim__5__1.rds"
#               )

file_vec <- batches_infos$file_name[1:num_batches]


for (i in 1:length(file_vec)) {
  this_sim <- readRDS(file_vec[[i]]) %>%
    as_tibble() %>%
    mutate_calibration_targets() %>%
    mutate(cc.dx.B = ifelse(is.nan(cc.dx.B), 0, cc.dx.B),
           cc.dx.H = ifelse(is.nan(cc.dx.H), 0, cc.dx.H),
           cc.dx.O = ifelse(is.nan(cc.dx.O), 0, cc.dx.O),
           cc.dx.W = ifelse(is.nan(cc.dx.W), 0, cc.dx.W),

           num_diagnosed.B = cc.dx.B*i.num.B,
           num_diagnosed.H = cc.dx.H*i.num.H,
           num_diagnosed.O = cc.dx.O*i.num.O,
           num_diagnosed.W = cc.dx.W*i.num.W,
           num_diagnosed = num_diagnosed.B + num_diagnosed.H + num_diagnosed.O
           + num_diagnosed.W,
           cc.dx = num_diagnosed/i.num,
           sim = file_vec[[i]])

  if (i == 1) {
    sim_targets <- this_sim
  } else {
    sim_targets <- bind_rows(sim_targets, this_sim)
  }

}


# Population Size

starting_n <- length(netstats$attr$age)

sim_targets %>% ggplot(aes(x = time, y = num, color = as.factor(sim))) +
  geom_line(alpha = .4) +
  ggtitle("Plot 1: Population Size") +
  theme_minimal() +
  geom_hline(aes(yintercept = 11612, color = "red"))


# Sam says we should be checking in this order:
# 1. Prep Usage (as this indirectly affects diagnosis given that PrEP users
# often get tested 4 times of year as a matter of course)
##### Right now our simulation doesn't have any PrEP use whatsoever. Do we want
##### to worry about this at present?



################################################################################
# 2. Diagnosis - Knowing one's HIV status affects whether or not they adopt ART,
# so we need to meet the diagnosis benchmarks BEFORE trying to calibrate ART and
# viral suppression

# Prop Diagnosed

sim_targets %>%
  ggplot(aes(x = time, y = cc.dx, color = as.factor(sim))) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 0.814, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 1: Proportion of HIV+ that are Diagnosed")


### Black
sim_targets %>%
  ggplot(aes(x = time, y = cc.dx.B, as.factor(sim))) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 0.804, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 2: Proportion of HIV+ that are Diagnosed (Black)")


### Hispanic
sim_targets %>%
  ggplot(aes(x = time, y = cc.dx.H, as.factor(sim))) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 0.799, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 3: Proportion of HIV+ that are Diagnosed (Hispanic)")


sim_targets %>%
  ggplot(aes(x = time, y = cc.dx.O, color = as.factor(sim))) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 0.826, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 4: Proportion of HIV+ that are Diagnosed (Other)")


sim_targets %>%
  ggplot(aes(x = time, y = cc.dx.W, as.factor(sim))) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 0.881, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 5: Proportion of HIV+ that are Diagnosed (White)")


################################################################################
# 3. Treatment (and Viral Suppression)
# Sam says that for the purposes of simplicity, we assume all folks on ART are
# in full viral suppression category. As far as calibration goes, we keep disengagement
# from treatment rates fixed and calibrate on reinitiation.
###### First we want proportions linked to care, then we do reinitiation
sim_targets %>% ggplot(aes(x = time, y = cc.linked1m.B, as.factor(sim))) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = .828, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 6: Proportion of HIV+ Nodes Linked to Care within One Month (Black)")


sim_targets %>% ggplot(aes(x = time, y = cc.linked1m.H, as.factor(sim))) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 0.867, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 6: Proportion of HIV+ Nodes Linked to Care within One Month (Hispanic)")

sim_targets %>% ggplot(aes(x = time, y = cc.linked1m.O, as.factor(sim))) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 0.875, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 6: Proportion of HIV+ Nodes Linked to Care within One Month (Other)")

sim_targets %>% ggplot(aes(x = time, y = cc.linked1m.W, as.factor(sim))) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 0.936, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 6: Proportion of HIV+ Nodes Linked to Care within One Month (White)")

########## Proportion viral suppressed
sim_targets %>% ggplot(aes(x = time, y = cc.vsupp.B, as.factor(sim))) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 0.571, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 6: Proportion of HIV+ Nodes with Viral Suppression (Black)")
sim_targets$cc.vsupp.H
sim_targets %>% ggplot(aes(x = time, y = cc.vsupp.H, color = as.factor(sim))) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 0.675, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 7: Proportion of HIV+ Nodes with Viral Suppression\n(Hispanic)")
 sim_targets$cc.vsupp.O
sim_targets %>% ggplot(aes(x = time, y = cc.vsupp.O, color = as.factor(sim))) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 0.586, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 8: Proportion of HIV+ Nodes with Viral Suppression\n(Other)")
sim_targets$cc.vsupp.W
sim_targets %>% ggplot(aes(x = time, y = cc.vsupp.H, color = as.factor(sim))) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 0.617, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 9: Proportion of HIV+ Nodes with Viral Suppression\n(White)")

################################################################################
# Exogenous Infections


mean_incid <-  sim_targets %>%
  mutate(total.incid.B = incid.B + exo.incid.B,
         total.incid.H = incid.H + exo.incid.H,
         total.incid.O = incid.O + exo.incid.O,
         total.incid.W = incid.W + exo.incid.W) %>%
  group_by(sim, time) %>%
  summarize(incid.B = mean(incid.B),
            exo.incid.B = mean(exo.incid.B),
            total.incid.B = mean(total.incid.B),

            incid.H = mean(incid.H),
            exo.incid.H = mean(exo.incid.H),
            total.incid.H = mean(total.incid.H),

            incid.O = mean(incid.O),
            exo.incid.O = mean(exo.incid.O),
            total.incid.O = mean(total.incid.O),

            incid.W = mean(incid.W),
            exo.incid.W = mean(exo.incid.W),
            total.incid.W = mean(total.incid.W)) %>%
  ungroup()


for (i in 1:nrow(mean_incid)) {
  this_row <- mean_incid[i,]
  past_year <- mean_incid %>%
    filter(time <= this_row$time & time > (this_row$time-52)) %>%
    filter(sim == this_row$sim)
  sums <- as.data.frame(t(colSums(past_year[,3:ncol(past_year)])))
  sums$sim <- this_row$sim
  sums$time <- this_row$time

  if (i == 1) {
    annual_incid <- sums
  } else {
    annual_incid <- dplyr::bind_rows(annual_incid, sums)
  }
}

target_plot(data = annual_incid,
            var = "exo.incid.B",
            group = "sim",
            benchmark = 32.57621,
            title = "Annual Exogenous Infections (Black)")

target_plot(data = annual_incid,
            var = "exo.incid.H",
            group = "sim",
            benchmark = .40*64.27536,
            title = "Annual Exogenous Infections (Hispanic)")

target_plot(data = annual_incid,
            var = "exo.incid.O",
            group = "sim",
            benchmark = .37*19.59955,
            title = "Annual Exogenous Infections (Other)")


target_plot(data = annual_incid,
            var = "exo.incid.W",
            group = "sim",
            benchmark = .37*19.59955,
            title = "Annual Exogenous Infections (White)")

### Endogenous Incidence

target_plot(data = annual_incid,
            var = "incid.B",
            group = "sim",
            benchmark = 116.3436*(1-.28),
            title = "Annual Endogenous Infections (Black)")

target_plot(data = annual_incid,
            var = "incid.H",
            group = "sim",
            benchmark = 64.27536*(1-.40),
            title = "Annual Endogenous Infections (Hispanic)")

target_plot(data = annual_incid,
            var = "incid.O",
            group = "sim",
            benchmark = 19.59955*(1-.37),
            title = "Annual Endogenous Infections (Other)")

target_plot(data = annual_incid,
            var = "incid.W",
            group = "sim",
            benchmark = 24.91874*(1-.44),
            title = "Annual Endogenous Infections (White)")


################################################################################
# 4. Incidence and Prevalence

# Incidence
sim_targets %>% ggplot(aes(x = time, y = i.prev, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = mean(netstats$attr$diag.status), color = "red")) +
  theme_minimal() +
  ggtitle("Plot 1: Prevalence, Infected")


sim_targets %>% ggplot(aes(x = time, y = i.prev.B, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = .354, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 2: Prevalence, Infected (Black)")

sim_targets %>% ggplot(aes(x = time, y = i.prev.H, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = .114, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 3: Prevalence, Infected (Hispanic)")

sim_targets %>% ggplot(aes(x = time, y = i.prev.O, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = .103, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 4: Prevalence, Infected (Other)")

sim_targets %>% ggplot(aes(x = time, y = i.prev.W, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = .106, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 2: Prevalence, Infected (White)")

sim_targets %>% ggplot(aes(x = time, y = incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 5.08, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 3: Weekly Incidence (Endogenous)")

sim_targets %>% dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(incid, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 5.08, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 3: Weekly Incidence (Endogenous)")

sim_targets %>% dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(exo.incid, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 5.08, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 3: Weekly Incidence (Exogenous)")

sim_targets %>%
  mutate(total.incid = incid + exo.incid) %>%
  dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(total.incid, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 5.08, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 3: Weekly Incidence (Combined)")


sim_targets %>%
  mutate(total.incid = incid + exo.incid) %>%
  summarize(mean(total.incid, na.rm = T))

# Black
sim_targets %>% dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(incid.B, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 3.02, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 3: Weekly Incidence (Black, Endogenous)")

sim_targets %>% dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(exo.incid.B, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 3.02, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 3: Weekly Incidence (Black, Exogenous)")

sim_targets %>%
  mutate(total.incid = incid.B + exo.incid.B) %>%
  dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(total.incid, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 3.02, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 6: Weekly Incidence (Black, Combined)")


# Hispanic
sim_targets %>% dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(incid.H, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 1.15, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 3: Weekly Incidence (Hispanic, Endogenous)")

sim_targets %>% dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(exo.incid.H, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 1.15, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 3: Weekly Incidence (Black, Exogenous)")

sim_targets %>%
  mutate(total.incid = incid.H + exo.incid.H) %>%
  dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(total.incid, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 1.15, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 7: Weekly Incidence (Hispanic, Combined)")


# Other
sim_targets %>% dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(incid.O, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = .31, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 3: Weekly Incidence (Other, Endogenous)")

sim_targets %>% dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(exo.incid.O, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = .31, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 3: Weekly Incidence (Other, Exogenous)")

sim_targets %>%
  mutate(total.incid = incid.O + exo.incid.O) %>%
  dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(total.incid, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = .31, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 8: Weekly Incidence (Other, Combined)")

# White
sim_targets %>% dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(incid.W, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = .6, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 3: Weekly Incidence (White, Endogenous)")

sim_targets %>% dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(exo.incid.W, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = .6, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 3: Weekly Incidence (White, Exogenous)")

sim_targets %>%
  mutate(total.incid = incid.W + exo.incid.W) %>%
  dplyr::group_by(sim, time) %>%
  dplyr::summarise(med_incid = mean(total.incid, na.rm = TRUE)) %>%
  ggplot(aes(x = time, y = med_incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = .6, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 9: Weekly Incidence (White, Combined)")



sim_targets %>% ggplot(aes(x = time, y = exo.incid, color = sim)) +
  geom_line(alpha = .4) +
  # geom_hline(aes(yintercept = 5.08, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 4: Weekly Incidence (Exogenous)")

sim_targets %>%
  mutate(total.incid = incid + exo.incid) %>%
  ggplot(aes(x = time, y = total.incid, color = sim)) +
  geom_line(alpha = .4) +
  geom_hline(aes(yintercept = 5.08, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 5: Weekly Incidence (Combined)")


# Percentage of infections that are exogenous
data_5 <- data_5 %>%
  mutate(prop.exo.incid.B = exo.incid.B/(incid.B + exo.incid.B),
         prop.exo.incid.H = exo.incid.H/(incid.H + exo.incid.H),
         prop.exo.incid.O = exo.incid.O/(incid.O + exo.incid.O),
         prop.exo.incid.W = exo.incid.W/(incid.W + exo.incid.W))

mean(data_5$prop.exo.incid.B, na.rm = T)
mean(data_5$prop.exo.incid.H, na.rm = T)
mean(data_5$prop.exo.incid.O, na.rm = T)
mean(data_5$prop.exo.incid.W, na.rm = T)


data.frame(race = netstats$attr$race,
           hiv.status = netstats$attr$diag.status) %>%
  dplyr::group_by(race) %>%
  dplyr::summarize(prev = mean(hiv.status)) %>%
  dplyr::ungroup()





# 1. Linked to care









