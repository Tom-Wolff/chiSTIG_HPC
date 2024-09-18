## Example interactive epidemic simulation run script with basic
## parameterization and all parameters defined in `param_msm`, with example of
## writing/debugging modules

# Libraries  -------------------------------------------------------------------
library("EpiModelHIV")
# library("chiSTIGmodules")

# library("chiSTIGmodules")
devtools::load_all("~/Desktop/chiSTIGmodules")

# Settings ---------------------------------------------------------------------
source("R/utils-0_project_settings.R")
source("./R/utils-targets.R")

# Necessary files
epistats <- readRDS("data/intermediate/estimates/epistats-local.rds")
epistats$age.breaks <- c(16, 21, 30)
# age limits probably needs to be maxed at 31
epistats$age.limits <- c(16, 30)

netstats_files <- list.files("data/intermediate/estimates/")
netstats_files <- netstats_files[stringr::str_detect(netstats_files, "^netstats")]
netstats_files <- netstats_files[stringr::str_detect(netstats_files, "000")]

netest_files <- list.files("data/intermediate/estimates/")
netest_files <- netest_files[stringr::str_detect(netest_files, "^basic_netest")]
netest_files <- netest_files[stringr::str_detect(netest_files, "000")]

sim_run_list <- list()

for (j in 5:9) {

  print(j)

# netstats <- readRDS("data/intermediate/estimates/netstats-local.rds")
# # `netstats` without venues target stats, in case this changes anything
# netstats <- readRDS("data/intermediate/estimates/netstats-novenues-local.rds")

# netstats <- readRDS("data/intermediate/estimates/netstats-level-local.rds")
netstats <- readRDS(paste0("data/intermediate/estimates/", netstats_files[[j]]))



##### `est` object for basic ERGMs (no venues or apps)
# est      <- readRDS("data/intermediate/estimates/basic_netest-local.rds")
##### `est` object for full ERGMs (venues and apps)
# est      <- readRDS("data/intermediate/estimates/full_netest-local.rds")
##### `est` object for test ERGMs
# est      <- readRDS("data/intermediate/estimates/test_netest-local.rds")

# # Necessary files for original template workflow (saved for reference)
# epistats <- readRDS("data/intermediate/estimates/template_epistats-local.rds")
# netstats <- readRDS("data/intermediate/estimates/template_netstats-local.rds")
# #### `est` object for template ERMGs
# est <- readRDS("data/intermediate/estimates/template_netest-local.rds")
# est <- readRDS("data/intermediate/estimates/basic_netest-level-local.rds")
est <- readRDS(paste0("data/intermediate/estimates/", netest_files[[j]]))


age_race <- as.matrix(table(netstats$attr$race, netstats$attr$age.grp))
age_race

# Expected age group counts for level age distribution
level_vals <- data.frame(under21 = round(rowSums(age_race)*(5/14)),
                         over21 = round(rowSums(age_race)*(9/14)))

for (i in 1:4) {

  race_positions <- which(netstats$attr$race == i)
  under21_pos <- which(netstats$attr$race == i & netstats$attr$age < 21)
  over21_pos <- which(netstats$attr$race == i & netstats$attr$age >= 21)

  if (length(under21_pos) < level_vals[i,1]) {
    netstats$attr$age[sample(over21_pos, size = (length(over21_pos) - level_vals[i,2]))] <- sample(16:20, size = ((length(over21_pos) - level_vals[i,2])), replace = TRUE) + sample(1:1000,  size = ((length(over21_pos) - level_vals[i,2])), replace = TRUE)/1000
  } else if (length(under21_pos) > level_vals[i,1]) {
    netstats$attr$age[sample(under21_pos, size = (length(under21_pos) - level_vals[i,1]))]<- sample(21:29, size = ((length(under21_pos) - level_vals[i,1])), replace = TRUE) + sample(1:1000,  size = ((length(under21_pos) - level_vals[i,1])), replace = TRUE)/1000
  } else {
    next
  }
}

netstats$attr$age.grp <- ifelse(netstats$attr$age < 21, 1, 2)

table(netstats$attr$race, netstats$attr$age.grp)


level_vals

# Is the aging out of the older initial nodes driving HIV extinction?
# Let's level out the age distribution and find out
netstats$attr$age <- sample(16:29, length(netstats$attr$age), replace = TRUE)
netstats$attr$age <- netstats$attr$age + sample(1:1000, length(netstats$attr$age), replace = TRUE)/1000
netstats$attr$age.grp <- ifelse(netstats$attr$age < 21, 1, 2)

data.frame(age.grp = netstats$attr$age.grp,
           age = netstats$attr$age)


prep_start <- 52 * 2
param <- EpiModel::param.net(
  data.frame.params = readr::read_csv("data/input/params_chistig_apr18.csv"),
  netstats          = netstats,
  epistats          = epistats,
  prep.start        = Inf,
  riskh.start       = Inf
)

# param$venue.treat <- 1
# param$app.treat <- 1

# param <- EpiModel::param.net(
#   data.frame.params = readr::read_csv("data/input/params.csv"),
#   netstats          = netstats,
#   epistats          = epistats,
#   prep.start        = prep_start,
#   riskh.start       = prep_start - 53
# )

# See full listing of parameters
# See ?param_msm for definitions
print(param)

# Initial conditions (default prevalence initialized in epistats)
# For models with bacterial STIs, these must be initialized here with non-zero
# values
init <- EpiModelHIV::init_msm(
  prev.ugc = 0,
  prev.rct = 0,
  prev.rgc = 0,
  prev.uct = 0
)


# Control settings
control <- EpiModelHIV::control_msm(
  #  nsteps = 200,
  nsteps = calibration_end + 520,
  # nsteps =  prep_start + 52 * 2,
  nsims = 10,
  ncores = 10,
  cumulative.edgelist = TRUE,
  truncate.el.cuml = Inf,
  raw.output = TRUE,
  save.nwstats = TRUE,
  .tracker.list       = calibration_trackers,
  initialize.FUN =                  chiSTIGmodules::initialize_msm_chi,
  aging.FUN =                       chiSTIGmodules::aging_msm_chi,
  departure.FUN =                   chiSTIGmodules::departure_msm_chi,
  arrival.FUN =                     chiSTIGmodules::arrival_msm_chi,
  # venues.FUN =                      chiSTIGmodules:::venues_msm_chi,
  partident.FUN =                   chiSTIGmodules::partident_msm_chi,
  hivtest.FUN =                     chiSTIGmodules::hivtest_msm_chi,
  hivtx.FUN =                       chiSTIGmodules::hivtx_msm_chi,
  hivprogress.FUN =                 chiSTIGmodules::hivprogress_msm_chi,
  hivvl.FUN =                       chiSTIGmodules::hivvl_msm_chi,
  resim_nets.FUN =                  chiSTIGmodules:::simnet_msm_chi,
  acts.FUN =                        chiSTIGmodules::acts_msm_chi,
  condoms.FUN =                     chiSTIGmodules::condoms_msm_chi,
  position.FUN =                    chiSTIGmodules::position_msm_chi,
  prep.FUN =                        chiSTIGmodules::prep_msm_chi,
  hivtrans.FUN =                    chiSTIGmodules::hivtrans_msm_chi,
  exotrans.FUN =                    chiSTIGmodules:::exotrans_msm_chi,
  stitrans.FUN =                    chiSTIGmodules::stitrans_msm_chi,
  stirecov.FUN =                    chiSTIGmodules::stirecov_msm_chi,
  stitx.FUN =                       chiSTIGmodules::stitx_msm_chi,
  prev.FUN =                        chiSTIGmodules::prevalence_msm_chi,
  cleanup.FUN =                     chiSTIGmodules::cleanup_msm_chi,
  module.order = c( "aging.FUN", "departure.FUN", "arrival.FUN", # "venues.FUN",
                    "partident.FUN", "hivtest.FUN", "hivtx.FUN", "hivprogress.FUN",
                    "hivvl.FUN", "resim_nets.FUN", "acts.FUN", "condoms.FUN",
                    "position.FUN", "prep.FUN", "hivtrans.FUN", "exotrans.FUN",
                    "stitrans.FUN", "stirecov.FUN", "stitx.FUN", "prev.FUN", "cleanup.FUN")
)


# # Read in custom modifications of functions for debugging
# module_edits <- list.files("./R/epimodel_edits")
# for (i in 1:length(module_edits)) {
#   source(paste("./R/epimodel_edits/", module_edits[[i]], sep = ""))
# }
#
# # Control settings
# control <- EpiModelHIV::control_msm(
#   nsteps = 250,
#   nsims = 1,
#   ncores = 1,
#   initialize.FUN =                  initialize_msm_chi,
#   aging.FUN =                       aging_msm_chi,
#   departure.FUN =                   departure_msm_chi,
#   arrival.FUN =                     arrival_msm_chi,
#   partident.FUN =                   partident_msm_chi,
#   hivtest.FUN =                     hivtest_msm_chi,
#   hivtx.FUN =                       hivtx_msm_chi,
#   hivprogress.FUN =                 hivprogress_msm_chi,
#   hivvl.FUN =                       hivvl_msm_chi,
#   resim_nets.FUN =                  simnet_msm_chi,
#   acts.FUN =                        acts_msm_chi,
#   condoms.FUN =                     condoms_msm_chi,
#   position.FUN =                    position_msm_chi,
#   prep.FUN =                        prep_msm_chi,
#   hivtrans.FUN =                    hivtrans_msm_chi,
#   stitrans.FUN =                    stitrans_msm_chi,
#   stirecov.FUN =                    stirecov_msm_chi,
#   stitx.FUN =                       stitx_msm_chi,
#   prev.FUN =                        prevalence_msm_chi,
#   cleanup.FUN =                     cleanup_msm_chi
# )

# Control settings (DEFAULT EXAMPLE)
# control <- control_msm(
#   nsteps = 3100,
#   nsims = 1,
#   ncores = 1,
#
#   .tracker.list       = calibration_trackers,
#
#   initialize.FUN = initialize_msm,
#   aging.FUN = aging_msm,
#   departure.FUN = departure_msm,
#   arrival.FUN = arrival_msm,
#   partident.FUN = partident_msm,
#   hivtest.FUN = hivtest_msm,
#   hivtx.FUN = hivtx_msm,
#   hivprogress.FUN = hivprogress_msm,
#   hivvl.FUN = hivvl_msm,
#   resim_nets.FUN = simnet_msm,
#   acts.FUN = acts_msm,
#   condoms.FUN = condoms_msm,
#   position.FUN = position_msm,
#   prep.FUN = prep_msm,
#   hivtrans.FUN = hivtrans_msm,
#   stitrans.FUN = stitrans_msm,
#   stirecov.FUN = stirecov_msm,
#   stitx.FUN = stitx_msm,
#   prev.FUN = prevalence_msm,
#   verbose.FUN = verbose.net,
#   cleanup.FUN = cleanup_msm
# )


# See listing of modules and other control settings
# Module function defaults defined in ?control_msm
print(control)


# Epidemic simulation
# debug(acts_msm_chi)
# debug(condoms_msm_chi)
# debug(position_msm_chi)
# debug(hivtrans_msm_chi)
sim <- EpiModel::netsim(est, param, init, control)

sim_run_list[[j]] <- sim

}

saveRDS(sim_run_list, "~/Desktop/sim_runs_aug26_2.rds")

sim_run_list <- readRDS("~/Desktop/sim_runs_aug26_1.rds")
sim_run_list2 <- readRDS("~/Desktop/sim_runs_aug26_2.rds")
sim_run_list2 <- sim_run_list2[5:9]






# sim_run_list <- readRDS("~/Desktop/sim_runs_aug15.rds")

list_names <- stringr::str_replace_all(netstats_files, "netstats-level-local_", "")
list_names <- stringr::str_replace_all(list_names, ".rds", "")

# names(sim_run_list) <- netest_files
# names(sim_run_list) <- as.character(c(.0014, .0018, .0015, .0020, .0030))
# names(sim_run_list1) <- list_names[1:]
names(sim_run_list) <- list_names[1:4]
names(sim_run_list2) <- list_names[5:9]

compile_epi <- function(x) {

for (i in 1:length(x)) {

  this_run <- as.data.frame(x[[i]]$epi) %>% dplyr::mutate(
    cc.dx.B         = i_dx__B / i__B,
    cc.dx.H         = i_dx__H / i__H,
    cc.dx.O         = i_dx__O / i__O,
    cc.dx.W         = i_dx__W / i__W,
    cc.linked1m.B   = linked1m__B / i_dx__B,
    cc.linked1m.H   = linked1m__H / i_dx__H,
    cc.linked1m.O   = linked1m__O / i_dx__O,
    cc.linked1m.W   = linked1m__W / i_dx__W,
    cc.vsupp.B      = i_sup__B / i_dx__B,
    cc.vsupp.H      = i_sup__H / i_dx__H,
    cc.vsupp.O      = i_sup__O / i_dx__O,
    cc.vsupp.W      = i_sup__W / i_dx__W,
    gc_s            = gc_s__B + gc_s__H + gc_s__O + gc_s__W,
    ir100.gc        = incid.gc / gc_s * 5200,
    ct_s            = ct_s__B + ct_s__H + ct_s__O + ct_s__W,
    ir100.ct        = incid.ct / ct_s * 5200,
    i.prev.dx.B     = i_dx__B / n__B,
    i.prev.dx.H     = i_dx__H / n__H,
    i.prev.dx.O     = i_dx__O / n__O,
    i.prev.dx.W     = i_dx__W / n__W,
    cc.prep.B       = s_prep__B / s_prep_elig__B,
    cc.prep.H       = s_prep__H / s_prep_elig__H,
    cc.prep.O       = s_prep__O / s_prep_elig__O,
    cc.prep.W       = s_prep__W / s_prep_elig__W,
    # Adding additional measures where denominator is all MSM
    cc.prep.B_all       = s_prep__B / n__B,
    cc.prep.H_all       = s_prep__H / n__H,
    cc.prep.O_all       = s_prep__O / n__O,
    cc.prep.W_all       = s_prep__W / n__W,

    prep_users      = s_prep__B + s_prep__H + s_prep__O + s_prep__W,
    prep_elig       = s_prep_elig__B + s_prep_elig__H + s_prep__O + s_prep_elig__W,
    cc.prep         = prep_users / prep_elig,
    prep_prop_ret1y = prep_ret1y / lag(prep_startat, 52),
    prep_prop_ret2y = prep_ret2y / lag(prep_startat, 104)
  )
  this_run$time <- 1:nrow(this_run)

  if (i == 1) {
    run_data <- this_run
  } else {
    run_data <- dplyr::bind_rows(run_data, this_run)
  }

}

  return(run_data)

}

epi_dfs1 <- lapply(sim_run_list, compile_epi)
epi_dfs2 <- lapply(sim_run_list2, compile_epi)

for (i in 1:length(epi_dfs1)) {
  epi_dfs1[[i]]$drate <- names(epi_dfs1)[[i]]
}

for (i in 1:length(epi_dfs2)) {
  epi_dfs2[[i]]$drate <- names(epi_dfs2)[[i]]
}

full_data1 <- dplyr::bind_rows(epi_dfs1)
full_data2 <- dplyr::bind_rows(epi_dfs2)

full_data <- dplyr::bind_rows(full_data1, full_data2)

write.csv(full_data, "~/Desktop/drate_check/aug26_runs.csv")

# Compile Edgelists
compile_el <- function(x) {

  for (i in 1:length(x)) {

    this_el <- x[[i]]$el.cuml
    for (j in 1:3) {
      this_el[[j]]$type <- j
    }

    # Merge and indicate which run this is
    this_el <- dplyr::bind_rows(this_el) %>% dplyr::mutate(run = i)

    if (i == 1) {
      run_data <- this_el
    } else {
      run_data <- dplyr::bind_rows(run_data, this_el)
    }

  }

  return(run_data)

}

el_dfs1 <- lapply(sim_run_list, compile_el)
el_dfs2 <- lapply(sim_run_list2, compile_el)

for (i in 1:length(el_dfs1)) {
  el_dfs1[[i]]$drate <- names(el_dfs1)[[i]]
}

for (i in 1:length(el_dfs2)) {
  el_dfs2[[i]]$drate <- names(el_dfs2)[[i]]
}

full_el1 <- dplyr::bind_rows(el_dfs1)
full_el2 <- dplyr::bind_rows(el_dfs2)

full_el <- dplyr::bind_rows(full_el1, full_el2)

write.csv(full_el, "~/Desktop/drate_check/aug26_edgelists.csv")


library(tidyverse)

full_data %>%
  ggplot(aes(x = time, y = n_edges_main, color = drate)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 1630.057344, linetype = "dotted") +
  geom_hline(yintercept = netstats$main$edges) +
  geom_hline(yintercept = 1172.186555, linetype = "dotted") +
  labs(title = "Number of Main Partnerships",
       x = "Time",
       y = NULL)

full_data %>%
  ggplot(aes(x = time, y = n_edges_casual, color = drate)) +
  geom_line(alpha = .1) +
  theme_classic() +
  geom_hline(yintercept = 2386.663237, linetype = "dotted") +
  geom_hline(yintercept = netstats$casl$edges) +
  geom_hline(yintercept = 1591.149932, linetype = "dotted") +
  labs(title = "Number of Casual Partnerships",
       x = "Time",
       y = NULL)


full_data %>%
  ggplot(aes(x = time, y = mean_deg_main, color = drate)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.200339524, linetype = "dotted") +
  geom_hline(yintercept = 0.239467091) +
  geom_hline(yintercept = 0.278594658, linetype = "dotted") +
  labs(title = "Mean Degree, Main Partnerships",
       x = "Time",
       y = NULL)


full_data %>%
  mutate(drate_cas = stringr::str_extract(drate, "\\d*$")) %>%
  ggplot(aes(x = time, y = mean_deg_cas, color = drate_cas)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.271944955, linetype = "dotted") +
  geom_hline(yintercept = 0.339925924) +
  geom_hline(yintercept = 0.407906894, linetype = "dotted") +
  labs(title = "Mean Degree, Casual Partnerships",
       x = "Time",
       y = NULL)



full_data %>%
  filter(drate == "00018_00014") %>%
  ggplot(aes(x = time, y = mean_deg_main)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.200339524, linetype = "dotted") +
  geom_hline(yintercept = 0.239467091) +
  geom_hline(yintercept = 0.278594658, linetype = "dotted") +
  labs(title = "Mean Degree, Main Partnerships (drate = .0018)",
       x = "Time",
       y = NULL)


full_data %>%
  filter(drate == "00018_00014") %>%
  ggplot(aes(x = time, y = mean_deg_cas)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.271944955, linetype = "dotted") +
  geom_hline(yintercept = 0.339925924) +
  geom_hline(yintercept = 0.407906894, linetype = "dotted") +
  labs(title = "Mean Degree, Casual Partnerships (drate = .0014)",
       x = "Time",
       y = NULL)






# 16 to 20
full_data %>%
  ggplot(aes(x = time, y = mean_deg_main.under21, color = drate)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.121093802, linetype = "dotted") +
  geom_hline(yintercept = 0.175437267) +
  geom_hline(yintercept = 0.229780732, linetype = "dotted") +
  labs(title = "Mean Degree, Main Partnerships (Under 21)",
       x = "Time",
       y = NULL)

full_data %>%
  ggplot(aes(x = time, y = mean_deg_cas.under21, color = drate)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.268377431, linetype = "dotted") +
  geom_hline(yintercept = 0.380896673) +
  geom_hline(yintercept = 0.493415915, linetype = "dotted") +
  labs(title = "Mean Degree, Casual Partnerships (Under 21)",
       x = "Time",
       y = NULL)

# 21plus
full_data %>%
  ggplot(aes(x = time, y = mean_deg_main.over21, color = drate)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.213707794, linetype = "dotted") +
  geom_hline(yintercept = 0.263479241) +
  geom_hline(yintercept = 0.313250688, linetype = "dotted") +
  labs(title = "Mean Degree, Main Partnerships (Over 21)",
       x = "Time",
       y = NULL)


full_data %>%
  ggplot(aes(x = time, y = mean_deg_cas.over21, color = drate)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.241136199, linetype = "dotted") +
  geom_hline(yintercept = 0.324561275) +
  geom_hline(yintercept = 0.407986352, linetype = "dotted") +
  labs(title = "Mean Degree, Casual Partnerships (Over 21)",
       x = "Time",
       y = NULL)


full_data %>%
  ggplot(aes(x = time, y = dep.gen, color = drate)) +
  geom_line(alpha = .4) +
  theme_classic() +
  # geom_hline(yintercept = 0.358558943, linetype = "dotted") +
  # geom_hline(yintercept = 0.531606308) +
  # geom_hline(yintercept = 0.704653674, linetype = "dotted") +
  labs(title = "Departures",
       x = "Time",
       y = NULL)

full_data %>%
  mutate(aids_deaths = dep.AIDS.on.tx + dep.AIDS.off.tx) %>%
  ggplot(aes(x = time, y = aids_deaths, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  # geom_hline(yintercept = 0.358558943, linetype = "dotted") +
  # geom_hline(yintercept = 0.531606308) +
  # geom_hline(yintercept = 0.704653674, linetype = "dotted") +
  labs(title = "Departures",
       x = "Time",
       y = NULL)


# 5-Year Age Bins
### Main
full_data %>%
  ggplot(aes(x = time, y = mean_deg_main.15to20, color = drate)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.121093802, linetype = "dotted") +
  geom_hline(yintercept = 0.175437267) +
  geom_hline(yintercept = 0.229780732, linetype = "dotted") +
  labs(title = "Mean Degree, Main Partnerships (15-20)",
       x = "Time",
       y = NULL)

full_data %>%
  ggplot(aes(x = time, y = mean_deg_main.21to25, color = drate)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.213707794, linetype = "dotted") +
  geom_hline(yintercept = 0.263479241) +
  geom_hline(yintercept = 0.313250688, linetype = "dotted") +
  labs(title = "Mean Degree, Main Partnerships (21-25)",
       x = "Time",
       y = NULL)

full_data %>%
  ggplot(aes(x = time, y = mean_deg_main.25to30, color = drate)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.213707794, linetype = "dotted") +
  geom_hline(yintercept = 0.263479241) +
  geom_hline(yintercept = 0.313250688, linetype = "dotted") +
  labs(title = "Mean Degree, Main Partnerships (25+)",
       x = "Time",
       y = NULL)


### Casual NEED TO RERUN
full_data %>%
  ggplot(aes(x = time, y = mean_deg_cas.15to20, color = drate)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.268377431, linetype = "dotted") +
  geom_hline(yintercept = 0.380896673) +
  geom_hline(yintercept = 0.493415915, linetype = "dotted") +
  labs(title = "Mean Degree, Casual Partnerships (15-20)",
       x = "Time",
       y = NULL)

full_data %>%
  ggplot(aes(x = time, y = mean_deg_cas.21to25, color = drate)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.241136199, linetype = "dotted") +
  geom_hline(yintercept = 0.324561275) +
  geom_hline(yintercept = 0.407986352, linetype = "dotted") +
  labs(title = "Mean Degree, Casual Partnerships (21-25)",
       x = "Time",
       y = NULL)

full_data %>%
  ggplot(aes(x = time, y = mean_deg_cas.25to30, color = drate)) +
  geom_line(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.241136199, linetype = "dotted") +
  geom_hline(yintercept = 0.324561275) +
  geom_hline(yintercept = 0.407986352, linetype = "dotted") +
  labs(title = "Mean Degree, Casual Partnerships (25+)",
       x = "Time",
       y = NULL)




# Proportion Diagnosed

full_data %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = cc.dx.B, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.546535643) +
  labs(title = "Proportion Diagnosed (Black)",
       x = "Time",
       y = NULL)

full_data %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = cc.dx.H, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.5431367893) +
  labs(title = "Proportion Diagnosed (Hispanic)",
       x = "Time",
       y = NULL)

full_data %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = cc.dx.O, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.5601310) +
  labs(title = "Proportion Diagnosed (Other)",
       x = "Time",
       y = NULL)

full_data %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = cc.dx.W, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.5988779867) +
  labs(title = "Proportion Diagnosed (White)",
       x = "Time",
       y = NULL)

# Linked to Care

full_data %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = cc.linked1m.B, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = .828) +
  labs(title = "Proportion Linked to Care (Black)",
       x = "Time",
       y = NULL)

full_data %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = cc.linked1m.H, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.867) +
  labs(title = "Proportion Linked to Care (Hispanic)",
       x = "Time",
       y = NULL)

full_data %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = cc.linked1m.O, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.875) +
  labs(title = "Proportion Linked to Care (Other)",
       x = "Time",
       y = NULL)

full_data %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = cc.linked1m.W, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.936) +
  labs(title = "Proportion Linked to Care (White)",
       x = "Time",
       y = NULL)

# Proportion with Suppression

full_data %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = cc.vsupp.B, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.571) +
  labs(title = "Proportion with Viral Suppression (Black)",
       x = "Time",
       y = NULL)


full_data %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = cc.vsupp.H, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.675) +
  labs(title = "Proportion with Viral Suppression (Hispanic)",
       x = "Time",
       y = NULL)

full_data %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = cc.vsupp.O, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.586) +
  labs(title = "Proportion with Viral Suppression (Other)",
       x = "Time",
       y = NULL)

full_data %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = cc.vsupp.W, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.617) +
  labs(title = "Proportion with Viral Suppression (White)",
       x = "Time",
       y = NULL)

# Total Incidence Rate

full_data %>%
  filter(drate == "00018_00014") %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = ir100.B, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 4.44, linetype = "dotted") +
  geom_hline(yintercept = 6.42) +
  geom_hline(yintercept = 9.30, linetype = "dotted") +
  labs(title = "Total Incidence Rate (Black)",
       x = "Time",
       y = NULL)


full_data %>%
  filter(drate == "00018_00014") %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = ir100.H, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 1.10, linetype = "dotted") +
  geom_hline(yintercept = 2.04) +
  geom_hline(yintercept = 3.79, linetype = "dotted") +
  labs(title = "Total Incidence Rate (Hispanic)",
       x = "Time",
       y = NULL)

full_data %>%
  filter(drate == "00018_00014") %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = ir100.O, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.55, linetype = "dotted") +
  geom_hline(yintercept = 1.71) +
  geom_hline(yintercept = 5.31, linetype = "dotted") +
  labs(title = "Total Incidence Rate (Other)",
       x = "Time",
       y = NULL)

full_data %>%
  filter(drate == "00018_00014") %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = ir100.W, color = drate)) +
  geom_smooth(alpha = .4) +
  theme_classic() +
  geom_hline(yintercept = 0.24, linetype = "dotted") +
  geom_hline(yintercept = 0.73) +
  geom_hline(yintercept = 2.26, linetype = "dotted") +
  labs(title = "Total Incidence Rate (White)",
       x = "Time",
       y = NULL)

# Prevalence


full_data %>%
  filter(drate == "00018_00014") %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = i.prev.B, color = drate)) +
  geom_smooth(alpha = .4) +
  geom_hline(yintercept = 0.354) +
  theme_classic() +
  labs(title = "Prevalence (Black)",
       x = "Time",
       y = NULL)


full_data %>%
  filter(drate == "00018_00014") %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = i.prev.H, color = drate)) +
  geom_smooth(alpha = .4) +
  geom_hline(yintercept = 0.114) +
  theme_classic() +
  labs(title = "Prevalence (Hispanic)",
       x = "Time",
       y = NULL)

full_data %>%
  filter(drate == "00018_00014") %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = i.prev.O, color = drate)) +
  geom_smooth(alpha = .4) +
  geom_hline(yintercept = 0.103) +
  theme_classic() +
  labs(title = "Prevalence (Other)",
       x = "Time",
       y = NULL)

full_data %>%
  filter(drate == "00018_00014") %>%
  filter(time > 3000) %>%
  ggplot(aes(x = time, y = i.prev.W, color = drate)) +
  geom_smooth(alpha = .4) +
  geom_hline(yintercept = 0.106) +
  theme_classic() +
  labs(title = "Prevalence (White)",
       x = "Time",
       y = NULL)



### Durations
full_el %>%
  filter(type ==1 & start > 3000 & start < (3640 - 115)) %>%
  mutate(duration = stop-start# ,
         # drate2 = stringr::str_extract(drate, "\\d*.rds$")
         ) %>%
  # filter(type == 1) %>%
  ggplot(aes(y = duration, x = drate)) +
  geom_boxplot() +
  theme_classic() +
  geom_hline(yintercept = 87) +
  labs(title = "Partnership Durations (Main)",
       x = "d.rate",
       y = "Duration (Weeks)")

full_el %>%
  filter(type == 2 & start > 3000 & start < (3640 - 72)) %>%
  mutate(duration = stop-start,
         drate2 = stringr::str_extract(drate, "\\d*.rds$")) %>%
  # filter(type == 1) %>%
  ggplot(aes(y = duration, x = drate)) +
  geom_boxplot() +
  theme_classic() +
  geom_hline(yintercept = 57) +
  labs(title = "Partnership Durations (Casual)",
       x = "d.rate",
       y = "Duration (Weeks)")

check <- full_el %>%
  filter(type ==1 & start > 3000 & start < (3640 - 172)) %>%
  mutate(duration = stop-start,
         drate2 = stringr::str_extract(drate, "\\d*.rds$")) %>%
  group_by(drate2, type, head, tail) %>%
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  ungroup() %>%
  slice(1)

full_el %>%
  mutate(duration = stop-start,
         drate2 = stringr::str_extract(drate, "\\d*.rds$")) %>%
  filter(drate2 == check$drate2 & type == 1 & head == check$head & tail == check$tail)
  # filter(type == 1) %>%
  ggplot(aes(y = duration, x = drate2)) +
  geom_boxplot() +
  theme_classic() +
  geom_hline(yintercept = 172) +
  labs(title = "Partnership Durations (Main)",
       x = "d.rate",
       y = "Duration (Weeks)")

full_el %>%
  mutate(duration = stop-start,
         drate2 = stringr::str_extract(drate, "\\d*.rds$")) %>%
  filter(type == 2) %>%
  ggplot(aes(y = duration, x = drate)) +
  geom_boxplot() +
  theme_classic() +
  geom_hline(yintercept = 115)





test2 <- as.data.frame(test[[1]]$epi)


sim_list <- as.list(sim)


plot(sim_list[[1]]$epi$mean_deg_main, main = "Main Partnerships")


plot(sim_list[[1]]$epi$n_edges_main/sim_list[[1]]$epi$num, main = "Main Partnerships")

library(tidyverse)

# GENERAL NUMBER OF EDGES

data.frame(y = sim_list[[1]]$epi$n_edges_main) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  geom_hline(yintercept = 1630.057344, linetype = "dotted") +
  geom_hline(yintercept = netstats$main$edges) +
  geom_hline(yintercept = 1172.186555, linetype = "dotted") +
  labs(title = "Number of Main Partnerships",
       x = "Time",
       y = NULL)

data.frame(y = sim_list[[1]]$epi$mean_deg_main) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  # geom_hline(yintercept = 1630.057344, linetype = "dotted") +
  # geom_hline(yintercept = netstats$main$edges) +
  # geom_hline(yintercept = 1172.186555, linetype = "dotted") +
  labs(title = "Mean Degree (Main Partnerships)",
       x = "Time",
       y = NULL)


data.frame(y = sim_list[[1]]$epi$mean_deg_main.under21) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  # geom_hline(yintercept = 1630.057344, linetype = "dotted") +
  # geom_hline(yintercept = netstats$main$edges) +
  # geom_hline(yintercept = 1172.186555, linetype = "dotted") +
  labs(title = "Mean Degree (Main Partnerships, Under 21)",
       x = "Time",
       y = NULL)

data.frame(y = sim_list[[1]]$epi$mean_deg_main.over21) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  # geom_hline(yintercept = 1630.057344, linetype = "dotted") +
  # geom_hline(yintercept = netstats$main$edges) +
  # geom_hline(yintercept = 1172.186555, linetype = "dotted") +
  labs(title = "Mean Degree (Main Partnerships, Over 21)",
       x = "Time",
       y = NULL)


data.frame(y = sim_list[[1]]$epi$mean_deg_cas.under21) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  # geom_hline(yintercept = 1630.057344, linetype = "dotted") +
  # geom_hline(yintercept = netstats$main$edges) +
  # geom_hline(yintercept = 1172.186555, linetype = "dotted") +
  labs(title = "Mean Degree (Casual Partnerships, Under 21)",
       x = "Time",
       y = NULL)

data.frame(y = sim_list[[1]]$epi$mean_deg_cas.over21) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  # geom_hline(yintercept = 1630.057344, linetype = "dotted") +
  # geom_hline(yintercept = netstats$main$edges) +
  # geom_hline(yintercept = 1172.186555, linetype = "dotted") +
  labs(title = "Mean Degree (Casual Partnerships, Over 21)",
       x = "Time",
       y = NULL)


data.frame(y = sim_list[[1]]$epi$n_edges_casual) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  geom_hline(yintercept = 2386.663237, linetype = "dotted") +
  geom_hline(yintercept = netstats$casl$edges) +
  geom_hline(yintercept = 1591.149932, linetype = "dotted") +
  labs(title = "Number of Casual Partnerships",
       x = "Time",
       y = NULL)

data.frame(y = sim_list[[1]]$epi$mean_deg_cas) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  # geom_hline(yintercept = 1630.057344, linetype = "dotted") +
  # geom_hline(yintercept = netstats$main$edges) +
  # geom_hline(yintercept = 1172.186555, linetype = "dotted") +
  labs(title = "Mean Degree (Casual Partnerships)",
       x = "Time",
       y = NULL)

data.frame(y = sim_list[[1]]$epi$n_edges_onetime) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  geom_hline(yintercept = 95.26738158, linetype = "dotted") +
  geom_hline(yintercept = netstats$inst$edges) +
  geom_hline(yintercept = 60.6132855832735, linetype = "dotted") +
  labs(title = "Number of One-Time Partnerships",
       x = "Time",
       y = NULL)

# POPULATION SIZE

data.frame(y = sim_list[[1]]$epi$num) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  geom_hline(yintercept = length(netstats$attr$numeric.id)) +
  labs(title = "Population Size",
       x = "Time",
       y = NULL)

# PROPORTION BY RACE

data.frame(y = sim_list[[1]]$epi$n__B/sim_list[[1]]$epi$num) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  geom_hline(yintercept = mean(netstats$attr$race == 1)) +
  labs(title = "Proportion of Population Black",
       x = "Time",
       y = "Proportion") +
  scale_y_continuous(breaks = seq(0, .4, .1),
                     limits = c(0, .4))

data.frame(y = sim_list[[1]]$epi$n__H/sim_list[[1]]$epi$num) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  geom_hline(yintercept = mean(netstats$attr$race == 2)) +
  labs(title = "Proportion of Population Hispanic",
       x = "Time",
       y = "Proportion") +
  scale_y_continuous(breaks = seq(0, .4, .1),
                     limits = c(0, .4))

data.frame(y = sim_list[[1]]$epi$n__O/sim_list[[1]]$epi$num) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  geom_hline(yintercept = mean(netstats$attr$race == 3)) +
  labs(title = "Proportion of Population Other",
       x = "Time",
       y = "Proportion") +
  scale_y_continuous(breaks = seq(0, .4, .1),
                     limits = c(0, .4))

data.frame(y = sim_list[[1]]$epi$n__W/sim_list[[1]]$epi$num) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  geom_hline(yintercept = mean(netstats$attr$race == 4)) +
  labs(title = "Proportion of Population White",
       x = "Time",
       y = "Proportion") +
  scale_y_continuous(breaks = seq(0, .4, .1),
                     limits = c(0, .4))

# PROPORTION BY AGE GROUP

data.frame(y = sim_list[[1]]$epi$num.16to20/sim_list[[1]]$epi$num) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
  theme_classic() +
  geom_hline(yintercept = mean(netstats$attr$age.grp == 1)) +
  labs(title = "Proportion of Population Aged 16-20",
       x = "Time",
       y = "Proportion") +
  scale_y_continuous(breaks = seq(0, .4, .1),
                     limits = c(0, .4))

data.frame(y = sim_list[[1]]$epi$num.21plus/sim_list[[1]]$epi$num) %>%
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x,
             y = y)) +
  geom_line() +
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



# NUMBER OF EDGES (AGE GROUP)
### Main
data.frame(under21 = sim_list[[1]]$epi$n_edges_main.under21,
           over21 = sim_list[[1]]$epi$n_edges_main.over21,
           # diffage = sim_list[[1]]$epi$n_edges_main.diffage,
           sameage = sim_list[[1]]$epi$n_edges_main.sameage) %>%
  mutate(x = 1:n()) %>%
  tidyr::pivot_longer(cols = under21:sameage) %>%
  ggplot(aes(x = x,
             y = value,
             color = name)) +
  geom_line() +
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
data.frame(under21 = sim_list[[1]]$epi$n_edges_casual.under21,
           over21 = sim_list[[1]]$epi$n_edges_casual.over21,
           # diffage = sim_list[[1]]$epi$n_edges_casual.diffage,
           sameage = sim_list[[1]]$epi$n_edges_casual.sameage) %>%
  mutate(x = 1:n()) %>%
  tidyr::pivot_longer(cols = under21:sameage) %>%
  ggplot(aes(x = x,
             y = value,
             color = name)) +
  geom_line() +
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


data.frame(under21 = sim_list[[1]]$epi$n_edges_onetime.under21,
           over21 = sim_list[[1]]$epi$n_edges_onetime.over21,
           # diffage = sim_list[[1]]$epi$n_edges_onetime.diffage,
           sameage = sim_list[[1]]$epi$n_edges_onetime.sameage) %>%
  mutate(x = 1:n()) %>%
  tidyr::pivot_longer(cols = under21:sameage) %>%
  ggplot(aes(x = x,
             y = value,
             color = name)) +
  geom_line() +
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



# T1 Edgelist for István
sim_list[[1]]$raw.records[[1]]


plot(sim_list[[1]]$epi$n_edges_casual, main = "Casual Partnerships")
plot(sim_list[[1]]$epi$n_edges_onetime)

plot(sim_list[[1]]$epi$num.B)
plot(sim_list[[1]]$epi$num.mainzero)
plot(sim_list[[1]]$epi$num.main1plus)

plot(sim_list[[1]]$epi$num.caslzero)
plot(sim_list[[1]]$epi$num.casl1)
plot(sim_list[[1]]$epi$num.casl2plus)

plot(sim_list[[1]]$epi$num.16to20)
plot(sim_list[[1]]$epi$num.21plus)



plot(sim_list[[1]]$epi$n_edges_main.under21, main = "Main Partnerships, Both Under 21")
plot(sim_list[[1]]$epi$n_edges_main.over21, main = "Main Partnerships, Both Over 21")
plot(sim_list[[1]]$epi$n_edges_main.diffage, main = "Main Partnerships, Cross-Age Group")


plot(
  sim_list[[1]]$epi$num.W,
  sim_list[[1]]$epi$n_edges_main,
  main = "Main Partnerships")

# Combine pertinent elements of `sim_list` to a single data frame
epi_wrangle <- function(x) {
  this_epi <- dplyr::bind_rows(x$epi)
  this_epi$time = 1:nrow(this_epi)
  return(this_epi)
}

el_wrangle <- function(x) {
  this_epi <- dplyr::bind_rows(x$object)
  # this_epi$time = 1:nrow(this_epi)
  return(this_epi)
}

full_data <- dplyr::bind_rows(lapply(sim_list, epi_wrangle))
full_el <- dplyr::bind_rows(lapply(sim_list[[1]]$raw.records, el_wrangle))
full_data$sim <- rep(1:20, each = 3640)

library(tidyverse)

target_plot(full_data, var = "num", group = "sim", benchmark = 11612)

target_plot(full_data, var = "num.16to20", group = "sim", benchmark = 3125)
target_plot(full_data, var = "num.21plus", group = "sim", benchmark = 8487)

target_plot(full_data, var = "n_edges_main", group = "sim")
target_plot(full_data, var = "n_edges_casual", group = "sim")

i = 1
target_plot(full_data, var = "num.mainzero", group = "sim",
            title = paste("Plot ", i, ": Number of Nodes with 0 Ties (Main Partnerships)", sep = ""))
i = 1 + i
target_plot(full_data, var = "num.main1plus", group = "sim",
            title = paste("Plot ", i, ": Number of Nodes with 1+ Ties (Main Partnerships)", sep = ""))
i = 1 + i
target_plot(full_data, var = "num.caslzero", group = "sim",
            title = paste("Plot ", i, ": Number of Nodes with 0 Ties (Casual Partnerships)", sep = ""))
i = 1 + i
target_plot(full_data, var = "num.casl1", group = "sim",
            title = paste("Plot ", i, ": Number of Nodes with 1 Tie (Casual Partnerships)", sep = ""))
i = 1 + i
target_plot(full_data, var = "num.casl2plus", group = "sim",
            title = paste("Plot ", i, ": Number of Nodes with 2+ Ties (Casual Partnerships)", sep = ""))
i = 1 + i
target_plot(full_data, var = "num.totzero", group = "sim",
            title = paste("Plot ", i, ": Number of Nodes with 0 Ties (All Partnerships)", sep = ""))
i = 1 + i
target_plot(full_data, var = "num.tot1", group = "sim",
            title = paste("Plot ", i, ": Number of Nodes with 1 Tie (All Partnerships)", sep = ""))
i = 1 + i
target_plot(full_data, var = "num.tot2", group = "sim",
            title = paste("Plot ", i, ": Number of Nodes with 2 Ties (All Partnerships)", sep = ""))
i = 1 + i
target_plot(full_data, var = "num.tot3plus", group = "sim",
            title = paste("Plot ", i, ": Number of Nodes with 3+ Ties (All Partnerships)", sep = ""))



# full_data <- dplyr::bind_rows(lapply(sim_list, function(x){x$epi}))



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
merged_data <- edges %>% dplyr::left_join(inf_events, by = c("head_uid", "tail_uid"))

# What percentage of cases in the infection act list successfully merged?
# Note that we use `nrow(inf_events2)` here for the accurate denominator
sum(merged_data$infection, na.rm = TRUE)/nrow(inf_events2)

merged_data %>% group_by(type) %>%
  summarize(sum(infection, na.rm = TRUE))

# Now extract just the partnerships that led to infection
inf_edges <- merged_data %>%
  filter(infection == 1)


for (i in 1:100) {

  inf_events <- full_test %>% dplyr::filter(type == 4) %>%
    # Filter time on infections to be last 520 weeks (we store past 600) for edgelist
    filter(time < (max(time) - i)) %>%
    # filter(time == 400) %>%
    select(head_uid, tail_uid, inf_time = time) %>%
    mutate(infection = 1)

  inf_events2 <- inf_events %>% dplyr::select(head_uid = tail_uid,
                                              tail_uid = head_uid,
                                              inf_time, infection)
  inf_events <- dplyr::bind_rows(inf_events, inf_events2)


  inf_events <- inf_events %>% select(head_uid, tail_uid, infection)

  edges <- full_test %>% dplyr::filter(type != 4)

  edge_test <- edges %>% dplyr::left_join(inf_events, by = c("head_uid", "tail_uid"))

  print(sum(edge_test$infection, na.rm = TRUE)/nrow(inf_events2))

}


test2 <- as.data.frame(dplyr::bind_rows(test[[1]]$raw.records))$object
colnames(test2) <- c()

test3 <- test2 %>% dplyr::filter(type == "hivtrans")

test_main <- test2 %>% dplyr::filter(type == "main") %>%
  dplyr::group_by(ego_id, alter_id) %>%
  summarize(count = n())

test3 <- test2 %>%
  dplyr::arrange(time) %>%
  dplyr::group_by(ego_id, alter_id, type) %>%
  dplyr::slice(1)

test4 <- test2 %>%
  dplyr::group_by(ego_id, alter_id) %>%
  summarize(count = n())

sim$netstats

test <- sim[[1]]$el.cuml[[1]]

test <- sim[[1]]$el[[1]]

check_this <- sim[[1]]

# Total Network Size
popsize <- check_this$epi$num
tail(popsize, n = 20)
# Number eligible for PrEP
prepElig <- check_this$epi$prepElig
tail(prepElig, n = 20)
# Number currently on PrEP
prepCurr <- check_this$epi$prepCurr
tail(prepCurr, n = 20)

# Proportion of eligible nodes on PrEP
tail(prepCurr/prepElig, n = 20)

# Proportion of entire population on PrEP
tail(prepCurr/popsize, n = 20)

# Proportion of entire population eligible for PrEP
tail(prepElig/popsize, n = 20)



saveRDS(sim, "./data/intermediate/control_sim.rds")

sim <- readRDS("./data/intermediate/test_sim.rds")

# Examine the model object output
print(sim)

# Plot outcomes
par(mar = c(3, 3, 2, 2), mgp = c(2, 1, 0))
plot(sim, y = "n_edges_main", main = "Main edges")
plot(sim, y = "n_edges_casual", main = "Casual edges")
plot(sim, y = "n_edges_onetime", main = "One-time edges")
plot(sim, y = "ir100", main = "Incidence")

par(mar = c(3, 3, 2, 2), mgp = c(2, 1, 0))
plot(sim, y = "incid", main = "Plot 2: Incidence")
plot(x = 1:250, y = sim$epi$incid$sim1)
abline(h = 5.08)
plot(sim, y = "incid.B", main = "Incidence (Black)")
plot(sim, y = "incid.H", main = "Incidence (Hispanic)")
plot(sim, y = "incid.O", main = "Incidence (Other)")
plot(sim, y = "incid.W", main = "Incidence (White)")


plot(x = 1:250, y = sim$epi$cc.dx$sim1)

plot(sim, y = "num_disc_acts")

plot(sim, y = "num", main = "Plot 1: Population Size")
plot(sim, y = "dep.gen", main = "Plot 2: General Departures")
plot(sim, y = "dep.AIDS.on.tx", main = "Plot 3: AIDS-Related Deaths (On Treatment)")
plot(sim, y = "dep.AIDS.off.tx", main = "Plot 4: AIDS-Related Deaths (No Treatment)")

plot(sim, y = "nNew", main = "Plot 3: New Arrivals")

# Convert to data frame
df <- as.data.frame(sim)
head(df)
tail(df)

death_plot <- df %>%
  mutate(gen_dep_rate = dep.gen / num,
         AIDStx_dep_rate = dep.AIDS.on.tx/num,
         AIDSnotx_dep_rate = dep.AIDS.off.tx/num) %>%
  ggplot(aes(x = time, y = dep.AIDS.on.tx)) +
  geom_line()
death_plot

## Run 2 simulations on 2 cores
## Note: this will not run generate a progress tracker in the console
control <- control_msm(
  nsteps = 250,
  nsims = 2,
  ncores = 2
)
sim <- netsim(est, param, init, control)

par(mfrow = c(2, 1))
plot(sim, y = "i.num", main = "Prevalence")
plot(sim, y = "ir100", main = "Incidence")

## Example debugging of HIV transmission module in debug mode
# Start by sourcing local version of EpiModelHIV
pkgload::load_all("~/Desktop/epimodelhiv_p")

# Rerun control settings (to source in local set of module functions)
# Note: debugging needs to run with 1 simulation on 1 core
control <- control_msm(
  nsteps = 250,
  nsims = 1,
  ncores = 1,
)

library(tidyverse)
sim_tib <- sim %>% as_tibble()


sim_targets <- sim_tib %>%
  mutate_calibration_targets()

sim_targets %>% ggplot(aes(x = time, y = num)) +
  geom_line() +
  ggtitle("Plot 1: Population Size") +
  theme_minimal() +
  geom_hline(aes(yintercept = length(netstats$attr$age), color = "red"))


sim_targets %>%
  mutate(dep.total = dep.gen + dep.AIDS.off.tx + dep.AIDS.on.tx) %>%
  ggplot(aes(x = time, y = dep.gen)) +
  geom_line() +
  ggtitle("Plot 2: Weekly Departures") +
  theme_minimal()

mean(sim_targets$dep.gen, na.rm = TRUE)

sim_targets %>% ggplot(aes(x = time, y = nNew)) +
  geom_line() +
  ggtitle("Plot 3: Weekly Arrivals") +
  theme_minimal()

mean(sim_targets$nNew, na.rm = TRUE)




sim_targets %>% ggplot(aes(x = time, y = n_edges_main)) +
  geom_line()  +
  ggtitle("Plot 4: Number of Edges (Main)") +
  theme_minimal() +
  geom_hline(aes(yintercept = 1163, color = "red")) +
  geom_hline(aes(yintercept = 1617, color = "blue"))
sim_targets %>% ggplot(aes(x = time, y = n_edges_casual)) +
  geom_line() +
  ggtitle("Plot 5: Number of Edges (Casual)") +
  theme_minimal() +
  geom_hline(aes(yintercept = 1578, color = "red")) +
  geom_hline(aes(yintercept = 2368, color = "blue"))
sim_targets %>% ggplot(aes(x = time, y = n_edges_onetime)) +
  geom_line() +
  ggtitle("Plot 6: Number of Edges (One-Time)") +
  theme_minimal() +
  geom_hline(aes(yintercept = 60, color = "red")) +
  geom_hline(aes(yintercept = 94, color = "blue"))



sim_targets %>% ggplot(aes(x = time, y = num_acts)) +
  geom_line() +
  ggtitle("Plot 6a: Number of Sexual Acts in Network") +
  theme_minimal()

sim_targets %>% ggplot(aes(x = time, y = num)) +
  geom_line() +
  geom_hline(aes(yintercept = length(netstats$attr$diag.status), color = "red"))

popsize_check <- sim_targets %>% select(time, num)
nodes <- data.frame(id = 1:length(netstats$attr$age),
                    age = netstats$attr$age) %>%
  filter(age > 29) %>%
  mutate(age = age + 1/52) %>%
  filter(age >= 30)



sim_targets %>% ggplot(aes(x = time, y = i.prev)) +
  geom_line() +
  geom_hline(aes(yintercept = mean(netstats$attr$diag.status), color = "red")) +
  theme_minimal() +
  ggtitle("Plot 7: Prevalence, Infected")

sim_targets %>% ggplot(aes(x = time, y = incid)) +
  geom_line() +
  geom_hline(aes(yintercept = 5.08, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 8: Weekly Incidence (Endogenous)")

sim_targets %>% ggplot(aes(x = time, y = exo.incid)) +
  geom_line() +
  # geom_hline(aes(yintercept = 5.08, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 8a: Weekly Incidence (Exogenous)")

sim_targets %>%
  mutate(total.incid = incid + exo.incid) %>%
  ggplot(aes(x = time, y = total.incid)) +
  geom_line() +
  geom_hline(aes(yintercept = 5.08, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 8b: Weekly Incidence (Combined)")





sim_targets %>% select(time, dep.HIV, incid) %>%
  tidyr::pivot_longer(cols = c("dep.HIV", "incid"),
                      names_to = "Var",
                      values_to = "Count") %>%
  ggplot(aes(x = time, y = Count, color = Var)) + geom_line() +
  theme_minimal()

sim_targets %>% select(time, dep.HIV, incid) %>%
  filter(time < 500) %>%
  ggplot(aes(x = incid, y = dep.HIV)) + geom_point() +
  geom_jitter()+
  theme_minimal()

sim_targets %>% ggplot(aes(x = time, y = cc.dx.B)) +
  geom_line()

sim_targets$cc.dx.B
sim_targets$cc.dx.H
sim_targets$cc.dx.O
sim_targets$cc.dx.W

sim_targets$cc.linked1m.B
sim_targets$cc.linked1m.H
sim_targets$cc.linked1m.O
sim_targets$cc.linked1m.W

sim_targets$cc.vsupp.B
sim_targets %>% ggplot(aes(x = time, y = cc.vsupp.B)) +
  geom_line() +
  geom_hline(aes(yintercept = 0.571, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 7: Proportion of HIV+ Nodes with Viral Suppression (Black)")
sim_targets$cc.vsupp.H
sim_targets %>% ggplot(aes(x = time, y = cc.vsupp.H)) +
  geom_line() +
  geom_hline(aes(yintercept = 0.675, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 8: Proportion of HIV+ Nodes with Viral Suppression\n(Hispanic)")
sim_targets$cc.vsupp.O
sim_targets %>% ggplot(aes(x = time, y = cc.vsupp.O)) +
  geom_line() +
  geom_hline(aes(yintercept = 0.586, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 9: Proportion of HIV+ Nodes with Viral Suppression\n(Other)")
sim_targets$cc.vsupp.W
sim_targets %>% ggplot(aes(x = time, y = cc.vsupp.H)) +
  geom_line() +
  geom_hline(aes(yintercept = 0.617, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 9: Proportion of HIV+ Nodes with Viral Suppression\n(White)")


mean(sim_tib$num_unprotected_acts/sim_tib$num_acts, na.rm = T)

mean(sim_tib$num_unprotected_disc_acts/sim_tib$num_disc_acts, na.rm = T)

sim_targets <- sim_targets %>%
  mutate(prop_unprotected_discordant = num_unprotected_disc_acts/num_disc_acts)

sim_targets %>% ggplot(aes(x = time, y = prop_unprotected_discordant)) +
  geom_line() +
  theme_minimal() +
  ggtitle("Plot 10: Proportion of HIV-discordant Acts that are Unprotected")

sim_targets %>% ggplot(aes(x = time, y = num_unprotected_disc_acts)) +
  geom_line() +
  theme_minimal() +
  ggtitle("Plot 11: Number of HIV-discordant Acts that are Unprotected")

sim_targets <- sim_targets %>%
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

sim_targets %>%
  ggplot(aes(x = time, y = cc.dx)) +
  geom_line() +
  geom_hline(aes(yintercept = 0.814, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 12: Proportion of HIV+ that are Diagnosed")

sim_targets %>%
  select(time, cc.dx, i.prev) %>%
  tidyr::pivot_longer(cols = c('cc.dx', 'i.prev'),
                      names_to = "Variable",
                      values_to = "Prop") %>%
  ggplot(aes(x = time, y = Prop, color = Variable)) +
  geom_line() +
  theme_minimal() +
  ggtitle("Plot 13: Proportion of HIV+ that are Diagnosed x Infected Prevalence")


sim_targets <- sim_targets %>%
  mutate(cc.vsupp.B = ifelse(is.nan(cc.vsupp.B), 0, cc.vsupp.B),
         cc.vsupp.H = ifelse(is.nan(cc.vsupp.H), 0, cc.vsupp.H),
         cc.vsupp.O = ifelse(is.nan(cc.vsupp.O), 0, cc.vsupp.O),
         cc.vsupp.W = ifelse(is.nan(cc.vsupp.W), 0, cc.vsupp.W),

         num_supp.B = num_diagnosed.B*cc.vsupp.B,
         num_supp.H = num_diagnosed.H*cc.vsupp.H,
         num_supp.O = num_diagnosed.O*cc.vsupp.O,
         num_supp.W = num_diagnosed.W*cc.vsupp.W,
         num_supp = num_supp.B + num_supp.H + num_supp.O
         + num_supp.W,
         cc.vsupp = num_supp/num_diagnosed)

sim_targets %>%
  ggplot(aes(x = time, y = cc.vsupp)) +
  geom_line() +
  geom_hline(aes(yintercept = 0.617, color = "red")) +
  theme_minimal() +
  ggtitle("Plot 14: Proportion of Diagnosed with Viral Suppression")

sim_targets %>%
  select(time, cc.dx, i.prev, cc.vsupp) %>%
  tidyr::pivot_longer(cols = c('cc.dx', 'i.prev', "cc.vsupp"),
                      names_to = "Variable",
                      values_to = "Prop") %>%
  ggplot(aes(x = time, y = Prop, color = Variable)) +
  geom_line() +
  theme_minimal() +
  ggtitle("Plot 15: Proportion of HIV+ that are Diagnosed x Infected Prevalence x \nProportion Viral Suppressed") +
  geom_hline(aes(yintercept = 0.617, color = "cc.vsupp Benchmark")) +
  geom_hline(aes(yintercept = .814, color = "cc.dx Benchmark"))


sim_targets %>%
  select(time, i.prev, i.prev.B, i.prev.H, i.prev.O, i.prev.W) %>%
  tidyr::pivot_longer(cols = c('i.prev', "i.prev.B", "i.prev.H", "i.prev.O", "i.prev.W"),
                      names_to = "Variable",
                      values_to = "Prop") %>%
  ggplot(aes(x = time, y = Prop, color = Variable)) +
  geom_line() +
  theme_minimal() +
  ggtitle("Plot 16: HIV Prevalence, by Race")

sim_targets %>%
  select(time, i.num, i.num.B, i.num.H, i.num.O, i.num.W) %>%
  tidyr::pivot_longer(cols = c('i.num', "i.num.B", "i.num.H", "i.num.O", "i.num.W"),
                      names_to = "Variable",
                      values_to = "Prop") %>%
  ggplot(aes(x = time, y = Prop, color = Variable)) +
  geom_line() +
  theme_minimal() +
  ggtitle("Plot 17: Number Infected, by Race")


sim_targets %>%
  select(time, n_edges_main.B, n_edges_main.H, n_edges_main.O, n_edges_main.W, n_edges_main.RHet, n_edges_main) %>%
  tidyr::pivot_longer(cols = c("n_edges_main.B", "n_edges_main.H", "n_edges_main.O", "n_edges_main.W", "n_edges_main.RHet", "n_edges_main"),
                      names_to = "Variable",
                      values_to = "Prop") %>%
  ggplot(aes(x = time, y = Prop, color = Variable)) +
  geom_line() +
  theme_minimal() +
  ggtitle("Plot 18: Number of Edges, Main Partnerships")

sim_targets %>%
  select(time, n_edges_casual.B, n_edges_casual.H, n_edges_casual.O, n_edges_casual.W, n_edges_casual.RHet, n_edges_casual) %>%
  tidyr::pivot_longer(cols = c("n_edges_casual.B", "n_edges_casual.H", "n_edges_casual.O", "n_edges_casual.W", "n_edges_casual.RHet", "n_edges_casual"),
                      names_to = "Variable",
                      values_to = "Prop") %>%
  ggplot(aes(x = time, y = Prop, color = Variable)) +
  geom_line() +
  theme_minimal() +
  ggtitle("Plot 19: Number of Edges, Casual Partnerships")

sim_targets %>%
  select(time, n_edges_onetime.B, n_edges_onetime.H, n_edges_onetime.O, n_edges_onetime.W, n_edges_onetime.RHet, n_edges_onetime) %>%
  tidyr::pivot_longer(cols = c("n_edges_onetime.B", "n_edges_onetime.H", "n_edges_onetime.O", "n_edges_onetime.W", "n_edges_onetime.RHet", "n_edges_onetime"),
                      names_to = "Variable",
                      values_to = "Prop") %>%
  ggplot(aes(x = time, y = Prop, color = Variable)) +
  geom_line() +
  theme_minimal() +
  ggtitle("Plot 20: Number of Edges, One-Time Partnerships")

sim_targets %>%
  mutate(n_edges_total = n_edges_main + n_edges_casual + n_edges_onetime) %>%
  select(time, n_edges_total.B, n_edges_total.H, n_edges_total.O, n_edges_total.W, n_edges_total.RHet, n_edges_total) %>%
  tidyr::pivot_longer(cols = c("n_edges_total.B", "n_edges_total.H", "n_edges_total.O", "n_edges_total.W", "n_edges_total.RHet", "n_edges_total"),
                      names_to = "Variable",
                      values_to = "Prop") %>%
  ggplot(aes(x = time, y = Prop, color = Variable)) +
  geom_line() +
  theme_minimal() +
  ggtitle("Plot 21: Number of Edges, Total Partnerships")


sim_targets <- sim_targets %>%
  mutate(n_edges_main.RHom = n_edges_main.B + n_edges_main.H + n_edges_main.O + n_edges_main.W,
         n_edges_casual.RHom = n_edges_casual.B + n_edges_casual.H + n_edges_casual.O + n_edges_casual.W,
         n_edges_onetime.RHom = n_edges_onetime.B + n_edges_onetime.H + n_edges_onetime.O + n_edges_onetime.W,
         n_edges_total.RHom = n_edges_total.B + n_edges_total.H + n_edges_total.O + n_edges_total.W,
         n_edges_total = n_edges_main + n_edges_casual + n_edges_onetime,

         prop_edges_main.RHom = n_edges_main.RHom/n_edges_main,
         prop_edges_casual.RHom = n_edges_casual.RHom/n_edges_casual,
         prop_edges_onetime.RHom = n_edges_onetime.RHom/n_edges_onetime,
         prop_edges_total.RHom = n_edges_total.RHom/n_edges_total)

sim_targets %>%
  select(time, prop_edges_main.RHom, prop_edges_casual.RHom, prop_edges_onetime.RHom, prop_edges_total.RHom) %>%
  tidyr::pivot_longer(cols = starts_with("prop"), names_to = "Variable", values_to = "Prop") %>%
  ggplot(aes(x = time, y = Prop, color = Variable)) +
  geom_line() +
  theme_minimal()


netstats$main$nodematch_race/netstats$main$edges
netstats$casl$nodematch_race/netstats$casl$edges
netstats$inst$nodematch_race/netstats$inst$edges



sim_targets$i.prev

# Run in debug mode, more details and examples here:
# https://github.com/EpiModel/EpiModeling/wiki/Writing-and-Debugging-EpiModel-Code
debug(hivtrans_msm)
sim <- netsim(est, param, init, control)
undebug(hivtrans_msm)
