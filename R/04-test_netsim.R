## Example interactive epidemic simulation run script with basic
## parameterization and all parameters defined in `param_msm`, with example of
## writing/debugging modules

# Libraries  -------------------------------------------------------------------
library("EpiModelHIV")
library("chiSTIGmodules")

# library("chiSTIGmodules")
# devtools::load_all("~/Desktop/chiSTIGmodules")

# Settings ---------------------------------------------------------------------
source("R/utils-0_project_settings.R")
source("./R/utils-targets.R")

# Necessary files
epistats <- readRDS("data/intermediate/estimates/epistats-local.rds")
epistats$age.breaks <- c(16, 20, 30)
# age limits probably needs to be maxed at 31
epistats$age.limits <- c(16, 30)

netstats <- readRDS("data/intermediate/estimates/netstats-local.rds")
# `netstats` without venues target stats, in case this changes anything
netstats <- readRDS("data/intermediate/estimates/netstats-novenues-local.rds")


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
est <- readRDS("data/intermediate/estimates/basic_netest-local.rds")

# Is the aging out of the older initial nodes driving HIV extinction?
# Let's level out the age distribution and find out
 netstats$attr$age <- sample(16:29, length(netstats$attr$age), replace = TRUE)
 netstats$attr$age <- netstats$attr$age + sample(1:1000, length(netstats$attr$age), replace = TRUE)/1000




prep_start <- 52 * 2
param <- EpiModel::param.net(
  data.frame.params = readr::read_csv("data/input/params_chistig_feb19.csv"),
  netstats          = netstats,
  epistats          = epistats,
  prep.start        = Inf,
  riskh.start       = Inf
)

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
  #  nsteps = 400,
    nsteps = calibration_end + 520,
   # nsteps =  prep_start + 52 * 2,
  nsims = 10,
  ncores = 10,
  cumulative.edgelist = TRUE,
  truncate.el.cuml = 400,
  raw.output = TRUE,
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
  module.order = c("aging.FUN", "departure.FUN", "arrival.FUN", # "venues.FUN",
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

sim_list <- as.list(sim)

# Combine pertinent elements of `sim_list` to a single data frame
full_data <- dplyr::bind_rows(lapply(sim_list[[1]]$raw.records, function(x){x$object}))

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
