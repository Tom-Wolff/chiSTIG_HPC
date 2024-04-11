# Libraries --------------------------------------------------------------------
library("tidyverse")
library("future.apply")
library("EpiModelHIV")

# Data ------------------------------------------------------------------------
#sim_targets <- readRDS("./data/intermediate/alldata_oct26.rds")
netstats <- readRDS("./data/intermediate/estimates/netstats-novenues-local.rds")

# Necessary files
context <- "local"
est_dir <- "blah"
prep_start = 52*2
source("./R/utils-chistig_basic_inputs.R") # generate `path_to_est`, `param` and `initchis`
# path_to_est <- "./data/intermediate/estimates/basic_netest-local.rds"
# path_to_est      <- "/Users/wms1212/Desktop/ChiSTIG_model/epimodel/data/intermediate/estimates/venue_only_netest-local.rds"
path_to_est <- "./data/intermediate/estimates/basic_netest-local.rds"
# Controls
source("./R/utils-targets.R")

#

load_diag <- function(this_dir, nsim = 1) {

  for (i in 1:nsim) {

    sim_dir <- paste("/data/intermediate/calibration/sim__", i, "__1.rds", sep = "")

    # this_targets <- readRDS(paste(this_dir, sim_dir, sep = ""))

    this_targets <- tibble::as_tibble(readRDS(paste(this_dir, sim_dir, sep = ""))) %>%
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
             sim = sim_dir)

    if (i == 1) {
      sim_targets <- this_targets
    } else {
      sim_targets <- dplyr::bind_rows(sim_targets, this_targets)
    }
  }

  return(sim_targets)
}

sim_targets <- load_diag(this_dir = getwd(), nsim = 4)

# Custom function for generating summary plots ---------------------------------
target_plot <- function(data, group, var, benchmark = NULL, title = NULL,
                        target_range = NULL) {

  # Create placeholder of variable we need
  data2 <- data[, c("time", group, var)]
  colnames(data2) <- c("time", "venues_treat", "this_var")

  if (is.null(target_range)) {

        this_plot <- data2 %>%
          group_by(venues_treat, time) %>%
          summarize(this_var = median(this_var)) %>%
          dplyr::ungroup() %>%
          ggplot(aes(x = time, y = this_var, color = as.factor(venues_treat))) +
          geom_line() +
          theme_minimal() +
          labs(y = var)

  } else {
        this_plot <- data2 %>%
          group_by(venues_treat, time) %>%
          summarize(this_var = median(this_var)) %>%
          dplyr::ungroup() %>%
          ggplot(aes(x = time, y = this_var, color = as.factor(venues_treat))) +
          geom_ribbon(aes(x = time, ymin = sort(target_range)[[1]], ymax = sort(target_range)[[2]]), fill = "grey", alpha = 0.15, linetype = 0) +
          geom_line() +
          theme_minimal() +
          labs(y = var)
  }

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


# Annualized Incidence Counts --------------------------------------------------

mean_incid <-  sim_targets %>%
  # mutate(sim = treat) %>%
  mutate(total.incid.B = incid.B + exo.incid.B,
         total.incid.H = incid.H + exo.incid.H,
         total.incid.O = incid.O + exo.incid.O,
         total.incid.W = incid.W + exo.incid.W) %>%
  group_by(sim, time) %>%
  summarize(incid.B = mean(incid.B),
            exo.incid.B = mean(exo.incid.B),
            total.incid.B = mean(total.incid.B),
            endo.ir100.B = mean(endo.ir100.B),
            exo.ir100.B = mean(exo.ir100.B),
            ir100.B = mean(ir100.B),

            incid.H = mean(incid.H),
            exo.incid.H = mean(exo.incid.H),
            total.incid.H = mean(total.incid.H),
            endo.ir100.H = mean(endo.ir100.H),
            exo.ir100.H = mean(exo.ir100.H),
            ir100.H = mean(ir100.H),

            incid.O = mean(incid.O),
            exo.incid.O = mean(exo.incid.O),
            total.incid.O = mean(total.incid.O),
            endo.ir100.O = mean(endo.ir100.O),
            exo.ir100.O = mean(exo.ir100.O),
            ir100.O = mean(ir100.O),

            incid.W = mean(incid.W),
            exo.incid.W = mean(exo.incid.W),
            total.incid.W = mean(total.incid.W),
            endo.ir100.W = mean(endo.ir100.W),
            exo.ir100.W = mean(exo.ir100.W),
            ir100.W = mean(ir100.W),



            ) %>%
  ungroup()


for (j in 1:nrow(mean_incid)) {
  this_row <- mean_incid[j,]
  past_year <- mean_incid %>%
    filter(time <= this_row$time & time > (this_row$time-52)) %>%
    filter(sim == this_row$sim)
  sums <- as.data.frame(t(colSums(past_year[,3:ncol(past_year)])))
  sums$sim <- this_row$sim
  sums$time <- this_row$time

  if (j == 1) {
    annual_incid <- sums
  } else {
    annual_incid <- dplyr::bind_rows(annual_incid, sums)
  }
}

for (j in 1:nrow(mean_incid)) {
  this_row <- mean_incid[j,]
  past_year <- mean_incid %>%
    filter(time <= this_row$time & time > (this_row$time-52)) %>%
    filter(sim == this_row$sim)
  means <- as.data.frame(t(colMeans(past_year[,3:ncol(past_year)])))
  means$sim <- this_row$sim
  means$time <- this_row$time

  if (j == 1) {
    annual_incid2 <- means
  } else {
    annual_incid2 <- dplyr::bind_rows(annual_incid2, means)
  }
}

# Filter for last year
# sim_targets <- sim_targets %>% filter(time > 3000)

# Settings ---------------------------------------------------------------------
source("./R/utils-0_project_settings.R")

# ------------------------------------------------------------------------------
# Necessary files
# source("R/utils-chistig_basic_inputs.R") # generate `path_to_est`, `param` and `init`
# path_to_est <- "./data/intermediate/estimates/basic_netest-local.rds"
# source("./R/utils-targets.R")

sim_targets <- sim_targets %>% dplyr::filter(time >= 3000)
annual_incid <- annual_incid %>% dplyr::filter(time >= 3000)
annual_incid2 <- annual_incid2 %>% dplyr::filter(time >= 3000)

library(tidyverse)
# Population Size --------------------------------------------------------------
i = 1

target_plot(data = sim_targets,
            group = "sim",
            var = "num",
            benchmark = 11612,
            title = paste("Plot ", i, ": Population Size", sep = ""))


# Proportion HIV+ Diagnosed ----------------------------------------------------

i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.dx",
            benchmark = 0.55333,
            title = paste("Plot ", i, ": Proportion of HIV+ that are Diagnosed", sep = ""))

i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.dx.B",
            benchmark = 0.546535643,
            title = paste("Plot ", i, ": Proportion of HIV+ that are Diagnosed (Black)", sep = ""))


i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.dx.H",
            benchmark = 0.5431367893,
            title = paste("Plot ", i, ": Proportion of HIV+ that are Diagnosed (Hispanic)", sep = ""))

i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.dx.O",
            benchmark = 0.5601310,
            title = paste("Plot ", i, ": Proportion of HIV+ that are Diagnosed (Other)", sep = ""))

i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.dx.W",
            benchmark = 0.5988779867,
            title = paste("Plot ", i, ": Proportion of HIV+ that are Diagnosed (White)", sep = ""))



# Proportion HIV+ Linked to Care in 1st Month ----------------------------------

i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.linked1m.B",
            benchmark = .828,
            title = paste("Plot ", i, ": Proportion of HIV+ Nodes Linked to Care within One Month (Black)", sep = ""))


i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.linked1m.H",
            benchmark = 0.867,
            title = paste("Plot ", i, ": Proportion of HIV+ Nodes Linked to Care within One Month (Hispanic)", sep = ""))

i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.linked1m.O",
            benchmark = 0.875,
            title = paste("Plot ", i, ": Proportion of HIV+ Nodes Linked to Care within One Month (Other)", sep = ""))

i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.linked1m.W",
            benchmark = 0.936,
            title = paste("Plot ", i, ": Proportion of HIV+ Nodes Linked to Care within One Month (White)", sep = ""))


# Proportion HIV+ Linked to Care in 1st Month ----------------------------------


i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.vsupp.B",
            benchmark = 0.571,
            title = paste("Plot ", i, ": Proportion of HIV+ Nodes with Viral Suppression (Black)", sep = ""))


i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.vsupp.H",
            benchmark = 0.675,
            title = paste("Plot ", i, ": Proportion of HIV+ Nodes with Viral Suppression (Hispanic)", sep = ""))

i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.vsupp.O",
            benchmark = 0.586,
            title = paste("Plot ", i, ": Proportion of HIV+ Nodes with Viral Suppression (Other)", sep = ""))

i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.vsupp.W",
            benchmark = 0.617,
            title = paste("Plot ", i, ": Proportion of HIV+ Nodes with Viral Suppression (White)", sep = ""))


# Proportion of Indicated MSM currently using PrEP -----------------------------

i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.prep.B",
            benchmark = 0.350,
            title = paste("Plot ", i, ": Proportion of Indicated MSM currently using PrEP (Black)", sep = ""))

i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.prep.H",
            benchmark = 0.386,
            title = paste("Plot ", i, ": Proportion of Indicated MSM currently using PrEP (Hispanic)", sep = ""))

i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.prep.O",
            benchmark = 0.357,
            title = paste("Plot ", i, ": Proportion of Indicated MSM currently using PrEP (Other)", sep = ""))

i <- i+1

target_plot(data = sim_targets,
            group = "sim",
            var = "cc.prep.W",
            benchmark = 0.368,
            title = paste("Plot ", i, ": Proportion of Indicated MSM currently using PrEP (White)", sep = ""))


# Exogenous Infections ---------------------------------------------------------
# i <- i+1
# target_plot(data = annual_incid,
#             var = "exo.incid.B",
#             group = "sim",
#             benchmark = (131/.75)*(.28),
#             title = paste("Plot ", i, ": Annual Exogenous Infections (Black)", sep = ""))
# i <- i+1
# target_plot(data = annual_incid,
#             var = "exo.incid.H",
#             group = "sim",
#             benchmark = (47/.75)*(.40),
#             title = paste("Plot ", i, ": Annual Exogenous Infections (Hispanic)", sep = ""))
# i <- i+1
# target_plot(data = annual_incid,
#             var = "exo.incid.O",
#             group = "sim",
#             benchmark = (14/.75)*(.37),
#             title = paste("Plot ", i, ": Annual Exogenous Infections (Other)", sep = ""))
# i <- i+1
# target_plot(data = annual_incid,
#             var = "exo.incid.W",
#             group = "sim",
#             benchmark = (19/.8)*(.44),
#             title = paste("Plot ", i, ": Annual Exogenous Infections (White)", sep = ""))

# Exogenous Incidence Rate -----------------------------------------------------
i <- i+1
target_plot(data = sim_targets,
            var = "exo.ir100.B",
            group = "sim",
            benchmark = mean(c(1.438, 1.798)),
            title = paste("Plot ", i, ": Exogenous Incidence Rate (Black)", sep = ""))
i <- i+1
target_plot(data = sim_targets,
            var = "exo.ir100.H",
            group = "sim",
            benchmark = mean(c(0.653, 0.816)),
            title = paste("Plot ", i, ": Exogenous Incidence Rate (Hispanic)", sep = ""))
i <- i+1
target_plot(data = sim_targets,
            var = "exo.incid.O",
            group = "sim",
            benchmark = mean(c(0.506, 0.633)),
            title = paste("Plot ", i, ": Exogenous Incidence Rate (Other)", sep = ""))
i <- i+1
target_plot(data = sim_targets,
            var = "exo.incid.W",
            group = "sim",
            benchmark = mean(c(0.257, 0.3212)),
            title = paste("Plot ", i, ": Exogenous Incidence Rate (White)", sep = ""))

# Annualized Exogenous Incidence Rate -----------------------------------------
i <- i+1
target_plot(data = annual_incid2,
            var = "exo.ir100.B",
            group = "sim",
            benchmark = mean(c(1.438, 1.798)),
            target_range = c(1.438, 1.798),
            title = paste("Plot ", i, ": Exogenous Incidence Rate (Black, Annualized)", sep = ""))

i <- i+1
target_plot(data = annual_incid2,
            var = "exo.ir100.H",
            group = "sim",
            benchmark = mean(c(0.653, 0.816)),
            target_range = c(0.653, 0.816),
            title = paste("Plot ", i, ": Exogenous Incidence Rate (Hispanic, Annualized)", sep = ""))
i <- i+1
target_plot(data = annual_incid2,
            var = "exo.ir100.O",
            group = "sim",
            benchmark = mean(c(0.506, 0.633)),
            target_range = c(0.506, 0.633),
            title = paste("Plot ", i, ": Exogenous Incidence Rate (Other, Annualized)", sep = ""))
i <- i+1
target_plot(data = annual_incid2,
            var = "exo.ir100.W",
            group = "sim",
            benchmark = mean(c(0.257, 0.3212)),
            target_range = c(0.257, 0.3212),
            title = paste("Plot ", i, ": Exogenous Incidence Rate (White, Annualized)", sep = ""))


# Endogenous Infections --------------------------------------------------------
# i <- i+1
# target_plot(data = annual_incid,
#             var = "incid.B",
#             group = "sim",
#             benchmark = (131/.75)*(1-.28),
#             title = paste("Plot ", i, ": Annual Endogenous Infections (Black)", sep = ""))
# i <- i+1
# target_plot(data = annual_incid,
#             var = "incid.H",
#             group = "sim",
#             benchmark = (47/.75)*(1-.40),
#             title = paste("Plot ", i, ": Annual Endogenous Infections (Hispanic)", sep = ""))
# i <- i+1
# target_plot(data = annual_incid,
#             var = "incid.O",
#             group = "sim",
#             benchmark = (14/.75)*(1-.37),
#             title = paste("Plot ", i, ": Annual Endogenous Infections (Other)", sep = ""))
# i <- i+1
# target_plot(data = annual_incid,
#             var = "incid.W",
#             group = "sim",
#             benchmark = (19/.8)*(1-.44),
#             title = paste("Plot ", i, ": Annual Endogenous Infections (White)", sep = ""))

# Endogenous Incidence Rate ----------------------------------------------------
i <- i+1
target_plot(data = sim_targets,
            var = "endo.ir100.B",
            group = "sim",
            benchmark = 6.42 - mean(c(1.438, 1.798)),
            title = paste("Plot ", i, ": Endogenous Incidence Rate (Black)", sep = ""))
i <- i+1
target_plot(data = sim_targets,
            var = "endo.ir100.H",
            group = "sim",
            benchmark = 2.04 - mean(c(0.653, 0.816)),
            title = paste("Plot ", i, ": Endogenous Incidence Rate (Hispanic)", sep = ""))
i <- i+1
target_plot(data = sim_targets,
            var = "endo.ir100.O",
            group = "sim",
            benchmark = 1.71 - mean(c(0.506, 0.633)),
            title = paste("Plot ", i, ": Endogenous Incidence Rate (Other)", sep = ""))
i <- i+1
target_plot(data = sim_targets,
            var = "endo.ir100.W",
            group = "sim",
            benchmark = 0.73 - mean(c(0.257, 0.3212)),
            title = paste("Plot ", i, ": Endogenous Incidence Rate (White)", sep = ""))

# Annualized Endogenous Incidence Rate -----------------------------------------
i <- i+1
target_plot(data = annual_incid2,
            var = "endo.ir100.B",
            group = "sim",
            benchmark = 6.42 - mean(c(1.438, 1.798)),
            target_range = c(4.44-1.438, 9.30-1.798),
            title = paste("Plot ", i, ": Endogenous Incidence Rate (Black, Annualized)", sep = ""))
i <- i+1
target_plot(data = annual_incid2,
            var = "endo.ir100.H",
            group = "sim",
            benchmark = 2.04 - mean(c(0.653, 0.816)),
            target_range = c(1.10-0.653, 3.79-0.816),
            title = paste("Plot ", i, ": Endogenous Incidence Rate (Hispanic, Annualized)", sep = ""))
i <- i+1
target_plot(data = annual_incid2,
            var = "endo.ir100.O",
            group = "sim",
            benchmark = 1.71 - mean(c(0.506, 0.633)),
            target_range = c(0.55-0.506, 5.31-0.633),
            title = paste("Plot ", i, ": Endogenous Incidence Rate (Other, Annualized)", sep = ""))
i <- i+1
target_plot(data = annual_incid2,
            var = "endo.ir100.W",
            group = "sim",
            benchmark = 0.73 - mean(c(0.257, 0.3212)),
            target_range = c(0.24, 2.26-0.3212),
            title = paste("Plot ", i, ": Endogenous Incidence Rate (White, Annualized)", sep = ""))


# Total Incidence Rate (Exogenous + Endogenous) --------------------------------
i <- i+1
target_plot(data = sim_targets,
            var = "ir100.B",
            group = "sim",
            benchmark = mean(6.42),
            title = paste("Plot ", i, ": Endogenous Incidence Rate (Black)", sep = ""))
i <- i+1
target_plot(data = sim_targets,
            var = "ir100.H",
            group = "sim",
            benchmark = mean(2.04),
            title = paste("Plot ", i, ": Endogenous Incidence Rate (Hispanic)", sep = ""))
i <- i+1
target_plot(data = sim_targets,
            var = "ir100.O",
            group = "sim",
            benchmark = mean(1.71),
            title = paste("Plot ", i, ": Endogenous Incidence Rate (Other)", sep = ""))

i <- i+1
target_plot(data = sim_targets,
            var = "ir100.W",
            group = "sim",
            benchmark = mean(0.73),
            title = paste("Plot ", i, ": Endogenous Incidence Rate (White)", sep = ""))

# Annualized Total Incidence Rate ----------------------------------------------
i <- i+1
target_plot(data = annual_incid2,
            var = "ir100.B",
            group = "sim",
            benchmark = 6.42,
            target_range = c(4.44, 9.30),
            title = paste("Plot ", i, ": Total Incidence Rate (Black, Annualized)", sep = ""))
i <- i+1
target_plot(data = annual_incid2,
            var = "ir100.H",
            group = "sim",
            benchmark = 2.04,
            target_range = c(1.10, 3.79),
            title = paste("Plot ", i, ": Total Incidence Rate (Hispanic, Annualized)", sep = ""))
i <- i+1
target_plot(data = annual_incid2,
            var = "ir100.O",
            group = "sim",
            benchmark = 1.71,
            target_range = c(0.55, 5.31),
            title = paste("Plot ", i, ": Total Incidence Rate (Other, Annualized)", sep = ""))
i <- i+1
target_plot(data = annual_incid2,
            var = "ir100.W",
            group = "sim",
            benchmark = 0.73,
            target_range = c(0.24, 2.26),
            title = paste("Plot ", i, ": Total Incidence Rate (White, Annualized)", sep = ""))


# Prevalence -------------------------------------------------------------------
se_prop <- function(p, n) {
  se <- sqrt((p*(1-p))/n)
  return(c((p-1.95*se), (p+1.95*se)))
}

i <- i+1
target_plot(data = sim_targets,
            var = "i.prev",
            group = "sim",
            benchmark = 0.165,
            target_range = se_prop(.165, 1015),
            title = paste("Plot ", i, ": Prevalence, Infected", sep = ""))

i <- i+1
target_plot(data = sim_targets,
            var = "i.prev.B",
            group = "sim",
            benchmark = 0.32,
            target_range = se_prop(.32, 244),
            title = paste("Plot ", i, ": Prevalence, Infected (Black)", sep = ""))

i <- i+1
target_plot(data = sim_targets,
            var = "i.prev.H",
            group = "sim",
            benchmark = 0.125,
            target_range = se_prop(.125, 304),
            title = paste("Plot ", i, ": Prevalence, Infected (Hispanic)", sep = ""))

i <- i+1
target_plot(data = sim_targets,
            var = "i.prev.O",
            group = "sim",
            benchmark = 0.122,
            target_range = se_prop(.122, 115),
            title = paste("Plot ", i, ": Prevalence, Infected (Other)", sep = ""))

i <- i+1
target_plot(data = sim_targets,
            var = "i.prev.W",
            group = "sim",
            benchmark = .02,
            target_range = se_prop(.02, 252),
            title = paste("Plot ", i, ": Prevalence, Infected (White)", sep = ""))

