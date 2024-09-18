# Libraries --------------------------------------------------------------------
library("tidyverse")
library("future.apply")
library("EpiModelHIV")

# Data ------------------------------------------------------------------------
#sim_targets <- readRDS("./data/intermediate/alldata_oct26.rds")
netstats <- readRDS("~/Desktop/chiSTIG_hpc/data/intermediate/estimates/netstats-level-local_00018_00014_alt.rds")

# Necessary files
context <- "local"
est_dir <- "blah"
prep_start = 52*2
source("~/Desktop/chiSTIG_hpc/R/utils-chistig_basic_inputs.R") # generate `path_to_est`, `param` and `initchis`
# path_to_est <- "./data/intermediate/estimates/basic_netest-local.rds"
# path_to_est      <- "/Users/wms1212/Desktop/ChiSTIG_model/epimodel/data/intermediate/estimates/venue_only_netest-local.rds"
path_to_est <- "~/Desktop/chiSTIG_hpc/data/intermediate/estimates/basic_netest-local.rds"
# Controls
source("~/Desktop/chiSTIG_hpc/R/utils-targets.R")


# Directory containing runs
this_dir <- "~/Desktop/sept12_runs"

files <- list.files(this_dir)

for (i in 1:length(files)) {
  this_targets <- tibble::as_tibble(readRDS(paste(this_dir, files[[i]], sep = "/"))) %>%
    mutate_calibration_targets() %>%
    mutate(i.prev.disp.BW = i.prev.B - i.prev.W,
           i.prev.disp.HW = i.prev.H - i.prev.W,
           cc.dx.B = ifelse(is.nan(cc.dx.B), 0, cc.dx.B),
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
           sim = i,
           treat = case_when(str_detect(files[[i]], "^apps_") ~ "Apps",
                                 str_detect(files[[i]], "^both_") ~ "Venues and Apps",
                                 str_detect(files[[i]], "^control_") ~ "Control",
                                 str_detect(files[[i]], "^venues_") ~ "Venues",
                                 TRUE ~ NA)#,
           #target_shift = stringr::str_detect(files[[i]], "plus")
           )

  if (i == 1) {
    sim_targets <- this_targets
  } else {
    sim_targets <- dplyr::bind_rows(sim_targets, this_targets)
  }
}


# Annualized Incidence Counts --------------------------------------------------

mean_incid <-  sim_targets %>%
  # mutate(sim = treat) %>%
  mutate(total.incid.B = incid.B + exo.incid.B,
         total.incid.H = incid.H + exo.incid.H,
         total.incid.O = incid.O + exo.incid.O,
         total.incid.W = incid.W + exo.incid.W) %>%
  group_by(treat, time) %>%
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
    filter(treat == this_row$treat)
  sums <- as.data.frame(t(colSums(past_year[,3:ncol(past_year)])))
  sums$treat <- this_row$treat
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
    filter(treat == this_row$treat)
  means <- as.data.frame(t(colMeans(past_year[,3:ncol(past_year)])))
  means$treat <- this_row$treat
  means$time <- this_row$time

  if (j == 1) {
    annual_incid2 <- means
  } else {
    annual_incid2 <- dplyr::bind_rows(annual_incid2, means)
  }
}

annual_incid2$total.incid.rate.dispar.BW <- annual_incid2$total.incid.B - annual_incid2$total.incid.W
annual_incid2$total.incid.rate.dispar.HW <- annual_incid2$total.incid.H - annual_incid2$total.incid.W

# sim_targets <- sim_targets %>% dplyr::filter(time >= 3000)
# annual_incid <- annual_incid %>% dplyr::filter(time >= 3000)
# annual_incid2 <- annual_incid2 %>% dplyr::filter(time >= 3000)

sim_targets %>%
  ggplot(aes(x = time, y = mean_deg_main, color = treat)) +
  geom_line(alpha = .4) +
  geom_hline(yintercept = 0.200339524, linetype = "dotted") +
  geom_hline(yintercept = 0.239467091) +
  geom_hline(yintercept = 0.278594658, linetype = "dotted") +
  labs(title = "Mean Degree, Main Partnerships",
       x = "Time",
       y = NULL) +
  theme_classic() #+
  # facet_grid(cols = vars(target_shift))


sim_targets %>%
  ggplot(aes(x = time, y = mean_deg_cas, color = treat)) +
  geom_line(alpha = .4) +
  geom_hline(yintercept = 0.271944955, linetype = "dotted") +
  geom_hline(yintercept = 0.339925924) +
  geom_hline(yintercept = 0.407906894, linetype = "dotted") +
  labs(title = "Mean Degree, Casual Partnerships",
       x = "Time",
       y = NULL) +
  theme_classic() # +
  # facet_grid(cols = vars(target_shift))


sim_targets %>%
  ggplot(aes(x = time, y = n_edges_onetime, color = treat)) +
  geom_line(alpha = .4) +
  # geom_hline(yintercept = 0.271944955, linetype = "dotted") +
  # geom_hline(yintercept = 0.339925924) +
  # geom_hline(yintercept = 0.407906894, linetype = "dotted") +
  labs(title = "Number of Edges, Onetime Partnerships",
       x = "Time",
       y = NULL) +
  theme_classic() #+
  # facet_grid(cols = vars(target_shift))
