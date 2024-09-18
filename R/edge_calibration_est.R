library("EpiModelHIV")
library("ARTnet")

context <- "local"

# Extract the row from `drate_mat` that will preside over this run
# Read in the arguments from the commandline
args <- commandArgs(trailingOnly = TRUE)
sbatch_run_num <- args[1]


# Read in Epistats (update directory on HPC)
epistats <- readRDS("./data/intermediate/estimates/epistats-local.rds")

# Read in Target stat calibration inputs
drate_mat <- read.csv("./data/intermediate/estimates/edge_target_calibration_vals.csv")

# drate_mat$fit_no <- drate_mat$fit_no + 2
# drate_mat$venues_main <- drate_mat$venues_main + c(-50, -100)

drate_mat$fit_no <- paste("sept17", drate_mat$fit_no, sep = "_")


library(tidyverse)

#egos <- read.csv("./data/input/egos_v4_1.csv") %>%
egos <- read.csv("./data/synthpop_gen/egos_v4_1_uniform_age_dist.csv") %>%
  dplyr::mutate(race_art = dplyr::case_when(race_ethnicity == "whiteNH" ~ "white",
                                            race_ethnicity == "blackNH" ~ "black",
                                            race_ethnicity == "hispanic" ~ "hispanic",
                                            race_ethnicity == "otherNH" ~ "other",
                                            TRUE ~ NA),
                race_art2 = dplyr::case_when(race_ethnicity == "whiteNH" ~ "4_white",
                                             race_ethnicity == "blackNH" ~ "1_black",
                                             race_ethnicity == "hispanic" ~ "2_hispanic",
                                             race_ethnicity == "otherNH" ~ "3_other",
                                             TRUE ~ NA),
                race = dplyr::case_when(race_ethnicity == "whiteNH" ~ 4,
                                        race_ethnicity == "blackNH" ~ 1,
                                        race_ethnicity == "hispanic" ~ 2,
                                        race_ethnicity == "otherNH" ~ 3,
                                        TRUE ~ NA),
                age.grp = agegroup,
                age = age,
                deg.main = init_ser_cat,
                deg.casl = init_cas_cat,
                deg.tot = init_pers_cat,

                # venues_all = venue_list)

                venues_all = venue_list_1week)


#
# rel <- read.csv("./from_chistig/egos_rel_full.csv") %>%
#   dplyr::mutate(deg.main = init_ser_cat,
#                 deg.casl = init_cas_cat,
#                 deg.tot = init_pers_cat)
# venues <- read.csv("./from_chistig/egos_venues_full.csv")
# apps <- read.csv("./from_chistig/egos_apps_full.csv")
# hiv <- read.csv("./from_chistig/egos_hiv_full.csv")
#
# egos <- egos %>% left_join(venues, by=c("numeric_id"="egoid")) %>%
#   left_join(apps, by = c("numeric_id"="egoid")) %>%
#   left_join(rel, by = c("numeric_id"="egoid")) %>%
#   left_join(hiv, by = c("numeric_id" = "egoid"))

# Let's just select the variables we need so we have a clean inventory
# Also, ARTnet uses periods to replace spaces whereas we use underscores.
# I think the underscores we use might cause naming problems with ERGM fitting
# down the line, so I'm changing variable names here accordingly.
egos <- egos %>%
  mutate(sqrt.age = sqrt(age),
         active.sex = 1,
         age.grp = ifelse(age.grp == "16to20", 1, 2),
         apps.all = app_list) %>%
  select(numeric_id, egoid,
         age, sqrt.age, agegroup, age.grp,
         race.ethnicity = race_ethnicity, race,
         deg.casl, deg.main, deg.tot,
         # risk.grp ?
         diag.status = hiv_status,
         venues.all = venues_all,
         apps.all,
         active.sex)


##################################
# ASMR data frame (for netstats) #
##################################


## Age-sex-specific mortality rates (B, H, W)
#  in 1-year age decrements starting with age 1
#  from CDC NCHS Underlying Cause of Death database (for 2020)
asmr.B <- c(0.00079, 0.00046, 0.00030, 0.00025, 0.00024, 0.00025, 0.00019,
            0.00019, 0.00021, 0.00020, 0.00026, 0.00026, 0.00038, 0.00056,
            0.00077, 0.00100, 0.00151, 0.00227, 0.00271, 0.00264, 0.00297,
            0.00302, 0.00315, 0.00319, 0.00322, 0.00319, 0.00336, 0.00337,
            0.00330, 0.00363, 0.00396, 0.00392, 0.00407, 0.00428, 0.00411,
            0.00453, 0.00485, 0.00486, 0.00533, 0.00513, 0.00575, 0.00580,
            0.00628, 0.00671, 0.00669, 0.00750, 0.00773, 0.00858, 0.00934,
            0.00947, 0.00999, 0.01141, 0.01216, 0.01360, 0.01432, 0.01517,
            0.01699, 0.01853, 0.02021, 0.02099, 0.02366, 0.02547, 0.02877,
            0.02979, 0.03104, 0.03467, 0.03653, 0.03941, 0.04114, 0.04320,
            0.04487, 0.04879, 0.05100, 0.05678, 0.05611, 0.06384, 0.06891,
            0.07399, 0.07682, 0.08209, 0.08938, 0.09737, 0.10400, 0.11336,
            0.16336)
asmr.H <- c(0.00032, 0.00021, 0.00018, 0.00011, 0.00011, 0.00009, 0.00010,
            0.00009, 0.00009, 0.00012, 0.00013, 0.00015, 0.00016, 0.00025,
            0.00036, 0.00058, 0.00076, 0.00106, 0.00125, 0.00134, 0.00145,
            0.00156, 0.00164, 0.00166, 0.00164, 0.00159, 0.00176, 0.00172,
            0.00201, 0.00198, 0.00192, 0.00191, 0.00202, 0.00204, 0.00219,
            0.00223, 0.00251, 0.00246, 0.00272, 0.00272, 0.00298, 0.00307,
            0.00321, 0.00351, 0.00367, 0.00391, 0.00442, 0.00484, 0.00512,
            0.00521, 0.00616, 0.00649, 0.00714, 0.00790, 0.00863, 0.00938,
            0.00992, 0.01094, 0.01222, 0.01217, 0.01464, 0.01483, 0.01630,
            0.01731, 0.01850, 0.02054, 0.02269, 0.02321, 0.02515, 0.02734,
            0.02937, 0.03064, 0.03349, 0.03670, 0.03980, 0.04387, 0.04724,
            0.05151, 0.05591, 0.05902, 0.06345, 0.07317, 0.07849, 0.08617,
            0.13436)
asmr.W <- c(0.00034, 0.00023, 0.00019, 0.00014, 0.00014, 0.00010, 0.00010,
            0.00009, 0.00009, 0.00012, 0.00014, 0.00015, 0.00022, 0.00028,
            0.00036, 0.00050, 0.00059, 0.00082, 0.00096, 0.00104, 0.00126,
            0.00128, 0.00134, 0.00144, 0.00153, 0.00163, 0.00172, 0.00186,
            0.00194, 0.00205, 0.00220, 0.00225, 0.00238, 0.00245, 0.00247,
            0.00264, 0.00274, 0.00280, 0.00306, 0.00312, 0.00324, 0.00329,
            0.00344, 0.00354, 0.00371, 0.00405, 0.00442, 0.00479, 0.00511,
            0.00547, 0.00599, 0.00653, 0.00706, 0.00768, 0.00827, 0.00922,
            0.00978, 0.01065, 0.01151, 0.01235, 0.01349, 0.01437, 0.01548,
            0.01664, 0.01730, 0.01879, 0.01986, 0.02140, 0.02263, 0.02419,
            0.02646, 0.02895, 0.03031, 0.03625, 0.03753, 0.04268, 0.04631,
            0.05235, 0.05724, 0.06251, 0.06934, 0.07589, 0.08669, 0.09582,
            0.16601)
asmr.O <- c(0.00034, 0.00023, 0.00019, 0.00014, 0.00014, 0.00010, 0.00010,
            0.00009, 0.00009, 0.00012, 0.00014, 0.00015, 0.00022, 0.00028,
            0.00036, 0.00050, 0.00059, 0.00082, 0.00096, 0.00104, 0.00126,
            0.00128, 0.00134, 0.00144, 0.00153, 0.00163, 0.00172, 0.00186,
            0.00194, 0.00205, 0.00220, 0.00225, 0.00238, 0.00245, 0.00247,
            0.00264, 0.00274, 0.00280, 0.00306, 0.00312, 0.00324, 0.00329,
            0.00344, 0.00354, 0.00371, 0.00405, 0.00442, 0.00479, 0.00511,
            0.00547, 0.00599, 0.00653, 0.00706, 0.00768, 0.00827, 0.00922,
            0.00978, 0.01065, 0.01151, 0.01235, 0.01349, 0.01437, 0.01548,
            0.01664, 0.01730, 0.01879, 0.01986, 0.02140, 0.02263, 0.02419,
            0.02646, 0.02895, 0.03031, 0.03625, 0.03753, 0.04268, 0.04631,
            0.05235, 0.05724, 0.06251, 0.06934, 0.07589, 0.08669, 0.09582,
            0.16601)

# transformed to rates by time unit
trans.asmr.H <- 1 - (1 - asmr.H)^(1 / (364 / epistats$time.unit))
trans.asmr.W <- 1 - (1 - asmr.W)^(1 / (364 / epistats$time.unit))
trans.asmr.B <- 1 - (1 - asmr.B)^(1 / (364 / epistats$time.unit))
trans.asmr.O <- 1 - (1 - asmr.O)^(1 / (364 / epistats$time.unit))

# Transformed rates, 85+ rate for ages 85 - 100
vec.asmr.B <- c(trans.asmr.B, rep(tail(trans.asmr.B, n = 1), 15))
vec.asmr.H <- c(trans.asmr.H, rep(tail(trans.asmr.H, n = 1), 15))
vec.asmr.W <- c(trans.asmr.W, rep(tail(trans.asmr.W, n = 1), 15))
vec.asmr.O <- c(trans.asmr.O, rep(tail(trans.asmr.O, n = 1), 15))

asmr <- data.frame(age = 1:100,
                   vec.asmr.B,
                   vec.asmr.H,
                   vec.asmr.W,
                   vec.asmr.O)


# Setting deterministic mortality prob = 1 at upper age limit
max.age <- epistats$age.limits[2]
asmr[asmr$age >= max.age, ] <- 1


#######################################
# Read in target stats (for netstats) #
#######################################

old_target_df <- read.csv("./from_chistig/target_values_full.csv")
target_df <- read.csv("./data/synthpop_gen/target_values_v4_1_uniform_age_dist.csv")

# Function for quickly extracting target values from dataframe
target_extract = function(df = target_df, term, model) {
  this_row <- which(df$X == term)
  this_col <- which(colnames(df) == paste("mean_", model, sep = ""))
  target_val <- df[this_row, this_col]
  return(target_val)
}



# Get Duration/dissolution coefficients based on ARTNet data
source("./R/dur_coefs.R")
dur_coefs <- out

sbatch_run_num <- 2
est_dir <- "./data/intermediate/estimates/"

drate_mat <- drate_mat[sbatch_run_num, ]


# Netstats
netstats <- list(
  # demog : list of demographic information for network
  demog = list(

    # num : network size (nodes)
    num = nrow(egos),

    # props : proportion of nodes in each racial/ethnic category
    # Looks like it's a data frame
    props = data.frame("White" = sum(egos$race == 4)/nrow(egos),
                       "Black" = sum(egos$race == 1)/nrow(egos),
                       "Hispanic" = sum(egos$race == 2)/nrow(egos),
                       "Other" = sum(egos$race == 3)/nrow(egos)),

    # num.B : proportion of nodes black
    num.B = sum(egos$race == 1)/nrow(egos),

    # num.H : proportion of nodes hispanic
    num.H = sum(egos$race == 2)/nrow(egos),

    # num.W : proportion of nodes white/other (adjust for our own categorization
    # schema)
    num.W = sum(egos$race == 4)/nrow(egos),

    # num.O : proportion of nodes other race?
    num.O = sum(egos$race == 3)/nrow(egos),

    # asmr : dataframe containing 100 rows (possible age range) with age-specific
    # mortality rates
    ##### age (`1:100`)
    ##### vec.asmr.B (something black)
    ##### vec.asmr.H (something hispanic)
    ##### vec.asmr.W (somethign white/other; adjust for our own categorization)

    asmr = asmr,

    # ages : vector of valid age values in simulation
    ages = epistats$age.limits[[1]]:epistats$age.limits[[2]],

    # age.breaks : vector of categorical age cutoffs
    age.breaks = epistats$age.breaks),

  # geog.lvl : character of geographic level
  geog.lvl = NULL,

  # race : logical if things should be broken down by race
  race = TRUE,

  # time.unit : numeric value indicating number of days in each network step
  time.unit = 7,

  # attr : list of node-level attributs
  attr = list(

    # Original ego identifiers
    numeric.id = egos$numeric_id,
    egoid = egos$egoid,

    # age : numeric vector of node ages
    age = egos$age,

    # sqrt.age : square root of ages
    sqrt.age = egos$sqrt.age,

    # age.grp : numeric designation of age group membership
    age.grp = egos$age.grp,

    # active.sex : 1/0 indicator of if node is sexually active
    active.sex = egos$active.sex,

    # race : numeric designation of racial/ethnic categorization
    race = egos$race,

    # deg.casl : degree in casual network
    deg.casl = egos$deg.casl,

    # deg.main : degree in main network
    deg.main = egos$deg.main,

    # deg.tot : total degree
    deg.tot = egos$deg.tot,

    # risk.grp : risk group (investigate for what this means practically)

    # role.class : Preference for sexual position during acts
    # Randomly assign to nodes according to proportions derived from RADAR by
    # Morgan et al. (2021)
    role.class = sample(0:2, size = nrow(egos), replace = TRUE, prob = c(73, 87, 475)/sum(c(73, 87, 475))),
    # For now give everyone "versatile"
    # role.class = rep(2, nrow(egos)),

    # diag.status : I believe this is HIV status
    diag.status = egos$diag.status,

    # venues_all
    venues.all = egos$venues.all,

    # apps_all
    apps.all = egos$apps.all



  ),

  # main : list of target stats and dissolution model for main partnerships
  main = list(

    # edges
    edges = target_extract(df = target_df,
                           term = "edges",
                           model = "main"),

    # concurrent
    concurrent = target_extract(df = target_df,
                                term = "concurrent",
                                model = "main"),

    # nodefactor race
    nodefactor_race = c(target_extract(term = "nodefactor.race_ethnicity.blackNH", model = "main"),
                        target_extract(term = "nodefactor.race_ethnicity.hispanic", model = "main"),
                        target_extract(term = "nodefactor.race_ethnicity.otherNH", model = "main"),
                        target_extract(term = "nodefactor.race_ethnicity.whiteNH", model = "main")),


    # nodematch race
    nodematch_race = target_extract(term = "nodematch.race_ethnicity", model = "main"),

    # nodematch black
    nodematch_race.1 = target_extract(term = "nodematch.race_ethnicity.blackNH", model = "main"),

    # nodematch_age.grp
    # nodematch_age.grp = target_extract(term = "nodematch.age", model = "main"),
    nodematch_age.grp = c(target_extract(term = "nodematch.age.16to20", model = "main"),
                          target_extract(term = "nodematch.age.21to29", model = "main")),

    # nodefactor_age.grp
    # nodefactor_age.grp = c(target_extract(term = "nodefactor.age.16to20", model = "main"),
    #                        target_extract(term = "nodefactor.age.21to29", model = "main")),

    # nodefactor_init_cas_cat
    nodefactor_deg.casl = c(target_extract(term = "nodefactor.init_cas_cat.0", model = "main"),
                            target_extract(term = "nodefactor.init_cas_cat.1", model = "main"),
                            target_extract(term = "nodefactor.init_cas_cat.2+", model = "main")),

    # fuzzynodematch_venues_all
    # fuzzynodematch_venues.all = target_extract(term = "fuzzynodematch.venues_all.TRUE", model = "main"),
    fuzzynodematch_venues.all = drate_mat$venues_main,

    # fuzzynodematch_apps_all
    fuzzynodematch_apps.all = drate_mat$apps_main,


    # fuzzynodematch_apps_dating
    fuzzynodematch_apps.dating = target_extract(term = "fuzzynodematch.apps_dating.TRUE", model = "main"),

    # dissolution model
    dissolution = dissolution_coefs(~offset(edges), duration = 87, d.rate = drate_mat$drate_main),
    diss.homog = dissolution_coefs(dissolution = ~offset(edges),
                                   duration = dur_coefs$main$durs.main.homog$mean.dur.adj,
                                   d.rate = drate_mat$drate_main),
    diss.byage = dissolution_coefs(dissolution = ~offset(edges) +
                                     offset(nodematch("age.grp", diff = TRUE)),
                                   duration = dur_coefs$main$durs.main.byage$mean.dur.adj,
                                   d.rate = drate_mat$drate_main)



  ),

  # casl : same deal as above but for casual network
  casl = list(

    edges = target_extract(df = target_df,
                           term = "edges",
                           model = "casual"),

    # concurrent
    concurrent = target_extract(df = target_df,
                                term = "concurrent",
                                model = "casual"),

    # nodefactor race
    nodefactor_race = c(target_extract(term = "nodefactor.race_ethnicity.blackNH", model = "casual"),
                        target_extract(term = "nodefactor.race_ethnicity.hispanic", model = "casual"),
                        target_extract(term = "nodefactor.race_ethnicity.otherNH", model = "casual"),
                        target_extract(term = "nodefactor.race_ethnicity.whiteNH", model = "casual")),


    # nodematch race
    nodematch_race = target_extract(term = "nodematch.race_ethnicity", model = "casual"),

    # nodematch black
    nodematch_race.1 = target_extract(term = "nodematch.race_ethnicity.blackNH", model = "casual"),

    # nodematch_age.grp
    # nodematch_age.grp = target_extract(term = "nodematch.age", model = "casual"),
    nodematch_age.grp = c(target_extract(term = "nodematch.age.16to20", model = "casual"),
                          target_extract(term = "nodematch.age.21to29", model = "casual")),


    # nodefactor_age.grp
    # nodefactor_age.grp = c(target_extract(term = "nodefactor.age.16to20", model = "casual"),
    #                       target_extract(term = "nodefactor.age.21to29", model = "casual")),

    # nodefactor_init_cas_cat
    nodefactor_deg.main = c(target_extract(term = "nodefactor.init_ser_cat.0", model = "casual"),
                            target_extract(term = "nodefactor.init_ser_cat.1+", model = "casual")),

    # fuzzynodematch_venues_all
    # fuzzynodematch_venues.all = target_extract(term = "fuzzynodematch.venues_all.TRUE", model = "casual"),
    fuzzynodematch_venues.all = drate_mat$venues_casual,

    # fuzzynodematch_apps_all
    fuzzynodematch_apps.all = drate_mat$apps_casual,

    # fuzzynodematch_apps_dating
    fuzzynodematch_apps.dating = target_extract(term = "fuzzynodematch.apps_dating.TRUE", model = "casual"),

    # dissolution model
    dissolution = dissolution_coefs(~offset(edges), duration = 57, d.rate = drate_mat$drate_cas),
    diss.homog = dissolution_coefs(dissolution = ~offset(edges),
                                   duration = dur_coefs$casl$durs.casl.homog$mean.dur.adj,
                                   d.rate = drate_mat$drate_cas),
    diss.byage = dissolution_coefs(dissolution = ~offset(edges) +
                                     offset(nodematch("age.grp", diff = TRUE)),
                                   duration = dur_coefs$casl$durs.casl.byage$mean.dur.adj,
                                   d.rate = drate_mat$drate_cas)


  ),

  # inst : more or less same as above but for one-off network
  inst = list(

    edges = target_extract(df = target_df,
                           term = "edges",
                           model = "one.time"),

    # nodefactor race
    nodefactor_race = c(target_extract(term = "nodefactor.race_ethnicity.blackNH", model = "one.time"),
                        target_extract(term = "nodefactor.race_ethnicity.hispanic", model = "one.time"),
                        target_extract(term = "nodefactor.race_ethnicity.otherNH", model = "one.time"),
                        target_extract(term = "nodefactor.race_ethnicity.whiteNH", model = "one.time")),


    # nodematch race
    nodematch_race = target_extract(term = "nodematch.race_ethnicity", model = "one.time"),

    # nodematch black
    nodematch_race.1 = target_extract(term = "nodematch.race_ethnicity.blackNH", model = "one.time"),

    # nodematch_age.grp
    # nodematch_age.grp = target_extract(term = "nodematch.age", model = "one.time"),
    nodematch_age.grp = c(target_extract(term = "nodematch.age.16to20", model = "one.time"),
                          target_extract(term = "nodematch.age.21to29", model = "one.time")),

    # nodefactor_age.grp
    # nodefactor_age.grp = c(target_extract(term = "nodefactor.age.16to20", model = "one.time"),
    #                        target_extract(term = "nodefactor.age.21to29", model = "one.time")),

    # nodefactor_init_pers_cat
    nodefactor_deg.tot = c(target_extract(term = "nodefactor.init_pers_cat.0", model = "one.time"),
                           target_extract(term = "nodefactor.init_pers_cat.1", model = "one.time"),
                           target_extract(term = "nodefactor.init_pers_cat.2", model = "one.time"),
                           target_extract(term = "nodefactor.init_pers_cat.3+", model = "one.time")),


    # fuzzynodematch_venues_all
    # fuzzynodematch_venues.all = target_extract(term = "fuzzynodematch.venues_all.TRUE", model = "one.time"),
    fuzzynodematch_venues.all = drate_mat$venues_onetime,

    # fuzzynodematch_apps_all
    fuzzynodematch_apps.all = drate_mat$apps_onetime,

    # fuzzynodematch_apps_dating
    fuzzynodematch_apps.dating = target_extract(term = "fuzzynodematch.apps_dating.TRUE", model = "one.time")

    # edges
    # nodefactor_race
    # nodematch_race
    # nodematch_race_diffF (Ask)
    # nodefactor_age.grp
    # nodematch_age.grp
    # absdiff_age
    # absdiff_sqrtage
    # nodefactor_deg.tot
    # concurrent
    # nodefactor_diag.status

  )
)


saveRDS(netstats, paste0(est_dir, "netstats-level-", context, "_", drate_mat$fit_no, ".rds"))


### Initialize network
nw <- network::network.initialize(n = nrow(egos),
                                  loops = FALSE,
                                  directed = FALSE)

attr_names <- names(netstats$attr)
attr_values <- netstats$attr

nw_main <- EpiModel::set_vertex_attribute(nw, attr_names, attr_values)
nw_casl <- nw_main
nw_inst <- nw_main


# 1. Have variable codings match ARTnet as much as possible
# 2. Make sure network attributes and target stats are similarly labeled
# 3. Create streamlined way of accessing target stats from CSV (this script has
# examples for ARTnet)

# build.epistats requires a  value for the geog.lvl argument. So we’ll need to come up with a placeholder value in whatever dataset we need
# and geog.cat value too

### Set vertex attributes
# nw <- network::set.vertex.attribute(nw, "race", egos$race)
# nw <- network::set.vertex.attribute(nw, "age", egos$age)
# nw <- network::set.vertex.attribute(nw, "age.grp", egos$age.grp)
#
# nw <- network::set.vertex.attribute(nw, "venues_all", egos$venues_all)
# # nw <- network::set.vertex.attribute(nw, "venues_dating", egos$venues_dating)
# # nw <- network::set.vertex.attribute(nw, "venues_nondating", egos$venues_nondating)
#
# nw <- network::set.vertex.attribute(nw, "apps_all", egos$apps_all)
# # nw <- network::set.vertex.attribute(nw, "apps_dating", egos$apps_dating)
# # nw <- network::set.vertex.attribute(nw, "apps_nondating", egos$apps_nondating)
#
# nw <- network::set.vertex.attribute(nw, "deg.casl", as.character(egos$deg.casl))
# nw <- network::set.vertex.attribute(nw, "deg.main", as.character(egos$deg.main))
#
# nw <- network::set.vertex.attribute(nw, "hiv", egos$hiv_status)
#
# nw_main <- nw
# nw_casl <- nw_main
# nw_inst <- nw_main


# A1. Main Model (No Apps or Venues) -------------------------------------------

### Specify target stats
# I figured that if I use the custom `target.stats.main` function to compile
# the vector of target stats, it'll be easier for us to know which value
# corresponds to which ERGM term:

print("A1 Main")

model_main <- ~ edges +
  # nodefactor("age.grp", levels = 1) +
  nodematch("age.grp", diff = TRUE) +
  # nodematch("age.grp", levels = -1) +
  concurrent +
  nodefactor("race", levels = -4) +
  nodematch("race") +
  # nodematch("race", diff = TRUE, levels = 1) +
  nodefactor("deg.casl", levels= -1)


target.stats.main <- c(
  edges = netstats$main$edges,
  # nodefactor_age.grp = netstats$main$nodefactor_age.grp[1],
  nodematch_age.grp = netstats$main$nodematch_age.grp,
  # nodematch_age.grp = netstats$main$nodematch_age.grp[-1],
  concurrent = netstats$main$concurrent,
  nodefactor_race = netstats$main$nodefactor_race[1:3],
  nodematch_race = netstats$main$nodematch_race,
  # nodematch_race.1 = netstats$main$nodematch_race.1,
  nodefactor_deg.casl = netstats$main$nodefactor_deg.casl[-1]
)
target.stats.main <- unname(target.stats.main)


# Node attribute names should match what's in epimodelHIV
# Age and race need to be updated
### And check order of categorical variables
# Change name for initial degree (`deg.main`)
# EpiModel supports up to 4 racial categories, but make sure that the values
# line up
# Look at the age matching module/part of EpiModelHIV to keep track of
# age category handling

# TARGETS ABOVE SHOULD FOLLOW THE ORDER OF THE FORMULA BELOW
fit_main <- netest(
  nw = nw_main,
  formation = model_main,
  target.stats = target.stats.main,
  coef.diss = netstats$main$diss.byage,
  set.control.ergm =
    control.ergm(
      parallel = 4,
      MCMC.interval = 10000,
      MCMLE.effectiveSize=NULL,
      MCMC.burnin = 1000,
      MCMC.samplesize = 20000,
      SAN.maxit = 20,
      SAN.nsteps.times = 10
    )
)

fit_main <- trim_netest(fit_main)

coef_df <- data.frame(treatment = "Basic",
                      model = "Main",
                      term = names(fit_main$coef.form),
                      estimate = fit_main$coef.form)


# A2. Casual Model (No Apps or Venues) -----------------------------------------

print("A2 Casual")

model_casl <- ~ edges +
  # nodefactor("age.grp", levels= 1) +
  nodematch("age.grp", diff = TRUE) +
  concurrent +
  nodefactor("race", levels=-4) +
  nodematch("race") +
  # nodematch("race", diff=TRUE, levels=1) +
  nodefactor("deg.main", levels=-1)


target.stats.casl <- c(
  edges =                           netstats$casl$edges,
  # nodefactor_age.grp =            netstats$casl$nodefactor_age.grp[1],
  nodematch_age.grp =             netstats$casl$nodematch_age.grp,
  concurrent =                      netstats$casl$concurrent,
  nodefactor_race =               netstats$casl$nodefactor_race[1:3],
  nodematch_race =                netstats$casl$nodematch_race,
  # nodematch_race.1 =              netstats$casl$nodematch_race.1,
  nodefactor_deg.main =           netstats$casl$nodefactor_deg.main[-1]
)
target.stats.casl <- unname(target.stats.casl)



# Node attribute names should match what's in epimodelHIV
# Age and race need to be updated
### And check order of categorical variables
# Change name for initial degree (`deg.main`)
# EpiModel supports up to 4 racial categories, but make sure that the values
# line up
# Look at the age matching module/part of EpiModelHIV to keep track of
# age category handling

# TARGETS ABOVE SHOULD FOLLOW THE ORDER OF THE FORMULA BELOW
fit_casl <- netest(
  nw = nw_casl,
  formation = model_casl,
  target.stats = target.stats.casl,
  coef.diss = netstats$casl$diss.byage,
  set.control.ergm =
    control.ergm(
      parallel = 4,
      MCMC.interval = 10000,
      MCMLE.effectiveSize=NULL,
      MCMC.burnin = 1000,
      MCMC.samplesize = 20000,
      SAN.maxit = 20,
      SAN.nsteps.times = 10
    )
)

fit_casl <- trim_netest(fit_casl)

casl_df <- data.frame(treatment = "Basic",
                      model = "Casual",
                      term = names(fit_casl$coef.form),
                      estimate = fit_casl$coef.form)

coef_df <- dplyr::bind_rows(coef_df, casl_df)



# A3. One-Off Model (No Apps or Venues) ----------------------------------------

print("A3 Inst")

model_inst <-  ~ edges +
  # nodefactor("age.grp", levels = 1) +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("race", levels=-4) +
  nodematch("race") +
  # nodematch("race", diff=TRUE, levels=1) +
  nodefactor("deg.tot", levels=-1)

target.stats.inst <- c(
  edges =                           netstats$inst$edges,
  # nodefactor_age.grp =            netstats$inst$nodefactor_age.grp[1],
  nodematch_age.grp =             netstats$inst$nodematch_age.grp,
  nodefactor_race =               netstats$inst$nodefactor_race[1:3],
  nodematch_race =                netstats$inst$nodematch_race,
  # nodematch_race.1 =              netstats$inst$nodematch_race.1,
  nodefactor_deg.tot =           netstats$inst$nodefactor_deg.tot[-1]
)
target.stats.inst <- unname(target.stats.inst)

fit_inst <- netest(
  nw = nw_inst,
  formation = model_inst,
  target.stats = target.stats.inst,
  coef.diss = dissolution_coefs(~offset(edges), duration = 1),
  set.control.ergm =
    control.ergm(
      parallel = 4,
      MCMC.interval = 10000,
      MCMLE.effectiveSize=NULL,
      MCMC.burnin = 1000,
      MCMC.samplesize = 20000,
      SAN.maxit = 20,
      SAN.nsteps.times = 10
    )
)

fit_inst <- trim_netest(fit_inst)

inst_df <- data.frame(treatment = "Basic",
                      model = "Onetime",
                      term = names(fit_inst$coef.form),
                      estimate = fit_inst$coef.form)

coef_df <- dplyr::bind_rows(coef_df, inst_df)



# A4. Save Data (No Apps or Venues)---------------------------------------------
out <- list(fit_main = fit_main, fit_casl = fit_casl, fit_inst = fit_inst)
saveRDS(out, paste0(est_dir, "basic_netest-level-allvenues", context, "_", drate_mat$fit_no, ".rds"))





# B1. Main Model (With Apps and Venues) ----------------------------------------

### Specify target stats
# I figured that if I use the custom `target.stats.main` function to compile
# the vector of target stats, it'll be easier for us to know which value
# corresponds to which ERGM term:

print("B1 Main")

model_main <- ~ edges +
  # nodefactor("age.grp", levels= 1) +
  nodematch("age.grp", diff = TRUE) +
  concurrent +
  nodefactor("race", levels = -4) +
  nodematch("race") +
  # nodematch("race", diff = TRUE, levels = 1) +

  nodefactor("deg.casl", levels= -1) +
  fuzzynodematch("venues.all", binary=TRUE) +
  fuzzynodematch("apps.all", binary = TRUE)
# fuzzynodematch("apps_nondating", binary=TRUE)


target.stats.main <- c(
  edges = netstats$main$edges,
  # nodefactor_age.grp = netstats$main$nodefactor_age.grp[1],
  nodematch_age.grp = netstats$main$nodematch_age.grp,
  concurrent = netstats$main$concurrent,
  nodefactor_race = netstats$main$nodefactor_race[1:3],
  nodematch_race = netstats$main$nodematch_race,
  # nodematch_race.1 = netstats$main$nodematch_race.1,

  nodefactor_deg.casl = netstats$main$nodefactor_deg.casl[-1],
  fuzzynodematch_venues.all = netstats$main$fuzzynodematch_venues.all,
  fuzzynodematch_apps.all = netstats$main$fuzzynodematch_apps.all
  # fuzzynodematch_apps.nondating = netstats$main$fuzzynodematch_apps.dating
)
target.stats.main <- unname(target.stats.main)


# Node attribute names should match what's in epimodelHIV
# Age and race need to be updated
### And check order of categorical variables
# Change name for initial degree (`deg.main`)
# EpiModel supports up to 4 racial categories, but make sure that the values
# line up
# Look at the age matching module/part of EpiModelHIV to keep track of
# age category handling

# TARGETS ABOVE SHOULD FOLLOW THE ORDER OF THE FORMULA BELOW
fit_main <- netest(
  nw = nw_main,
  formation = model_main,
  target.stats = target.stats.main,
  coef.diss = netstats$main$diss.byage,
  set.control.ergm =
    control.ergm(
      parallel = 4,
      MCMC.interval = 10000,
      MCMLE.effectiveSize=NULL,
      MCMC.burnin = 1000,
      MCMC.samplesize = 20000,
      SAN.maxit = 20,
      SAN.nsteps.times = 10
    )
)

fit_main <- trim_netest(fit_main)

main_df <- data.frame(treatment = "Venues and Apps",
                      model = "Main",
                      term = names(fit_main$coef.form),
                      estimate = fit_main$coef.form)

# coef_df <- dplyr::bind_rows(coef_df, main_df)

# B2. Casual Model (With Apps and Venues) --------------------------------------

print("B2 Casual")

model_casl <- ~ edges +
  # nodefactor("age.grp", levels= 1) +
  nodematch("age.grp", diff = TRUE) +
  concurrent +
  nodefactor("race", levels=-4) +
  nodematch("race") +
  # nodematch("race", diff=TRUE, levels=1) +

  nodefactor("deg.main", levels=-1) +
  fuzzynodematch("venues.all", binary=TRUE) +
  fuzzynodematch("apps.all", binary = TRUE)
#fuzzynodematch("apps_nondating", binary=TRUE)


target.stats.casl <- c(
  edges =                           netstats$casl$edges,
  #  nodefactor_age.grp =            netstats$casl$nodefactor_age.grp[1],
  nodematch_age.grp =             netstats$casl$nodematch_age.grp,
  concurrent =                      netstats$casl$concurrent,
  nodefactor_race =               netstats$casl$nodefactor_race[1:3],
  nodematch_race =                netstats$casl$nodematch_race,
  # nodematch_race.1 =              netstats$casl$nodematch_race.1,

  nodefactor_deg.main =           netstats$casl$nodefactor_deg.main[-1],
  fuzzynodematch_venues.all =     netstats$casl$fuzzynodematch_venues.all,
  fuzzynodematch_apps.all =       netstats$casl$fuzzynodematch_apps.all
  # fuzzynodematch_apps.nondating = netstats$casl$fuzzynodematch_apps.dating
)
target.stats.casl <- unname(target.stats.casl)



# Node attribute names should match what's in epimodelHIV
# Age and race need to be updated
### And check order of categorical variables
# Change name for initial degree (`deg.main`)
# EpiModel supports up to 4 racial categories, but make sure that the values
# line up
# Look at the age matching module/part of EpiModelHIV to keep track of
# age category handling

# TARGETS ABOVE SHOULD FOLLOW THE ORDER OF THE FORMULA BELOW
fit_casl <- netest(
  nw = nw_casl,
  formation = model_casl,
  target.stats = target.stats.casl,
  coef.diss = netstats$casl$diss.byage,
  set.control.ergm =
    control.ergm(
      parallel = 4,
      MCMC.interval = 10000,
      MCMLE.effectiveSize=NULL,
      MCMC.burnin = 1000,
      MCMC.samplesize = 20000,
      SAN.maxit = 20,
      SAN.nsteps.times = 10
    )
)

fit_casl <- trim_netest(fit_casl)

casl_df <- data.frame(treatment = "Venues and Apps",
                      model = "Casual",
                      term = names(fit_casl$coef.form),
                      estimate = fit_casl$coef.form)

# coef_df <- dplyr::bind_rows(coef_df, casl_df)


# B3. One-Off Model (With Apps and Venues) -------------------------------------

print("B3 Inst")

model_inst <-  ~ edges +
  # nodefactor("age.grp", levels=1) +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("race", levels=-4) +
  nodematch("race") +
  # nodematch("race", diff=TRUE, levels=1) +

  nodefactor("deg.tot", levels=-1) +
  fuzzynodematch("venues.all", binary=TRUE) +
  fuzzynodematch("apps.all", binary = TRUE)
#fuzzynodematch("apps_nondating", binary=TRUE)

target.stats.inst <- c(
  edges =                           netstats$inst$edges,
  # nodefactor_age.grp =            netstats$inst$nodefactor_age.grp[1],
  nodematch_age.grp =             netstats$inst$nodematch_age.grp,
  nodefactor_race =               netstats$inst$nodefactor_race[1:3],
  nodematch_race =                netstats$inst$nodematch_race,
  # nodematch_race.1 =              netstats$inst$nodematch_race.1,

  nodefactor_deg.tot =           netstats$inst$nodefactor_deg.tot[-1],
  fuzzynodematch_venues.all =       netstats$inst$fuzzynodematch_venues.all,
  fuzzynodematch_apps.all = netstats$inst$fuzzynodematch_apps.all
  # fuzzynodematch_apps.nondating = netstats$inst$fuzzynodematch_apps.dating
)
target.stats.inst <- unname(target.stats.inst)

fit_inst <- netest(
  nw = nw_inst,
  formation = model_inst,
  target.stats = target.stats.inst,
  coef.diss = dissolution_coefs(~offset(edges), duration = 1),
  set.control.ergm =
    control.ergm(
      parallel = 4,
      MCMC.interval = 10000,
      MCMLE.effectiveSize=NULL,
      MCMC.burnin = 1000,
      MCMC.samplesize = 20000,
      SAN.maxit = 20,
      SAN.nsteps.times = 10
    )
)

fit_inst <- trim_netest(fit_inst)

inst_df <- data.frame(treatment = "Venues and Apps",
                      model = "Onetime",
                      term = names(fit_inst$coef.form),
                      estimate = fit_inst$coef.form)

# coef_df <- dplyr::bind_rows(coef_df, inst_df)

# B4. Save Data (With Apps and Venues) -----------------------------------------
out <- list(fit_main = fit_main, fit_casl = fit_casl, fit_inst = fit_inst)
saveRDS(out, paste0(est_dir, "venues_apps_netest-level-allvenues", context, "_", drate_mat$fit_no, ".rds"))

################################

# C1. Main Model (Apps Only) ----------------------------------------

### Specify target stats
# I figured that if I use the custom `target.stats.main` function to compile
# the vector of target stats, it'll be easier for us to know which value
# corresponds to which ERGM term:

print("C1 Main")

model_main <- ~ edges +
  # nodefactor("age.grp", levels= 1) +
  nodematch("age.grp", diff = TRUE) +
  concurrent +
  nodefactor("race", levels = -4) +
  nodematch("race") +
  # nodematch("race", diff = TRUE, levels = 1) +

  nodefactor("deg.casl", levels= -1) +
  fuzzynodematch("apps.all", binary = TRUE)
# fuzzynodematch("apps_nondating", binary=TRUE)


target.stats.main <- c(
  edges = netstats$main$edges,
  # nodefactor_age.grp = netstats$main$nodefactor_age.grp[1],
  nodematch_age.grp = netstats$main$nodematch_age.grp,
  concurrent = netstats$main$concurrent,
  nodefactor_race = netstats$main$nodefactor_race[1:3],
  nodematch_race = netstats$main$nodematch_race,
  # nodematch_race.1 = netstats$main$nodematch_race.1,
  nodefactor_deg.casl = netstats$main$nodefactor_deg.casl[-1],
  fuzzynodematch_apps.all = netstats$main$fuzzynodematch_apps.all
  # fuzzynodematch_apps.nondating = netstats$main$fuzzynodematch_apps.dating
)
target.stats.main <- unname(target.stats.main)


# Node attribute names should match what's in epimodelHIV
# Age and race need to be updated
### And check order of categorical variables
# Change name for initial degree (`deg.main`)
# EpiModel supports up to 4 racial categories, but make sure that the values
# line up
# Look at the age matching module/part of EpiModelHIV to keep track of
# age category handling

# TARGETS ABOVE SHOULD FOLLOW THE ORDER OF THE FORMULA BELOW
fit_main <- netest(
  nw = nw_main,
  formation = model_main,
  target.stats = target.stats.main,
  coef.diss = netstats$main$diss.byage,
  set.control.ergm =
    control.ergm(
      parallel = 4,
      MCMC.interval = 10000,
      MCMLE.effectiveSize=NULL,
      MCMC.burnin = 1000,
      MCMC.samplesize = 20000,
      SAN.maxit = 20,
      SAN.nsteps.times = 10
    )
)

fit_main <- trim_netest(fit_main)

main_df <- data.frame(treatment = "Apps Only",
                      model = "Main",
                      term = names(fit_main$coef.form),
                      estimate = fit_main$coef.form)

# coef_df <- dplyr::bind_rows(coef_df, main_df)



# C2. Casual Model (Apps Only) --------------------------------------

print("C2 Casual")

model_casl <- ~ edges +
  # nodefactor("age.grp", levels= 1) +
  nodematch("age.grp", diff = TRUE) +
  concurrent +
  nodefactor("race", levels=-4) +
  nodematch("race") +
  # nodematch("race", diff=TRUE, levels=1) +

  nodefactor("deg.main", levels=-1) +
  fuzzynodematch("apps.all", binary = TRUE)
#fuzzynodematch("apps_nondating", binary=TRUE)


target.stats.casl <- c(
  edges =                           netstats$casl$edges,
  # nodefactor_age.grp =            netstats$casl$nodefactor_age.grp[1],
  nodematch_age.grp =             netstats$casl$nodematch_age.grp,
  concurrent =                      netstats$casl$concurrent,
  nodefactor_race =               netstats$casl$nodefactor_race[1:3],
  nodematch_race =                netstats$casl$nodematch_race,
  # nodematch_race.1 =              netstats$casl$nodematch_race.1,

  nodefactor_deg.main =           netstats$casl$nodefactor_deg.main[-1],
  fuzzynodematch_apps.all =       netstats$casl$fuzzynodematch_apps.all
  # fuzzynodematch_apps.nondating = netstats$casl$fuzzynodematch_apps.dating
)
target.stats.casl <- unname(target.stats.casl)



# Node attribute names should match what's in epimodelHIV
# Age and race need to be updated
### And check order of categorical variables
# Change name for initial degree (`deg.main`)
# EpiModel supports up to 4 racial categories, but make sure that the values
# line up
# Look at the age matching module/part of EpiModelHIV to keep track of
# age category handling

# TARGETS ABOVE SHOULD FOLLOW THE ORDER OF THE FORMULA BELOW
fit_casl <- netest(
  nw = nw_casl,
  formation = model_casl,
  target.stats = target.stats.casl,
  coef.diss = netstats$casl$diss.byage,
  set.control.ergm =
    control.ergm(
      parallel = 4,
      MCMC.interval = 10000,
      MCMLE.effectiveSize=NULL,
      MCMC.burnin = 1000,
      MCMC.samplesize = 20000,
      SAN.maxit = 20,
      SAN.nsteps.times = 10
    )
)

fit_casl <- trim_netest(fit_casl)

casl_df <- data.frame(treatment = "Apps Only",
                      model = "Casual",
                      term = names(fit_casl$coef.form),
                      estimate = fit_casl$coef.form)

# coef_df <- dplyr::bind_rows(coef_df, casl_df)

# C3. One-Off Model (Apps Only) -------------------------------------

print("C3 Inst")

model_inst <-  ~ edges +
  # nodefactor("age.grp", levels=1) +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("race", levels=-4) +
  nodematch("race") +
  # nodematch("race", diff=TRUE, levels=1) +

  nodefactor("deg.tot", levels=-1) +
  fuzzynodematch("apps.all", binary = TRUE)
#fuzzynodematch("apps_nondating", binary=TRUE)

target.stats.inst <- c(
  edges =                           netstats$inst$edges,
  # nodefactor_age.grp =            netstats$inst$nodefactor_age.grp[1],
  nodematch_age.grp =             netstats$inst$nodematch_age.grp,
  nodefactor_race =               netstats$inst$nodefactor_race[1:3],
  nodematch_race =                netstats$inst$nodematch_race,
  # nodematch_race.1 =              netstats$inst$nodematch_race.1,

  nodefactor_deg.tot =           netstats$inst$nodefactor_deg.tot[-1],
  fuzzynodematch_apps.all = netstats$inst$fuzzynodematch_apps.all
  # fuzzynodematch_apps.nondating = netstats$inst$fuzzynodematch_apps.dating
)
target.stats.inst <- unname(target.stats.inst)

fit_inst <- netest(
  nw = nw_inst,
  formation = model_inst,
  target.stats = target.stats.inst,
  coef.diss = dissolution_coefs(~offset(edges), duration = 1),
  set.control.ergm =
    control.ergm(
      parallel = 4,
      MCMC.interval = 10000,
      MCMLE.effectiveSize=NULL,
      MCMC.burnin = 1000,
      MCMC.samplesize = 20000,
      SAN.maxit = 20,
      SAN.nsteps.times = 10
    )
)

fit_inst <- trim_netest(fit_inst)

inst_df <- data.frame(treatment = "Apps Only",
                      model = "Onetime",
                      term = names(fit_inst$coef.form),
                      estimate = fit_inst$coef.form)

# coef_df <- dplyr::bind_rows(coef_df, inst_df)

# C4. Save Data (With Apps and Venues) -----------------------------------------
out <- list(fit_main = fit_main, fit_casl = fit_casl, fit_inst = fit_inst)
saveRDS(out, paste0(est_dir, "apps_only_netest-level-allvenues", context, "_", drate_mat$fit_no, ".rds"))

#################################

# D1. Main Model (Venues Only) ----------------------------------------

### Specify target stats
# I figured that if I use the custom `target.stats.main` function to compile
# the vector of target stats, it'll be easier for us to know which value
# corresponds to which ERGM term:

print("D1 Main")

model_main <- ~ edges +
  # nodefactor("age.grp", levels= 1) +
  nodematch("age.grp", diff = TRUE) +
  concurrent +
  nodefactor("race", levels = -4) +
  nodematch("race") +
  # nodematch("race", diff = TRUE, levels = 1) +

  nodefactor("deg.casl", levels= -1) +
  fuzzynodematch("venues.all", binary=TRUE)
# fuzzynodematch("apps_nondating", binary=TRUE)


target.stats.main <- c(
  edges = netstats$main$edges,
  # nodefactor_age.grp = netstats$main$nodefactor_age.grp[1],
  nodematch_age.grp = netstats$main$nodematch_age.grp,
  concurrent = netstats$main$concurrent,
  nodefactor_race = netstats$main$nodefactor_race[1:3],
  nodematch_race = netstats$main$nodematch_race,
  # nodematch_race.1 = netstats$main$nodematch_race.1,

  nodefactor_deg.casl = netstats$main$nodefactor_deg.casl[-1],
  fuzzynodematch_venues.all = netstats$main$fuzzynodematch_venues.all
  # fuzzynodematch_apps.nondating = netstats$main$fuzzynodematch_apps.dating
)
target.stats.main <- unname(target.stats.main)


# Node attribute names should match what's in epimodelHIV
# Age and race need to be updated
### And check order of categorical variables
# Change name for initial degree (`deg.main`)
# EpiModel supports up to 4 racial categories, but make sure that the values
# line up
# Look at the age matching module/part of EpiModelHIV to keep track of
# age category handling

# TARGETS ABOVE SHOULD FOLLOW THE ORDER OF THE FORMULA BELOW
fit_main <- netest(
  nw = nw_main,
  formation = model_main,
  target.stats = target.stats.main,
  coef.diss = netstats$main$diss.byage,
  set.control.ergm =
    control.ergm(
      parallel = 4,
      MCMC.interval = 10000,
      MCMLE.effectiveSize=NULL,
      MCMC.burnin = 1000,
      MCMC.samplesize = 20000,
      SAN.maxit = 20,
      SAN.nsteps.times = 10
    )
)

fit_main <- trim_netest(fit_main)

main_df <- data.frame(treatment = "Venues Only",
                      model = "Main",
                      term = names(fit_main$coef.form),
                      estimate = fit_main$coef.form)

# coef_df <- dplyr::bind_rows(coef_df, main_df)

# D2. Casual Model (Venues Only) --------------------------------------

print("D2 Casual")

model_casl <- ~ edges +
  # nodefactor("age.grp", levels= 1) +
  nodematch("age.grp", diff = TRUE) +
  concurrent +
  nodefactor("race", levels=-4) +
  nodematch("race") +
  # nodematch("race", diff=TRUE, levels=1) +

  nodefactor("deg.main", levels=-1) +
  fuzzynodematch("venues.all", binary=TRUE)
#fuzzynodematch("apps_nondating", binary=TRUE)


target.stats.casl <- c(
  edges =                           netstats$casl$edges,
  # nodefactor_age.grp =            netstats$casl$nodefactor_age.grp[1],
  nodematch_age.grp =             netstats$casl$nodematch_age.grp,
  concurrent =                      netstats$casl$concurrent,
  nodefactor_race =               netstats$casl$nodefactor_race[1:3],
  nodematch_race =                netstats$casl$nodematch_race,
  # nodematch_race.1 =              netstats$casl$nodematch_race.1,
  nodefactor_deg.main =           netstats$casl$nodefactor_deg.main[-1],
  fuzzynodematch_venues.all =     netstats$casl$fuzzynodematch_venues.all
  # fuzzynodematch_apps.nondating = netstats$casl$fuzzynodematch_apps.dating
)
target.stats.casl <- unname(target.stats.casl)



# Node attribute names should match what's in epimodelHIV
# Age and race need to be updated
### And check order of categorical variables
# Change name for initial degree (`deg.main`)
# EpiModel supports up to 4 racial categories, but make sure that the values
# line up
# Look at the age matching module/part of EpiModelHIV to keep track of
# age category handling

# TARGETS ABOVE SHOULD FOLLOW THE ORDER OF THE FORMULA BELOW
fit_casl <- netest(
  nw = nw_casl,
  formation = model_casl,
  target.stats = target.stats.casl,
  coef.diss = netstats$casl$diss.byage,
  set.control.ergm =
    control.ergm(
      parallel = 4,
      MCMC.interval = 10000,
      MCMLE.effectiveSize=NULL,
      MCMC.burnin = 1000,
      MCMC.samplesize = 20000,
      SAN.maxit = 20,
      SAN.nsteps.times = 10
    )
)

fit_casl <- trim_netest(fit_casl)

casl_df <- data.frame(treatment = "Venues Only",
                      model = "Casual",
                      term = names(fit_casl$coef.form),
                      estimate = fit_casl$coef.form)

# coef_df <- dplyr::bind_rows(coef_df, casl_df)



# D3. One-Off Model (Venues Only) -------------------------------------

print("D3 Inst")

model_inst <-  ~ edges +
  # nodefactor("age.grp", levels=1) +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("race", levels=-4) +
  nodematch("race") +
  # nodematch("race", diff=TRUE, levels=1) +

  nodefactor("deg.tot", levels=-1) +
  fuzzynodematch("venues.all", binary=TRUE)
#fuzzynodematch("apps_nondating", binary=TRUE)

target.stats.inst <- c(
  edges =                           netstats$inst$edges,
  # nodefactor_age.grp =            netstats$inst$nodefactor_age.grp[1],
  nodematch_age.grp =             netstats$inst$nodematch_age.grp,
  nodefactor_race =               netstats$inst$nodefactor_race[1:3],
  nodematch_race =                netstats$inst$nodematch_race,
  # nodematch_race.1 =              netstats$inst$nodematch_race.1,

  nodefactor_deg.tot =           netstats$inst$nodefactor_deg.tot[-1],
  fuzzynodematch_venues.all =       netstats$inst$fuzzynodematch_venues.all
  # fuzzynodematch_apps.nondating = netstats$inst$fuzzynodematch_apps.dating
)
target.stats.inst <- unname(target.stats.inst)

fit_inst <- netest(
  nw = nw_inst,
  formation = model_inst,
  target.stats = target.stats.inst,
  coef.diss = dissolution_coefs(~offset(edges), duration = 1),
  set.control.ergm =
    control.ergm(
      parallel = 4,
      MCMC.interval = 10000,
      MCMLE.effectiveSize=NULL,
      MCMC.burnin = 1000,
      MCMC.samplesize = 20000,
      SAN.maxit = 20,
      SAN.nsteps.times = 10
    )
)

fit_inst <- trim_netest(fit_inst)

inst_df <- data.frame(treatment = "Venues Only",
                      model = "Onetime",
                      term = names(fit_inst$coef.form),
                      estimate = fit_inst$coef.form)

# coef_df <- dplyr::bind_rows(coef_df, inst_df)


# D4. Save Data (Venues Only) -----------------------------------------
out <- list(fit_main = fit_main, fit_casl = fit_casl, fit_inst = fit_inst)
saveRDS(out, paste0(est_dir, "venue_only_netest-level-allvenues", context, "_", drate_mat$fit_no, ".rds"))

