## Example interactive epidemic simulation run script with basic
## parameterization and all parameters defined in `param_msm`, with example of
## writing/debugging modules

set.seed(15)

# Define which "Treatment" we're running here
#### 1 - "Control" Simulation (No apps, no venues)
#### 2 - Venues, no apps
#### 3 - Apps, no venues
#### 4 - Venues and Apps
treatment <- 1


# Libraries  -------------------------------------------------------------------
library("chiSTIGmodules")
library("EpiModelHIV")

### 0. Set up python and R environments ###
# load reticulate
library(reticulate)
reticulate::use_python("/projects/p31728/pythonenvs/env1/bin/python")
#
#
#### ChiSTIG model prelim ------------------------------------------------------
python_funs <- environment(reticulate::source_python(paste0(this_dir,"chistig_reticulate_python_script.py")))

# necessary emppop files
empop_demo_file = paste0(this_dir, "data/input/empop_egoid_demo_def.csv")
empop_vmatrix_file = paste0(this_dir, "data/input/empop_ego_venue_matrix.csv")
empop_amatrix_file = paste0(this_dir, "data/input/empop_ego_app_matrix_binary.csv")
empop_relbin_file = paste0(this_dir, "data/input/ser_rel_data.csv")

# necessary synthpop files
efile <- paste0(this_dir, "data/input/egos_v4_0.csv")
vfile <- paste0(this_dir, "data/input/venue-attendance_v4_0.csv")

# necessary app and venue files
atypefile <- paste0(this_dir, "data/input/empop_appid_type_def.csv")
vtypefile <- paste0(this_dir, "data/input/empop_venueid_type_def.csv")

# format the empirical population datasets
empopdfs <- python_funs$init_empop(empop_demo_file, empop_vmatrix_file, empop_amatrix_file, empop_relbin_file)
empopegos2venues <- python_funs$create_emppop_venue_attendance_dict(empopdfs$vmatrix)
empopegos2demos <- python_funs$create_emppop_buckets(empopdfs$demo)
empopegos2apps <- python_funs$create_emppop_appuse_dict(empopdfs$amatrix)
empopegos2rel <- python_funs$create_rel_binary_emppop_buckets(empopegos2demos, empopdfs$rel)

# setting up the mapping of venues to venue type
venues2types <- python_funs$categorize_venues(vtypefile)

# setting up initial synthetic population and attendance/appuse objects
segos2venues <- python_funs$obtain_egos_vdict(vfile)
segos2apps <- python_funs$categorize_apps_segos(efile, atypefile)
sego2apps_norel <- python_funs$assign_apps_norel(efile, segos2apps)

# set a variable to toggle whether or not appuse is assigned via relationship status
appuse_by_rel_status <- FALSE


# # Settings ---------------------------------------------------------------------
source(paste0(this_dir, "R/utils-0_project_settings.R"))
source(paste0(this_dir, "R/utils-targets.R"))
#
# Necessary files
if (treatment == 1) {
  epistats <- readRDS("./data/intermediate/estimates/epistats-local.rds")
  netstats <- readRDS("./data/intermediate/estimates/netstats-novenues-local.rds")
  est <- readRDS("./data/intermediate/basic_netest-local.rds")
} else if (treatment == 2) {
  epistats <- readRDS("./data/intermediate/estimates/epistats-local.rds")
  netstats <- readRDS("./data/intermediate/estimates/netstats-novenues-local.rds")
  est <- readRDS("./data/intermediate/venue_only_netest-local.rds")
} else if (treatment == 3) {
  epistats <- readRDS("./data/intermediate/estimates/epistats-local.rds")
  netstats <- readRDS("./data/intermediate/estimates/netstats-novenues-local.rds")
  est <- readRDS("./data/intermediate/apps_only_netest-local.rds")
} else {
  epistats <- readRDS("./data/intermediate/estimates/epistats-local.rds")
  netstats <- readRDS("./data/intermediate/estimates/netstats-novenues-local.rds")
  est <- readRDS("./data/intermediate/venues_apps_netest-local.rds")
}

epistats$age.breaks <- c(16, 20, 30)
epistats$age.limits <- c(16, 30)

netstats$attr$age <- sample(16:29, length(netstats$attr$age), replace = TRUE)
netstats$attr$age <- netstats$attr$age + sample(1:1000, length(netstats$attr$age), replace = TRUE)/1000


param <- EpiModel::param.net(
  data.frame.params = readr::read_csv(paste0(this_dir, "data/input/params_chistig_feb19.csv")),
  netstats          = netstats,
  epistats          = epistats,
  prep.start        = Inf,
  riskh.start       = Inf
)
print(param)


# Initial conditions
init <- EpiModelHIV::init_msm()

if (treatment == 1) {

  # Control settings
  control <- control_msm(
    nsteps = 52 * 70,
    nsims  = 1,
    ncores = 1,
    .tracker.list       = calibration_trackers,

    initialize.FUN =              chiSTIGmodules::initialize_msm_chi,
    aging.FUN =                   chiSTIGmodules::aging_msm_chi,
    departure.FUN =               chiSTIGmodules::departure_msm_chi,
    arrival.FUN =                 chiSTIGmodules::arrival_msm_chi,
    # venues.FUN =                  chiSTIGmodules:::venues_msm_chi,
    partident.FUN =               chiSTIGmodules::partident_msm_chi,
    hivtest.FUN =                 chiSTIGmodules::hivtest_msm_chi,
    hivtx.FUN =                   chiSTIGmodules::hivtx_msm_chi,
    hivprogress.FUN =             chiSTIGmodules::hivprogress_msm_chi,
    hivvl.FUN =                   chiSTIGmodules::hivvl_msm_chi,
    resim_nets.FUN =              chiSTIGmodules::simnet_msm_chi,
    acts.FUN =                    chiSTIGmodules::acts_msm_chi,
    condoms.FUN =                 chiSTIGmodules::condoms_msm_chi,
    position.FUN =                chiSTIGmodules::position_msm_chi,
    prep.FUN =                    chiSTIGmodules::prep_msm_chi,
    hivtrans.FUN =                chiSTIGmodules::hivtrans_msm_chi,
    exotrans.FUN =                chiSTIGmodules:::exotrans_msm_chi,
    stitrans.FUN =                chiSTIGmodules::stitrans_msm_chi,
    stirecov.FUN =                chiSTIGmodules::stirecov_msm_chi,
    stitx.FUN =                   chiSTIGmodules::stitx_msm_chi,
    prev.FUN =                    chiSTIGmodules::prevalence_msm_chi,
    cleanup.FUN =                 chiSTIGmodules::cleanup_msm_chi,

    module.order = c("aging.FUN", "departure.FUN", "arrival.FUN", # "venues.FUN",
                     "partident.FUN", "hivtest.FUN", "hivtx.FUN", "hivprogress.FUN",
                     "hivvl.FUN", "resim_nets.FUN", "acts.FUN", "condoms.FUN",
                     "position.FUN", "prep.FUN", "hivtrans.FUN", "exotrans.FUN",
                     "stitrans.FUN", "stirecov.FUN", "stitx.FUN", "prev.FUN", "cleanup.FUN")
  )

} else {

  # Control settings
  control <- control_msm(
    nsteps = 52 * 70,
    nsims  = 1,
    ncores = 1,
    .tracker.list       = calibration_trackers,

    initialize.FUN =              chiSTIGmodules::initialize_msm_chi,
    aging.FUN =                   chiSTIGmodules::aging_msm_chi,
    departure.FUN =               chiSTIGmodules::departure_msm_chi,
    arrival.FUN =                 chiSTIGmodules::arrival_msm_chi,
    venues.FUN =                  chiSTIGmodules:::venues_msm_chi,
    partident.FUN =               chiSTIGmodules::partident_msm_chi,
    hivtest.FUN =                 chiSTIGmodules::hivtest_msm_chi,
    hivtx.FUN =                   chiSTIGmodules::hivtx_msm_chi,
    hivprogress.FUN =             chiSTIGmodules::hivprogress_msm_chi,
    hivvl.FUN =                   chiSTIGmodules::hivvl_msm_chi,
    resim_nets.FUN =              chiSTIGmodules::simnet_msm_chi,
    acts.FUN =                    chiSTIGmodules::acts_msm_chi,
    condoms.FUN =                 chiSTIGmodules::condoms_msm_chi,
    position.FUN =                chiSTIGmodules::position_msm_chi,
    prep.FUN =                    chiSTIGmodules::prep_msm_chi,
    hivtrans.FUN =                chiSTIGmodules::hivtrans_msm_chi,
    exotrans.FUN =                chiSTIGmodules:::exotrans_msm_chi,
    stitrans.FUN =                chiSTIGmodules::stitrans_msm_chi,
    stirecov.FUN =                chiSTIGmodules::stirecov_msm_chi,
    stitx.FUN =                   chiSTIGmodules::stitx_msm_chi,
    prev.FUN =                    chiSTIGmodules::prevalence_msm_chi,
    cleanup.FUN =                 chiSTIGmodules::cleanup_msm_chi,

    module.order = c("aging.FUN", "departure.FUN", "arrival.FUN", "venues.FUN",
                     "partident.FUN", "hivtest.FUN", "hivtx.FUN", "hivprogress.FUN",
                     "hivvl.FUN", "resim_nets.FUN", "acts.FUN", "condoms.FUN",
                     "position.FUN", "prep.FUN", "hivtrans.FUN", "exotrans.FUN",
                     "stitrans.FUN", "stirecov.FUN", "stitx.FUN", "prev.FUN", "cleanup.FUN")
  )

}


#
start_time <- Sys.time()
# Epidemic simulation
sim <- netsim(est, param, init, control)
end_time <- Sys.time()

saveRDS(sim, paste0(this_dir, "data/intermediate/treatment", treatment, "_tom_test.rds"))
