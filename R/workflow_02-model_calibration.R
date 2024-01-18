##
## Epidemic Model Parameter Calibration, HPC setup
##

# Libraries --------------------------------------------------------------------
library("slurmworkflow")
library("EpiModelHPC")
library("EpiModelHIV")

# Settings ---------------------------------------------------------------------
source("./R/utils-0_project_settings.R")
context <- "hpc"
max_cores <- 1

source("./R/utils-chistig_basic_inputs.R") # make `path_to_est`, `param` and `init`
source("./R/utils-hpc_configs.R") # creates `hpc_configs`

# ------------------------------------------------------------------------------

# Workflow creation
wf <- create_workflow(
  wf_name = "model_calibration",
  default_sbatch_opts = hpc_configs$default_sbatch_opts
)

# Update RENV on the HPC
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_renv_restore(
    git_branch = current_git_branch,
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = hpc_configs$renv_sbatch_opts
)

# Controls
source("./R/utils-targets.R")
control <- control_msm(
  nsteps              = calibration_end,
  nsims               = 1,
  ncores              = 1,
  cumulative.edgelist = TRUE,
  truncate.el.cuml    = 0,
  .tracker.list       = calibration_trackers,
  # .checkpoint.dir     = "./temp/cp_calib",
  # .checkpoint.clear   = TRUE,
  # .checkpoint.steps   = 15 * 52,
  verbose             = FALSE,

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

# insert test values here
n_scenarios <- 2
scenarios_df <- tibble(
  .scenario.id = as.character(seq_len(n_scenarios)),
  .at = 1,
  ugc.prob = seq(0.3225, 0.3275, length.out = n_scenarios), # best 0.325
  rgc.prob = plogis(qlogis(ugc.prob) + log(1.25)),
  uct.prob = seq(0.29, 0.294, length.out = n_scenarios), # best 0.291
  rct.prob = plogis(qlogis(uct.prob) + log(1.25))
)
scenarios_list <- EpiModel::create_scenario_list(scenarios_df)

wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_netsim_scenarios(
    path_to_est, param, init, control,
    scenarios_list = NULL, #cenarios_list,
    output_dir = "./data/intermediate/calibration",
    # libraries = c(#"EpiModelHIV",
    #               "slurmworkflow",
    #               "EpiModelHPC",
    #               "chiSTIGmodules"
    #               ),
    n_rep = 2,
    n_cores = max_cores,
    save_pattern = "simple",
    max_array_size = 999,
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "mail-type" = "FAIL,TIME_LIMIT",
    "cpus-per-task" = max_cores,
    "time" = "04:00:00",
    "mem" = "5G"
  )
)

# # Process calibrations
# #
# # produce a data frame with the calibration targets for each scenario
# wf <- add_workflow_step(
#   wf_summary = wf,
#   step_tmpl = step_tmpl_do_call_script(
#     r_script = "./R/11-calibration_process.R",
#     args = list(
#       context = "hpc",
#       ncores = 15
#     ),
#     setup_lines = hpc_configs$r_loader
#   ),
#   sbatch_opts = list(
#     "cpus-per-task" = 15,
#     "time" = "04:00:00",
#     "mem-per-cpu" = "4G",
#     "mail-type" = "END"
#   )
# )

# Send the workflow folder to the <HPC> and run it
#
# $ scp -r ./workflows/model_calibration <HPC>:<project_dir>/workflows/
#
# on the HPC:
# $ ./workflows/model_calibration/start_workflow.sh

# Once the worfklow is finished download the data from the HPC
#
# $ scp -r <HPC>:<project_dir>/data/intermediate/calibration/assessments.rds ./data/intermediate/calibration/
#
# and analyse them locally using: "./R/12-calibration_eval.R"
