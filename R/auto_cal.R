model <- function(proposal) {
  # Load all required elements
  library(EpiModelHIV)
  library(dplyr)

  epistats <- readRDS("data/input/epistats.rds")
  netstats <- readRDS("data/input/netstats.rds")
  est      <- readRDS("data/input/netest.rds")

  param <- param.net(
    data.frame.params = read.csv("data/input/params.csv"),
    netstats          = netstats,
    epistats          = epistats
  )

  init <- init_msm()

  control <- control_msm(
    nsteps = 52 * 60,
    nsims  = 1,
    ncores = 1
  )

  # Proposal to scenario -------------------------------------------------------
  scenario <- EpiModelHPC::swfcalib_proposal_to_scenario(proposal)
  param_sc <- EpiModel::use_scenario(param, scenario)

  # Run the simulation ---------------------------------------------------------
  sim <- netsim(est, param_sc, init, control)

  # Process the results  -------------------------------------------------------
  results <- as_tibble(sim) |>
    mutate_calibration_targets() |>
    filter(time >= max(time) - 52) |>
    select(
      cc.dx.B, cc.dx.H, cc.dx.W,
      cc.linked1m.B, cc.linked1m.H, cc.linked1m.W,
      i.prev.dx.B, i.prev.dx.H, i.prev.dx.W,
    ) |>
    summarise(across( everything(), ~ mean(.x, na.rm = TRUE)))

  # Return the one row `tibble`
  return(results)
}



config = list(
  simulator = model,
  root_directory = "data/calib",
  n_sims = n_sims,
  max_iteration = 100,
  default_proposal = dplyr::tibble(
    # Testing
    hiv.test.rate_1   = 0.0038,
    hiv.test.rate_2   = 0.004192807,
    hiv.test.rate_3   = .0055,
    hiv.test.rate_4   = 0.0046625,

    # Linked to care
    tx.init.rate_1    = 0.3622703,
    tx.init.rate_2    = 0.39,
    tx.init.rate_3    = 0.4202679,
    tx.init.rate_4    = 0.52,

    # Treatment discontinuation
    tx.halt.full.or_1 = 0.9,
    tx.halt.full.or_2 = 0.63,
    tx.halt.full.or_3 = 1.45,
    tx.halt.full.or_4 = 1.25,

    # trans.scale
    hiv.trans.scale_1 = 21,
    hiv.trans.scale_2 = 2.7,
    hiv.trans.scale_3 = 3.3,
    hiv.trans.scale_4 = 1.6
  )
)

# Complete Config

n_sims  <- 400

calib_object <- list(
  config = list(
    simulator = model,
    default_proposal = dplyr::tibble(
    # hiv.test.rate
      hiv.test.rate_1   = 0.0038,
      hiv.test.rate_2   = 0.004192807,
      hiv.test.rate_3   = .0055,
      hiv.test.rate_4   = 0.0046625,
    # tx.init
    tx.init.rate_1    = 0.3622703,
    tx.init.rate_2    = 0.39,
    tx.init.rate_3    = 0.4202679,
    tx.init.rate_4    = 0.52,
    # tx.halt.full.or
      tx.halt.full.or_1 = 0.9,
      tx.halt.full.or_2 = 0.63,
      tx.halt.full.or_3 = 1.45,
      tx.halt.full.or_4 = 1.25,
    # trans.scale
      hiv.trans.scale_1 = 21,
      hiv.trans.scale_2 = 2.7,
      hiv.trans.scale_3 = 3.3,
      hiv.trans.scale_4 = 1.6
    ),
    root_directory = "data/calib",
    max_iteration = 100,
    n_sims = n_sims
  ),
  waves = list(
    wave1 = list(
      job1 = list(
        targets = "cc.dx.B",
        targets_val = 0.5465356428,
        params = c("hiv.test.rate_1"), # target: 0.00385
        initial_proposals = dplyr::tibble(
          hiv.test.rate_1 = seq(0.002, 0.006, length.out = n_sims), # ADD THESE
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
      ),
      job2 = list(
        targets = "cc.dx.H",
        targets_val = 0.5431367893,
        params = c("hiv.test.rate_2"), # target: 0.0038
        initial_proposals = dplyr::tibble(
          hiv.test.rate_2 = seq(0.002, 0.006, length.out = n_sims), # ADD THESE
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
      ),

      job3 = list(
        targets = "cc.dx.O",
        targets_val = 0.5988779867,
        params = c("hiv.test.rate_3"), # target: 0.0038
        initial_proposals = dplyr::tibble(
          hiv.test.rate_2 = seq(0.002, 0.006, length.out = n_sims), # ADD THESE
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
      ),

      job4 = list(
        targets = "cc.dx.W",
        targets_val = 0.5614905982,
        params = c("hiv.test.rate_4"), # target: 0.0069
        initial_proposals = dplyr::tibble(
          hiv.test.rate_3 = seq(0.004, 0.008, length.out = n_sims),
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
      ),

      # Linked to care
      job5 = list(
        targets = paste0("cc.linked1m.", c("B", "H", "O", "W")),
        targets_val = c(.828, 0.867, 0.875, 0.936),
        params = paste0("tx.init.rate_", 1:4),
        initial_proposals = dplyr::tibble(
          tx.init.rate_1 = sample(seq(0.1, 0.6, length.out = n_sims)),
          tx.init.rate_2 = sample(tx.init.rate_1),
          tx.init.rate_3 = sample(tx.init.rate_1),
          tx.init.rate_4 = sample(tx.init.rate_1)
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 3)
      ),

      # ART Discontinuation
      job5 = list(
        targets = paste0("cc.vsupp.", c("B", "H", "O", "W")),
        targets_val = c(0.571, 0.675, 0.586, 0.617),
        params = paste0("tx.init.rate_", 1:4),
        initial_proposals = dplyr::tibble(
          tx.init.rate_1 = sample(seq(0.2, 0.6, length.out = n_sims)),
          tx.init.rate_2 = sample(tx.init.rate_1),
          tx.init.rate_3 = sample(tx.init.rate_1),
          tx.init.rate_4 = sample(tx.init.rate_1)
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 3)
      )
    ),
    wave2 = list(
      job1 = list(
        targets = paste0("i.prev.dx.", c("B", "H", "O", "W")),
        targets_val = c(0.5398688, 0.1417571, 0.1303336, 0.0681581),
        params = paste0("hiv.trans.scale_", 1:4),
        initial_proposals = dplyr::tibble(
          hiv.trans.scale_1 = sample(seq(18, 24, length.out = n_sims)),
          hiv.trans.scale_2 = sample(seq(2, 4, length.out = n_sims)),
          hiv.trans.scale_3 = sample(seq(2, 4, length.out = n_sims)),
          hiv.trans.scale_4 = sample(seq(1, 2, length.out = n_sims))
        ),
        make_next_proposals =
          swfcalib::make_proposer_se_range(n_sims, retain_prop = 0.3),
        get_result = swfcalib::determ_end_thresh(
          thresholds = rep(0.02, 3),
          n_enough = 100
        )
      )
    )
  )
)
################################################################################
# Running a Calibration System
#### Calibrate population at the same time as testing since it's independent from tests
#### Exo trans could probably go at the beginning if you need it
#### Always put in first wave unless you really

library("slurmworkflow")

# Define the `model` function
model <- function(proposal) {
  # Load all required elements
  library(EpiModelHIV)
  library(dplyr)

  epistats <- readRDS("data/input/epistats.rds")
  netstats <- readRDS("data/input/netstats.rds")
  est      <- readRDS("/Users/wms1212/Desktop/ChiSTIG_model/epimodel/data/intermediate/estimates/basic_netest-local.rds")
  # Is the aging out of the older initial nodes driving HIV extinction?
  # Let's level out the age distribution and find out
  netstats$attr$age <- sample(16:29, length(netstats$attr$age), replace = TRUE)
  netstats$attr$age <- netstats$attr$age + sample(1:1000, length(netstats$attr$age), replace = TRUE)/1000

  prep_start <- 52 * 2
  param <- param.net(
    data.frame.params   = read.csv("/Users/wms1212/Desktop/ChiSTIG_model/epimodel/data/input/params_chistig_nov30.csv"),
    netstats            = netstats,
    epistats            = epistats,
    prep.start          = prep_start,
    riskh.start         = prep_start - 53
  )

  init <- init_msm()

  control <- control_msm(
    nsteps = 52 * 60,
    nsims  = 1,
    ncores = 1
  )

  # Proposal to scenario -------------------------------------------------------
  scenario <- EpiModelHPC::swfcalib_proposal_to_scenario(proposal)
  param_sc <- EpiModel::use_scenario(param, scenario)

  # Run the simulation ---------------------------------------------------------
  sim <- netsim(est, param_sc, init, control)

  # Process the results  -------------------------------------------------------
  results <- as_tibble(sim) |>
    mutate_calibration_targets() |>
    filter(time >= max(time) - 52) |>
    select(
      cc.dx.B, cc.dx.H, cc.dx.W,
      cc.linked1m.B, cc.linked1m.H, cc.linked1m.W,
      i.prev.dx.B, i.prev.dx.H, i.prev.dx.W,
    ) |>
    summarise(across( everything(), ~ mean(.x, na.rm = TRUE)))

  # Return the one row `tibble`
  return(results)
}

# Create the `calib_object`
n_sims  <- 400

calib_object <- list(
  config = list(
    simulator = model,
    default_proposal = dplyr::tibble(
      # hiv.test.rate
      hiv.test.rate_1   = 0.0038,
      hiv.test.rate_2   = 0.004192807,
      hiv.test.rate_3   = .0055,
      hiv.test.rate_4   = 0.0046625,
      # tx.init
      tx.init.rate_1    = 0.3622703,
      tx.init.rate_2    = 0.39,
      tx.init.rate_3    = 0.4202679,
      tx.init.rate_4    = 0.52,
      # tx.halt.full.or
      tx.halt.full.or_1 = 0.9,
      tx.halt.full.or_2 = 0.63,
      tx.halt.full.or_3 = 1.45,
      tx.halt.full.or_4 = 1.25,
      # trans.scale
      hiv.trans.scale_1 = 21,
      hiv.trans.scale_2 = 2.7,
      hiv.trans.scale_3 = 3.3,
      hiv.trans.scale_4 = 1.6
    ),
    root_directory = "data/calib",
    max_iteration = 100,
    n_sims = n_sims
  ),
  waves = list(
    wave1 = list(
      job1 = list(
        targets = "cc.dx.B",
        targets_val = 0.5465356428,
        params = c("hiv.test.rate_1"), # target: 0.00385
        initial_proposals = dplyr::tibble(
          hiv.test.rate_1 = seq(0.002, 0.006, length.out = n_sims), # ADD THESE
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
      ),
      job2 = list(
        targets = "cc.dx.H",
        targets_val = 0.5431367893,
        params = c("hiv.test.rate_2"), # target: 0.0038
        initial_proposals = dplyr::tibble(
          hiv.test.rate_2 = seq(0.002, 0.006, length.out = n_sims), # ADD THESE
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
      ),

      job3 = list(
        targets = "cc.dx.O",
        targets_val = 0.5988779867,
        params = c("hiv.test.rate_3"), # target: 0.0038
        initial_proposals = dplyr::tibble(
          hiv.test.rate_2 = seq(0.002, 0.006, length.out = n_sims), # ADD THESE
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
      ),

      job4 = list(
        targets = "cc.dx.W",
        targets_val = 0.5614905982,
        params = c("hiv.test.rate_4"), # target: 0.0069
        initial_proposals = dplyr::tibble(
          hiv.test.rate_3 = seq(0.002, 0.006, length.out = n_sims),
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 5)
      ),

      # Linked to care
      job5 = list(
        targets = paste0("cc.linked1m.", c("B", "H", "O", "W")),
        targets_val = c(.828, 0.867, 0.875, 0.936),
        params = paste0("tx.init.rate_", 1:4),
        initial_proposals = dplyr::tibble(
          tx.init.rate_1 = sample(seq(0.1, 0.6, length.out = n_sims)),
          tx.init.rate_2 = sample(tx.init.rate_1),
          tx.init.rate_3 = sample(tx.init.rate_1),
          tx.init.rate_4 = sample(tx.init.rate_1)
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 3)
      ),

      # ART Discontinuation
      job5 = list(
        targets = paste0("cc.vsupp.", c("B", "H", "O", "W")),
        targets_val = c(0.571, 0.675, 0.586, 0.617),
        params = paste0("tx.halt.full.or_", 1:4),
        initial_proposals = dplyr::tibble(
          tx.halt.full.or_1 = sample(seq(0.2, 0.6, length.out = n_sims)),
          tx.halt.full.or_2 = sample(tx.halt.full.or_1),
          tx.halt.full.or_3 = sample(tx.halt.full.or_1),
          tx.halt.full.or_4 = sample(tx.halt.full.or_1)
        ),
        make_next_proposals = swfcalib::make_shrink_proposer(n_sims, shrink = 2),
        get_result = swfcalib::determ_poly_end(0.001, poly_n = 3)
      )
    ),
    wave2 = list(
      job1 = list(
        targets = paste0("i.prev.dx.", c("B", "H", "O", "W")),
        targets_val = c(0.5398688, 0.1417571, 0.1303336, 0.0681581),
        params = paste0("hiv.trans.scale_", 1:4),
        initial_proposals = dplyr::tibble(
          hiv.trans.scale_1 = sample(seq(18, 24, length.out = n_sims)),
          hiv.trans.scale_2 = sample(seq(2, 4, length.out = n_sims)),
          hiv.trans.scale_3 = sample(seq(2, 4, length.out = n_sims)),
          hiv.trans.scale_4 = sample(seq(1, 2, length.out = n_sims))
        ),
        make_next_proposals =
          swfcalib::make_proposer_se_range(n_sims, retain_prop = 0.3),
        get_result = swfcalib::determ_end_thresh(
          thresholds = rep(0.02, 3),
          n_enough = 100
        )
      )
    )
  )
)


# Workflow ---------------------------------------------------------------------

# Use preconfigured HPC settings
hpc_configs <- EpiModelHPC::swf_configs_rsph(
  partition = "epimodel",
  r_version = "4.3",
  mail_user = "<your_user>@emory.edu"
)

wf <- create_workflow(
  wf_name = "Vignette_auto_calib",
  default_sbatch_opts = hpc_configs$default_sbatch_opts
)

# Calibration step 1
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_do_call(
    what = swfcalib::calibration_step1,
    args = list(
      n_cores = 8,
      calib_object = calib_object
    ),
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "cpus-per-task" = 8,
    "time" = "00:20:00",
    "mem-per-cpu" = "4G",
    "mail-type" = "FAIL"
  )
)

# Calibration step 2
batch_size <- 8
batch_numbers <- swfcalib:::get_batch_numbers(calib_object, batch_size)
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_map(
    FUN = swfcalib::calibration_step2,
    batch_num = batch_numbers,
    setup_lines = hpc_configs$r_loader,
    max_array_size = 500,
    MoreArgs = list(
      n_cores = batch_size,
      n_batches = max(batch_numbers),
      calib_object = calib_object
    )
  ),
  sbatch_opts = list(
    "cpus-per-task" = batch_size,
    "time" = "05:00:00",
    "mem-per-cpu" = "5G",
    "mail-type" = "FAIL"
  )
)

# Calibration step 3
wf <- add_workflow_step(
  wf_summary = wf,
  step_tmpl = step_tmpl_do_call(
    what = swfcalib::calibration_step3,
    args = list(
      calib_object = calib_object
    ),
    setup_lines = hpc_configs$r_loader
  ),
  sbatch_opts = list(
    "cpus-per-task" = 1,
    "time" = "00:20:00",
    "mem-per-cpu" = "8G",
    "mail-type" = "END"
  )
)
