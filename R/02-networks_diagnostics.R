#
## 02. Network Model Diagnostics
##

# Libraries  -------------------------------------------------------------------
library("EpiModelHIV")

# Settings ---------------------------------------------------------------------
context <- if (!exists("context")) "local" else context
source("R/utils-0_project_settings.R")

if (context == "local") {
  ncores <- 2
  nsims  <- 10
} else if (context == "hpc") {
  # hpc configs not passed through the workflow
} else {
  stop("The `context` variable must be set to either 'local' or 'hpc'")
}

nsteps <- 500

# Main -------------------------------------------------------------------------
est <- readRDS(paste0(est_dir, "basic_netest-level-", context, ".rds"))
fit_main <- est[["fit_main"]]

# model_main_dx <- ~edges +
#   nodematch("age.grp", diff = TRUE) +
#   nodefactor("age.grp", levels = TRUE) +
#   nodematch("race", diff = TRUE) +
#   nodefactor("race", levels = TRUE) +
#   nodefactor("deg.casl", levels = TRUE) +
#   degrange(from = 3) +
#   concurrent +
#   nodematch("role.class", diff = TRUE) +
#   degree(0:3)

model_main_dx <- ~ edges +
  concurrent +
  nodefactor("race", levels = -4) +
  nodematch("race") +
  # nodematch("race", diff = TRUE, levels = 1) +
  nodefactor("age.grp", levels= 1) +
  nodematch("age.grp") +
  nodefactor("deg.casl", levels= -1)# +
  # fuzzynodematch("venues.all", binary=TRUE) +
  # fuzzynodematch("apps.all", binary = TRUE)

dx_main <- netdx(
  fit_main,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  nwstats.formula = model_main_dx,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5),
  set.control.tergm = control.simulate.formula.tergm(MCMC.burnin.min = 2e5)
)

dx_main_static <- netdx(
  fit_main,
  dynamic = FALSE,
  nsims = 10000,
  nwstats.formula = model_main_dx,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5)
)

dx <- list(dx_main = dx_main, dx_main_static = dx_main_static)
saveRDS(dx, paste0(diag_dir, "netdx-main-", context, ".rds"))
rm(dx, dx_main, dx_main_static)

# Casual -----------------------------------------------------------------------
fit_casl <- est[["fit_casl"]]

model_casl_dx <- ~edges +
  concurrent +
  fuzzynodematch("venues.all", binary=TRUE)

dx_casl <- netdx(
  fit_casl,
  nsims = nsims,
  ncores = ncores,
  nsteps = nsteps,
  nwstats.formula = model_casl_dx,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5),
  set.control.tergm = control.simulate.formula.tergm(MCMC.burnin.min = 2e5)
)

dx_casl_static <- netdx(
  fit_casl,
  dynamic = FALSE,
  nsims = 10000,
  nwstats.formula = model_casl_dx,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5)
)

dx <- list(dx_casl = dx_casl, dx_casl_static = dx_casl_static)
saveRDS(dx, paste0(diag_dir, "netdx-casl-", context, ".rds"))
rm(dx, dx_casl, dx_casl_static)

# One-Off ----------------------------------------------------------------------
fit_inst <- est[["fit_inst"]]

model_inst_dx <- ~edges +
  fuzzynodematch("venues.all", binary=TRUE)

dx_inst <- netdx(
  fit_inst,
  nsims = 10000,
  dynamic = FALSE,
  nwstats.formula = model_inst_dx,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5)
)

dx <- list(dx_inst = dx_inst)
saveRDS(dx, paste0(diag_dir, "netdx-inst-", context, ".rds"))
