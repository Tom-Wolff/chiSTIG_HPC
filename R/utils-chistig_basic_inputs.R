
# if (!context %in% c("local", "hpc")) {
#   stop("The `context` variable must be set to either 'local' or 'hpc'")
# }

context <- "hpc"

epistats <- readRDS("./data/intermediate/estimates/epistats-local.rds")
netstats <- readRDS("./data/intermediate/estimates/netstats-local.rds")
# netstats <- readRDS("data/intermediate/estimates/netstats-local.rds")
# Is the aging out of the older initial nodes driving HIV extinction?
# Let's level out the age distribution and find out
netstats$attr$age <- sample(16:29, length(netstats$attr$age), replace = TRUE)
netstats$attr$age <- netstats$attr$age + sample(1:1000, length(netstats$attr$age), replace = TRUE)/1000

##### `est` object for basic ERGMs (no venues or apps)
path_to_est <- "./data/intermediate/estimates/basic_netest-local.rds"
path_to_restart <- paste0(est_dir, "restart-", context, ".rds")

# `netsim` Parameters
prep_start <- 52 * 2
param <- EpiModel::param.net(
  data.frame.params = readr::read_csv("data/input/params_chistig_feb19.csv"),
  netstats          = netstats,
  epistats          = epistats,
  prep.start        = Inf,
  riskh.start       = Inf
)

# Initial conditions (default prevalence initialized in epistats)
# For models without bacterial STIs, these must be initialized here
# with non-zero values
# init <- init_msm(
#   prev.ugc = 0.1,
#   prev.rct = 0.1,
#   prev.rgc = 0.1,
#   prev.uct = 0.1
# )
init <- EpiModelHIV::init_msm()

