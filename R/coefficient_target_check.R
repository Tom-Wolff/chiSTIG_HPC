#####################################################################################
#                                                                                   #
#    C H E C K I N G   T A R G E T   S T A T S   A N D   C O E F F I C I E N T S    #
#                                                                                   #
#####################################################################################

# The files needed for inspection are located on the HPC at `/projects/p32153/ChiSTIG_model/data/input/`
# File names are as follows:

### ERGM Fit and Target Stats - Control: `basic_netest-local.rds`
### ERGM Fit and Target Stats - Venues Only: `venue_only_netest-local.rds`
### ERGM Fit and Target Stats - Apps Only: `apps_only_netest-local.rds`
### ERGM Fit and Target Stats - Venues and Apps: `venues_apps_netest-local.rds`

### `netstats` Object Containing Attributes of Initial Population - `netstats-local.rds`

# Once you've downloaded these files, specify the diretory you're storing them in here:
this_dir <- "~/Desktop/chiSTIG_hpc/data/intermediate/estimates/"

###################################################################
#    T A R G E T   S T A T S   A N D   C O E F F I C I E N T S    #
###################################################################

# Control Model
control <- readRDS(paste(this_dir, "basic_netest-local.rds", sep = ""))
### Target Stats
control$fit_main$target.stats
control$fit_casl$target.stats
control$fit_inst$target.stats
### ERGM Coefficients
control$fit_main$coef.form
control$fit_casl$coef.form
control$fit_inst$coef.form

# Venues Only
venues <- readRDS(paste(this_dir, "venue_only_netest-local.rds", sep = ""))
### Target Stats
venues$fit_main$target.stats
venues$fit_casl$target.stats
venues$fit_inst$target.stats
### ERGM Coefficients
venues$fit_main$coef.form
venues$fit_casl$coef.form
venues$fit_inst$coef.form

# Apps Only
apps <- readRDS(paste(this_dir, "apps_only_netest-local.rds", sep = ""))
### Target Stats
apps$fit_main$target.stats
apps$fit_casl$target.stats
apps$fit_inst$target.stats
### ERGM Coefficients
apps$fit_main$coef.form
apps$fit_casl$coef.form
apps$fit_inst$coef.form

# Venues and Apps
both <- readRDS(paste(this_dir, "venues_apps_netest-local.rds", sep = ""))
### Target Stats
both$fit_main$target.stats
both$fit_casl$target.stats
both$fit_inst$target.stats
### ERGM Coefficients
both$fit_main$coef.form
both$fit_casl$coef.form
both$fit_inst$coef.form


#########################################
#    A G E   D I S T R I B U T I O N    #
#########################################

# Read in `netstats` object
netstats <- readRDS(paste(this_dir, "netstats-local.rds", sep = ""))

# View age distribution that used to fit ERGMs
hist(netstats$attr$age)

# This is how I level the age distribution before feeding it into the simulation
netstats$attr$age <- sample(16:29, length(netstats$attr$age), replace = TRUE)
netstats$attr$age <- netstats$attr$age + sample(1:1000, length(netstats$attr$age), replace = TRUE)/1000

# View updated age distribution
hist(netstats$attr$age)
