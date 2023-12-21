### 0. Set up python and R environments ###
# load reticulate
library(reticulate)
# point to python binary
use_python("/opt/homebrew/bin/python3.10")
# load python library/script being used
source_python("chistig_reticulate_python_script.py")


### 1. Set location of necessary files, and setup the apps and venues variables ###
chistig_synthpop_repo <- "/Users/rimersara/dev/repos/ChiSTIG/ChiSTIG_synthpop/"
# venue attendance matrix file
afile <- paste(chistig_synthpop_repo,"v4/venue-attendance_v4_0.csv", sep="")
# synthetic egos file
efile <- paste(chistig_synthpop_repo,"v4/egos_v4_0.csv", sep="")
# app category file
atypefile <- paste(chistig_synthpop_repo, "data/raw/v4/empop_appid_type_def.csv", sep="")
# venue category file
vtypefile <- paste(chistig_synthpop_repo, "data/raw/v4/empop_venueid_type_def.csv", sep="")

#format a dataframe of apps and egos
segos2apps <- categorize_apps(efile, atypefile)

# setup dataframe of venuetypes
vtypes <- categorize_venues(vtypefile)


### 2. Attend venues and use apps ##

# Pass the attendance matrix file and the dataframe of venue types categorized above.
# This function will Have ALL of the synthetic egos their venues for 7 days according to their corresponding
# attendance frequencies.
# Returns a dataframe where first column is egoid (only the numeric value), and columns of venues attended,
# including venues_all, venues_dating, and venues_nondating.
egos2venues <- attend_venues(afile, vtypes)

# Create and pass a dataframe of two columns, where first column is egoid of the egos,
# and the second column is the current relationship status of the corresponding ego,
# where 0 is "not in a serious relationship" and 1 is "in a serious relationship."
# Returns a dataframe where first column is egoid (only the numeric value), and columns of apps used,
# including apps_all, apps_dating, and apps_nondating
#NOTE: the following creates an example dataframe of relationship status to pass to the `use_apps` function
rel_status_df <- create_test_rel_status(efile)
egos2apps <- use_apps(rel_status_df, segos2apps)

