###########################################
#                                         #
#    R E G R E S S I O N   M O D E L S    #
#                                         #
###########################################

###################
#    S E T U P    #
###################

source("./R/utils-targets.R")

# Function for processing data
load_diag <- function(dir, start_time) {

  # Get all file names
  files <- list.files(dir)
  # Extract which treatment was applied from filename
  treatment <- unlist(stringr::str_extract_all(files, "^[a-z]*"))
  # Get full directory path to each file
  files <- paste(dir, files, sep = "")

  # Loop over each file and process outcome measures
  for (i in 1:length(files)) {

    this_targets <- tibble::as_tibble(readRDS(files[[i]])) %>%
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
             treatment = treatment[[i]]) %>%
      filter(time > start_time)

    if (i == 1) {
      sim_targets <- this_targets
    } else {
      sim_targets <- dplyr::bind_rows(sim_targets, this_targets)
    }
  }

  # Create dummy variables for each treatment
  sim_targets <- sim_targets %>%
    mutate(control = ifelse(treatment == "control", TRUE, FALSE),
           venues = ifelse(treatment == "venues", TRUE, FALSE),
           apps = ifelse(treatment == "apps", TRUE, FALSE),
           both = ifelse(treatment == "both", TRUE, FALSE))

  return(sim_targets)
}

###########################
#    L O A D   D A T A    #
###########################

output_data <- load_diag(dir = "./data/quest_output/", start_time = 3120)


#####################
#    M O D E L S    #
#####################

# Annualized Incidence by Race
# HIV Prevalence by Race
# Disparities in annualized incidence and HIV Prevalence

# Example regression model for comparison
summary(lm(i.prev.B ~ venues + apps + both + time
           + venues*time + apps*time + both*time, data = output_data))

# Will want to compare model fits to see if interactions are worth including


#####################################
#    V I S U A L I Z A T I O N S    #
#####################################


output_data %>%
  group_by(treatment, time) %>%
  summarize(i.prev.B = median(i.prev.B)) %>%
  ungroup() %>%
  ggplot(aes(x = time, y = i.prev.B, color = as.factor(treatment))) +
  geom_line() +
  theme_minimal()



#########################################
#                                       #
#    H I G H - R I S K   V E N U E S    #
#                                       #
#########################################

apps_test <- readRDS("./data/quest_output/apps_simtest_5mar2024.rds")
View(bind_rows(apps_test$el.cuml))
