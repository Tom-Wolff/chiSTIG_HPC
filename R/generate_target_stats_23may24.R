###########################################################################
##### WRITE OUT THE APP TYPE FILE TO SEND TO ARGONNE
###########################################################################

if (!require(tidyverse))  install.packages("tidyverse")
if (!require(survey)) install.packages("survey")
library(tidyverse); library(survey)

# Path to your Northwestern Northwestern FSM chiSTIG data shared drive folder - change this
# based on your local computer
file_loc <- "/Volumes/fsmresfiles/MSS/Birkett_Lab/Projects/chiSTIG/Data/"

# read in ego data. Pasting together your file_loc assigned above with the filepath
# within the chiSTIG folder structure to hopefully avoid needing to change all of the filepaths
# throughout the script
egodat <- read.csv(paste0(file_loc,
                          "Internal/chiSTIG/Data Sent to Argonne/Target stats data/ego_level_target_stats_data.csv"))

# Read in synthetic population
synthpop <- read.csv("./data/synthpop_gen/synthpop-demo-counts_v3_tom.csv")

alter_data <- read.csv(paste0(file_loc,
                              "Internal/chiSTIG/Data Sent to Argonne/Target stats data/alter_level_target_stats_data.csv"))

##################################################################################
##################################################################################
##################################################################################
############# Generate weights to account for difference in distribution of the
############# empirical population compared to the synthetic population
##################################################################################
##################################################################################
##################################################################################

# code synthetic population to match empirical data
synthpop <- synthpop %>%
    rename(RaceEth4 = race_ethnicity,ego_agecat = age, total = population) %>%
    mutate(RaceEth4 = recode(RaceEth4,"blackNH" = "BlackNH",
                             "whiteNH" = "WhiteNH", "hispanic" = "Latinx",
                             "otherNH" = "OtherNH"),
           ego_agecat = recode(ego_agecat, "21to29" = "21-29", "16to20" = "16-20"))

# bin empirical population by race/ethnicity and age
emppop <- egodat %>% group_by(RaceEth4,ego_agecat) %>%
    summarise(total_emp = n())

# merge the synthetic population and empirical population bins
weightdat <- merge(emppop,synthpop)

# weight1 = the relative size of the synthpop in that demographic group compared to the
# size of that demographic group in the empirical population
weightdat$weight1 <- weightdat$total / weightdat$total_emp
# weight2 = the overall relative size of the empirical popualtion to the size of the
# synthetic population
weightdat$weight2 <- sum(weightdat$total_emp) / sum(weightdat$total)
# combining weights 1 and 2 into weights
weightdat$weights <- weightdat$weight1 * weightdat$weight2
# "newpop" is the empirical population weighted to have the same relative makeup by
# 8 demographic categories as are in the synthetic population (can use newpop to check
# that the relative proportions are right!)
weightdat$newpop <- weightdat$weights * weightdat$total_emp

# merge the weights we created back into the ego dataset
egodat <- merge(egodat, select(weightdat,ego_agecat,RaceEth4,weights))

# merge the weights into the alter data as well
alter_data <- merge(alter_data,select(egodat,egoid,weights), all.x=TRUE,
                    by=c("egoid"))


##################################################################################
##################################################################################
##################################################################################
########## EXPLANATION OF TARGET STATS CALCULATION
########## The following section where the target stats are calculated are structured
########## in a not-very-nice way, basically because things don't fit nicely into
########## a table without some weird syntax
########## What happens, over and over, is calling svymean or svyciprop from the survey
########## package to calculate either mean degree (svymean) or proportions (svyciprop)
########## accounting for the weights we've generated above.  The weights account for the
########## the differences in the relative population sizes by age category & race/ethnicity
########## in the empirical data rather than in the synthetic population.  This is necessary
########## to reconcile between the stratified mean degree & proportion estimates and the overall
########## mean estimates when inputting target stats into the netest function (i.e.
########## internal consistency within our target stats).
########## svymean & svyciprop take an equation and a survey design object. The survey design object
########## is created using svydesign() where you tell the ego's ID, the weights, and the dataset.
########## The equation is "~" followed by the outcome variable we are interested in
########## getting the mean or survey-weighted proportion for.
########## For instance, if we were to get the mean overall degree, we'd do:
########## svymean(~ num_partners, svydesign(id = ~ egoid, weights = ~ weights, data = egodat), na.rm = TRUE)
########## (the na.rm=TRUE just specifies that we want to ignore NA's when calculating our mean,
########## if na.rm=FALSE, the mean is NA if there are any NA values)
########## To get the 95% CIs on svymean and svyciprop we nest those two commands in confint()
########## Below there are many different specifications of the svydesign for different
########## subsets of the data, since we are frequently interested in mean degree
########## given that you have a specific characteristic etc
##################################################################################
##################################################################################
##################################################################################

# row names for the different stats we want to output
ind_names <- c("Mean degree","Mean degree BlackNH","Mean degree WhiteNH","Mean degree Latinx",
               "Mean degree OtherNH","Mean degree age 16-20","Mean degree age 21-29",
               "Prop concurrent","Prop same race","Prop same race | BlackNH","Prop same race | WhiteNH",
               "Prop same race | Latinx","Prop same race | OtherNH",
               "Prop same age cat","Prop same age cat | 16-20","Prop same age cat | 21-29",
               "Prop met in person (any)",
               "Prop met in Bar/Club", "Prop met in Bathhouse", "Prop met other in-person",
               "Prop met online hookup site","Prop met online soc network","Prop met elsewhere",
               "Mean degree | casual degree 0", "Mean degree | casual degree 1", "Mean degree | casual degree 2+",
               "Mean degree | main degree 0", "Mean degree | main degree 1 or 2",
               "Mean 1-time | persistent degree 0", "Mean 1-time | persistent degree 1",
               "Mean 1-time | persistent degree 2", "Mean 1-time | persistent degree 3+")

# data frame to hold our target stats - currently empty
target_stats_out <- data.frame(ind = ind_names,
                               all_p = NA, all_p_lb = NA, all_p_ub = NA,
                               main_p = NA, main_p_lb = NA, main_p_ub = NA,
                               cas_p = NA, cas_p_lb = NA, cas_p_ub = NA,
                               ot_p = NA, ot_p_lb = NA, ot_p_ub = NA,
                               tp_all_p = NA, tp_all_p_lb = NA, tp_all_p_ub = NA,
                               tp_main_p = NA, tp_main_p_lb = NA, tp_main_p_ub = NA,
                               tp_cas_p = NA, tp_cas_p_lb = NA, tp_cas_p_ub = NA,
                               tp_ot_p = NA, tp_ot_p_lb = NA, tp_ot_p_ub = NA)

# using a for loop to get through the different variables we care about
# so specifying which indices we'll be wanting to insert things into with nbottoms and ntops
nbottoms <- seq(2,by=3,length.out=8)
ntops <- seq(4,by=3,length.out=8)

# list of our variables that we want to use as the equation for svymean() within our for loop
# variables defined here:
# num_partners = number of ongoing anal sex partners 3 months prior to interview (regardless of type of partner)
# num_ser_partners = number of ongoing serious anal sex partners 3 months prior to interview
# num_cas_partners = number of ongoing casual anal sex partners 3 months prior to interview
# num_one_time = number of one-time anal sex partners in the ENTIRE 6 MONTHS PRIOR TO INTERVIEW
# tp_num_partners = number of ongoing anal sex partners 3 months prior to interview who are in the
#                   "target population" - i.e. excluding partners who aren't cismen or transwomen and
#                   who aren't in our target ages (16-29) - regardless of partner type
# tp_num_ser_partners = number of ongoing serious anal sex partners 3 months prior to interview in the "target population"
# tp_num_cas_partners = number of ongoing casual anal sex partners 3 months prior to interview in the "target population"
# tp_num_one_time = number of one-time anal sex partners in the ENTIRE 6 MONTHS PRIOR TO INTERVIEW in the "target population"

eqs <- c(as.formula("~ num_partners"),as.formula("~ num_ser_partners"),as.formula("~ num_cas_partners"),
         as.formula("~ num_one_time"), as.formula("~ tp_num_partners"), as.formula("~ tp_num_ser_partners"),
         as.formula("~ tp_num_cas_partners"), as.formula("~ tp_num_one_time"))

# for loop to insert
for(i in 1:length(eqs)) {
    # which columns do we want to fill?
    ind1 <- nbottoms[i]
    ind2 <- ntops[i]
    # svydesign for entire population
    des1 <- svydesign(id=~egoid, weights=~weights, data=egodat)
    # output mean degree + 95% lower and upper bound
    target_stats_out[target_stats_out$ind=="Mean degree",ind1:ind2] <- c(svymean(eqs[[i]],des1,na.rm=TRUE),
                                                                         confint(svymean(eqs[[i]],des1,na.rm=TRUE)))
    # svydesign for population restricted to those who are black non-hispanic... etc below
    des2 <- svydesign(id=~egoid, weights=~weights, data=filter(egodat,RaceEth4=="BlackNH"))
    target_stats_out[target_stats_out$ind=="Mean degree BlackNH",ind1:ind2] <- c(svymean(eqs[[i]],des2,na.rm=TRUE),
                                                                         confint(svymean(eqs[[i]],des2,na.rm=TRUE)))
    des3 <- svydesign(id=~egoid, weights=~weights, data=filter(egodat,RaceEth4=="WhiteNH"))
    target_stats_out[target_stats_out$ind=="Mean degree WhiteNH",ind1:ind2] <- c(svymean(eqs[[i]],des3,na.rm=TRUE),
                                                                                 confint(svymean(eqs[[i]],des3,na.rm=TRUE)))
    des4 <- svydesign(id=~egoid, weights=~weights, data=filter(egodat,RaceEth4=="Latinx"))
    target_stats_out[target_stats_out$ind=="Mean degree Latinx",ind1:ind2] <- c(svymean(eqs[[i]],des4,na.rm=TRUE),
                                                                                 confint(svymean(eqs[[i]],des4,na.rm=TRUE)))
    des5 <- svydesign(id=~egoid, weights=~weights, data=filter(egodat,RaceEth4=="OtherNH"))
    target_stats_out[target_stats_out$ind=="Mean degree OtherNH",ind1:ind2] <- c(svymean(eqs[[i]],des5,na.rm=TRUE),
                                                                                 confint(svymean(eqs[[i]],des5,na.rm=TRUE)))
    des6 <- svydesign(id=~egoid, weights=~weights, data=filter(egodat,ego_agecat=="16-20"))
    target_stats_out[target_stats_out$ind=="Mean degree age 16-20",ind1:ind2] <- c(svymean(eqs[[i]],des6,na.rm=TRUE),
                                                                                 confint(svymean(eqs[[i]],des6,na.rm=TRUE)))
    des7 <- svydesign(id=~egoid, weights=~weights, data=filter(egodat,ego_agecat=="21-29"))
    target_stats_out[target_stats_out$ind=="Mean degree age 21-29",ind1:ind2] <- c(svymean(eqs[[i]],des7,na.rm=TRUE),
                                                                                 confint(svymean(eqs[[i]],des7,na.rm=TRUE)))
}

# Repeat the above, but for cross-partnership type degree for main partnerships
# getting out the mean main degree given that you have 0, 1 or 2+ casual partners
nbottoms <- c(5,17)
ntops <- nbottoms + 2

eqs <- c(as.formula("~ num_ser_partners"), as.formula("~ tp_num_ser_partners"))

# here we're specifying the svydesign outside of the forloop - but stratifying to different
# levels of casual partnerships
des_list <- list(svydesign(id=~egoid, weights=~weights, data=filter(egodat,cas_part_cat=="0")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,cas_part_cat=="1")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,cas_part_cat=="2+")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,tp_cas_part_cat=="0")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,tp_cas_part_cat=="1")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,tp_cas_part_cat=="2+")))

for(i in 1:length(eqs)) {
    ind1 <- nbottoms[i]
    ind2 <- ntops[i]
    target_stats_out[target_stats_out$ind=="Mean degree | casual degree 0",ind1:ind2] <- c(svymean(eqs[[i]],des_list[[(i-1)*3+1]],na.rm=TRUE),
                                                                                 confint(svymean(eqs[[i]],des_list[[(i-1)*3+1]],na.rm=TRUE)))
    target_stats_out[target_stats_out$ind=="Mean degree | casual degree 1",ind1:ind2] <- c(svymean(eqs[[i]],des_list[[(i-1)*3+2]],na.rm=TRUE),
                                                                                           confint(svymean(eqs[[i]],des_list[[(i-1)*3+2]],na.rm=TRUE)))
    target_stats_out[target_stats_out$ind=="Mean degree | casual degree 2+",ind1:ind2] <- c(svymean(eqs[[i]],des_list[[(i-1)*3+3]],na.rm=TRUE),
                                                                                           confint(svymean(eqs[[i]],des_list[[(i-1)*3+3]],na.rm=TRUE)))
}

# Repeat the above but for cross-partnership type degree for casual partners
# getting out the mean casual degree given that you have 0 or "1 or 2" main partners
nbottoms <- c(8,20)
ntops <- nbottoms + 2

eqs <- c(as.formula("~ num_cas_partners"), as.formula("~ tp_num_cas_partners"))

des_list <- list(svydesign(id=~egoid, weights=~weights, data=filter(egodat,ser_part_cat=="0")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,ser_part_cat=="1 or 2")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,tp_ser_part_cat=="0")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,tp_ser_part_cat=="1 or 2")))

for(i in 1:length(eqs)) {
    ind1 <- nbottoms[i]
    ind2 <- ntops[i]
    target_stats_out[target_stats_out$ind=="Mean degree | main degree 0",ind1:ind2] <- c(svymean(eqs[[i]],des_list[[(i-1)*2+1]],na.rm=TRUE),
                                                                                           confint(svymean(eqs[[i]],des_list[[(i-1)*2+1]],na.rm=TRUE)))
    target_stats_out[target_stats_out$ind=="Mean degree | main degree 1 or 2",ind1:ind2] <- c(svymean(eqs[[i]],des_list[[(i-1)*2+2]],na.rm=TRUE),
                                                                                           confint(svymean(eqs[[i]],des_list[[(i-1)*2+2]],na.rm=TRUE)))
}

# Repeat the above for cross-partnership type degree for one-time partners
# getting out the mean number of one-time partners over the 6 month period given that
# you have 0, 1, 2 or 3+ persistent partners (persistent partners = main or casual)
nbottoms <- c(11,23)
ntops <- nbottoms + 2

eqs <- c(as.formula("~ num_one_time"), as.formula("~ tp_num_one_time"))

des_list <- list(svydesign(id=~egoid, weights=~weights, data=filter(egodat,pers_part_cat=="0")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,pers_part_cat=="1")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,pers_part_cat=="2")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,pers_part_cat=="3+")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,tp_pers_part_cat=="0")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,tp_pers_part_cat=="1")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,tp_pers_part_cat=="2")),
                 svydesign(id=~egoid, weights=~weights, data=filter(egodat,tp_pers_part_cat=="3+")))

for(i in 1:length(eqs)) {
    ind1 <- nbottoms[i]
    ind2 <- ntops[i]
    target_stats_out[target_stats_out$ind=="Mean 1-time | persistent degree 0",ind1:ind2] <- c(svymean(eqs[[i]],des_list[[(i-1)*4+1]],na.rm=TRUE),
                                                                                           confint(svymean(eqs[[i]],des_list[[(i-1)*4+1]],na.rm=TRUE)))
    target_stats_out[target_stats_out$ind=="Mean 1-time | persistent degree 1",ind1:ind2] <- c(svymean(eqs[[i]],des_list[[(i-1)*4+2]],na.rm=TRUE),
                                                                                           confint(svymean(eqs[[i]],des_list[[(i-1)*4+2]],na.rm=TRUE)))
    target_stats_out[target_stats_out$ind=="Mean 1-time | persistent degree 2",ind1:ind2] <- c(svymean(eqs[[i]],des_list[[(i-1)*4+3]],na.rm=TRUE),
                                                                                           confint(svymean(eqs[[i]],des_list[[(i-1)*4+3]],na.rm=TRUE)))
    target_stats_out[target_stats_out$ind=="Mean 1-time | persistent degree 3+",ind1:ind2] <- c(svymean(eqs[[i]],des_list[[(i-1)*4+4]],na.rm=TRUE),
                                                                                            confint(svymean(eqs[[i]],des_list[[(i-1)*4+4]],na.rm=TRUE)))
}

# do some quick clean up since confint on svymean can have a lower bound <0 (treating it as
# a continuous variable) - set any LBs<0 to be 0
target_stats_out$main_p_lb[target_stats_out$main_p_lb<0] <- 0
target_stats_out$tp_main_p_lb[target_stats_out$tp_main_p_lb<0] <- 0
target_stats_out$cas_p_lb[target_stats_out$cas_p_lb<0] <- 0
target_stats_out$tp_cas_p_lb[target_stats_out$tp_cas_p_lb<0] <- 0
target_stats_out$tp_ot_p_lb[target_stats_out$tp_ot_p_lb<0] <- 0

# Repeat the above, but for concurrent partnerships
nbottoms <- c(2,5,8,14,17,20)
ntops <- nbottoms + 2

# These variables are all 0/1 variables to say whether the individual had 2+ ongoing partnerships
# at 3 months prior to interview, by different types (all, ser - serious, or cas - casual) and whether or not
# in the target population (tp_ indicates this is partners within the target population)
eqs <- c(as.formula("~ all_p_conc"), as.formula("~ ser_p_conc"), as.formula("~ cas_p_conc"),
         as.formula("~ tp_all_p_conc"), as.formula("~ tp_ser_p_conc"), as.formula("~ tp_cas_p_conc"))

for(i in 1:length(eqs)) {
    ind1 <- nbottoms[i]
    ind2 <- ntops[i]
    des1 <- svydesign(id=~egoid, weights=~weights, data=egodat)
    target_stats_out[target_stats_out$ind=="Prop concurrent",ind1:ind2] <- c(svyciprop(eqs[[i]],des1,na.rm=TRUE),
                                                                             confint(svyciprop(eqs[[i]],des1,na.rm=TRUE)))
}

# output from alter_data (previous variables came from the ego data):
# proportion of partnerships met in-person
# proportion of partnerships met in a sex club/bathhouse
# proportion of partnerships met in a bar/club
# proportion of partnerships met elsewhere in person
# proportion of partnerships met online via sex seeking app
# proportion of partnerships met online via social networking app
# proportion of partnerships met "Somewhere else"
# proportion of partnerships in same age group
# proportion of partnerships in same race/ethnicity
#
# These use all partners in the past 6 months

nbottoms <- seq(2,by=3,length.out=8)
ntops <- seq(4,by=3,length.out=8)

# here specifying the subsets of the alter_data dataset that we'll want to use in our svydesign() calls below
datasets <- list(# all anal sex partners
                 filter(alter_data,dyad_edge.analSex=="Anal Sex"),
                 # all serious anal sex partners
                 filter(alter_data,dyad_edge.analSex=="Anal Sex" & seriousRel==1),
                 # all casual anal sex partners
                 filter(alter_data,dyad_edge.analSex=="Anal Sex" & casual==1),
                 # all one-time anal sex partners
                 filter(alter_data,dyad_edge.analSex=="Anal Sex" & onetime==1),
                 # all anal sex partners where the partner is in the target population
                 filter(alter_data,dyad_edge.analSex=="Anal Sex" & partner_in_model==1),
                 # all serious anal sex partners where the partner is int he target population
                 filter(alter_data,dyad_edge.analSex=="Anal Sex" & seriousRel==1 & partner_in_model==1),
                 # all casual anal sex partners where the partner is in the target population
                 filter(alter_data,dyad_edge.analSex=="Anal Sex" & casual==1 & partner_in_model==1),
                 # all one-time partners anal sex partners where the partner is in the target population
                 filter(alter_data,dyad_edge.analSex=="Anal Sex" & onetime==1 & partner_in_model==1))

for(i in 1:length(datasets)) {
    ind1 <- nbottoms[i]
    ind2 <- ntops[i]
    des1 <- svydesign(id = ~egoid, weights = ~ weights, data=datasets[[i]])

    target_stats_out[target_stats_out$ind=="Prop same race",ind1:ind2] <- c(svyciprop(~samerace,des1,na.rm=TRUE),
                                                                      confint(svyciprop(~samerace,des1,na.rm=TRUE)))
    des2 <- svydesign(id = ~egoid, weights = ~ weights, data=filter(datasets[[i]],RaceEth4_ego=="BlackNH"))
    target_stats_out[target_stats_out$ind=="Prop same race | BlackNH",ind1:ind2] <- c(svyciprop(~samerace,des2,na.rm=TRUE),
                                                                            confint(svyciprop(~samerace,des2,na.rm=TRUE)))
    des3 <- svydesign(id = ~egoid, weights = ~ weights, data=filter(datasets[[i]],RaceEth4_ego=="WhiteNH"))
    target_stats_out[target_stats_out$ind=="Prop same race | WhiteNH",ind1:ind2] <- c(svyciprop(~samerace,des3,na.rm=TRUE),
                                                                            confint(svyciprop(~samerace,des3,na.rm=TRUE)))
    des4 <- svydesign(id = ~egoid, weights = ~ weights, data=filter(datasets[[i]],RaceEth4_ego=="Latinx"))
    target_stats_out[target_stats_out$ind=="Prop same race | Latinx",ind1:ind2] <- c(svyciprop(~samerace,des4,na.rm=TRUE),
                                                                            confint(svyciprop(~samerace,des4,na.rm=TRUE)))
    des5 <- svydesign(id = ~egoid, weights = ~ weights, data=filter(datasets[[i]],RaceEth4_ego=="OtherNH"))
    target_stats_out[target_stats_out$ind=="Prop same race | OtherNH",ind1:ind2] <- c(svyciprop(~samerace,des5,na.rm=TRUE),
                                                                            confint(svyciprop(~samerace,des5,na.rm=TRUE)))
    target_stats_out[target_stats_out$ind=="Prop same age cat",ind1:ind2] <- c(svyciprop(~sameage,des1,na.rm=TRUE),
                                                                         confint(svyciprop(~sameage,des1,na.rm=TRUE)))
    des6 <- svydesign(id = ~egoid, weights = ~ weights, data=filter(datasets[[i]],ego_agecat=="16-20"))
    target_stats_out[target_stats_out$ind=="Prop same age cat | 16-20",ind1:ind2] <- c(svyciprop(~sameage,des6,na.rm=TRUE),
                                                                                     confint(svyciprop(~sameage,des6,na.rm=TRUE)))
    des7 <- svydesign(id = ~egoid, weights = ~ weights, data=filter(datasets[[i]],ego_agecat=="21-29"))
    target_stats_out[target_stats_out$ind=="Prop same age cat | 21-29",ind1:ind2] <- c(svyciprop(~sameage,des7,na.rm=TRUE),
                                                                                      confint(svyciprop(~sameage,des7,na.rm=TRUE)))
    target_stats_out[target_stats_out$ind=="Prop met in person (any)",ind1:ind2] <- c(svyciprop(~met_any_inperson,des1,na.rm=TRUE),
                                                                                confint(svyciprop(~met_any_inperson,des1,na.rm=TRUE)))
    target_stats_out[target_stats_out$ind=="Prop met in Bar/Club",ind1:ind2] <- c(svyciprop(~met_barclub,des1,na.rm=TRUE),
                                                                            confint(svyciprop(~met_barclub,des1,na.rm=TRUE)))
    target_stats_out[target_stats_out$ind=="Prop met in Bathhouse",ind1:ind2] <- c(svyciprop(~met_bathhouse,des1,na.rm=TRUE),
                                                                             confint(svyciprop(~met_bathhouse,des1,na.rm=TRUE)))
    target_stats_out[target_stats_out$ind=="Prop met other in-person",ind1:ind2] <- c(svyciprop(~met_other_inperson,des1,na.rm=TRUE),
                                                                                confint(svyciprop(~met_other_inperson,des1,na.rm=TRUE)))
    target_stats_out[target_stats_out$ind=="Prop met online hookup site",ind1:ind2] <- c(svyciprop(~met_online_sex,des1,na.rm=TRUE),
                                                                                   confint(svyciprop(~met_online_sex,des1,na.rm=TRUE)))
    target_stats_out[target_stats_out$ind=="Prop met online soc network",ind1:ind2] <- c(svyciprop(~met_online_soc,des1,na.rm=TRUE),
                                                                                   confint(svyciprop(~met_online_soc,des1,na.rm=TRUE)))
    target_stats_out[target_stats_out$ind=="Prop met elsewhere",ind1:ind2] <- c(svyciprop(~met_other,des1,na.rm=TRUE),
                                                                          confint(svyciprop(~met_other,des1,na.rm=TRUE)))
}

target_stats_out <- target_stats_out %>%
    dplyr::rename("All partners" = all_p, "All partners, lower bound" = all_p_lb, "All partners, upper bound" = all_p_ub,
                  "Main partners" = main_p, "Main partners, lower bound" = main_p_lb, "Main partners, upper bound" = main_p_ub,
                  "Casual partners" = cas_p, "Casual partners, lower bound" = cas_p_lb, "Casual partners, upper bound" = cas_p_ub,
                  "One-time partners" = ot_p, "One-time partners, lower bound" = ot_p_lb, "One-time partners, upper bound" = ot_p_ub,
                  "All partners, modeled pop" = tp_all_p, "All partners, modeled pop, lower bound" = tp_all_p_lb, "All partners, modeled pop, upper bound" = tp_all_p_ub,
                  "Main partners, modeled pop" = tp_main_p, "Main partners, modeled pop, lower bound" = tp_main_p_lb, "Main partners, modeled pop, upper bound" = tp_main_p_ub,
                  "Casual partners, modeled pop" = tp_cas_p, "Casual partners, modeled pop, lower bound" = tp_cas_p_lb, "Casual partners, modeled pop, upper bound" = tp_cas_p_ub,
                  "One-time partners, modeled pop" = tp_ot_p, "One-time partners, modeled pop, lower bound" = tp_ot_p_lb, "One-time partners, modeled pop, upper bound" = tp_ot_p_ub)


write.csv(target_stats_out,
          paste0(file_loc,"Internal/chiSTIG/Data Sent to Argonne/Target stats data/target_stats_v8.csv"),
          row.names=FALSE)


############### Get target stats by multiplying to subpopulation sizes


names <- target_stats_out$ind

multiplier <- c(sum(synthpop$total)/2,
                sum(synthpop[synthpop$RaceEth4 == "BlackNH", "total"])/2,
                sum(synthpop[synthpop$RaceEth4 == "WhiteNH", "total"])/2,
                sum(synthpop[synthpop$RaceEth4 == "Latinx", "total"])/2,
                sum(synthpop[synthpop$RaceEth4 == "OtherNH", "total"])/2,
                sum(synthpop[synthpop$ego_agecat == "16-20", "total"])/2,
                sum(synthpop[synthpop$ego_agecat == "21-29", "total"])/2
)

target_stats_out$`Main partners`[1:length(multiplier)] * multiplier

target_stats_out$`Casual partners`[1:length(multiplier)] * multiplier

target_stats_out$`One-time partners`[1:length(multiplier)] * multiplier
