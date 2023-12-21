netstats <- readRDS("~/Desktop/ChiSTIG_model/epimodel/data/intermediate/estimates/netstats-novenues-local.rds")

# Known prevalence numbers (from CDPH 2019)
white_known <- 209
black_known <- 1173
hisp_known <- 382
other_known <- 126

black_exp <- black_known/.75
black_exp

hisp_exp <- hisp_known/.75
hisp_exp

other_exp <- other_known/.75
other_exp

white_exp <- white_known/.8
white_exp

total_exp <- white_exp + black_exp + hisp_exp + other_exp

synthpop_size <- length(netstats$attr$race)
synthpop_size

prev_exp <- total_exp/synthpop_size
prev_exp

black_prev_exp <- black_exp/sum(netstats$attr$race == 1)
black_prev_exp

hisp_prev_exp <- hisp_exp/sum(netstats$attr$race == 2)
hisp_prev_exp

other_prev_exp <- other_exp/sum(netstats$attr$race == 3)
other_prev_exp

white_prev_exp <- white_exp/sum(netstats$attr$race == 4)
white_prev_exp
