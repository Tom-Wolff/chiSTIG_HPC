
#####################################
# Load libraries 
#####################################
library(EpiModel); library(tidyverse)

#####################################
# Set seed
#####################################
set.seed(4)

#####################################
# Load and specify data
#####################################
egos <- read.csv("egos_v3_3.csv")

venues <- read.csv("egos_venues.csv")
egos <- left_join(egos,venues,by=c("numeric_id"="egoid"))

apps <- read.csv("egos_apps.csv")
egos <- left_join(egos,apps,by=c("numeric_id"="egoid"))

rel <- read.csv("egos_rel.csv")
egos <- left_join(egos,rel,by=c("numeric_id"="egoid"))


#####################################
# Initialize network and set vertex attributes 
#####################################
nw <- network.initialize(n = nrow(egos), 
			loops = FALSE, 
			directed = FALSE)

# Set vertex attributes
nw <- set.vertex.attribute(nw, "race_ethnicity", egos$race_ethnicity)
nw <- set.vertex.attribute(nw, "age", egos$agerange)

nw <- set.vertex.attribute(nw, "venues_all", egos$venues_all)
nw <- set.vertex.attribute(nw, "venues_dating", egos$venues_dating)
nw <- set.vertex.attribute(nw, "venues_nondating", egos$venues_nondating)

nw <- set.vertex.attribute(nw, "apps_all", egos$apps_all)
nw <- set.vertex.attribute(nw, "apps_dating", egos$apps_dating)
nw <- set.vertex.attribute(nw, "apps_nondating", egos$apps_nondating)

nw <- set.vertex.attribute(nw, "init_ser_cat", as.character(egos$init_ser_cat))


#####################################
# Specify target stats  
#####################################
target.stats <- c(
	1973.6099165043618,
	729.9160170974485,
	969.8622937528471,
	1105.1213023742796,
	614.6964285714287,
	933.7690082948353,
	362.57994164913316,
	1190.3021032077875,
	1429.1260872140865,
	178.04335253707455,
	285.3968620867489,
	145.40249058632142
)


#####################################
# Specify model formation  
#####################################
est <- netest(
		nw = nw,
		formation =
			~ edges + 
			concurrent + 
			nodefactor("race_ethnicity", levels=-LARGEST) + 
			nodematch("race_ethnicity") + 
			nodematch("race_ethnicity", diff=TRUE, levels=c("blackNH")) + 
			nodefactor("age", levels=-LARGEST) + 
			nodematch("age") + 
			nodefactor("init_ser_cat", levels=-LARGEST) + 
			fuzzynodematch("venues_all", binary=TRUE) + 
			fuzzynodematch("apps_nondating", binary=TRUE)
		,
		target.stats = target.stats,
		coef.diss = dissolution_coefs(~offset(edges), duration = 76),
		set.control.ergm = 
			control.ergm(
				parallel = 4,
				MCMC.interval = 10000,
				MCMLE.effectiveSize=NULL,
				MCMC.burnin = 1000,
				MCMC.samplesize = 20000,
				SAN.maxit = 20,
				SAN.nsteps.times = 10
			)                                          
	)
summary(est)




dx <- netdx(
		est, 
		dynamic = TRUE,
		nsims = 100, 
		ncores = 4, 
		nsteps = 1,    
		set.control.ergm = control.simulate.formula(MCMC.burnin = 5e5),
		set.control.tergm = control.simulate.formula.tergm(MCMC.burnin.min = 3e5),
		nwstats.formula = #used to monitor in each of the simulated networks against the est object 
			~ edges + 
			concurrent +
			nodefactor("race_ethnicity", levels=-LARGEST)  +  
			nodematch("race_ethnicity", diff=FALSE) + 
			nodematch("race_ethnicity", diff=TRUE) + 
			nodefactor("age", levels=-LARGEST) + 
			nodematch("age", diff=FALSE) +
			nodematch("age", diff=TRUE) +
 			nodefactor("init_ser_cat", levels=-LARGEST) + 
			fuzzynodematch("venues_all", binary=TRUE) +  						
			fuzzynodematch("venues_dating", binary=TRUE) + 
			fuzzynodematch("venues_nondating", binary=TRUE) + 
			fuzzynodematch("apps_all", binary=TRUE) +  
			fuzzynodematch("apps_dating", binary=TRUE) +
			fuzzynodematch("apps_nondating", binary=TRUE)   
		,
		keep.tedgelist = TRUE,
		keep.tnetwork = TRUE
	)
dx
