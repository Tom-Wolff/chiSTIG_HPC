build_epistats_chistig <- function(
                           d,
                           l,

                           # Arguments for RADAR data

                           geog.lvl = NULL,
                           geog.cat = NULL,
                           race = TRUE,
                           age.limits = c(16, 29),
                           age.breaks = c(16, 20, 29),
                           age.sexual.cessation = NULL,
                           init.hiv.prev = NULL,
                           time.unit = 7,
                           browser = FALSE) {

  # Fix global binding check errors
  duration.time <- anal.acts.time <- anal.acts.time.cp <- NULL

  if (browser == TRUE) {
    browser()
  }


  ## Data ##
  # d <- ARTnet.wide
  # l <- ARTnet.long

  out <- list()

  geog_names <- c("city", "county", "state", "region", "division")
  if (!is.null(geog.lvl)) {
    if (!(geog.lvl %in% geog_names)) {
      stop("Selected geographic level must be one of: city, county, state, region or division")
    }
  }

  # Data Processing ---------------------------------------------------------

  # Geography
  if (length(geog.lvl) > 1) {
    stop("Only one geographical level may be chosen at a time.")
  }

  if (!is.null(geog.lvl)) {
    if (geog.lvl == "city") {
      if (sum(geog.cat %in% unique(d$city)) == 0) {
        stop("None of the city names found in the data")
      }
      l <- suppressMessages(left_join(l, d[, c("AMIS_ID", "city2")]))
      l$geogYN <- ifelse(l[, "city2"] %in% geog.cat, 1, 0)
      l$geog <- l$city2
      d$geogYN <- ifelse(d[, "city2"] %in% geog.cat, 1, 0)
      d$geog <- d$city2
    }

    if (geog.lvl == "county") {
      if (sum(geog.cat %in% unique(d$COUNTYFIPS)) == 0) {
        stop("None of the county FIPS codes found in the data")
      }
      l <- suppressMessages(left_join(l, d[, c("AMIS_ID", "COUNTYFIPS")]))
      l$geogYN <- ifelse(l[, "COUNTYFIPS"] %in% geog.cat, 1, 0)
      l$geog <- l$COUNTYFIPS
      d$geogYN <- ifelse(d[, "COUNTYFIPS"] %in% geog.cat, 1, 0)
      d$geog <- d$COUNTYFIPS
    }

    if (geog.lvl == "state") {
      if (sum(geog.cat %in% unique(d$State)) == 0) {
        stop("None of the states found in the data")
      }
      l <- suppressMessages(left_join(l, d[, c("AMIS_ID", "State")]))
      l$geogYN <- ifelse(l[, "State"] %in% geog.cat, 1, 0)
      l$geog <- l$State
      d$geogYN <- ifelse(d[, "State"] %in% geog.cat, 1, 0)
      d$geog <- d$State
    }

    if (geog.lvl == "division") {
      if (sum(geog.cat %in% unique(d$DIVCODE)) == 0) {
        stop("None of the census division codes found in the data")
      }
      l <- suppressMessages(left_join(l, d[, c("AMIS_ID", "DIVCODE")]))
      l$geogYN <- ifelse(l[, "DIVCODE"] %in% geog.cat, 1, 0)
      l$geog <- l$DIVCODE
      d$geogYN <- ifelse(d[, "DIVCODE"] %in% geog.cat, 1, 0)
      d$geog <- d$DIVCODE
    }

    if (geog.lvl == "region") {
      if (sum(geog.cat %in% unique(d$REGCODE)) == 0) {
        stop("None of the census region codes found in the data")
      }
      l <- suppressMessages(left_join(l, d[, c("AMIS_ID", "REGCODE")]))
      l$geogYN <- ifelse(l[, "REGCODE"] %in% geog.cat, 1, 0)
      l$geog <- l$REGCODE
      d$geogYN <- ifelse(d[, "REGCODE"] %in% geog.cat, 1, 0)
      d$geog <- d$REGCODE
    }
  }


  ## Age Processing ##
  if (length(age.limits) != 2 || age.limits[1] > age.limits[2]) {
    stop("age.limits must be a vector of length 2, where age.limits[2] > age.limits[1]")
  }

  # Warning if age range is out of allowed range
  flag.ll <- age.limits[1] >= 15 & age.limits[1] <= 100
  flag.ul <- age.limits[2] >= 15 & age.limits[2] <= 100
  flag.lim <- flag.ll * flag.ul
  if (flag.lim == 0) {
    stop("Age range specified in `age.limits` must be >= 15 and <= 100")
  }

  # Warning if age breaks fall outside age limits
  flag.bks <- prod(age.breaks < age.limits[2] & age.breaks >= age.limits[1])
  if (flag.bks == 0) {
    stop("Age breaks must be between specified age limits")
  }

  # Set default age.sexual.cessation and error if > ARTnet data
  if (is.null(age.sexual.cessation)) {
    age.sexual.cessation <- age.limits[2]
  }
  if (age.sexual.cessation > 66) {
    stop("Maximum allowed age of sexual cessation is 66, corresponding to the upper age eligilibity
         criteria of 65 (inclusive) in ARTnet")
  }

  # Composite age.breaks are now union of age.limits, age.breaks, and age.sexual.cessation
  age.breaks <- unique(sort(c(age.limits[1], age.breaks, age.sexual.cessation, age.limits[2])))

  # p_age_imp initialization for lintr
  p_age_imp <- NULL

  # Subset datasets by lower age limit and age.sexual.cessation
  # Now applies to both index (respondents) and partners for long dataset
  l <- subset(l, age >= age.limits[1] & age < age.sexual.cessation &
                p_age_imp >= age.limits[1] & p_age_imp < age.sexual.cessation)
  d <- subset(d, age >= age.limits[1] & age < age.sexual.cessation)

  if (age.limits[2] > age.sexual.cessation) {
    sex.cess.mod <- TRUE
  } else {
    sex.cess.mod <- FALSE
  }

  # Calculate combine age of index and partners
  l$comb.age <- l$age + l$p_age_imp
  l$diff.age <- abs(l$age - l$p_age_imp)


  ## Race ethnicity ##
  if (race == TRUE) {
    d$race.cat4 <- rep(NA, nrow(d))
    d$race.cat4[d$race.cat == "blackNH"] <- 1
    d$race.cat4[d$race.cat == "hispanic"] <- 2
    d$race.cat4[d$race.cat == "whiteNH"] <- 4
    d$race.cat4[d$race.cat == "other"] <- 3

    l$race.cat4 <- rep(NA, nrow(l))
    l$race.cat4[l$race.cat == "blackNH"] <- 1
    l$race.cat4[l$race.cat == "hispanic"] <- 2
    l$race.cat4[l$race.cat == "whiteNH"] <- 4
    l$race.cat4[l$race.cat == "other"] <- 3

    l$p_race.cat4 <- rep(NA, nrow(l))
    l$p_race.cat4[l$p_race.cat == "blackNH"] <- 1
    l$p_race.cat4[l$p_race.cat == "hispanic"] <- 2
    l$p_race.cat4[l$p_race.cat == "whiteNH"] <- 4
    l$p_race.cat4[l$p_race.cat == "other"] <- 3

    # redistribute NAs in proportion to non-missing partner races
    probs <- prop.table(table(l$race.cat4, l$p_race.cat4), 1)

    imp_black <- which(is.na(l$p_race.cat4) & l$race.cat4 == 1)
    l$p_race.cat4[imp_black] <- sample(1:4, length(imp_black), TRUE, probs[1, ])

    imp_hisp <- which(is.na(l$p_race.cat4) & l$race.cat4 == 2)
    l$p_race.cat4[imp_hisp] <- sample(1:4, length(imp_hisp), TRUE, probs[2, ])

    imp_white <- which(is.na(l$p_race.cat4) & l$race.cat4 == 4)
    l$p_race.cat4[imp_white] <- sample(1:4, length(imp_white), TRUE, probs[4, ])

    imp_other <- which(is.na(l$p_race.cat4) & l$race.cat4 == 3)
    l$p_race.cat4[imp_white] <- sample(1:4, length(imp_white), TRUE, probs[3, ])

    # WE NEED TO CONFER ON WHAT TO DO HERE RE: DYADIC RACE/ETHNICITY COMBOS
    l$race.combo <- rep(NA, nrow(l))
    l$race.combo[l$race.cat4 == 1 & l$p_race.cat4 == 1] <- 1
    l$race.combo[l$race.cat4 == 1 & l$p_race.cat4 %in% 2:4] <- 2
    l$race.combo[l$race.cat4 == 2 & l$p_race.cat4 %in% c(1, 3, 4)] <- 3
    l$race.combo[l$race.cat4 == 2 & l$p_race.cat4 == 2] <- 4
    l$race.combo[l$race.cat4 == 3 & l$p_race.cat4 %in% 1:2] <- 5
    l$race.combo[l$race.cat4 == 3 & l$p_race.cat4 == 3] <- 6
    # New additions
    l$race.combo[l$race.cat4 == 4 & l$p_race.cat4 == 4] <- 7
    l$race.combo[l$race.cat4 == 4 & l$p_race.cat4 %in% 1:3] <- 8

    l <- select(l, -c(race.cat3, p_race.cat3))
  }

  ## HIV diagnosed status of index and partners ##
  l$p_hiv2 <- ifelse(l$p_hiv == 1, 1, 0)

  hiv.combo <- rep(NA, nrow(l))
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 0] <- 1
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 1] <- 2
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 0] <- 3
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 1] <- 3
  hiv.combo[l$hiv2 == 0 & l$p_hiv == 2] <- 4
  hiv.combo[l$hiv2 == 1 & l$p_hiv == 2] <- 5
  l$hiv.concord.pos <- ifelse(hiv.combo == 2, 1, 0)

  ## PrEP ##
  d$prep <- ifelse(d$artnetPREP_CURRENT == 0 | is.na(d$artnetPREP_CURRENT), 0, 1)

  ##### CHANGE AMIS_ID TO WHATEVER ELSE HERE
  dlim <- select(d, c(AMIS_ID, survey.year, prep))
  l <- left_join(l, dlim, by = "AMIS_ID")

  ## Time unit processing ##

  # Set time.unit limits from 1 to 30
  if (time.unit < 1 || time.unit > 30) {
    stop("time.unit must be between 1 and 30")
  }

  # Scale time-based ARTnet data by time.unit
  l$duration.time <- l$duration * 7 / time.unit
  l$anal.acts.time <- l$anal.acts.week * time.unit / 7
  l$anal.acts.time.cp <- l$anal.acts.week.cp * time.unit / 7


  # Act Rates ---------------------------------------------------------------

  # acts/per week/per partnership for main and casual partnerships

  # NOTE: THIS IS WHERE WE WANT TO MODEL FROM RADAR DATA

  # # Pull Data
  # if (race == TRUE) {
  #   if (is.null(geog.lvl)) {
  #     la <- select(l, ptype, duration.time, comb.age,
  #                  race.combo, RAI, IAI, hiv.concord.pos, prep,
  #                  acts = anal.acts.time, cp.acts = anal.acts.time.cp) %>%
  #       filter(ptype %in% 1:2) %>%
  #       filter(RAI == 1 | IAI == 1)
  #     la <- select(la, -c(RAI, IAI))
  #   } else {
  #     la <- select(l, ptype, duration.time, comb.age, geogYN = geogYN,
  #                  race.combo, RAI, IAI, hiv.concord.pos, prep,
  #                  acts = anal.acts.time, cp.acts = anal.acts.time.cp) %>%
  #       filter(ptype %in% 1:2) %>%
  #       filter(RAI == 1 | IAI == 1)
  #     la <- select(la, -c(RAI, IAI))
  #   }
  # }  else {
  #   if (is.null(geog.lvl)) {
  #     la <- select(l, ptype, duration.time, comb.age,
  #                  RAI, IAI, hiv.concord.pos, prep,
  #                  acts = anal.acts.time, cp.acts = anal.acts.time.cp) %>%
  #       filter(ptype %in% 1:2) %>%
  #       filter(RAI == 1 | IAI == 1)
  #     la <- select(la, -c(RAI, IAI))
  #   } else {
  #     la <- select(l, ptype, duration.time, comb.age, geogYN = geogYN,
  #                  RAI, IAI, hiv.concord.pos, prep,
  #                  acts = anal.acts.time, cp.acts = anal.acts.time.cp) %>%
  #       filter(ptype %in% 1:2) %>%
  #       filter(RAI == 1 | IAI == 1)
  #     la <- select(la, -c(RAI, IAI))
  #   }
  # }
  #
  # # Poisson Model
  # if (race == TRUE) {
  #   if (is.null(geog.lvl)) {
  #     acts.mod <- glm(floor(acts * 364 / time.unit) ~ duration.time + I(duration.time^2) +
  #                       as.factor(race.combo) + as.factor(ptype) + duration.time *
  #                       as.factor(ptype) + comb.age + I(comb.age^2) + hiv.concord.pos,
  #                     family = poisson(), data = la)
  #   } else {
  #     acts.mod <- glm(floor(acts * 364 / time.unit) ~ duration.time + I(duration.time^2) +
  #                       as.factor(race.combo) + as.factor(ptype) +
  #                       duration.time * as.factor(ptype) + comb.age + I(comb.age^2) +
  #                       hiv.concord.pos + geogYN,
  #                     family = poisson(), data = la)
  #   }
  # }  else {
  #   if (is.null(geog.lvl)) {
  #     acts.mod <- glm(floor(acts * 364 / time.unit) ~ duration.time + I(duration.time^2) +
  #                       as.factor(ptype) + duration.time * as.factor(ptype) + comb.age +
  #                       I(comb.age^2) + hiv.concord.pos,
  #                     family = poisson(), data = la)
  #   } else {
  #     acts.mod <- glm(floor(acts * 364 / time.unit) ~ duration.time + I(duration.time^2) +
  #                       as.factor(ptype) + duration.time * as.factor(ptype) + comb.age +
  #                       I(comb.age^2) + hiv.concord.pos + geogYN,
  #                     family = poisson(), data = la)
  #   }
  # }
  #
  # # Condom Use // Main Casual -----------------------------------------------
  #
  # la$prob.cond <- la$cp.acts / la$acts
  # la$any.cond <- ifelse(la$prob.cond > 0, 1, 0)
  # la$never.cond <- ifelse(la$prob.cond == 0, 1, 0)
  #
  # if (race == TRUE) {
  #   if (is.null(geog.lvl)) {
  #     cond.mc.mod <- glm(any.cond ~ duration.time + I(duration.time^2) + as.factor(race.combo) +
  #                          as.factor(ptype) + duration.time * as.factor(ptype) + comb.age +
  #                          I(comb.age^2) + hiv.concord.pos + prep,
  #                        family = binomial(), data = la)
  #   } else {
  #     cond.mc.mod <- glm(any.cond ~ duration.time + I(duration.time^2) + as.factor(race.combo) +
  #                          as.factor(ptype) + duration.time * as.factor(ptype) + comb.age +
  #                          I(comb.age^2) + hiv.concord.pos + prep + geogYN,
  #                        family = binomial(), data = la)
  #   }
  # }  else {
  #   if (is.null(geog.lvl)) {
  #     cond.mc.mod <- glm(any.cond ~ duration.time + I(duration.time^2) +
  #                          as.factor(ptype) + duration.time * as.factor(ptype) + comb.age +
  #                          I(comb.age^2) + hiv.concord.pos + prep,
  #                        family = binomial(), data = la)
  #   } else {
  #     cond.mc.mod <- glm(any.cond ~ duration.time + I(duration.time^2) +
  #                          as.factor(ptype) + duration.time * as.factor(ptype) + comb.age +
  #                          I(comb.age ^ 2) + hiv.concord.pos + prep + geogYN,
  #                        family = binomial(), data = la)
  #   }
  # }
  #
  # # Condom Use // Inst ------------------------------------------------------
  # if (race == TRUE) {
  #   if (is.null(geog.lvl)) {
  #     lb <- select(l, ptype, comb.age,
  #                  race.combo, hiv.concord.pos, prep,
  #                  RAI, IAI, RECUAI, INSUAI) %>%
  #       filter(ptype == 3) %>%
  #       filter(RAI == 1 | IAI == 1)
  #   } else {
  #     lb <- select(l, ptype, comb.age, geogYN = geogYN,
  #                  race.combo, hiv.concord.pos, prep,
  #                  RAI, IAI, RECUAI, INSUAI) %>%
  #       filter(ptype == 3) %>%
  #       filter(RAI == 1 | IAI == 1)
  #   }
  # } else {
  #   if (is.null(geog.lvl)) {
  #     lb <- select(l, ptype, comb.age,
  #                  hiv.concord.pos, prep,
  #                  RAI, IAI, RECUAI, INSUAI) %>%
  #       filter(ptype == 3) %>%
  #       filter(RAI == 1 | IAI == 1)
  #   } else {
  #     lb <- select(l, ptype, comb.age, geogYN = geogYN,
  #                  hiv.concord.pos, prep,
  #                  RAI, IAI, RECUAI, INSUAI) %>%
  #       filter(ptype == 3) %>%
  #       filter(RAI == 1 | IAI == 1)
  #   }
  # }
  #
  # lb$prob.cond <- rep(NA, nrow(lb))
  # lb$prob.cond[lb$RAI == 1 & lb$IAI == 0] <- lb$RECUAI[lb$RAI == 1 & lb$IAI == 0] /
  #   lb$RAI[lb$RAI == 1 & lb$IAI == 0]
  # lb$prob.cond[lb$RAI == 0 & lb$IAI == 1] <- lb$INSUAI[lb$RAI == 0 & lb$IAI == 1] /
  #   lb$IAI[lb$RAI == 0 & lb$IAI == 1]
  # lb$prob.cond[lb$RAI == 1 & lb$IAI == 1] <- (lb$RECUAI[lb$RAI == 1 & lb$IAI == 1] +
  #                                               lb$INSUAI[lb$RAI == 1 & lb$IAI == 1]) /
  #   (lb$RAI[lb$RAI == 1 & lb$IAI == 1] + lb$IAI[lb$RAI == 1 & lb$IAI == 1])
  # lb$prob.cond[which(lb$prob.cond == 0.5)] <- 0
  # lb$prob.cond[which(lb$prob.cond %in% c(88, 99, 44))] <- NA
  # lb <- select(lb, -c(RAI, IAI, RECUAI, INSUAI))
  #
  # if (race == TRUE) {
  #   if (is.null(geog.lvl)) {
  #     cond.oo.mod <- glm(prob.cond ~ as.factor(race.combo) +
  #                          comb.age + I(comb.age^2) +
  #                          hiv.concord.pos + prep,
  #                        family = binomial(), data = lb)
  #   } else {
  #     cond.oo.mod <- glm(prob.cond ~ as.factor(race.combo) +
  #                          comb.age + I(comb.age^2) +
  #                          hiv.concord.pos + prep + geogYN,
  #                        family = binomial(), data = lb)
  #   }
  # } else {
  #   if (is.null(geog.lvl)) {
  #     cond.oo.mod <- glm(prob.cond ~ comb.age + I(comb.age^2) +
  #                          hiv.concord.pos + prep,
  #                        family = binomial(), data = lb)
  #   } else {
  #     cond.oo.mod <- glm(prob.cond ~ comb.age + I(comb.age^2) +
  #                          hiv.concord.pos + prep + geogYN,
  #                        family = binomial(), data = lb)
  #   }
  # }
  #
  # Init HIV Status ---------------------------------------------------------
  if (is.null(init.hiv.prev)) {
    if (race == TRUE) {
      if (is.null(geog.lvl)) {
        d1 <- select(d, race.cat4, age, hiv2)

        hiv.mod <- glm(hiv2 ~ age + as.factor(race.cat4),
                       data = d1, family = binomial())
      } else {
        d1 <- select(d, race.cat3, geogYN, age, hiv2)
        hiv.mod <- glm(hiv2 ~ age + geogYN + as.factor(race.cat3) + geogYN * as.factor(race.cat3),
                       data = d1, family = binomial())
      }
    } else {
      if (is.null(geog.lvl)) {
        d1 <- select(d, age, hiv2)

        hiv.mod <- glm(hiv2 ~ age,
                       data = d1, family = binomial())
      } else {
        d1 <- select(d, geogYN, age, hiv2)

        hiv.mod <- glm(hiv2 ~ age + geogYN,
                       data = d1, family = binomial())
      }
    }
    # Output
    out$hiv.mod <- hiv.mod
  } else {
    if (length(init.hiv.prev) != 4) {
      stop("Input parameter init.prev.hiv must be a vector of size four")
    }
    if (prod(init.hiv.prev < 1) == 0  || prod(init.hiv.prev > 0) == 0) {
      stop("All elements of init.hiv.prev must be between 0 and 1 non-inclusive")
    }
  }

  # Save Out File -----------------------------------------------------------

  if (!is.null(geog.lvl)) {
    out$geogYN.l <- l$geogYN
    out$geogYN.d <- d$geogYN
    out$geog.cat  <- geog.cat
  }

  out$geog.lvl <- geog.lvl
  out$race <- race
  out$acts.mod <- acts.mod
  out$cond.mc.mod <- cond.mc.mod
  out$cond.oo.mod <- cond.oo.mod
  out$geog.l <- as.character(l$geog)
  out$geog.d <- as.character(d$geog)
  out$age.limits <- age.limits
  out$age.breaks <- age.breaks
  out$age.grps <- length(age.breaks) - 1
  out$age.sexual.cessation <- age.sexual.cessation
  out$sex.cess.mod <- sex.cess.mod
  out$init.hiv.prev <- init.hiv.prev
  out$time.unit <- time.unit
  return(out)
}
