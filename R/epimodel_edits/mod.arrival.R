
#' @title Arrivals Module
#'
#' @description Module function for arrivals into the sexually active population.
#'
#' @inheritParams aging_msm
#'
#' @details
#' New population members are added based on expected numbers of entries,
#' stochastically determined with draws from Poisson distributions. For each new
#' entry, a set of attributes is added for that node, and the nodes are added onto
#' the network objects.
#'
#' @return
#' This function updates the `attr` list with new attributes for each new population member.
#'
#' @export
#'
arrival_msm_chi <- function(dat, at) {

  print("arrival")
  assign(x = "arrival_dat", value = dat, .GlobalEnv)
  assign(x = "arrival_at", value = at, .GlobalEnv)


  ## Input
  # Attributes

  # Parameters
  a.rate   <- get_param(dat, "a.rate")

  ## Process
  num <- get_epi(dat, "num", at = 1)

  nNew <- rpois(1, a.rate * num)
print("Number of new nodes added")
print(nNew)

  ## Update Attr
  if (nNew > 0) {
    dat <- setNewAttr_msm(dat, at, nNew)
  }

  # Update Networks
  if (nNew > 0) {
    for (i in 1:3) {
      dat[["el"]][[i]] <- add_vertices(dat[["el"]][[i]], nNew)
    }
  }

  ## Output
  dat <- set_epi(dat, "nNew", at, nNew)

  return(dat)
}


setNewAttr_msm <- function(dat, at, nNew) {

  dat <- append_core_attr(dat, at, nNew)

  ## Inputs
  # Attributes
  active <- get_attr(dat, "active")

  # Parameters
  netstats             <- get_param(dat, "netstats")
  tt.partial.supp.prob <- get_param(dat, "tt.partial.supp.prob")
  tt.full.supp.prob    <- get_param(dat, "tt.full.supp.prob")
  tt.durable.supp.prob <- get_param(dat, "tt.durable.supp.prob")
  hiv.circ.prob        <- get_param(dat, "hiv.circ.prob")
  arrival.age          <- get_param(dat, "arrival.age")


  # Set all attributes NA by default
  attrToInit <- names(which(vapply(dat[["attr"]], length, 0) < length(active)))
  for (attr_name in attrToInit) {
    attr_val <- get_attr(dat, attr_name)
    attr_val <- c(attr_val, rep(NA, nNew))
    dat <- set_attr(dat, attr_name, attr_val, override.length.check = TRUE)
  }
  newIds <- which(get_attr(dat, "entrTime") == at)

  race.dist <- prop.table(table(netstats[["attr"]][["race"]]))

  print("Sampling Race")
  race <- sample(sort(unique(get_attr(dat, "race"))), nNew, TRUE, race.dist)
print(sort(unique(get_attr(dat, "race"))))
print("Races sampled")
print(race)
  dat <- set_attr(dat, "race", race, posit_ids = newIds)

  dat <- set_attr(dat, "age", rep(arrival.age, nNew), posit_ids = newIds)
  age.breaks <- netstats[["demog"]][["age.breaks"]]
  attr_age.grp <- cut(
    get_attr(dat, "age", posit_ids = newIds),
    age.breaks,
    labels = FALSE,
    right = FALSE
  )
  dat <- set_attr(dat, "age.grp", attr_age.grp, posit_ids = newIds)

  # Disease status and related
  dat <- set_attr(dat, "status", rep(0, nNew), posit_ids = newIds)
  dat <- set_attr(dat, "diag.status", rep(0, nNew), posit_ids = newIds)
  dat <- set_attr(dat, "rGC", 0, posit_ids = newIds)
  dat <- set_attr(dat, "uGC", 0, posit_ids = newIds)
  dat <- set_attr(dat, "rCT", 0, posit_ids = newIds)
  dat <- set_attr(dat, "uCT", 0, posit_ids = newIds)
  dat <- set_attr(dat, "rGC.timesInf", 0, posit_ids = newIds)
  dat <- set_attr(dat, "uGC.timesInf", 0, posit_ids = newIds)
  dat <- set_attr(dat, "rCT.timesInf", 0, posit_ids = newIds)
  dat <- set_attr(dat, "uCT.timesInf", 0, posit_ids = newIds)

  rates <- get_param(dat, "hiv.test.late.prob")[race]
  late.tester <- runif(length(rates)) < rates
  dat <- set_attr(dat, "late.tester", late.tester, posit_ids = newIds)

  race.new <- get_attr(dat, "race", posit_ids = newIds)
print(race.new)
  races <- sort(unique(race.new))
print("Races before tt.traj")
print(races)
  tt.traj <- rep(NA, nNew)
  for (i in races) {
    ids.race <- which(get_attr(dat, "race")[newIds] == i)
    print("tt.traj sample")
    print(ids.race)
    print(c(tt.partial.supp.prob[i], tt.full.supp.prob[i], tt.durable.supp.prob[i]))

    tt.traj[ids.race] <- sample(
      1:3, length(ids.race), TRUE,
      c(tt.partial.supp.prob[i], tt.full.supp.prob[i], tt.durable.supp.prob[i])
    )

  }
  dat <- set_attr(dat, "tt.traj", tt.traj, posit_ids = newIds)

  # Circumcision
  circ <- rep(NA, nNew)
  for (i in races) {
    ids.race <- which(race.new == i)
    circ[ids.race] <- runif(length(ids.race)) < hiv.circ.prob[i]
  }
  dat <- set_attr(dat, "circ", circ, posit_ids = newIds)

  # Role
  ns <- netstats[["attr"]]
  role.class <- rep(NA, nNew)
  for (i in races) {
    ids.race <- which(race.new == i)
    rc.probs <- prop.table(table(ns[["role.class"]][ns[["race"]] == i]))
print("role.class sample")
    role.class[ids.race] <- sample(0:2, length(ids.race), TRUE, rc.probs)
    # If we want to assign all to versatile just to get things going
    # role.class[ids.race] <- rep(2, length(ids.race))
  }
  dat <- set_attr(dat, "role.class", role.class, posit_ids = newIds)

  ins.quot <- rep(NA, nNew)
  role.class.new <- get_attr(dat, "role.class", posit_ids = newIds)
  ins.quot[role.class.new == 0]  <- 1
  ins.quot[role.class.new == 1]  <- 0
  ins.quot[role.class.new == 2]  <- runif(sum(role.class.new == 2))

  dat <- set_attr(dat, "ins.quot", ins.quot, posit_ids = newIds)

  # Degree
  dat <- set_attr(dat, "deg.main", 0, posit_ids = newIds)
  dat <- set_attr(dat, "deg.casl", 0, posit_ids = newIds)
  dat <- set_attr(dat, "deg.tot", 0, posit_ids = newIds)

  # One-off risk group
  print("one-off risk group sample")
  dat <- set_attr(dat, "risk.grp", sample(1:5, nNew, TRUE), posit_ids = newIds)

  # PrEP
  dat <- set_attr(dat, "prepStat", 0, posit_ids = newIds)

  # Venue attendance (PLACEHOLDER)
  dat <- set_attr(dat, "venues.all", paste(newIds, "s", sep = ""), posit_ids = newIds)

  # App use (TBD)


  return(dat)


}
