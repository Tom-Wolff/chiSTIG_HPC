#' @title Cleanup Module
#'
#' @description Module function for cleaning up the Master data list object at the end of a
#' simulation step.
#'
#' @inheritParams aging_msm
#'
#' @return
#' This function returns `dat` after having cleaned up the `temp` list.
#'
#' @export
#'
cleanup_msm_chi <- function(dat, at) {
  print("cleanup")
  dat$temp <- list()
  return(dat)
}
