swf_configs_quest <- function(partition = "short",
                              account = "p32153",
                              r_version = "4.3.0",
                              conda_proj = NULL,
                              mail_user = "tom.wolff@northwestern.edu") {

  hpc_configs <- list()
  hpc_configs[["default_sbatch_opts"]] <-  list(
    "partition" = partition,
    "account" = account,
    "mail-type" = "FAIL"
  )

  if (!is.null(mail_user)) {
    hpc_configs[["default_sbatch_opts"]][["mail-user"]] <- mail_user
  }

  hpc_configs[["renv_sbatch_opts"]] <- swf_renv_sbatch_opts()

  hpc_configs[["r_loader"]] <- c(
    "module purge all",
    paste0("module load R/", r_version),
    "module load git"
  )

  if (!is.null(conda_proj)) {
    hpc_configs[["r_loader"]] <- c(hpc_configs[["r_loader"]], c(
      "module load python-miniconda3",
      paste0("source activate /projects/", conda_proj ,"/pythonenvs/env1")
    ))
  }

  return(hpc_configs)
}
