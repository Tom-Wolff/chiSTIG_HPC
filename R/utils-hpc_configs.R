# Must be sourced **AFTER** "./R/utils-0_project_settings.R"

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

hpc_configs[["renv_sbatch_opts"]] <- EpiModelHPC:::swf_renv_sbatch_opts()

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

######

hpc_configs <- swf_configs_quest(
  partition = "short", # ASK ABOUT THIS
  r_version = "4.3.0",
  mail_user = "tom.wolff@northwestern.edu"
)
#
# hpc_configs <- EpiModelHPC::swf_configs_hyak(
#   hpc = "mox",
#   partition = "ckpt",
#   r_version = "4.1.2",
#   mail_user = mail_user
# )

# hpc_configs <- EpiModelHPC::swf_configs_hyak(
#   hpc = "klone",
#   partition = "ckpt",
#   r_version = "4.1.1",
#   mail_user = mail_user
# )
#
# hpc_configs <- EpiModelHPC::swf_configs_rsph(
#   partition = c("epimodel", "preemptable")[1],
#   r_version = "4.3.0",
#   mail_user = mail_user
#   )
