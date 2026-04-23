dryrun <- FALSE
ask <- TRUE

cat(
  "Hello!\n\n",
  "This script will install required R packages for the\n",
  "'Analysing Ecological Data with Detection Error' workshop.\n\n"
)

if (ask) {
  cat("Do you want to proceed?")
  if (utils::menu(c("Yes", "No")) != 1L) {
    cat("\nNo packages installed.\n\n")
  } else {
    cat("\nHang tight ...\n\n* Checking R version\n\n")
    if (getRversion() < "4.0") {
      stop(paste0(
        "R version ",
        getRversion(),
        " detected.",
        " Please upgrade to R >= 4.0"
      ))
    } else {
      cat(
        "  - R version",
        as.character(getRversion()),
        "detected ... OK\n"
      )
    }

    cat("\n* Installing R packages from CRAN\n\n")
    pkgs <- c(
      c(
        "bSims",
        "corrplot",
        "detect",
        "dismo",
        "Distance",
        "dplyr",
        "gbm",
        "GGally",
        "ggplot2",
        "intrval",
        "knitr",
        "lme4",
        "MASS",
        "mefa4",
        "mgcv",
        "MuMIn",
        "opticut",
        "partykit",
        "plotrix",
        "pROC",
        "pscl",
        "quarto",
        "rconfig",
        "reactable",
        "remotes",
        "ResourceSelection",
        "rgl",
        "rmarkdown",
        "sf",
        "shiny",
        "skimr",
        "survival",
        "terra",
        "unmarked",
        "utils",
        "visreg"
      )
    )
    to_inst <- setdiff(pkgs, rownames(installed.packages()))
    if (length(to_inst) > 0L) {
      if (dryrun) {
        for (i in to_inst) {
          cat("  - Installing R package ", i, " ... OK\n")
        }
      } else {
        try(install.packages(
          to_inst,
          repos = c(
            "CRAN" = getOption("repos"),
            "INLA" = "https://inla.r-inla-download.org/R/stable"
          ),
          dependencies = TRUE
        ))
      }
    }

    cat("\n* Installing R packages from GitHub\n\n")
    gh_pkgs <- c(
      # bSims = "psolymos/bSims", # use dev
      # detect = "psolymos/detect", # use dev
      lhreg = "psolymos/lhreg"
    ) # not on CRAN
    for (i in gh_pkgs) {
      if (dryrun) {
        cat("  - Installing R package ", i, " ... OK\n")
      } else {
        try(remotes::install_github(
          i,
          upgrade = "never",
          force = FALSE
        ))
      }
    }

    still_missing <- setdiff(
      c(pkgs, names(gh_pkgs)),
      rownames(installed.packages())
    )
    if (length(still_missing)) {
      cat(
        "\nThe following packages could not be installed:\n",
        paste("\t-", still_missing, collapse = "\n"),
        "\n"
      )
    } else {
      cat("\nYou are all set! Have fun at the workshop.\n")
    }
  }
}
