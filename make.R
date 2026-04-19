# Render qmd files
# Rscript --vanilla make.R
# Rscript --vanilla make.R --day 1
# Rscript --vanilla make.R --session 01
# Rscript --vanilla make.R --session 01 --format beamer

if (!requireNamespace("rconfig", quietly = TRUE)) {
  install.packages("rconfig")
}
CONFIG <- rconfig::rconfig()
str(CONFIG)

FORMAT <- rconfig::value(CONFIG$format, "all")

if (!is.null(CONFIG$day) && !is.null(CONFIG$session)) {
  stop("Cannot specify both day and session")
}

if (!is.null(CONFIG$day)) {
  DAY <- rconfig::value(CONFIG$day, "0")
  dir <- switch(
    as.character(DAY[1]),
    "0" = ".",
    "1" = "day-01",
    "2" = "day-02",
    "3" = "day-03",
    "4" = "day-04",
    stop("Invalid day argument")
  )
  fl <- list.files(dir, recursive = TRUE)
  fl <- fl[grep("\\.qmd$", fl)]
}

if (!is.null(CONFIG$session)) {
  SESSION <- rconfig::value(CONFIG$session, "0")
  dir <- "."
  fl <- list.files(dir, recursive = TRUE)
  fl <- fl[grep("\\.qmd$", fl)]
  fl <- fl[grep(SESSION, fl)]
}

OK <- logical(length(fl))
for (i in fl) {
  r <- try(quarto::quarto_render(file.path(dir, i), output_format = FORMAT))
  OK[i == fl] <- !inherits(r, "try-error")
}

if (any(!OK)) {
  warning("Some files failed to render:")
  cat(fl[!OK], sep = "\n")
  cat("\n")
} else {
  message("All files rendered successfully.")
}
