# Convenience loader for the Sliced Wasserstein Regression codebase.

swr_current_file <- function() {
  frames <- sys.frames()
  ofiles <- vapply(frames, function(frame) {
    file <- frame$ofile
    if (is.null(file)) NA_character_ else file
  }, character(1))
  ofiles <- ofiles[!is.na(ofiles)]

  if (length(ofiles) > 0) {
    return(normalizePath(ofiles[length(ofiles)], mustWork = FALSE))
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[startsWith(args, "--file=")][1])
  if (!is.na(file_arg)) {
    return(normalizePath(file_arg, mustWork = FALSE))
  }

  NA_character_
}

.swr_loader_file <- swr_current_file()
.swr_r_dir <- if (!is.na(.swr_loader_file)) {
  dirname(.swr_loader_file)
} else {
  candidates <- c(file.path(getwd(), "R"), getwd())
  candidates <- candidates[file.exists(file.path(candidates, "utils.R"))]
  if (length(candidates) == 0) {
    stop("Cannot locate the R source directory. Source R/load_swr.R from the project root.")
  }
  normalizePath(candidates[1], mustWork = TRUE)
}

swr_root <- normalizePath(dirname(.swr_r_dir), mustWork = TRUE)

swr_path <- function(...) {
  file.path(swr_root, ...)
}

swr_source <- local({
  r_dir <- .swr_r_dir
  function(file) {
    source(file.path(r_dir, file), chdir = TRUE)
  }
})

swr_source("utils.R")
swr_source("SWW.R")
swr_source("SAW.R")
swr_source("FM.R")
swr_source("bayes_bivar_baseline.R")

rm(.swr_loader_file, .swr_r_dir)
