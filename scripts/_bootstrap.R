# Bootstrap paths for simulation scripts.

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

swr_find_root <- function() {
  script_file <- swr_current_file()
  starts <- unique(c(
    if (!is.na(script_file)) dirname(script_file),
    getwd(),
    dirname(getwd())
  ))

  for (start in starts) {
    candidate <- normalizePath(start, mustWork = FALSE)
    while (dirname(candidate) != candidate) {
      if (file.exists(file.path(candidate, "R", "load_swr.R")) &&
          file.exists(file.path(candidate, "scripts"))) {
        return(normalizePath(candidate, mustWork = TRUE))
      }
      candidate <- dirname(candidate)
    }
  }

  stop("Cannot locate the project root. Run from the repo or keep scripts/_bootstrap.R in place.")
}

PROJECT_ROOT <- swr_find_root()
source(file.path(PROJECT_ROOT, "R", "load_swr.R"))

DATA_DIR <- swr_path("data")
FIGURE_DIR <- swr_path("figures")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

rm(swr_current_file, swr_find_root)
