
##TODO: document all these functions

#' reduce output file from AnnotFA
#'
#'\code{reduce_output} creates a smaller file with filename prefixed by "less_"
#'containing the essential parts of the output of AnnotFA.
#'
#' @param filename path to RDS file containing the output from AnnotFA,
#' which is usually of the form output_<num_iters>.RDS.

#' @export
reduce_output <- function(filename) {
  stopifnot(file.exists(filename))

  fname <- basename(filename)
  directory <- dirname(filename)

  full_data <- readRDS(filename)

  reduced_data <- list(eta = full_data$vars$eta,
            ELBO = full_data$vFE,
            params = full_data$params,
            geneloadings = full_data$vars$WS$mom1,
            PIPs = full_data$vars$WS$gamma,
            cellscores = full_data$vars$A$mom1,
            Beta = full_data$vars$Beta$mom1)

  saveRDS(reduced_data, paste0(directory, "/less_", fname))
}

#' reduce all output files in a directory.
#'
#' For each AnnotFA output file in the provided directory,
#' \code{reduce_all_output_files} creates a smaller file with filename
#' prefixed by "less_" containing the essential parts of the output of
#' AnnotFA.
#'
#' @param directory path to directory in which files are to be reduced.
#'
#' @export
reduce_all_output_files <- function(directory) {

  stopifnot(dir.exists(directory))

  file_list <- list.files(directory)

  already_reduced_files <- grepl("^less_", file_list) # files with less_ prefix

  files_to_skip <- paste0("less_", file_list) %in% file_list # previously reduced

  output_files <- grepl("output_[0-9]+.RDS$", file_list)

  #includes all output files not already reduced
  files_to_reduce <- output_files & !files_to_skip & !already_reduced_files

  new_file_list <- file_list[files_to_reduce]
  if (length(new_file_list) > 0) {
    for (file in new_file_list) {
      reduce_output(paste0(directory, "/", file))
    }
  }
  return(0)
}

check_less_file_complete <- function(file) {
  # minimal at this point, checks the file exists, and
  #checks the names in the file are correct.  Could/should also check content.
  cat("Checking file", file, "\n")
  if (!file.exists(file)) {
    return(FALSE)
  }else{
    tmp <- readRDS(file)
  }
  if (class(tmp) != "list") {
    return(FALSE)
  }
  expected_names <- c("eta",
                      "ELBO",
                      "params",
                      "geneloadings",
                      "PIPs",
                      "cellscores",
                      "Beta")
  if (all(expected_names %in% names(tmp))) {
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#' remove all except last output file.
#'
#' remove all AnnotFA output files from a
#' directory except the last.
#'
#' @param directory path to directory.

#'@export
rm_output_files_except_last <- function(directory) {
  ## checks the directory for files with name "output_2000.RDS" for example, and
  ## removes main output file when the reduced file "less_output_2000.RDS" exists.
  ## Will not remove the file corresponding to the largest iteration (which
  ## can still be used for a checkpoint).

  stopifnot(dir.exists(directory))

  file_list <- list.files(directory)

  already_reduced_files <- grepl("^less_output_", file_list) #includes all files with less_ prefix

  files_to_skip <- paste0("less_", file_list) %in% file_list # identifies files which have previously been reduced

  output_files <- grepl("output_[0-9]+.RDS$", file_list)

  files_to_consider <- output_files & files_to_skip & !already_reduced_files #includes all possible files from outputs

  files_to_remove <- file_list[files_to_consider]

  itsrds <- lapply(base::strsplit(files_to_remove, "_"),  function(x) x[length(x)])
  iteration <- as.numeric(sapply(strsplit(unlist(itsrds), ".", fixed = TRUE), function(x) x[1]))

  ## do not reduce the maximal iteration file...
  cat("Will not remove", files_to_remove[which.max(iteration)], "\n")
  files_to_remove <- files_to_remove[-which.max(iteration)]
  cat("removing files ", files_to_remove, "\n")
  if (length(files_to_remove) > 0) {
    stopifnot(all(paste0("less_", files_to_remove) %in% file_list))
  }
  # final check that any file being removed has a reduced version present...
  #could check this more thoroughly..
  if (length(files_to_remove) > 0) {
    for (file in files_to_remove) {
      stopifnot(check_less_file_complete(paste0(directory, "/less_", file)))
      cmd <- paste0("rm ", directory, "/", file)
      cat("Running command", cmd, "\n")
      system(cmd)
    }
  }
}

#' reduce all output and remove
#'
#' reduces all output files in directory and
#'  then removes all except the last.  Assumes the output file will be
#'   labelled as "output_<iteration>.RDS" (the default from AnnotFA)
#'   and the reduced files will be
#'   labelled as "less_output_<iteration>.RDS".
#'
#'  @param directory : path to directory.

#' @export
reduce_all_and_remove <- function(directory) {
  stopifnot(dir.exists(directory))
  reduce_all_output_files(directory)
  rm_output_files_except_last(directory)
}
