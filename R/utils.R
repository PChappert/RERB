#' safely switch to lapply from mclapply if parallel is not found.
#'
#' \code{safe_mclapply}
#'
#' @keywords internal
#' 

safe_mclapply <- function(X, FUN, ..., mc.cores = 1) {
  if (requireNamespace("parallel", quietly = TRUE)) {
    return(parallel::mclapply(X, FUN, ..., mc.cores = mc.cores))
  } else {
    message("Package 'parallel' not found. Falling back to lapply().")
    return(lapply(X, FUN, ...))
  }
}

#' safely abort kaleido call if missing
#'
#' \code{safe_kaleido}
#'
#' @keywords internal
#' 

safe_kaleido <- function(scope = "auto") {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    message("Optional: 'plotly' package not available — skipping export.")
    return(invisible(NULL))
  }
  
  result <- tryCatch(
    plotly::kaleido(scope = scope),
    error = function(e) {
      message("plotly::kaleido() failed: ", e$message)
      return(invisible(NULL))
    }
  )
  
  return(result)
}

#### Function to capture safely all outputs from a function and send it to a log file ####
#' safely combine sink() and on.exit() so that even if function (expr) fails, the sink is properly closed.
#'
#' \code{time_and_log}
#'
#' @param expr        function
#' @param verbose     whether to write to console
#' @param time        whether to return elapsed time
#' @param log_file    log file name
#' @param log_title   title for written section in log file
#' @param open_mode   passed to file, "a" or "at" = append [default], "w" or "wt" = write (will erase prior text).
#' 
#' @keywords internal

time_and_log <- function(expr,
                         verbose = TRUE, 
                         time = TRUE,
                         log_file = NULL, 
                         log_title = NULL,
                         open_mode = "a") {
  
  error_occurred <- FALSE
  error_msg <- NULL
  start <- Sys.time()
  result <- NULL
  log_con <- NULL
  raw_txt <- deparse(substitute(expr))
  if (length(raw_txt) > 2 && raw_txt[1] == "{" && tail(raw_txt, 1) == "}") {
    # Strip brackets and combine
    expr_txt <- paste(raw_txt[-c(1, length(raw_txt))], collapse = "; ")
  } else {
    expr_txt <- paste(raw_txt, collapse = "; ")
  }
  
  # Helper to strip ANSI codes
  strip_ansi <- function(text) {
    gsub("\033\\[[0-9;]*m", "", text)
  }
  
  if (!is.null(log_file)) {
    log_file <- normalizePath(log_file, mustWork = FALSE)
    log_con <- file(log_file, open = open_mode)
    sink(log_con)
    sink(log_con, type = "message")
    if(!is.null(log_title)){
      cat(sprintf("[%s] =============================\n", format(start, "%Y-%m-%d %H:%M:%S")))
      cat(sprintf(paste0("[%s] ", log_title, "\n"), format(start, "%Y-%m-%d %H:%M:%S")))
    }
    cat(sprintf("[%s] =============================\n", format(start, "%Y-%m-%d %H:%M:%S")))
    #cat(sprintf("[%s] Started running: %s\n", format(start, "%Y-%m-%d %H:%M:%S"), expr_txt))
  }
  
  tryCatch(
    {
      withCallingHandlers(
        eval(substitute(expr), envir = parent.frame()), 
        warning = function(w) {
          msg <- strip_ansi(w$message)
          cat(sprintf("[%s] WARNING: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          msg <- strip_ansi(m$message)
          cat(sprintf("[%s] MESSAGE: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
          invokeRestart("muffleMessage")
        }
      )
    },
    error = function(e) {
      error_occurred <<- TRUE
      error_msg <<- strip_ansi(e$message)
      cat(sprintf("[%s] ERROR: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), error_msg))
    }
  )
  
  elapsed <- difftime(Sys.time(), start, units = "secs")
  
  if (!is.null(log_file)) {
    if(time){
      cat(sprintf("[%s] Elapsed: %.2f seconds\n", 
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  as.numeric(elapsed)))
    }
    #cat(sprintf("[%s] finished %s\n", 
    #            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    #            if (error_occurred) "with error" else "successfully"))
    cat(sprintf("[%s] =============================\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    sink(type = "message")
    sink()
    close(log_con)
  }
  
  if (verbose && !error_occurred && time) {
    cat(paste0("Elapsed: ", sprintf("%.2f %s", elapsed, units(elapsed)), "\n"))
  }
  
  if (error_occurred && !is.null(log_file)) {
    message(sprintf("Error occurred. Please check the log file at '%s' for details.", 
                    log_file))
  }
  
  invisible(result)  # return result quietly
}
