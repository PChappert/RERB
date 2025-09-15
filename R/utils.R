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

#### Function to check column types prior to using dplyr::bind_rows ####
#' safely bind data table when column of identical names have different types (e.g. only NA or missing values in one of the two (d_call for Light chains for exemple).
#'
#' \code{safe_bind_rows}
#'
#' @param list a list of data frames
#' 
#' @keywords internal

safe_bind_rows <- function(list,
                           airr_col_types = c(
                             "sequence_id" = "character",
                             "sequence" = "character",
                             "rev_comp" = "logical",
                             "productive" = "logical",
                             "v_call" = "character",
                             "d_call" = "character",
                             "j_call" = "character",
                             "c_call" = "character",
                             "sequence_alignment" = "character",
                             "germline_alignment" = "character",
                             "junction" = "character",
                             "junction_aa" = "character",
                             "v_cigar" = "character",
                             "d_cigar" = "character",
                             "j_cigar" = "character",
                             "stop_codon" = "logical",
                             "vj_in_frame" = "logical",
                             "locus" = "character",
                             "junction_length" = "double",
                             "np1_length" = "double",
                             "np2_length" = "double",
                             "v_sequence_start" = "double",
                             "v_sequence_end" = "double",
                             "v_germline_start" = "double",
                             "v_germline_end" = "double",
                             "d_sequence_start" = "double",
                             "d_sequence_end" = "double",
                             "d_germline_start" = "double",
                             "d_germline_end" = "double",
                             "j_sequence_start" = "double",
                             "j_sequence_end" = "double",
                             "j_germline_start" = "double",
                             "j_germline_end" = "double",
                             "v_score" = "double",
                             "v_identity" = "double",
                             "v_support" = "double",
                             "d_score" = "double",
                             "d_identity" = "double",
                             "d_support" = "double",
                             "j_score" = "double",
                             "j_identity" = "double",
                             "j_support" = "double",
                             "fwr1" = "character",
                             "fwr2" = "character",
                             "fwr3" = "character",
                             "fwr4" = "character",
                             "cdr1" = "character",
                             "cdr2" = "character",
                             "cdr3" = "character",
                             "primers" = "character",
                             "sequence_length" = "double",
                             "pct_under_30QC_in_trimmed"= "double",
                             "QC_passed" = "logical",
                             "v_germline_length" = "double",
                             "d_germline_length" = "double",
                             "j_germline_length" = "double",
                             "mu_freq_cdr_r" = "double",
                             "mu_freq_cdr_s" = "double",
                             "mu_freq_fwr_r" = "double",
                             "mu_freq_fwr_s" = "double",
                             "mu_freq" = "double",
                             "mu_count_cdr_r" = "double",
                             "mu_count_cdr_s" = "double",
                             "mu_count_fwr_r" = "double",
                             "mu_count_fwr_s" = "double",
                             "mu_count" = "double",
                             "c_call_igblast" = "character",
                             "c_call_alignmant_score" = "double",
                             "c_call_pct_match" = "double",
                             "c_call_alignmant_length" = "double"
                           )){
  
  #first make sure all AIRR columns are of the right type:
  coerce_cols <- function(df, type_map) {
    common_cols <- intersect(names(type_map), names(df))
    
    df %>%
      mutate(across(all_of(common_cols), ~ switch(
        type_map[cur_column()],
        character = as.character(.),
        logical   = as.logical(.),
        double    = as.double(.),
        integer   = as.integer(.),
        factor    = as.factor(.),
        .         # default: leave as is
      )))
  }
  list <- purrr::map(list, coerce_cols, type_map = airr_col_types)
  
  #then find remaining issues:
  col_types <- list %>% 
    purrr::map(~ sapply(.x, typeof))
  
  types_long <- col_types %>%
    dplyr::bind_rows(.id = "df_id") %>%
    tidyr::pivot_longer(-df_id, names_to = "colname", values_to = "type")
  
  cols_with_conflict <- types_long %>%
    dplyr::group_by(colname) %>%
    dplyr::summarise(
      n_types = n_distinct(type), 
      types_present = toString(unique(type)),
      type_to_coerce = {
        if ("character" %in% unique(type)){
          "character"
        } else {
          other_types <- setdiff(unique(type), c("NULL"))
          if(other_types[1] %in% c("double", "integer", "logical")){
            other_types[1]
          } else {"character"}
        }
      },
      .groups = "drop") %>%
    dplyr::filter(n_types > 1)
  
  coerce_lookup <- cols_with_conflict %>%
    dplyr::filter(!is.na(type_to_coerce)) %>%
    dplyr::select(colname, type_to_coerce) %>%
    tibble::deframe() 
  
  # all columns should be one of character, double, logical and eventually integer.
  coercers <- list(
    character = as.character,
    integer = as.integer,
    double = as.double,
    logical = as.logical
  )
  
  list <- list %>%
    purrr::map(
      function(df) {
        for (col in cols_with_conflict$colname){
          coercer <- coercers[[coerce_lookup[[col]]]]
          if(col %in% colnames(df)){
            df[[col]] <- coercer(df[[col]])
          }
        }
        return(df)
      })
  
  dataframe <- dplyr::bind_rows(list)
  
  #recheck all AIRR columns to ensure proper formatting (usefull if some AIRR columns are absent from some data table in the original list):
  dataframe <- coerce_cols(dataframe, type_map = airr_col_types)
  
  return(dataframe)
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
  
  elapsed <- difftime(Sys.time(), start, units = "mins")
  
  if (!is.null(log_file)) {
    if(time){
      cat(sprintf("[%s] Elapsed: %.2f mins\n", 
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
