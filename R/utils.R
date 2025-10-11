
#' Pipe operator
#'
#' See \code{dplyr::\link[dplyr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom dplyr %>%
#' @usage lhs %>% rhs
NULL

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


#### Function to quickly print an annotated table to the plot window usingb grid ####
#' adds formatted title and subtitle to keep track of what is being printed
#'
#' \code{plot_recap_gtable}
#'
#' @param df a dataframe containing at least the two columns defined in factors
#' @param factors names of columns to use, by default columns 1 and 2
#' @param title title to use
#' @param title_gp gpar for title
#' @param subtitle whether to add a subtitle resuming the comparison made ('factor1 (rows) versus factor2 (cols'))
#' @param subtitle_gp gpar for subtitle
#' @param padding padding to use
#' @param useNA one of 'ifany' [default], 'no' or "always', passed to table() 
#' @param print_table whether to print final table grob object
#' @param return_table whether to return the final table grob object
#' 
#' @export
#' 
plot_recap_gtable <- function(df, 
                              factors = colnames(df)[1:2], 
                              title = "", 
                              title_gp = grid::gpar(fontsize = 14, fontface = "bold"), 
                              subtitle = TRUE, 
                              subtitle_gp = grid::gpar(fontsize = 9, fontface = "italic"), 
                              padding = 5, 
                              useNA = c("ifany", "no", "always"),
                              print_table = TRUE,
                              return_table = FALSE) {
  
  useNA <- match.arg(useNA)
  padding <- unit(padding,"mm")
  
  table <- gridExtra::tableGrob(as.matrix(table(df[[factors[1]]], df[[factors[2]]], useNA=useNA)))
  
  title.list <- list()
  if(!is.null(title)){
    title <- grid::textGrob(title, gp = title_gp, just = "center", hjust = 0.5)
    table <- gtable::gtable_add_rows(table, heights = grid::grobHeight(title) + padding, pos = 0)
    title.list <- c(title.list, list(title))
  }
  if(subtitle){
    subtitle <- grid::textGrob(paste0(factors[1], " (rows) by ", factors[2], " (cols)"), x=0.5, hjust=0.5, gp= subtitle_gp)
    table <- gtable::gtable_add_rows(table, heights = grid::grobHeight(subtitle) + padding, pos = 0)
    title.list <- c(title.list, list(subtitle))
  }
  if(length(title.list)==2){
    table <- gtable::gtable_add_grob(table, title.list, t=c(1,2), b=c(1,2), l=c(1,1), r=ncol(table), clip = "off")
  } 
  if(length(title.list)==1){
    table <- gtable::gtable_add_grob(table, title.list[[1]], t=1, b=1, l=1, r=ncol(table), clip = "off")
  } 
  
  if (print_table){
    grid::grid.newpage()
    grid::grid.draw(table)
    invisible(table)
  }

  if (return_table){
    return(table)
  }
}

#### Function to quickly export an annotated table() call to excel ####
#' adds a new sheet to an existing opened excel workbook
#'
#' \code{excel_recap_table}
#'
#' @param wb which workbook to write into
#' @param df data frame containing at least the two columns defined in factors to use for table generation
#' @param factors which columns to use for table (factors[1] = rows; factors[2] = cols), by default colnames 1 and 2 of df
#' @param add_header whether to add recap header with info on data and factors used
#' @param writeData whether to write to a workbook sheet
#' @param sheet name of the sheet to write too. Will overwrite any data in it.
#' @param useNA whether to include NA in count table and or in freq table, if useNA="counts (default), only non-NA values will be used to calculate frequencies but NA counts will be included in the counts table.
#' @param add_freq whether to add a freq table (row wise) next to the count table
#' @param colNames whether to add rownames in the final excel sheet
#' @param rowNmaes whether to add colnames in the final excel sheet
#' 
#' @return create a new sheet in the provided Workbook, with raw table outputs, absolute numbers (row-wise and column-wise), and frequencies (row-wise))
#' 
#' @export
#' 
excel_recap_table <- function(wb, 
                              df, 
                              factors = colnames(df)[1:2], 
                              writeData = TRUE, 
                              sheet = NULL, 
                              add_header = TRUE, 
                              useNA=c("counts", "no", "all"), 
                              add_freq = TRUE, 
                              colNames = TRUE, rowNames = TRUE){
  
  cmd <- deparse(substitute(df))
  useNA <- match.arg(useNA)
  if(useNA %in% c("counts", "all")){
    useNA_counts <- "ifany"
  } else { useNA_counts <- "no" }
  
  tb <- as.data.frame.matrix(table(df[[factors[1]]], df[[factors[2]]], useNA = useNA_counts))
  # Clean up any NA row/column names
  rownames(tb)[rownames(tb) == "NA."] <- "missing"
  colnames(tb)[is.na(colnames(tb))] <- "missing"
  tb$total_events <- rowSums(tb)
  if("missing" %in% colnames(tb)){
    tb <- tb %>% dplyr::mutate(total_events_woNA = total_events - missing)
  } else {
    tb <- tb %>% dplyr::mutate(total_events_woNA = total_events)
  }
  tb <- rbind(tb, total_events = colSums(tb))
  tb <- rbind(tb, total_events_woNA = colSums(tb[!rownames(tb) %in% c("total_events", "missing"),]))
  
  n_rows_counts <- nrow(tb)
  n_cols_counts <- ncol(tb)
  
  if(add_freq){
    if(useNA == "all"){
      tb_counts <- tb %>% dplyr::select(-any_of(c("total_events", "total_events_woNA")))
    } else {
      tb_counts <- tb %>% dplyr::select(-any_of(c("total_events", "total_events_woNA", "missing")))
    }
    
    tb_freq <- tb_counts
    tb_freq[,] <- NA  # initialize
    for (r in 1:(nrow(tb_counts))) {
      if(useNA == "all"){
        row_total <- tb$total_events[r]
      } else {
        row_total <- tb$total_events_woNA[r]
      }
      if (row_total > 0) {
        tb_freq[r, ] <- round(tb_counts[r, ] / row_total, 4)*100
      }
    }
    colnames(tb_freq) <- paste0(colnames(tb_freq), "_freq")
    tb <- cbind(tb, tb_freq)
  }
  
  # Write to excel workbook
  if(writeData){
    if(!is.null(sheet)){
      openxlsx::addWorksheet(wb, sheet)
      totalStyle <- openxlsx::createStyle(fgFill = "lightblue", textDecoration = "bold")
      #freqStyle <- createStyle(fgFill = "#FFD580")
      total_events_rows <- which(grepl("total_events", rownames(df)))
      total_events_cols <- which(grepl("total_events", colnames(df)))
      
      if(add_header){
        header_text <- data.frame(
          note = c(paste0(factors[1]," (rows) vs ", factors[2]," (columns)"), cmd)
        )
        openxlsx::writeData(wb, sheet, header_text, startRow = 1, colNames = FALSE)
        openxlsx::writeData(wb, sheet, x = tb, startRow = 4, colNames = colNames, rowNames = rowNames)
        openxlsx::addStyle(wb, sheet, totalStyle,
                           rows = total_events_rows + 4, #including header
                           cols = 1:(length(colnames(tb)) + 1), 
                           gridExpand = TRUE)
        openxlsx::addStyle(wb, sheet, totalStyle,
                           rows = 4:(4 + n_rows_counts),
                           cols = total_events_cols + 1,
                           gridExpand = TRUE)
      } else {
        openxlsx::writeData(wb, sheet, x = tb, colNames = colNames, rowNames = rowNames)
        openxlsx::addStyle(wb, sheet, totalStyle,
                           rows = total_events_rows + 1,   
                           cols = 1:(length(colnames(tb)) + 1), 
                           gridExpand = TRUE)
        openxlsx::addStyle(wb, sheet, totalStyle,
                           rows = 1:(1 + n_rows_counts),
                           cols = total_events_cols + 1,
                           gridExpand = TRUE)
      }
    }else{
      warning("missing worksheet name, returning results instead")
      return(tb)
    }
  } else {
    return(tb)
  }
}

