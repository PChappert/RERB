
#### Create FlowJO-compatible barcodes from 10X barcodes and back ####
#'useful to import 10X ADT counts (or other) into FlowJo for analysis and fine gating of cells of interest for future import into FlowJo.
#'input can be "10X" or "FlowJo"
#'10X: will replace the "cell_id" column by 1 orig.idents (if !NULL) and 4 numbered id columns
#'FlowJo: will recreate the original 10X barcode based on the orig.idents (if !NULL) and numbered id columns
#'only needs as an input a vector with all orig.idents.

#TODO add compatibility to BD Rhapsody
FlowJoBarcode <- function(db, input = "10X", orig.idents = NULL, end = "-1"){
  if(input == "10X"){
    db$numbered_cell_id <- db$cell_id
    if(!is.null(orig.idents)){
      idents <- orig.idents[1]
      for(i in (2:length(orig.idents))){
        idents <- paste0(idents, "|", orig.idents[i])
      }
      db$orig.idents <- stringr::str_extract(db$cell_id, idents)
      db$orig.idents <- match(db$orig.idents, orig.idents)
      for(i in seq_along(orig.idents)){
        db$numbered_cell_id <- gsub(paste0(orig.idents[i], "_"), "", db$numbered_cell_id)
      }
    }
    if(!is.null(end)){
      db$numbered_cell_id <- stringr::str_sub(db$numbered_cell_id, start = 0, end = -(nchar(end)+1))
    }
    db$numbered_cell_id <- gsub("A", 1, db$numbered_cell_id)
    db$numbered_cell_id <- gsub("C", 2, db$numbered_cell_id)
    db$numbered_cell_id <- gsub("G", 3, db$numbered_cell_id)
    db$numbered_cell_id <- gsub("T", 4, db$numbered_cell_id)
    #cut new "numbered" cell_id in four to accomodate limitations in FlowJo
    db$numbered_cell_id_1 <- stringr::str_sub(db$numbered_cell_id, start = 1, end = 4)
    db$numbered_cell_id_2 <- stringr::str_sub(db$numbered_cell_id, start = 5, end = 8)
    db$numbered_cell_id_3 <- stringr::str_sub(db$numbered_cell_id, start = 9, end = 12)
    db$numbered_cell_id_4 <- stringr::str_sub(db$numbered_cell_id, start = 13, end = 16)
    return(db[,!colnames(db) %in% c("cell_id", "numbered_cell_id")])
  }
  if(input == "FlowJo"){
    db$cell_id <- paste0(db$numbered_cell_id_1, db$numbered_cell_id_2, db$numbered_cell_id_3, db$numbered_cell_id_4)
    db$cell_id <- gsub("1", "A", db$cell_id)
    db$cell_id <- gsub("2", "C", db$cell_id)
    db$cell_id <- gsub("3", "G", db$cell_id)
    db$cell_id <- gsub("4", "T", db$cell_id)
    if(!is.null(orig.idents)){
      db$orig.idents <- orig.idents[db$orig.idents]
      db$cell_id <- paste0(db$orig.idents, "_", db$cell_id)
    }
    if(!is.null(end)){
      db$cell_id <- paste0(db$cell_id, end)
    }
    return(db)
  }
}


#### Function to batch import export csv generated upon FlowJo gating of INX data ####
#' Import INX FlowJo gating results and merge it to an existing AIRR formatted data frame [optional]
#' 
#' \code{ImportFJGates} import a list of FlowJo export files and merge it 
#' @param FJ_files          a data frame with the following minimal columns : orig.ident = "sample_id", directory = "directory", filename = "filename", for plasticity, all three colnames can be modified and entered as separate arguments
#' @param db                an AIRR formatted data frame containing heavy and light chain sequences to which the imported FlowJo export data will be appended. if NULL, returning only the merged FlowJo export data
#' @param orig.ident        column in FJ_files to use to find the sample orig.ident [default = "sample_id"], will be used to recreate the beginning of the cell_id (should contain all info other than the plate_id and well_id)
#' @param directory         column in FJ_files to use to find the file directory
#' @param filename          column in FJ_files to use to find the full name of each file to be imported
#' @param import.columns    which columns to import from all imported files [default = "all, minus the ones removed (listed in the rm.columns parameter)]
#' @param rm.columns        which columns to remove from all imported files
#' @param overwrite.columns whether to overwrite columns in the imported sanger files with data provided in the sanger_dir recap file [default = FALSE].
#' @param plate_format      format of the plate used for sorting [default = "96"], should be one listed in Idx_to_character
#' @param Idx_to_character  list of dataframes (one per plate format) used to convert IdxRow column in the FlowJo export file to the full well_id.
#' @param export_missing_seq whether to export a table with all cell_id without corresponding sequence in the provided db.
#'
#' @return if a AIRR formatted data frame is provided, the same data frame with additional columns corresponding to the exported information from the INX files processed in FlowJo.
#' Else, a data frame with the merged data from the imported FlowJo export files.
#'
#' @export

ImportFJGates <- function(FJ_files, db = NULL, 
                          orig.ident = "sample_id",
                          directory = "directory",
                          filename = "filename",
                          import.columns = "all",
                          rm.columns = c("plate_id", "sort_id", "sort_order", "orig.ident", "donor_id", "time_point", "Time"),
                          overwrite.columns = TRUE,
                          plate_format = "96",
                          Idx_to_character = list("96" = data.frame(IdxRow = c(8, 7, 6, 5, 4, 3, 2, 1),
                                                                    Row = c("A", "B", "C", "D", "E", "F", "G", "H")),
                                                  "384" = data.frame(IdxRow = c(16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1),
                                                                     Row = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"))),
                          export_missing_seq = FALSE){
  
  library(dplyr)
  
  filenames <- ifelse(!stringr::str_ends(FJ_files[[directory]],"/"), paste0(FJ_files[[directory]], "/", FJ_files[[filename]]),
                      paste0(FJ_files[[directory]], FJ_files[[filename]]))
  
  if(any(duplicated(filenames))){
    message("Duplicated sample_id in the provided export FlowJO template: ", FJ_files[duplicated(FJ__files[[orig.ident]]), orig.ident])
    message("Will only import the first to avoid having issues with duplicated cell_ids")
    FJ_files <- FJ_files[!duplicated(FJ_files[[orig.ident]]), ]
    filenames <- ifelse(!stringr::str_ends(FJ_files[[directory]],"/"), paste0(FJ_files[[directory]], "/", FJ_files[[filename]]),
                        paste0(FJ_files[[directory]], FJ_files[[filename]]))
  }
  
  add_columns <- setdiff(colnames(FJ_files), c("sample_id", "directory", "filename"))
  
  file.list <- lapply(filenames, FUN=function(file){
    data <- readr::read_csv(file, show_col_types = FALSE)
    if(nrow(data)>0){
      colnames(data) <- gsub("-", "_", colnames(data))
      colnames(data) <- gsub(" ", "_", colnames(data))
      colnames(data) <- gsub(":", "_", colnames(data)) #'simplify parameter:dye labels in SONY or BD outputs
      colnames(data) <- gsub("TIME", "Time", colnames(data)) #'harmonize between SONY and BD
      data <- data %>%
        dplyr::mutate(across(everything(), as.numeric)) #'circumvent some remaining issues with NAs...
      data$orig.ident <- FJ_files[match(file, filenames), orig.ident]
      for(col in add_columns){
        if(overwrite.columns){
          data[[col]] <- FJ_files[match(file, filenames), col]
        } else {
          if(!col %in% colnames(data)){data[[col]] <- FJ_files[match(file, filenames), col]}
        }
      }
      return(data)
    } else {return(NULL)}
  })
  FJ_db <- dplyr::bind_rows(file.list)
  
  #recreate cell_ids
  FJ_db <- FJ_db %>%
    dplyr::mutate(
      cell_id = paste0(orig.ident, "_P", plate_id, "_", Idx_to_character[[plate_format]][match(FJ_db$IdxRow, Idx_to_character[[plate_format]]$IdxRow), "Row"], FJ_db$IdxCol)
    ) 
  
  #resolve cases of duplicated cell_id
  #'we use here FSC_A and Time available in both SONY and BD index sort data
  FJ_db <- FJ_db %>%
    dplyr::group_by(cell_id) %>%
    dplyr::mutate(
      FSC_unique = n_distinct(.data[["FSC_A"]]), #same time but different FSC could be a doublet
      time_unique = n_distinct(.data[["Time"]]) #different time can only be an issue in excel recaps/INXs...
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      warning_flag = ifelse(time_unique > 1 | FSC_unique > 1, TRUE, FALSE)
    ) 
  
  flagged_FJ_db <- FJ_db %>%
    dplyr::filter(warning_flag) 
  
  if(nrow(flagged_FJ_db)>0){
    flagged_FJ_db <- flagged_FJ_db %>%
      dplyr::group_by(cell_id) %>%
      dplyr::summarise(
        across(
          everything(),
          ~ first(.), #'only keep the first occurence and send a warning
          .names = "{.col}"
        ), 
        .groups = "drop"
      )
    warning("Multiple values for FSC_A and Time found for the following cell_id: ", paste(flagged_FJ_db$cell_id, collapse = ", "), " only the first occurence was kept for now")
  }
    
  resolved_FJ_db <- FJ_db %>%
    dplyr::filter(!warning_flag) %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(
      across(
        everything(),
        ~ ifelse(n_distinct(.) == 1, first(.), paste(unique(sort(.)), collapse = "|")),
        .names = "{.col}"
      ), 
      .groups = "drop"
    ) %>%
    dplyr::bind_rows(flagged_FJ_db)
  
  if(import.columns[1] == "all"){
    import.columns <- c(setdiff(colnames(FJ_db), rm.columns))
  }
  
  resolved_FJ_db <- resolved_FJ_db %>%
    dplyr::select(unique(c("cell_id", import.columns,  "orig.ident"))) 
  
  
  
  if(!is.null(db)){
    preexisting_columns <- setdiff(intersect(colnames(db), colnames(resolved_FJ_db)), "cell_id") 
    
    if(export_missing_seq){
      missing_seqs <- setdiff(resolved_FJ_db$cell_id, db$cell_id)
      missing_seqs_FJ_db <- resolved_FJ_db %>%
        dplyr::filter(cell_id %in% missing_seqs)
      readr::write_tsv(missing_seqs_FJ_db, file = "FJ_cells_missing_VDJ.tsv")
    }
    
    db <- db %>%
      dplyr::select(-all_of(preexisting_columns)) %>%
      dplyr::left_join(resolved_FJ_db, by = "cell_id")
    
    return(db)
  } else {
    return(resolved_FJ_db)
  }
}







