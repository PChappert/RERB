
#### Create FlowJO-compatible barcodes from 10X barcodes and back ####
#' creates FlowJO compatible barcodes and back
#'
#' \code{FlowJoBarcode} creates FlowJO compatible barcodes and back
#' @param db        an AIRR formatted data frame containing heavy and light chain sequences to which the imported FlowJo export data will be appended. if NULL, returning only the merged FlowJo export data
#' @param input     10X or FlowJo (TODO add compatibility to BD Rhapsody)
#' @param orig.idents  column in FJ_files to use to find the file directory
#' @param end          optional value to add or remove after the actual barcode [default for 10X: -1]
#'
#' @export
#'
#' @importFrom stringr str_sub
#'
#' @details
#' useful to import 10X ADT counts (or other) into FlowJo for analysis and fine gating of cells of interest for future import into FlowJo.
#' input can be "10X" or "FlowJo"
#' 10X: will replace the "cell_id" column by 1 orig.idents (if !NULL) and 4 numbered id columns
#' FlowJo: will recreate the original 10X barcode based on the orig.idents (if !NULL) and numbered id columns
#' only needs as an input a vector with all orig.idents.

FlowJoBarcode <- function(db,
                          input = "10X",
                          orig.idents = NULL,
                          end = "-1"){

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
#' Imports INX FlowJo gating results and merge it to an existing AIRR formatted data frame [optional]
#'
#' \code{importFJGates} imports a list of FlowJo export files and merge it
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
#'
#' @import dplyr
#' @import readr
#' @importFrom stringr str_ends

importFJGates <- function(FJ_files, db = NULL,
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

  #library(dplyr)

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
      colnames(data) <- gsub(":", "_", colnames(data)) #simplify parameter:dye labels in SONY or BD outputs
      colnames(data) <- gsub("TIME", "Time", colnames(data)) #'harmonize between SONY and BD
      data <- data %>%
        dplyr::mutate(across(everything(), as.numeric)) #circumvent some remaining issues with NAs...
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
  #we use here FSC_A and Time available in both SONY and BD index sort data
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
          ~ first(.), #only keep the first occurence and send a warning
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


#### Function to plot a 96w or 384w plate with highlighted wells ####
#'
#' \code{plotPlate} Plot a single 96 or 384 well plate with highlighted wells
#' @param db                  a data frame containing at least a column with well_ids.
#' @param highlighted_wells   name of column containing well_id to plot or a vector containing well-ids (set db = NULL in this case).
#' @param color.by             name of column containing group labels for highlighting wells or a vector of group labels the same size as highlighted_wells (set db = NULL in this case).
#' @param fill_colors         a named vector with one color per group in color.by
#' @param labels              name of column containing text labels for highlighting wells or a vector of text labels the same size as highlighted_wells (set db = NULL in this case) [default = highlighted_wells].
#' @param main_title          name of the graph title 
#' @param legend_title        name of the legend title
#' @param plate_type          type of plate to plot, should be one of "96w" or "384w" [default = "96w"]
#' @param well_size           size of well in the final plot, if set to NULL, will be automatically define based on plate_type
#' @param fixed_size          whether to fix the size of the inside panel (requires the egg package)
#' @param panel_width         panel width for fixed panel size
#' @param panel_height        panel height for fixed panel size
#' @param return_plot         whether to return ggplot objects else simply plotted
#'
#' @return a ggplot object
#'    
#' @import dplyr
#' @import tibble
#'
#' @export

plotPlate <- function(db = NULL,
                      highlighted_wells = NULL,
                      color.by = NULL, 
                      fill_colors = NULL,
                      plot_labels = TRUE,
                      labels = NULL,
                      main_title = "",
                      legend_title = NULL,
                      plate_type = c("96w", "384w"),
                      well_size = NULL,
                      fixed_size = TRUE,
                      panel_width = NULL,
                      panel_height = NULL,
                      return_plot = FALSE) {
                            
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Optional: 'ggplot2' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(ggplot2))
  
  if (!requireNamespace("scales", quietly = TRUE)) {
    message("Optional: 'scales' not installed — skipping plot.")
    return(invisible(NULL))
  }
  
  plate_type <- match.arg(plate_type)
  
  if(is.null(db)){
    if(!is.null(highlighted_wells)){
      if(is.null(labels)){
        labels <- highlighted_wells
      }
      db <- tibble::tibble(well_id = highlighted_wells, fill = color.by, label = labels)
      color.by <- "fill"
    } else {
      message("if not providing a full dataframe, highlighted_wells should be a vector of well-ids to highlight - skipping plot")
      return(invisible(NULL))
    }
  } else {
    if(is.null(highlighted_wells)){
      highlighted_wells <- "well_id"
    }
    if(!highlighted_wells %in% colnames(db)){
      message("if providing a full dataframe, highlighted_wells should be a column containing the well-ids to highlight - skipping plot")
      return(invisible(NULL))
    }
    if(!highlighted_wells == "well_id"){
      if("well_id" %in% colnames(db)){
        db <- db %>%
          dplyr::mutate(
            original_well_id = well_id
          )
        if(!is.null(labels) && labels == "well_id") {
          labels == "original_well_id"
        }
      }
      db <- db %>%
        dplyr::mutate(
          well_id = !!rlang::sym(highlighted_wells)
          )
    }
    if(is.null(labels)){
      labels <- "well_id"
    }
    db <- db %>%
      dplyr::mutate(
        label = !!rlang::sym(labels)
      )
  }
  
  # Define row and column names and create full plate layout for empty plate:
  plate_format <- list("96w" = c(8, 12), "384w" = c(16, 24))
  if(is.null(well_size)){
    well_sizes <- list("96w" = 8, "384w" = 8)
    well_size <- well_sizes[[plate_type]]
  }
  if(is.null(panel_width)){
    panel_widths <- list("96w" = 3.5, "384w" = 7)
    panel_width <- panel_widths[[plate_type]]
  }
  if(is.null(panel_height)){
    panel_heights <- list("96w" = 2.5, "384w" = 4.7)
    panel_height <- panel_heights[[plate_type]]
  }
  
  rows <- LETTERS[1:plate_format[[plate_type]][1]]
  cols <- 1:plate_format[[plate_type]][2]
  
  plate <- expand.grid(col = cols, row = rows)
  plate$well_id <- paste0(plate$row, plate$col)
  
  # Add info for wells to highlight:
  plate <- plate %>%
    left_join(
      db,
      by = "well_id"
    )
  
  # Add flag and color for highlighted wells:
  if(is.null(color.by)){
    if(is.null(legend_title)){
      legend_title <- "highlighted"
    } 
    plate[[legend_title]] <- plate$well_id %in% db$well_id
    if(is.null(fill_colors)){
      fill_colors <- c("TRUE" = "red", "FALSE" = "white")
    } else {
      fill_colors <- c("TRUE" = fill_colors[1], "FALSE" = "white")
    }
  } else {
    if(is.null(legend_title)){
      legend_title <- color.by
    } else {
      plate <- plate %>%
        dplyr::mutate(
          !!rlang::sym(legend_title) := !!rlang::sym(color.by)
        )
    }
    if(is.null(fill_colors)){
      fill_colors <- scales::hue_pal()(length(unique(db[[color.by]])))
      names(fill_colors) <- sample(unique(db[[color.by]]))
    }
  }
  
  if(plot_labels){
    g <- ggplot(plate, aes(x = col, y = row)) +
      geom_point(aes(fill = !!rlang::sym(legend_title)), shape = 21, size = well_size, color = "black") +
      geom_text(data = dplyr::filter(plate, !is.na(label)), aes(label = label), size = 3, color = "black") +
      scale_y_discrete(limits = rev(rows), position = "left") +
      scale_x_continuous(breaks = cols, position = "top") +
      scale_fill_manual(values = fill_colors, na.value = "white") +
      coord_fixed(ratio = 1, xlim = c(0.5, (plate_format[[plate_type]][2]+0.5)), ylim = c(0.5, (plate_format[[plate_type]][1]+0.5)), expand = FALSE, clip = "off") +
      theme_minimal(base_size = 12) +
      labs(title = main_title, x = NULL, y = NULL) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 16)
      )
  } else {
    g <- ggplot(plate, aes(x = col, y = row)) +
      geom_point(aes(fill = !!rlang::sym(legend_title)), shape = 21, size = well_size, color = "black") +
      scale_y_discrete(limits = rev(rows), position = "left") +
      scale_x_continuous(breaks = cols, position = "top") +
      scale_fill_manual(values = fill_colors, na.value = "white") +
      coord_fixed(ratio = 1, xlim = c(0.5, (plate_format[[plate_type]][2]+0.5)), ylim = c(0.5, (plate_format[[plate_type]][1]+0.5)), expand = FALSE, clip = "off") +
      theme_minimal(base_size = 12) +
      labs(title = main_title, x = NULL, y = NULL) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 16)
      )
  }
  if(fixed_size){
    if (requireNamespace("egg", quietly = TRUE)) {
      g <- egg::set_panel_size(
        g,
        width  = grid::unit(panel_width, "in"),
        height = grid::unit(panel_height, "in")
      )
    } else {
      message("Optional: 'egg' not installed — skipping fixed_size argument.")
    }
  }
  if(return_plot) {
    return(g)
  } else {
    plot(g)
  }
}

#### Function to plot a recap of 96w or 384w plates with highlighted wells ####
#'
#' \code{plotPlateRecap} Plot a single 96 or 384 well plate with highlighted wells
#' @param db                  a data frame containing at least a column with well_ids.
#' @param well_id             name of column containing well_id informations.
#' @param split.by            name of column(s) to use to split the dataset in individual plates.
#' @param return_plot         whether to return individual ggplot objects as a list
#' @param save_plot           whether to save the ggplot objects
#' @param plot_filename       name of the saved plot file
#' @param save_as             one of "pdf" or "png" [default = "pdf]
#' @param nrow                number of rows per page for saving
#' @param ncol                number of columns per page for saving
#' @param excel_report        whether to output a full excel recap with well infos and plate plots
#' @param excel_filename      name of the saved excel recap file
#' @param ...                 arguments to pass to plotPlate()
#'
#' @return a list of ggplot object
#'
#' @import dplyr
#' @import tibble
#' @import purrr
#'
#' @export

plotPlateRecap <- function(db,
                            well_id = "well_id",
                            split.by = NULL,
                            return_plot = FALSE,
                            save_plot = TRUE,
                            nrow = 2,
                            ncol = 2,
                            plot_filename = "multi_page_plate_plot",
                            save_as = c("pdf", "png"),
                            excel_report = FALSE,
                            excel_filename = "recap",
                            ...){
                       
  #define well_id column
  if(!well_id == "well_id"){
    if("well_id" %in% colnames(db)){
      db <- db %>%
        dplyr::mutate(
          original_well_id = well_id,
          well_id = !!rlang::sym(well_id)
        )
      if("well_id" %in% split.by){
        split.by <- c(split.by[!split.by == "well_id"], "original_well_id")
      }
    } else {
      db <- db %>%
        dplyr::mutate(
          well_id = !!rlang::sym(well_id)
        )
    }
  }
  #order properly based on well_id (A1, A2,...A10, A11,...) and group wells according to split.by argument (allowing for multiple plates)
  grouped <- db %>%
    dplyr::mutate(
      col = gsub("[0-9]", "", well_id),
      row = as.integer(gsub("[A-Z]", "", well_id))
    ) %>%
    dplyr::arrange(!!!rlang::syms(split.by), factor(col, levels = LETTERS[1:16]), row) %>%
    dplyr::select(-col, -row) %>%
    dplyr::group_by(!!!rlang::syms(split.by)) 
  
  split_plates <- grouped %>%
    dplyr::group_split(.keep = TRUE)
  
  group_names <- grouped %>%
    dplyr::group_keys() %>%
    tidyr::unite("group_label", everything(), sep = "_") %>%
    dplyr::pull(group_label)
  
  names(split_plates) <- group_names
  
  #plot.list <- purrr::imap(split_data,
  #                         ~ plotPlate(.x, highlighted_wells = "well_id", main_title = .y, return_plot = TRUE, ...))
  
  plot.list <- purrr::imap(split_plates, function(df, plate_name) {
    plotPlate(df,
              highlighted_wells = "well_id",
              labels = "well_id",
              main_title = plate_name,
              return_plot = TRUE,
              ...
    )
  })
  
  if(save_plot){
    if (requireNamespace("patchwork", quietly = TRUE)) {
      chunked_plots <- split(plot.list, ceiling(seq_along(plot.list) / (nrow * ncol)))
      
      save_as <- match.arg(save_as)
      if (save_as == "pdf") {
        pdf(paste0(plot_filename, ".pdf"), width = 11.7, height = 8.3)  # A4 landscape
        purrr::walk(chunked_plots, function(plot_page) {
          combined <- patchwork::wrap_plots(plot_page, nrow = nrow, ncol = ncol)
          print(combined)
        })
        dev.off()
      }
      if (save_as == "png") {
        png(paste0(plot_filename, ".png"), width = 11.7, height = 8.3)  # A4 landscape
        purrr::walk(chunked_plots, function(plot_page) {
          combined <- patchwork::wrap_plots(plot_page, nrow = nrow, ncol = ncol)
          print(combined)
        })
        dev.off()
      }
      
    } else {
      message("Optional: 'patchwork' not installed — skipping saving plot.")
      }
  }
  
  if(excel_report){
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      wb <- openxlsx::createWorkbook()
      temp_dir <- tempdir()  # save png plot to temp folder
      s <- purrr::imap(plot.list, function(p, plate_name) {
        # Save the plot to a PNG file
        plot_file <- file.path(temp_dir, paste0(plate_name, ".png"))
        ggsave(plot_file, p, width = 6, height = 4, dpi = 300)
        
        # Add a worksheet and insert image
        openxlsx::addWorksheet(wb, plate_name)
        openxlsx::insertImage(wb, sheet = plate_name, file = plot_file, width = 6, height = 4, startRow = 1, startCol = 1)
        
        # Optionally also write the data for that plate
        df <- split_plates[[plate_name]]
        openxlsx::writeData(wb, sheet = plate_name, x = df, startRow = 20, startCol = 1)
      })
      openxlsx::saveWorkbook(wb, file = paste0(excel_filename, ".xlsx"), overwrite = TRUE)
    } else {
      message("Optional: 'openxlsx' not installed — skipping excel recap.")
    }
  }
  
  if(return_plot) {
    return(plot.list)
  }
}


#### Function to create consolidated plate(s) based on selected well-ids from a list of plates ####
#'
#' \code{consolidatePlates} create consolidated plate(s)
#' @param db                  a data frame containing at least a column with well_ids.
#' @param well_id             name of column containing well_id informations.
#' @param split.by            name of column(s) to use to split the dataset in individual plates.
#' @param plate_type          type of plate to plot, should be one of "96w" or "384w" [default = "96w"]
#' @param fill.by             whether to fill the consolidated plate(s) by row (A1 -> A12 then B1 -> B12...) or by column (A1 -> H1 then A2 -> H2...) [default = "row"]
#' @param empty_wells         vector of well-ids to keep as empty wells in the final plate(s) [Eurofins expects "H12" wells to be empty for example]
#' @param return_plot         whether to return individual ggplot objects as a list
#' @param save_plot           whether to save the ggplot objects
#' @param nrow                number of rows per page for saving
#' @param ncol                number of columns per page for saving
#' @param plot_filename       name of the saved plot file
#' @param save_as             one of "pdf" or "png" [default = "pdf]
#' @param excel_report        whether to output a full excel recap with well infos and plate plots
#' @param excel_filename      name of the saved excel recap file
#' @param ...                 arguments to pass to plotPlate()
#'
#' @return new template(s) colored and labelled based on plate and well of origin respectively and a full excel report with all info from well of origin.
#'
#' @import dplyr
#' @import tibble
#' @import purrr
#'
#' @export

consolidatePlates <- function(db,
                              well_id = "well_id",
                              split.by = NULL,
                              plate_type = c("96w", "384w"),
                              fill.by = c("row", "col"),
                              empty_wells = NULL,
                              return_plot = FALSE,
                              save_plot = TRUE,
                              nrow = 2,
                              ncol = 2,
                              plot_filename = "multi_page_plate_plot",
                              save_as = c("pdf", "png"),
                              excel_report = FALSE,
                              excel_filename = "recap",
                              ...){
  
  plate_type <- match.arg(plate_type)
  fill.by <- match.arg(fill.by)
  
  plate_format <- list("96w" = c(8, 12), "384w" = c(16, 24))
  #if(is.null(well_size)){
  #  well_size <- list("96w" = 8, "384w" = 4)
  #  well_size <- well_size[[plate_type]]
  #}
  rows <- LETTERS[1:plate_format[[plate_type]][1]]
  cols <- 1:plate_format[[plate_type]][2]
  
  if(fill.by == "row"){
    new_plate_wells <- new_plate_wells <- expand.grid(row = rows, col = cols) %>%
      dplyr::arrange(row, col) %>%
      dplyr::mutate(well = paste0(row, col)) %>%
      dplyr::pull(well)
  }
  if(fill.by == "col"){
    new_plate_wells <- as.vector(outer(rows, cols, paste0))
  }
  new_plate_wells <- new_plate_wells[!new_plate_wells %in% empty_wells] 
  
  db <- db %>%
    dplyr::mutate(
      col = gsub("[0-9]", "", !!rlang::sym(well_id)),
      row = as.integer(gsub("[A-Z]", "", !!rlang::sym(well_id)))
    ) %>%
    dplyr::arrange(sort_id, plate_id, factor(col, levels = LETTERS[1:16]), row) 
  
  if(!is.null(split.by)){
    db <- db %>%
      dplyr::mutate(group_label = paste(!!!rlang::syms(split.by), sep = "_"))
  } else {
    db$group_label <- "original_plate"
  }
  
  n <- nrow(db)
  db_recap <- db %>%
    dplyr::mutate(
      original_plate_id = plate_id,
      original_well_id = well_id,
      new_plate_id = paste0("P", ceiling(row_number() / (96-length(empty_wells)))),
      new_well_id = rep(new_plate_wells, length.out = n)
    ) %>%
    dplyr::select(-col, -row, -well_id, -plate_id)
  
  group_colors <- scales::hue_pal()(length(unique(db$group_label)))
  names(group_colors) <- sample(unique(db$group_label)) #using here sample() to randomise the color order on the plate for better visualisation
  
  # Split data by new_plate_num
  split_plates <- split(db_recap, db_recap$new_plate_id)
  
  # Generate one plot per new plate
  new_plate_plots <- purrr::imap(split_plates, function(df, plate_num) {
    plotPlate(df,
              highlighted_wells = "new_well_id",
              color.by = "group_label",
              fill_colors = group_colors,
              labels = "original_well_id",
              legend_title = "plate of origin",
              main_title = paste0("LC_sequencing_", plate_num),
              return_plot = TRUE,
              ...
    )
  })
  
  if(save_plot){
    if (requireNamespace("patchwork", quietly = TRUE)) {
      chunked_plots <- split(new_plate_plots, ceiling(seq_along(new_plate_plots) / (nrow * ncol)))
      
      save_as <- match.arg(save_as)
      if (save_as == "pdf") {
        pdf(paste0(plot_filename, ".pdf"), width = 11.7, height = 8.3)  # A4 landscape
        purrr::walk(chunked_plots, function(plot_page) {
          combined <- patchwork::wrap_plots(plot_page, nrow = nrow, ncol = ncol)
          print(combined)
        })
        dev.off()
      }
      if (save_as == "png") {
        png(paste0(plot_filename, ".png"), width = 11.7, height = 8.3)  # A4 landscape
        purrr::walk(chunked_plots, function(plot_page) {
          combined <- patchwork::wrap_plots(plot_page, nrow = nrow, ncol = ncol)
          print(combined)
        })
        dev.off()
      }
    } else {
      message("Optional: 'patchwork' not installed — skipping saving plot.")
    }
  }
  
  if(excel_report){
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      wb <- openxlsx::createWorkbook()
      temp_dir <- tempdir()  # save png plot to temp folder
      s <- purrr::imap(new_plate_plots, function(p, plate_num) {
        # Save the plot to a PNG file
        plot_file <- file.path(temp_dir, paste0("LC_seq_", plate_num, ".png"))
        ggsave(plot_file, p, width = 6, height = 4, dpi = 300)
        
        # Add a worksheet and insert image
        sheet_name <- paste0("LC_seq_", plate_num)
        openxlsx::addWorksheet(wb, sheet_name)
        openxlsx::insertImage(wb, sheet = sheet_name, file = plot_file, width = 6, height = 4, startRow = 1, startCol = 1)
        
        # Optionally also write the data for that plate
        df <- split_plates[[as.character(plate_num)]]
        openxlsx::writeData(wb, sheet = sheet_name, x = df, startRow = 20, startCol = 1)
      })
      openxlsx::addWorksheet(wb, "Recap Table")
      openxlsx::writeData(wb, sheet = "Recap Table", x = db_recap)
      openxlsx::saveWorkbook(wb, file = paste0(excel_filename, ".xlsx"), overwrite = TRUE)
    } else {
      message("Optional: 'openxlsx' not installed — skipping excel recap.")
    }
  }
  
  if(return_plot) {
    return(new_plate_plots)
  }
}
