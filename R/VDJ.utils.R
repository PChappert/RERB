
#### Function to directly process plate-based .Ab1 file ####
#' unzip Eurofins-based results and run all steps from QC to final AIRR formatted datatable
#'
#' \code{Ab1toAIRR}
#'
#' @param files   a vector of full file paths to Ab1.zip files
#' @param primers primers used (one of IgG, IgM, IgL, IgK or Mix)
#' @param igblast whether to run igblast
#' @param igblast_dir location of the igblast database
#' @param igblast_method whether to run runAssignGenes or runIgblastn
#' @param seq_type one of Ig or TCR (passed to igblast)
#' @param organism one of human, mouse,
#' @param reference_dir location of the imgt gap database
#' @param update_c_call whether to update c_call using runBlastnC
#' @param SHM     whether to calculate mutations in V gene
#' @param full_seq_aa whether to reconstruct full AA sequences for each sequence
#' @param verbose whether to write to the console
#' @param return_db whether to return the final productive sequences data table
#' @param save    passed to runAb1QC, whether to save png or html plots
#' @param reticulate_py_env if you want to force a precise python environment to use for plotly. [default: NULL]
#' @param nproc   number of processors to use, if NULL, will be automatically set based on available processors.
#'
#'
#' @return A data table containing all sequences passing QC and igblast.
#' All QC analysis plots, intermediate tsv files with sequences passing and failing igblast and an excel worksheet containing sheets with AIRR formatted data table for sequences failing initial QC, failing igblast, non-productive and passing all three steps.
#'
#' @export
#'
#' @import dplyr
#' @import fs
#' @import stringr
#' @import shazam
#' @importFrom dowser createGermlines
#' @importFrom dowser readIMGT

Ab1toAIRR <- function(files,
                      primers = c("IgG", "IgM", "IgL", "IgK", "Mix"),
                      trim_cutoff = list("IgG" = 300, "IgM" = 275, "IgK" = 400, "IgL"= 275, "Mix" = 275, "Undefined" = 275),
                      QC = TRUE,
                      igblast = TRUE,
                      seq_type = c("Ig", "TCR"),
                      organism = c("human", "mouse", "rabbit", "rat", "rhesus_monkey"),
                      igblast_log = TRUE,
                      igblast_dir = "~/share/igblast/",
                      imgt_dir = "~/share/germlines/imgt/",
                      igblast_method = c("AssignGenes", "runIgblastn"),
                      update_c_call = TRUE,
                      SHM = TRUE,
                      full_seq_aa = TRUE,
                      update_info = FALSE,
                      verbose = TRUE,
                      return_db = FALSE,
                      save = c("png", "html"),
                      reticulate_py_env = "~/Library/r-miniconda-arm64/envs/r-reticulate/bin/python",
                      nproc = NULL){

  primers <- match.arg(primers)
  save <- match.arg(save)
  seq_type <- match.arg(seq_type)
  organism <- match.arg(organism)
  igblast_method <- match.arg(igblast_method)

  if(is.null(nproc)){
    if (requireNamespace("parallel", quietly = TRUE)) {
      cores <- parallel::detectCores()
      nproc <- max((cores[1]-2), 1) #not to overload computer
    } else {
      message("Optional: 'parallel' not installed — defaulting to nproc = 1")
      nproc <- 1
    }
  }

  if(!igblast){
    update_c_call = FALSE
    SHM = FALSE
    full_seq_aa = FALSE
  }
  if(update_info){
    QC = FALSE
    igblast = FALSE
    update_c_call = FALSE
    SHM = FALSE
    full_seq_aa = FALSE
  }

  # Check each filepath, remove missing one and send a warning
  file_exists <- file.exists(files)

  missing_files <- files[!file_exists]

  files <- files[file_exists]
  if(length(files) == 0){
    stop("no file found, check file(s) path(s)")
  }
  if(length(missing_files)>0){
    warning("the following files could not be found: ", paste0(missing_files, collapse = "; "))
  }

  airr.list <- lapply(files, function(filepath){
    start <- Sys.time()
    # Get filename without extension
    filename <- fs::path_ext_remove(fs::path_file(filepath))
    file_folder <- fs::path_dir(filepath)
    # Replace pattern in filename to create output folder name
    outfilename <- stringr::str_replace(filename, "_SCF_SEQ_ABI", "")
    outfolder <- paste0(file_folder, "/", outfilename)
    # Create log file path
    log_file <- paste0(outfolder, "/", outfilename, ".log")
    
    # Create output folder
    if(!dir.exists(outfolder)){
      dir.create(outfolder)
      # allows using temp_folder to unzip the data from eurofins, while outputting data in a seperate folder
    }
    
    if(update_info){
      #reopen already QCed and aligned db
      #check first for full-pass db:
      filename_pass <- paste0(outfolder, "/", outfilename, "_igblast-pass_c-call-pass_germ-pass_shm-pass_full-pass.tsv.gz")
      if(file.exists(filename_pass)){
        if(verbose){cat("opening", filename_pass, "\n")}
        VDJ_db <- readr::read_tsv(filename_pass, show_col_types = FALSE)
      } else {
        #check first for shm-pass db:
        filename_pass <- paste0(outfolder, "/", outfilename, "_igblast-pass_c-call-pass_germ-pass_shm-pass.tsv.gz")
        if(file.exists(filename_pass)){
          if(verbose){cat("opening", filename_pass, "\n")}
          VDJ_db <- readr::read_tsv(filename_pass, show_col_types = FALSE)
        } else {
          #check first for full-pass db:
          filename_pass <- paste0(outfolder, "/", outfilename, "_igblast-pass_c-call-pass_germ-pass.tsv.gz")
          if(file.exists(filename_pass)){
            if(verbose){cat("opening", filename_pass, "\n")}
            VDJ_db <- readr::read_tsv(filename_pass, show_col_types = FALSE)
          } else {
            #check first for full-pass db:
            filename_pass <- paste0(outfolder, "/", outfilename, "_igblast-pass_c-call-pass.tsv.gz")
            if(file.exists(filename_pass)){
              if(verbose){cat("opening", filename_pass, "\n")}
              VDJ_db <- readr::read_tsv(filename_pass, show_col_types = FALSE)
            } else {
              filename_pass <- paste0(outfolder, "/", outfilename, "_igblast-pass.tsv.gz")
              if(file.exists(filename_pass)){
                if(verbose){cat("opening", filename_pass, "\n")}
                VDJ_db <- readr::read_tsv(filename_pass, show_col_types = FALSE)
              } else {
                filename_pass <- paste0(outfolder, "/", outfilename, "_QC-pass.tsv")
                if(file.exists(filename_pass)){
                  stop("no aligned file found in : ", outfolder, " run the partial Ab1toAIRR pipeline with QC = FALSE option")
                } else {
                  stop("no aligned or QCed file found in : ", outfolder, " run the full Ab1toAIRR pipeline")
                }
              }
            }
          }
        }
      }
      filename_igblast_fail <- paste0(outfolder, "/", outfilename, "_igblast-fail.tsv.gz")
      if(file.exists(filename_igblast_fail)){
        failed_VDJ_db <- readr::read_tsv(filename_igblast_fail, show_col_types = FALSE)
      } else {
        if(verbose){cat("file: ",filename_igblast_fail, " not found")}
        failed_VDJ_db <- NULL
      }
      filename_fail <- paste0(outfolder, "/", outfilename, "_QC-fail.tsv")
      if(file.exists(filename_fail)){
        QC_failed <- readr::read_tsv(filename_fail, show_col_types = FALSE)
      } else {
        if(verbose){cat("file: ",filename_fail, " not found")}
        QC_failed <- NULL
      }
    }
    
    if(QC){
      time_and_log({
        cat("filename: ", filename, "\n")
        cat("organism: ", organism, "\n")
        cat("seq_type: ", seq_type, "\n")
        cat("primers: ", primers, "\n")
      }, verbose = FALSE, time = FALSE, log_file = log_file, log_title = "Ab1toAIRR", open_mode = "wt")
      
      # Step1: unzip to 'temp_folder' and redirect output to log file
      if(verbose){cat("unzipping", filepath, "\n")}
      time_and_log({
        cat(paste0("unzipping ", filepath))
        unzip(filepath, exdir = paste0(outfolder, "/temp_folder"))
      }, verbose = FALSE, log_file = log_file, log_title = "Unzipping .Ab1.zip file", open_mode = "a")

      # Step2: run QC:
      if(verbose){cat("starting QC", "\n")}
      time_and_log({
        QC_results <- runAb1QC(Ab1_folder = paste0(outfolder, "/temp_folder"),
                               outfolder = outfolder,
                               outfilename = outfilename,
                               primers = primers,
                               trim_cutoff = trim_cutoff,
                               save = save,
                               reticulate_py_env = reticulate_py_env,
                               nproc = nproc)
      }, verbose = FALSE, log_file = log_file, log_title = "Running QC", open_mode = "a")

      VDJ_db <- QC_results[["pass"]]
      QC_failed <- QC_results[["fail"]]
      unlink(paste0(outfolder, "/temp_folder"), recursive = TRUE) #remove temp folder
      
    } else {
      if(!update_info){
        time_and_log({
          cat("filename: ", filename, "\n")
          cat("organism: ", organism, "\n")
          cat("seq_type: ", seq_type, "\n")
          cat("primers: ", primers, "\n")
        }, verbose = FALSE, time = FALSE, log_file = log_file, log_title = "New run without QC", open_mode = "a")
        
        #reopen already QCed Ab1 files
        filename_pass <- paste0(outfolder, "/", outfilename, "_QC-pass.tsv")
        if(file.exists(filename_pass)){
          if(verbose){cat("opening", filename_pass, "\n")}
          VDJ_db <- readr::read_tsv(filename_pass, show_col_types = FALSE)
        } else {
          if(verbose){cat("file: ",filename_pass, " not found")}
          VDJ_db <- NULL
        }
        filename_fail <- paste0(outfolder, "/", outfilename, "_QC-fail.tsv")
        if(file.exists(filename_fail)){
          QC_failed <- readr::read_tsv(filename_fail, show_col_types = FALSE)
        } else {
          if(verbose){cat("file: ",filename_fail, " not found")}
          QC_failed <- NULL
        }
      }
    }
    
    if(nrow(VDJ_db)>0){
      # Step3: run Igblast:
      if(igblast){
        if(verbose){cat("running Igblast using", igblast_method, "\n")}
        if(igblast_method == "AssignGenes"){
          time_and_log({
            igblast_results <- runAssignGenes(VDJ_db,
                                              sequence = "sequence",
                                              sequence_id = "well_id",
                                              organism = organism,
                                              seq_type = seq_type,
                                              igblast_dir = igblast_dir,
                                              imgt_dir = imgt_dir,
                                              log = igblast_log)
          }, verbose = FALSE, log_file = log_file, log_title = igblast_method, open_mode = "a")
        }

        if(igblast_method == "runIgblastn"){
          time_and_log({
            igblast_results <- runIgblastn(VDJ_db,
                                           sequence = "sequence",
                                           sequence_id = "well_id",
                                           organism = organism,
                                           seq_type = seq_type,
                                           igblast_dir = igblast_dir,
                                           imgt_dir = imgt_dir)
          }, verbose = FALSE, log_file = log_file, log_title = igblast_method, open_mode = "a")
        }

        VDJ_db <- igblast_results[["pass"]]
        VDJ_db$c_call_igblast <- VDJ_db$c_call

        failed_VDJ_db <- igblast_results[["fail"]]

        VDJ_db <- dplyr::relocate(VDJ_db, c_call, .after = j_call)

        if(!any(is.null(failed_VDJ_db), nrow(failed_VDJ_db)==0)){
          outfilename_fail <- paste0(outfilename, "_igblast-fail")
          readr::write_tsv(failed_VDJ_db, file = paste0(outfolder, "/", outfilename_fail, ".tsv.gz"))
        }

        outfilename <- paste0(outfilename, "_igblast-pass")
        if(!(update_c_call|SHM|full_seq_aa)){#we only save this file if no other analysis is performed
          readr::write_tsv(VDJ_db, file = paste0(outfolder, "/", outfilename, ".tsv.gz"))
        }
      }
      # Step4: run BlastnC to update c_call:
      if(update_c_call){
        time_and_log({
          VDJ_db <- runBlastnC(VDJ_db,
                               igblast_dir = igblast_dir)
        }, verbose = FALSE, log_file = log_file, open_mode = "a")
        outfilename <- paste0(outfilename, "_c-call-pass")
        if(!(SHM|full_seq_aa)){#we only save this file if no other analysis is performed
          readr::write_tsv(VDJ_db, file = paste0(outfolder, "/", outfilename, ".tsv.gz"))
        }
      }

      # Step5: calculate mutation loads:
      if(SHM){
        if(verbose){cat("adding germline alignments and observed mutations \n")}
        time_and_log({
          #run first createGermlines() to add proper germline_d_mask collumn
          reference_dir <- paste0(imgt_dir, organism, "/vdj/")
          reference <- dowser::readIMGT(reference_dir)
          VDJ_db <- dowser::createGermlines(VDJ_db,
                                            references = reference,
                                            trim_lengths = TRUE,
                                            force_trim = TRUE,
                                            nproc = nproc,
                                            amino_acid = FALSE,
                                            id = "sequence_id",
                                            clone = "sequence_id",
                                            na.rm = FALSE)
          
          # then run shazam observedMutations()
          suppressMessages(library(shazam)) #cannot make it work without loading the entire package as shazam calls multiple objects from the package
          VDJ_db <- observedMutations(VDJ_db, sequenceColumn= "sequence_alignment", germlineColumn="germline_alignment_d_mask", regionDefinition=IMGT_V, frequency=FALSE, nproc=nproc)
          VDJ_db <- observedMutations(VDJ_db, sequenceColumn= "sequence_alignment", germlineColumn="germline_alignment_d_mask", regionDefinition=IMGT_V, frequency=FALSE, combine=TRUE, nproc=nproc)
          VDJ_db <- observedMutations(VDJ_db, sequenceColumn= "sequence_alignment", germlineColumn="germline_alignment_d_mask", regionDefinition=IMGT_V, frequency=TRUE, nproc=nproc)
          VDJ_db <- observedMutations(VDJ_db, sequenceColumn= "sequence_alignment", germlineColumn="germline_alignment_d_mask", regionDefinition=IMGT_V, frequency=TRUE, combine=TRUE, nproc=nproc)
          
        }, verbose = FALSE, log_file = log_file, log_title = "adding germline alignments and observed mutations", open_mode = "a")
        outfilename <- paste0(outfilename, "_germ-pass_shm-pass")
        if(!full_seq_aa){#we only save this file if no other analysis is performed
          readr::write_tsv(VDJ_db, file = paste0(outfolder, "/", outfilename, ".tsv.gz"))
        }
      }
      # Step6: reconstruct full VDJ AA sequence:
      if(full_seq_aa){
        time_and_log({
          VDJ_db <- reconstructFullVDJ(VDJ_db)
        }, verbose = FALSE, log_file = log_file, log_title = "full VDJ reconstruction", open_mode = "a")
        outfilename <- paste0(outfilename, "_full-pass")
        readr::write_tsv(VDJ_db, file = paste0(outfolder, "/", outfilename, ".tsv.gz"))
      }

      # Step7: add additional attributes:
      if(file.exists(paste0(outfolder, "_info.xlsx"))){
        if(verbose){cat("importing additional info from ", paste0(outfolder, "_info.xlsx"), "\n","\n")}
        additional_info <- openxlsx::read.xlsx(paste0(outfolder, "_info.xlsx"), 1)

        AddInfo <- function(db){
          AIRR_collumns <- colnames(db)[!colnames(db) == "well_id"]
          info_collumns <- setdiff(colnames(additional_info), AIRR_collumns)
          redundant_columns <- intersect(colnames(additional_info), AIRR_collumns) #any column in the platename_info.xlsx file that share a name with existing output columns from igblast will not be taken into account except if a primer column is provided
          if ("primers" %in% colnames(additional_info)){
            db$primers <- NULL
            redundant_columns <- redundant_columns[!redundant_columns == "primers"]
          }
          if (length(redundant_columns)>0){
            warning("the following columns from ", outfolder, "_info.xlsx were renamed as info_'old_colname' to avoid duplicated colnames: ", paste(redundant_columns, collapse = ", "))
            db <- dplyr::rename_with(db, .fn = ~ paste0("info_", .), .cols = all_of(redundant_columns))
            redundant_columns <- paste0("info_", redundant_columns)
            info_collumns <- c(info_collumns, redundant_columns)
          }
          #additional_info.xlsx file should have a "well_id" column
          if("well_id" %in% info_collumns){
            db <- dplyr::left_join(db, additional_info, by=join_by(well_id))
            db <- db[,c(info_collumns, AIRR_collumns)]
          } else {message(paste0("missing well_id collumn in ", outfolder, "_info.xlsx"))}
          return(db)
        }
        VDJ_db <- AddInfo(VDJ_db)
        if(!any(is.null(QC_failed), nrow(QC_failed)==0)){QC_failed <- AddInfo(QC_failed)}
        if(!any(is.null(failed_VDJ_db), nrow(failed_VDJ_db)==0)){failed_VDJ_db <- AddInfo(failed_VDJ_db)}
      }

      #remove unproductive sequences:
      VDJ_db_nonprod <- VDJ_db %>%
        dplyr::filter(!productive)
      VDJ_db <- VDJ_db %>%
        dplyr::filter(productive)

      # define style
      if (requireNamespace("openxlsx", quietly = TRUE)) {
        length_flag_style <- openxlsx::createStyle(fgFill="bisque")
        QC_flag_style <- openxlsx::createStyle(fgFill="bisque", fontColour = "red")
        cols <- which(colnames(VDJ_db)%in%colnames(VDJ_db))
        
        if(!primers %in% c("IgG", "IgM", "IgK", "IgL", "Mix")){
          message("incorrect primers provided, defaulting to Mix")
          primers <- "Mix"
        }
        v_sequence_end_cutoff <- c(270, 270, 270, 270, 270, 270)
        names(v_sequence_end_cutoff) <- c("IgG", "IgM", "IgK", "IgL", "Mix")
        
        length_cutoff <- c(350, 300, 400, 300, 300, 300)
        names(v_sequence_end_cutoff) <- c("IgG", "IgM", "IgK", "IgL", "Mix")
        
        OUT <- openxlsx::createWorkbook()
        openxlsx::addWorksheet(OUT, "productives")
        openxlsx::writeData(OUT, sheet = "productives", x = VDJ_db, colNames = TRUE, rowNames = FALSE)
        length_issue <- which(VDJ_db$v_sequence_end<=v_sequence_end_cutoff[primers]|"missing too much bp in VH sequence" %in% VDJ_db$comments|VDJ_db$sequence_length < length_cutoff[primers])
        openxlsx::addStyle(OUT, sheet="productives", style=length_flag_style, rows=length_issue+1, cols=cols, gridExpand=TRUE) # "+1" for header line
        low_QC <- which(VDJ_db$pct_under_30QC_in_trimmed>=10)
        openxlsx::addStyle(OUT, sheet="productives", style=QC_flag_style, rows=low_QC+1, cols=cols, gridExpand=TRUE) # "+1" for header line
        nb_prod <- nrow(VDJ_db)
        #!second style added will override the first
        nb_unprod <- nrow(VDJ_db_nonprod)
        openxlsx::addWorksheet(OUT, "non_productives")
        openxlsx::writeData(OUT, sheet = "non_productives", x = VDJ_db_nonprod, colNames = TRUE, rowNames = FALSE)
        
        nb_igblast_failed <- nrow(failed_VDJ_db)
        openxlsx::addWorksheet(OUT, "igblast_failed")
        openxlsx::writeData(OUT, sheet = "igblast_failed", x = failed_VDJ_db, colNames = TRUE, rowNames = FALSE)
        
        nb_QC_failed <- nrow(QC_failed)
        openxlsx::addWorksheet(OUT, "QC_failed")
        openxlsx::writeData(OUT, sheet = "QC_failed", x = QC_failed, colNames = TRUE, rowNames = FALSE)
        
        openxlsx::saveWorkbook(OUT, file = paste0(outfolder, "/", stringr::str_replace(filename, "_SCF_SEQ_ABI", ""), "_full_recap.xlsx"), overwrite = TRUE)
        
      } else {
        message("Optional: 'openxlsx' not installed — simply saving as tsv files")
        readr::write_tsv(VDJ_db, file = paste0(outfolder, "/", stringr::str_replace(filename, "_SCF_SEQ_ABI", ""), "_full_recap-prod.tsv"))
        readr::write_tsv(VDJ_db_nonprod, file = paste0(outfolder, "/", stringr::str_replace(filename, "_SCF_SEQ_ABI", ""), "_full_recap-non-prod.tsv"))
      }
      
      time_and_log({
        cat("nb total sequences: ", nb_prod+nb_unprod+nb_igblast_failed+nb_QC_failed, "\n")
        cat(paste0("nb failing initial QC: ",nb_QC_failed), "\n")
        cat(paste0("nb failing igblast: ",nb_igblast_failed), "\n")
        cat(paste0("nb non productive: ",nb_unprod), "\n")
        cat(paste0("nb productive: ",nb_prod), "\n")
        print(table(VDJ_db$c_call, useNA= "ifany"))
        cat("final output folder: ", outfolder, "\n")
      }, verbose = FALSE, time = FALSE, log_file = log_file, log_title = "Final recap", open_mode = "a")
      
      if(verbose & !update_info){
        cat("nb total sequences: ", nb_prod+nb_unprod+nb_igblast_failed+nb_QC_failed, "\n")
        cat(paste0("nb failing initial QC: ",nb_QC_failed), "\n")
        cat(paste0("nb failing igblast: ",nb_igblast_failed), "\n")
        cat(paste0("nb non productive: ",nb_unprod), "\n")
        cat(paste0("nb productive: ",nb_prod), "\n")
        print(table(VDJ_db$c_call, useNA= "ifany"))
        cat("final output folder: ", outfolder, "\n")
      }

      end <- Sys.time()
      step <- paste0("elapsed:", sprintf("%.2f %s", end-start, units(difftime(end, start))))
      if(verbose){cat(step, "\n","\n")}

    } else {
      if(!any(is.null(QC_failed), nrow(QC_failed)==0)){
        if(verbose){
          cat("nb total sequences: ", nrow(QC_failed), "\n")
          cat(paste0("nb failing initial QC: ",nb_QC_failed), "\n")
          cat("final output folder: ", outfolder, "\n")
        }
        time_and_log({
          cat("nb total sequences: ", nrow(QC_failed), "\n")
          cat(paste0("nb failing initial QC: ",nb_QC_failed), "\n")
          cat("final output folder: ", outfolder, "\n")
        }, verbose = FALSE, time = FALSE, log_file = log_file, log_title = "Final recap", open_mode = "a")
      }
      end <- Sys.time()
      step <- paste0("elapsed:", sprintf("%.2f %s", end-start, units(difftime(end, start))))
      if(verbose){cat(step, "\n","\n")}
    }
    return(VDJ_db)
  })
  names(airr.list) <- files
  if(return_db){
    return(airr.list)
  }
}

#### Function to QC plate-based .Ab1 Eurofins sequencing file ####
#' unzip Eurofins-based results and run all steps from QC to final AIRR formatted datatable
#'
#' \code{runAb1QC}
#'
#' @param Ab1_folder  path to a folder containing .Ab1 files
#' @param outfolder   folder to output files to
#' @param save        whether to save png or html plots
#' @param reticulate_py_env if you want to force a precise python environment to use for plotly. [should be "~/Library/r-miniconda-arm64/envs/r-reticulate/bin/python" if following installation instructions on MACOS 10.15.4]
#' @param primers     primers used for PCR
#' @param trim_cutoff cutoff thresholds to use depending on primers
#' @param nproc   number of processors to use, if NULL, will be automatically set based on available processors.
#'
#' @keywords internal
#'
#' @import sangeranalyseR
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings subseq
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings replaceAt

runAb1QC <- function(Ab1_folder,
                     outfolder = NULL,
                     outfilename = NULL,
                     save = c("png", "html"),
                     reticulate_py_env = NULL,
                     primers = c("IgG", "IgM", "IgL", "IgK", "Mix"),
                     trim_cutoff = list("IgG" = 300, "IgM" = 275, "IgK" = 400, "IgL"= 275, "Mix" = 275, "Undefined" = 275),
                     nproc = NULL){

  suppressMessages(library(sangeranalyseR))
  #suppressMessages(library(ggplot2))
  #suppressMessages(library(ggrepel))

  if(!dir.exists(Ab1_folder)){
    stop("folder : ", Ab1_folder," does not exist")
  }

  if(is.null(outfolder)){
    if(verbose){cat("output folder not defined; exporting to the provided data folder: ", Ab1_folder)}
    outfolder <- Ab1_folder
  }
  if(is.null(outfilename)){
    if(verbose){cat("output filename not defined; exporting to the provided data folder: ", Ab1_folder)}
    outfilename <- Ab1_folder
  }

  if(!dir.exists(outfolder)){
    dir.create(outfolder)
  }

  QC_folder <- paste0(outfolder, "/QC_files")
  if(!dir.exists(QC_folder)){
    dir.create(QC_folder)
  }

  save <- match.arg(save)
  primers <- match.arg(primers)

  if(!primers %in% c("IgG", "IgM", "IgK", "IgL", "Mix")){
    message("incorrect primers provided, defaulting to minimal parameters for trimming cutoff: 275bp")
    primers <- "Undefined"
  }

  if(is.null(nproc)){
    if (requireNamespace("parallel", quietly = TRUE)) {
      cores <- parallel::detectCores()
      nproc <- max((cores[1]-2), 1) #not to overload computer
    } else {
      message("Optional: 'parallel' not installed — defaulting to nproc = 1")
      nproc <- 1
    }
  }

  filenames <- list.files(Ab1_folder, pattern="*.ab1", full.names=TRUE)

  well_ids <- unlist(lapply(filenames, function(filename){
    well_id <- stringr::str_sub(filename,-7,-5)
    #based on Eurofins numbering scheme, well_id is always the last 3 characters before ".ab1".
    return(well_id)
  }))

  ##open Ab1 files using sangeranalyseR to trimm sequence:
  #M1 Trimming Method = Mott's modified trimming algorithm
  #M1 Trimming Cutoff: the cutoff at which you consider a base to be bad. This works on a logarithmic scale, such that if you want to consider a score of 10 as bad, you set cutoff to 0.1; for 20 set it at 0.01; for 30 set it at 0.001; for 40 set it at 0.0001; and so on. Contiguous runs of bases below this quality will be removed from the start and end of the sequence. Given the high quality reads expected of most modern ABI sequencers, the defualt is 0.0001.
  message("Importing and trimming all sequences using sangeranalyseR")
  seqs <- safe_mclapply(filenames, FUN=function(filename){
    capture.output({
      suppressMessages({
        seq <- sangeranalyseR::SangerRead(readFeature = "Reverse Read",
                                          readFileName          = filename,
                                          geneticCode           = GENETIC_CODE,
                                          TrimmingMethod        = "M1",
                                          M1TrimmingCutoff      = 0.001,
                                          baseNumPerRow         = 100,
                                          heightPerRow          = 200,
                                          signalRatioCutoff     = 0.33,
                                          showTrimmed           = TRUE)
      })
    }, type = "output")
    #generateReport(seq, outputDir = Ab1_folder)
    return(seq)
  }, mc.cores = nproc)
  names(seqs) <- well_ids

  ##open Ab1 files using sanqerseqR to plot chromatograms:
  if (requireNamespace("sangerseqR", quietly = TRUE)) {
    suppressMessages(library(sangerseqR))
    message("Plotting chromatograms using sangerseqR")
    seqs2 <- safe_mclapply(filenames, sangerseqR::readsangerseq, mc.cores = nproc)
    names(seqs2) <- well_ids
    
    if(save == "html"){
      for(j in seq_along(seqs2)){
        #message("outputting QC graphs for ", well_ids[j])
        seq2final <- sangerseqR::makeBaseCalls(seqs2[[j]], ratio = 0.33)
        QC_file <- paste0(QC_folder, "/", gsub(paste0(Ab1_folder, "/"), "", gsub(".ab1", "", filenames[j])))
        pdf(paste0(QC_file, "_chromatogram.pdf"))
        sangerseqR::chromatogram(seq2final, width = 80, height = 3, trim5 = 0, trim3 = 0,showcalls = "both")
        dev.off()
        htmlwidgets::saveWidget(sangeranalyseR::qualityBasePlot(seqs[[j]]), file = paste0(QC_file, "_QC_plot.html"))
      }
    }
    if(save == "png"){
      #save_image() from the plotly package can also be used to save only one image but much slower:
      #in all cases, the kaleido package should be installed in python as follow and run through the reticulate R package
      #install.packages('reticulate')
      #reticulate::install_miniconda()
      #reticulate::conda_install('r-reticulate', 'python-kaleido')
      #reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
      #reticulate::use_miniconda('r-reticulate')
      #reticulate::py_config()
      #!!when updating RStudio, reticulate python version might need to be reset to r-reticulate...
      #or add: RETICULATE_PYTHON="/Users/yourname/Library/r-miniconda-arm64/envs/r-reticulate/bin/python" to .Renviron
      
      #alternatively you can also force reticulate to use the right python environment by using the reticulate_py_env argument:
      if(!is.null(reticulate_py_env)){
        Sys.setenv(RETICULATE_PYTHON = reticulate_py_env)
        reticulate::use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
      }
         
      print(reticulate::py_config())
      
      scope <- safe_kaleido()
      
      if (!is.null(scope) && is.function(scope$transform)) {
        for (j in seq_along(seqs2)) {
          #message("outputting QC graphs for ", well_ids[j])
          QC_file <- paste0(QC_folder, "/", gsub(paste0(Ab1_folder, "/"), "", gsub(".ab1", "", filenames[j])))
          
          p <- sangeranalyseR::qualityBasePlot(seqs[[j]])
          scope$transform(p, paste0(QC_file, "_QC_plot.png"))
          
          seq2final <- sangerseqR::makeBaseCalls(seqs2[[j]], ratio = 0.33)
          sangerseqR::chromatogram(seq2final, width = 80, height = 3, trim5 = 0, trim3 = 0, showcalls = "both", filename = paste0(QC_file, "_chromatogram.pdf"))
        }
      } else {
        message("Skipping png QC plots step: 'plotly::kaleido()' failed or is not usable.")
      }
      
      rm(scope) 
      gc()
    }
  } else {
    message("Optional: 'sangerseqR' not installed — skipping QC plots")
  }
  
  primary_seqs <- lapply(seqs, FUN=function(seq){
    return(Biostrings::reverseComplement(Biostrings::subseq(seq@primarySeq, start = seq@QualityReport@trimmedStartPos, end = seq@QualityReport@trimmedFinishPos)))
  })
  names(primary_seqs) <- well_ids
  primary_seqs <- Biostrings::DNAStringSet(primary_seqs)
  #Biostrings::writeXStringSet(primary_seqs, paste0(out_folder, "/", out_folder, "_trimmed_sequences.fasta"))

  secondary_seqs <- lapply(seqs, FUN=function(seq){
    return(Biostrings::reverseComplement(Biostrings::subseq(seq@secondarySeq, start = seq@QualityReport@trimmedStartPos, end = seq@QualityReport@trimmedFinishPos)))
  })
  names(secondary_seqs) <- well_ids
  secondary_seqs <- Biostrings::DNAStringSet(secondary_seqs)

  nb_seqs <- length(primary_seqs)
  seq_df <- data.frame(well_id=character(nb_seqs),
                       sequence=character(nb_seqs),
                       sequence_length=numeric(nb_seqs),
                       pct_under_30QC_in_trimmed=numeric(nb_seqs),
                       low10QC_original_seq=character(nb_seqs),
                       low30QC_original_seq=character(nb_seqs),
                       low10QC_alternate_calls=character(nb_seqs),
                       low30QC_alternate_calls=character(nb_seqs),
                       primary_sequence=character(nb_seqs),
                       raw_sequence=character(nb_seqs),
                       stringsAsFactors=FALSE)

  for(i in 1:nb_seqs){
    #primary seq
    seq_df$well_id[i] <- names(primary_seqs)[[i]]
    sequence <- primary_seqs[[i]]
    seq_df$sequence[i] <- as.character(sequence)
    seq_df$primary_sequence[i] <- as.character(sequence)
    seq_df$sequence_length[i] <- primary_seqs[[i]]@length
    #bases with low Phred Quality Score
    phred <- seqs[[i]]@QualityReport@qualityPhredScores
    names(phred) <- c(1:length(phred))
    low_quality_bases_30 <- as.numeric(names(phred[phred<30]))
    low_quality_bases_30 <- low_quality_bases_30[low_quality_bases_30>seqs[[i]]@QualityReport@trimmedStartPos & low_quality_bases_30<seqs[[i]]@QualityReport@trimmedFinishPos]
    seq_df$pct_under_30QC_in_trimmed[i] <- length(low_quality_bases_30)/primary_seqs[[i]]@length*100

    low_quality_bases_10 <- as.numeric(names(phred[phred<10]))
    low_quality_bases_10 <- low_quality_bases_10[low_quality_bases_10>seqs[[i]]@QualityReport@trimmedStartPos & low_quality_bases_10<seqs[[i]]@QualityReport@trimmedFinishPos]
    seq_df$low10QC_original_seq[i] <- ifelse(length(low_quality_bases_10)>0, stringr::str_c(low_quality_bases_30, collapse=';') , NA)
    if(length(low_quality_bases_10>0)){
      #update position based to reverse complemented sequence
      reverted_low_quality_bases_10 <- seqs[[i]]@QualityReport@trimmedFinishPos - low_quality_bases_10 + 1

      #extract differential base calls
      DiffBaseCalls_10 <- lapply(reverted_low_quality_bases_10, function(position){
        primary_call <- as.character(Biostrings::subseq(primary_seqs[[i]], start = position, end = position))
        secondary_call <- as.character(Biostrings::subseq(secondary_seqs[[i]], start = position, end = position))
        if(primary_call != secondary_call){
          return(paste0(primary_call, position, secondary_call))
        } else {return(NA)}
      })
      DiffBaseCalls_10_positions <- lapply(reverted_low_quality_bases_10, function(position){
        primary_call <- as.character(Biostrings::subseq(primary_seqs[[i]], start = position, end = position))
        secondary_call <- as.character(Biostrings::subseq(secondary_seqs[[i]], start = position, end = position))
        if(primary_call != secondary_call){
          return(position)
        } else {return(NA)}
      })

      DiffBaseCalls_10 <- na.omit(unlist(DiffBaseCalls_10))
      DiffBaseCalls_10_positions <- na.omit(unlist(DiffBaseCalls_10_positions))
      if(length(DiffBaseCalls_10)>0){
        #store positions and alternate calls
        seq_df$low10QC_alternate_calls[i] <- stringr::str_c(rev(DiffBaseCalls_10), collapse=';')
        #replace base in sequence with an N
        for(pos in seq_along(DiffBaseCalls_10_positions)){
          position <- DiffBaseCalls_10_positions[pos]
          at <- IRanges(position:position, width=1)
          sequence <- Biostrings::replaceAt(sequence, at, "N")
        }
        seq_df$sequence[i] <- as.character(sequence) #update sequence
      } else {seq_df$low10QC_alternate_calls[i] <- NA}
    } else {seq_df$low10QC_alternate_calls[i] <- NA}

    intermediate_quality_bases <- setdiff(low_quality_bases_30, low_quality_bases_10)
    seq_df$low30QC_original_seq[i] <- ifelse(length(intermediate_quality_bases)>0, stringr::str_c(intermediate_quality_bases, collapse=';') , NA)
    if(length(intermediate_quality_bases>0)){
      #update position based to reverse complemented sequence
      reverted_intermediate_quality_bases <- seqs[[i]]@QualityReport@trimmedFinishPos - intermediate_quality_bases + 1
      #extract differential base calls
      DiffBaseCalls_30 <- na.omit(unlist(lapply(reverted_intermediate_quality_bases, function(position){
        primary_call <- as.character(Biostrings::subseq(primary_seqs[[i]], start = position, end = position))
        secondary_call <- as.character(Biostrings::subseq(secondary_seqs[[i]], start = position, end = position))
        if(primary_call != secondary_call){
          return(paste0(primary_call, position, secondary_call))
        } else {return(NA)}
      })))
      if(length(DiffBaseCalls_30)>0){
        seq_df$low30QC_alternate_calls[i] <- stringr::str_c(rev(DiffBaseCalls_30), collapse=';')
      } else {seq_df$low30QC_alternate_calls[i] <- NA}
    } else {seq_df$low30QC_alternate_calls[i] <- NA}
    seq_df$raw_sequence[i]<-as.character(seqs[[i]]@primarySeq)
  }

  seq_df$primers <- primers
  seq_df$QC_passed <- seq_df$sequence_length>trim_cutoff[[primers]] & seq_df$pct_under_30QC_in_trimmed<20
  options(warn=-1)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot(seq_df, aes(x=sequence_length, y=pct_under_30QC_in_trimmed, color=QC_passed, label = well_id)) +
      geom_point(size=3) +
      geom_hline(yintercept=20) +
      geom_vline(xintercept=trim_cutoff[[primers]]) +
      ggtitle(paste0("QC plot: ", outfolder)) +
      theme_classic()
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + 
        ggrepel::geom_text_repel(
          # Repel away from the left edge, not from the right.
          xlim = c(NA, Inf),
          # Do not repel from top or bottom edges.
          ylim = c(-Inf, Inf)
        )
    } else {message("Optional: 'ggrepel' not installed — skipping text part on QC plot")}
    ggsave(p, filename=paste0(outfolder, "/QCplot_", outfilename, ".pdf"))
  } else {message("Optional: 'ggplot2' not installed — skipping QC plot")}
  options(warn=0)

  filtered_seq_df <- seq_df[seq_df$QC_passed,]
  if(nrow(filtered_seq_df)>0){
    readr::write_tsv(filtered_seq_df, file = paste0(outfolder, "/", outfilename, "_QC-pass.tsv"))
  }

  failed_seq_df <- seq_df[!seq_df$QC_passed,]
  if(nrow(failed_seq_df)>0){
    readr::write_tsv(failed_seq_df, file = paste0(outfolder, "/", outfilename, "_QC-fail.tsv"))
  }
  return(list("pass" = filtered_seq_df, "fail"= failed_seq_df))
}

#### Function to reformat input data from 10x or BD ####
#' rename columns to harmonize with the immcantation pipeline and between technologies, also removes a few useless ones.
#'
#' \code{reformatVDJinput}
#'
#' @param db    a data frame directly imported from SevenBridges (BD) or CellRanger (10X).
#' @param tech  name of the scRNA-seq technology used (one of 10X or BD)
#' @param cell_id.  name of the column containing cell_ids
#' @param locus     name of the column containing locus informations
#' @param heavy     which vlue to use in locus to identify heavy chains [default: IGH for BCR]
#' @param productive  name of the column containing contig productiveness informations
#' @param complete_vdj  name of the column containing contig completness informations
#' @param sequence_id   name of the column containing sequence_ids
#' @param umi_count name of the column containing umi_count informations
#' @param consensus_count name of the column containing consensus_count informations
#' @param junction  name of the column containing identified junctions (nucleotide format)
#' @param junction_aa name of the column containing identified junctions (amino-acid format)
#' @param sequence  name of the column containing full sequences
#' @param v_call  name of the column containing v_call informations
#' @param d_call  name of the column containing d_call informations
#' @param j_call  name of the column containing j_call informations
#' @param c_call  name of the column containing c_call informations
#' @param remove_columns  name of collumns to remove from original 10X or BD data frames
#'
#' @return a reformatted dataframe with same name for both 10X or BD technologies. Re-running igblast is however advised to get a full AIRR formatted dataframe.
#'
#' @details
#' the following columns are unique to the SevenBridges pipeline from BD:
#' "cell_type_experimental"; "high_quality_cell"; "complete_vdj"; "dominant"; "putative_cell"
#' "sequence_length"; "sequence_aa"; "sequence_aa_length" ; "sequence_alignment_length"
#' "sequence_alignment_aa"; "sequence_alignment_aa_length"; "cdr3_length"; "fwr1_aa";"fwr2_aa"; "fwr3_aa"; "fwr4_aa"; "cdr1_aa"; "cdr2_aa"; "cdr3_aa"
#' "germline_alignment_aa"; "v_germline_alignment"; "v_germline_alignment_aa"; "d_germline_alignment"; "d_germline_alignment_aa"; "j_germline_alignment";"j_germline_alignment_aa"
#' some will be added or modified at the CreateGermline() step later after clone inference (current germline alignment is no gapped so not usable in observedMutations)
#' https://bd-rhapsody-bioinfo-docs.genomics.bd.com/outputs/outputs_vdj_contigs.htm
#'
#' the following columns are named differently or unique to the 10x pipeline from BD:
#' barcode, is_cell,	contig_id,	high_confidence,	length,	chain,	v_gene,	d_gene,	j_gene,	c_gene,
#' full_length,	fwr1,	fwr1_nt,	cdr1,	cdr1_nt,	fwr2,	fwr2_nt,	cdr2,	cdr2_nt,	fwr3,	fwr3_nt,	cdr3,	cdr3_nt,	fwr4,	fwr4_nt,
#' reads,	umis,	raw_clonotype_id,	raw_consensus_id,	exact_subclonotype_id,
#'
#' @import dplyr
#'
#' @export

reformatVDJinput <- function(db, tech,
                             cell_id = "cell_id", locus = "locus", heavy = "IGH", productive = "productive", complete_vdj = "complete_vdj",
                             sequence_id = "sequence_id", umi_count = "umi_count", consensus_count = "consensus_count",
                             junction = "junction", junction_aa = "junction_aa", sequence = "sequence",
                             v_call = "v_call", d_call = "d_call", j_call = "j_call", c_call = "c_call",
                             remove_columns = list("10X" = c("high_confidence", "raw_clonotype_id",	"raw_consensus_id",	"exact_subclonotype_id",
                                                             "fwr1_aa", "fwr2_aa", "fwr3_aa", "fwr4_aa", "cdr1_aa", "cdr2_aa", "cdr3_aa"),
                                                   "BD" = c("cell_type_experimental", "high_quality_cell", "high_quality_cell_tcr_bcr", "sequence_aa", "sequence_aa_length" , "sequence_alignment_length",
                                                            "cdr3_length", "fwr1_aa", "fwr2_aa", "fwr3_aa", "fwr4_aa", "cdr1_aa", "cdr2_aa", "cdr3_aa",
                                                            "germline_alignment_aa", "v_germline_alignment", "v_germline_alignment_aa",
                                                            "d_germline_alignment", "d_germline_alignment_aa", "j_germline_alignment", "j_germline_alignment_aa"))){

  suppressMessages(library(dplyr))

  if(tech == "10X"){
    #rename a few column in the 10X contig_annotations file:
    db <- db %>%
      dplyr::rename(!!rlang::sym(v_call) := v_gene,
                    !!rlang::sym(d_call) := d_gene,
                    !!rlang::sym(j_call) := j_gene,
                    !!rlang::sym(c_call) := c_gene,
                    !!rlang::sym(cell_id) := barcode,
                    !!rlang::sym(sequence_id) := contig_id,
                    putative_cell = is_cell,
                    !!rlang::sym(locus) := chain,
                    !!rlang::sym(complete_vdj) := full_length,
                    sequence_length = length,
                    !!rlang::sym(consensus_count) := reads,
                    !!rlang::sym(umi_count) := umis,
                    fwr1 = fwr1_nt, fwr1_aa = fwr1,
                    fwr2 = fwr2_nt, fwr2_aa = fwr2,
                    fwr3 = fwr3_nt, fwr3_aa = fwr3,
                    fwr4 = fwr4_nt, fwr4_aa = fwr4,
                    cdr1 = cdr1_nt, cdr1_aa = cdr1,
                    cdr2 = cdr2_nt, cdr2_aa = cdr2,
                    cdr3 = cdr3_nt, cdr3_aa = cdr3)

    db[[junction]] <- db$cdr3
    db[[junction_aa]] <- db$cdr3_aa
  }
  for(i in c(v_call, d_call, j_call, c_call, junction, junction_aa)){
    db[[paste0(i, "_", tech)]] <- db[[i]]
  }
  keep_columns <- setdiff(colnames(db), remove_columns[[tech]])
  db <- db %>%
    dplyr::select(all_of(keep_columns))

  return(db)
}


#### Function to assign seq to bin based on umi_count and a multimodel fitting ####
#' Resolve VDJ chains multiplets, eventually using prior clustering knowledge
#'
#' \code{resolveMultiContigs} resolve cases of multiple heavy or light chains within a cell
#'
#' @param data vector of values.
#' @param plot whether to export plot
#' @param title title of plot
#' @param file filename for plot
#'
#' @return   a list containing three slots:
#' "bins": the allocated bin for each initial value
#' "plot": the histogram + density plot with antimode-based bins highlighted (Kernel Density Estimation)
#' "top75%": the top bins corresponding to >75% of all values.
#'
#' @import dplyr
#'
#' @keywords internal

binContigs <- function(data,
                       plot = TRUE,
                       title = "Histogram + Density Plot with Antimode-Based Bins",
                       file = "DensityPlot.pdf"){

  #suppressMessages(library(dplyr))

  # Compute Kernel Density Estimation (KDE)
  dens <- density(data)
  # Function to find local minima (antimodes)
  which.minima <- function(x) {
    which(diff(sign(diff(x))) == 2) + 1
  }
  # Find antimodes (local minima of density), bin edges and assign each data point to a bin
  antimodes <- dens$x[which.minima(dens$y)]
  bin_edges <- c(min(data), antimodes, max(data))
  bins <- cut(data, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)

  # Create a frequency table
  bin_relative_freq <- (table(bins) / length(data)) * 100
  bin_freq_df <- data.frame(bin = as.numeric(names(bin_relative_freq)),
                            relative_frequency = as.numeric(bin_relative_freq))
  # Select bins covering at least 75% of the highest values
  bin_freq_df <- bin_freq_df %>% dplyr::arrange(desc(bin))
  bin_freq_df$cumulative_frequency <- cumsum(bin_freq_df$relative_frequency)
  selected_bins <- bin_freq_df[bin_freq_df$cumulative_frequency <= 75, ]
  # Ensure we include the first bin that crosses 75%
  if (nrow(selected_bins) > 0){
    if (max(selected_bins$cumulative_frequency) < 75){
      next_bin <- bin_freq_df[nrow(selected_bins) + 1, ]
      selected_bins <- rbind(selected_bins, next_bin)
    }
  } else { selected_bins <- bin_freq_df[1, ] }

  if(plot){
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      suppressMessages(library(ggplot2))
      # Convert density estimation and antimodes to a data frame for ggplot
      density_df <- data.frame(x = dens$x, y = dens$y)
      antimodes_df <- data.frame(x = antimodes, y = dens$y[which.minima(dens$y)])
      
      p  <- ggplot2::ggplot() +
        ggplot2::geom_histogram(data = data.frame(x = data), aes(x = x, y = after_stat(density)),
                                bins = 30, fill = "gray", alpha = 0.5, color = "black") +
        ggplot2::geom_line(data = density_df, aes(x = x, y = y), color = "blue", linewidth = 1) +
        ggplot2::geom_point(data = antimodes_df, aes(x = x, y = y), color = "red", size = 3) +
        ggplot2::geom_vline(xintercept = bin_edges, linetype = "dashed", color = "gray") +
        ggplot2::geom_vline(xintercept = bin_edges[selected_bins$Bin[length(selected_bins$Bin)]], linetype = "solid", color = "green", linewidth = 1) +
        ggplot2::labs(title = title,
                      x = "Value", y = "Density") +
        ggplot2::theme_minimal()
      
      if(grepl(".pdf", file)|grepl(".png", file)){
        ggplot2::ggsave(file, plot = p, width = 8, height = 5, dpi = 300)
      }
      return(list("bins" = bins, "plot" = p, "top75%" = selected_bins))
    } else {
      message("Optional: 'ggplot2' not installed — skipping plot.")
      return(list("bins" = bins, "top75%" = selected_bins))
    }
  }
  return(list("bins" = bins, "top75%" = selected_bins))
}

#### Function to select the best heavy or light chain for each cell_id ####
#' Resolve VDJ chains multiplets, eventually using prior clustering knowledge
#'
#' \code{resolveMultiContigs} resolve cases of multiple heavy or light chains within a cell
#'
#' @param db            a dataframe containing heavy and light chain sequences, cell_ids, umi_counts, consensus counts (reads), productive and complete_vdj (or frw/cdr infos) columns.
#' @param cell_id       name of the column containing cell identifier.
#' @param seq_type      type of VDJ sequence ("Ig" or "TCR" to match igblastb requirements)
#' @param resolve_chain chain(s) to filter, should be one of "heavy" or "light", will return one chain per cell for the selected chain(s).
#' @param resolve_multi_CDR3 whether to pool contigs with similar v_call, j_call and junction_aa in a given cell (cell_id)
#' @param use_clone     whether to use prior clustering knowledge, a clone_id column should be present.
#' @param clone_id      name of the column containing cell identifier.#'
#' @param split.by      name of the column in the dataframe to use to split the dataset prior to filtering heavy chain
#' @param output        whether to output graphs with umi_counts for dominant versus second IGH VDJ contig and the recap excel workbook. If set to FALSE, only the corrected database is returned.
#' @param output_folder name of the folder in which graph for light chain clustering will be saved.
#' @param analysis_name name to use for outputs prefixes.
#' @param consensus_count   name of the column containing the number of reads for this contig (usually called consensus_count)
#' @param umi_count     name of the column containing the number of unique molecules (UMI) for this contig. Previously called "duplicate_count" in an earlier AIRR standard
#' @param productive    name of the column containing the info whether a given sequence is productive.
#' @param sequence_id   name of the column containing sequence identifier.
#' @param junction_aa   name of the column containing identified junction in amino-acid format.
#' @param v_call        name of the column containing V-segment allele assignments. All
#'                      entries in this column should be identical to the gene level.
#' @param j_call        name of the column containing J-segment allele assignments. All
#'                      entries in this column should be identical to the gene level.
#' @param nproc         number of processor(s) to use (passed to scoper::hierarchicalClones()).
#' @param second_columns which columns info to keep for secondary contigs of interest
#' @param locus         name of the column containing locus informations
#' @param c_call        name of the column containing c_call informations
#' @param junction      name of the column containing junction informations (nucleotide format)
#' @param dominant      name of the column containing BD SevenBridges contig dominance informations
#' @param complete_vdj  name of the column containing contig completness informations
#'
#' @return   a data frame containing the filtered chain for all cells as well as all initial columns and rows not related to that chain (or pair of chains) in the provided dataframe.
#' the following columns are added: unfiltered_bin, filtered_bin (outputs of binContigs pre and post filtering), pct_of_reads (% of total read in a given cell for each contig).
#' If use_clone = TRUE and expanded heavy chain clones are found, will run scoper::hierarchicalClones() and return a light_group column corresponding to light chain clustering as well as a graph for minimum distance between and maximal distance inside light clones.
#' All columns listed in the second_columns argument are also added for cells with two or more dominant light chain detected with a "second_" prefix.
#'
#' @details
#' it is recommended to perform such filtering at the level of each independent sample (split.by argument) as it will attempt to infer multimodal bins based on sample-intrinsic counts distribution.
#' requires the following columns in the provided database: "cell_id", "locus", "umi_count", "productive", "junction", "junction_aa", "consensus_count", "sequence", "v_call", "j_call", "c_call".
#' 1. will first pool contigs with similar CDR3, remove non-productive (if at least one productive if present in the cell) and filter for contigs representing at least 20% of all reads in a given cell;
#' 2. if use_clone = TRUE, will first cluster IGK/IGL contigs from expanded clones using scoper::HierarchicalClones() and a threshold of 0.1 and cases of cells with multiple light chain calls will then be resolved within each heavy chain-based clones (provided clone_id), using the following rules:
#'    2.1 clonally related light chain contigs within one cell are pulled and the best quality one (highest umi_count and completeness) is kept;
#'    2.2 clonally related light chain contigs between cells in a heavy chain-based clone are selected in priority for all of these cells;
#'    2.3 if no clonally related light chain is found, noncontig is selected at this step.
#' 3. to break remaining ties, the contigs are then sorted in order of: productivity, umi_count bin, full length, highest molecule count, highest read count.
#' dependencies: dplyr, ggplot2, openxlsx, alakazam, dowser
#'
#' [note] technically, dowser::resolveLightChains can deal with multiple light chains but choice of dominant light chain doesn't take into account umi_counts (selection seem purely based on alphabetical ordering of light chains...).
#' we use her a similar approach to preprocess the best light chain(s) for a given cell_id. The preprocess db can then be runned through dowser::resolveLightChains for splitting of heavy chain only-based clones using light chain info.
#'
#' clonal grouping should be done beforehand (scoper::HierarchicalClones()or SpectralClones()), without splitting by light chain (only_heavy = TRUE and split_light = FALSE).
#' scoper::HierarchicalClones() also provides light chain selection but unlike dowser::resolveLightChains(), it doesn't work with cells missing a light chain.
#'
#' @import dplyr
#' @import scoper
#'
#' @export

resolveMultiContigs <- function(db,
                                split.by = NULL,
                                seq_type = c("Ig", "TCR"),
                                resolve_chain = c("heavy", "light"),
                                assay = "assay",
                                resolve_multi_CDR3 = TRUE,
                                use_clone = FALSE,
                                analysis_name = "All_sequences",
                                output = TRUE,
                                output_folder = "VDJ_QC/",
                                second_columns = c("sequence_id", "locus", "umi_count", "consensus_count", "sequence", "v_call", "d_call", "j_call", "c_call", "junction", "junction_aa", "productive", "complete_vdj"),
                                cell_id = "cell_id", 
                                sequence_id = "sequence_id", 
                                locus = "locus", 
                                consensus_count = "consensus_count", 
                                umi_count = "umi_count", 
                                v_call = "v_call", 
                                j_call = "j_call", 
                                c_call = "c_call", 
                                junction = "junction", 
                                junction_aa = "junction_aa",
                                dominant = "dominant", 
                                productive = "productive", 
                                complete_vdj = "complete_vdj", 
                                clone_id = "clone_id", 
                                nproc = 1){

  #suppressMessages(library(dplyr))

  if(!is.null(output_folder)){
    if(!stringr::str_ends(output_folder, "/")){output_folder = paste0(output_folder, "/")}
    if(!dir.exists(output_folder)){
      dir.create(output_folder)
    }
  }
  
  seq_type <- match.arg(seq_type)
  resolve_chain <- match.arg(resolve_chain)
  
  if(seq_type == "Ig"){
    if(resolve_chain == "heavy"){
      chain <- "IGH"
    }
    if(resolve_chain == "light"){
      chain <- c("IGL", "IGK")
    }
  }
  if(seq_type == "TCR"){
    if(resolve_chain == "heavy"){
      chain <- c("TRB", "TRD")
    }
    if(resolve_chain == "light"){
      chain <- c("TRA", "TRG")
    }
  }
  
  #log_file <- paste0(output_folder, Sys.time(), "_", analysis_name, "_", paste(chain, collapse = "-"), "_filtration_logfile.txt")
  #log_connection <- file(log_file, open = "a")  # a for apending or w for erasing onto previous log

  required_columns <- c(cell_id, second_columns)
  if(use_clone){required_columns <- c(required_columns, clone_id)}
  if(!complete_vdj %in% colnames(db)){
    required_columns <- required_columns[required_columns != complete_vdj]
    required_columns <- c(required_collumns, "fwr1", "fwr2", "fwr3", "fwr4", "cdr1", "cdr2", "cdr3")
  }
  missing_columns <- setdiff(required_columns, colnames(db))
  if(length(missing_columns)>0) {stop(paste0("missing the following collumns: ", missing_columns))}

  if(!complete_vdj %in% colnames(db)){ # should have been calculated upon HC filtering (filterHC())
    db[[complete_vdj]] <- !(nchar(db$fwr1)<3|is.na(db$fwr1)|nchar(db$cdr1)<3|is.na(db$cdr1)|nchar(db$fwr2)<3|is.na(db$fwr2)|nchar(db$cdr2)<3|is.na(db$cdr2)|nchar(db$fwr3)<3|is.na(db$fwr3)|nchar(db$cdr3)<3|is.na(db$cdr3)|is.na(db$fwr4)|nchar(db$fwr4)<3)
  }

  if(is.null(split.by)){
    split <- NULL
    db <- db %>%
      dplyr::mutate(
        split.by = "all_sequences"
      )
    split.by = "split.by"
  }

  db_to_filter <- db %>% dplyr::filter(!!rlang::sym(locus) %in% chain)
  db_not_to_filter <- db %>% dplyr::filter(!(!!rlang::sym(locus) %in% chain))

  if(nrow(db_to_filter)<1){stop(paste0("no ", paste(chain, sep = "-")," contig provided"))}

  if(any(duplicated(db_to_filter[[sequence_id]]))){stop("Sequence IDs in provided dataframe must be unique!")}

  if(use_clone & !all(chain %in% c("IGL", "IGK"))){warning("using prior clustering knowledge only developped for light chain filtration, won't be used for heavy chains")}

  ini_seq_nb_for_filtering <- nrow(db_to_filter)

  #first resolve cases of multiple contigs in a cell with identical CDR3s
  if(resolve_multi_CDR3){
    db_filtered <- db_to_filter %>%
      dplyr::group_by(!!rlang::sym(cell_id), !!rlang::sym(locus), !!rlang::sym(v_call), !!rlang::sym(j_call), !!rlang::sym(junction_aa)) %>%
      dplyr::mutate(
        total_umi_count = sum(!!rlang::sym(umi_count)),
        total_consensus_count = sum(!!rlang::sym(consensus_count)),
        c_call = paste(unique(na.omit(!!rlang::sym(c_call))), collapse = "/") #update c_call to reflect the multiple potential c_call found for the same CDR3.
      ) %>%
      dplyr::arrange(desc(!!rlang::sym(productive)), desc(!!rlang::sym(complete_vdj)), desc(!!rlang::sym(umi_count)), desc(!!rlang::sym(consensus_count))) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(
        !!rlang::sym(umi_count) := total_umi_count, #update counts
        !!rlang::sym(consensus_count) := total_consensus_count, #update total counts
      ) %>%
      dplyr::select(-c(total_umi_count, total_consensus_count)) %>%
      dplyr::ungroup()
  }

  # then remove non productive if more than one productive present
  db_filtered <- db_filtered %>%
    dplyr::group_by(!!rlang::sym(cell_id)) %>%
    # Count the number of productive sequences in each cell_id group
    dplyr::mutate(
      productive_seqs = sum(productive)
    ) %>%
    # Filter: Keep productive = FALSE only if there is at least one productive sequence in the group
    dplyr::filter(!(!!rlang::sym(productive) == FALSE & productive_seqs >= 1)) %>%
    dplyr::select(-productive_seqs) %>%
    dplyr::ungroup()

  # bin contigs
  db_filtered <- db_filtered %>%
    dplyr::group_by(!!rlang::sym(split.by)) %>%
    dplyr::mutate(
      unfiltered_bin = binContigs(log10(umi_count),
                                  plot = output,
                                  title = paste0(dplyr::first(!!rlang::sym(split.by)), " - ", paste(chain, collapse = "-"), " Density plot"),
                                  file = paste0(output_folder, analysis_name, "_", dplyr::first(!!rlang::sym(split.by)), "_", paste(chain, collapse = "-"),"_density_plot.pdf"))[["bins"]],
    ) %>%
    ungroup()

  # filter on frequency of reads per cell
  db_filtered <- db_filtered %>%
    dplyr::group_by(!!rlang::sym(cell_id)) %>%
    dplyr::mutate(
      pct_of_reads = consensus_count/sum(consensus_count)*100,
    ) %>%
    dplyr::filter(pct_of_reads>=20) %>%
    ungroup()
  
  #TODO add preferred selection of TRA if TRB is the dominant contig in a T cell? (same for TRG/TRD)

  # if use_clone = TRUE, look for shared seq inside clone (to be used for light chains but could work for heavy chain also - to be investigated):
  if(use_clone & resolve_chain == "light"){
    # we then check for cells with no clone_id:
    if(any(is.na(db_filtered[[clone_id]]))){
      no_clone_id <- dplyr::filter(db_filtered, is.na(!!rlang::sym(clone_id)))
      missing_clone_id <- length(unique(no_clone_id$cell_id))
      db_filtered <- dplyr::filter(db_filtered, !is.na(!!rlang::sym(clone_id)))
      no_clone_log_message <- paste0(missing_clone_id, " cells without clone_id, processed as non expanded clones.\n")
      cat(no_clone_log_message)
    } else {missing_clone_id <- 0}

    # we further separate singlets from expanded clones:
    db_filtered <- db_filtered %>%
      dplyr::group_by(!!rlang::sym(clone_id)) %>%
      dplyr::mutate(expanded_clone = n_distinct(!!rlang::sym(cell_id)) > 1) %>%
      dplyr::ungroup()

    exp_db <- dplyr::filter(db_filtered, expanded_clone)
    singlets_db <- dplyr::filter(db_filtered, !expanded_clone)
    if(missing_clone_id>0){
      singlets_db <- dplyr::bind_rows(singlets_db, no_clone_id)
    }

    # for expanded clones we then perform clonal clustering of light chains, select the best light chain(s) for each clone and then for each cell.
    if(nrow(exp_db)>0){
      #1. clonal clustering of all light chain contigs using a stringent threshold because light chains are less diverses and with overall smaller cdr3s...:
      cloned_exp_db <- exp_db
      cloned_exp_db[[locus]] <- "IGH" #'we trick scoper::defineClonesScoper in believing these are actually heavy chains...
      if("light_clone" %in% colnames(cloned_exp_db)){
        cloned_exp_db$light_clone <- NULL
      }
      
      if(seq_type == "Ig"){
        threshold <- 0.1
      }
      if(seq_type == "TCR"){
        threshold <- 0
      }
      message(paste0("clustering ", paste(chain, collapse = "-")," contigs"))
      suppressMessages(
        cloned_exp_db_analysis <- scoper::hierarchicalClones(cloned_exp_db,
                                                             threshold = threshold,
                                                             method = "nt",
                                                             normalize = "len",
                                                             linkage = "single",
                                                             cell_id = NULL,
                                                             locus = locus,
                                                             only_heavy = TRUE,
                                                             split_light = FALSE,
                                                             junction = junction,
                                                             v_call = v_call,
                                                             j_call = j_call,
                                                             clone = "light_clone",
                                                             first = FALSE,
                                                             cdr3 = FALSE,
                                                             mod3 = FALSE,
                                                             max_n = 0,
                                                             nproc = nproc,
                                                             verbose = FALSE,
                                                             log = NULL,
                                                             summarize_clones = TRUE)
      )
      if(output){
        pdf(file = paste0(output_folder, analysis_name, "_", paste(chain, collapse = "-"), "_clustering_0.2.pdf"))
        scoper::plot(cloned_exp_db_analysis, binwidth=0.02)
        dev.off()
      }
      cloned_exp_db <- scoper::as.data.frame(cloned_exp_db_analysis)
      # importing light_clone info into exp_db (with locus column still intact)
      exp_db <- dplyr::left_join(exp_db, cloned_exp_db[,c(sequence_id, "light_clone")], by = join_by(!!rlang::sym(sequence_id)))

      # 3. identify shared light chains between clones
      # step 1: identify shared light_clone_id across more than one cell_id
      shared_seqs <- exp_db %>%
        dplyr::group_by(!!rlang::sym(clone_id), light_clone) %>%
        dplyr::filter(n_distinct(!!rlang::sym(cell_id)) > 1) %>%  # Only keep l_clone_id shared across more than one cell_id
        dplyr::ungroup() %>%
        dplyr::distinct()

      # step 2: in case a cell is sharing multiple light chains within a given clone, keep only the sequence_id with the most represented light_clone_id
      shared_seqs_filtered <- shared_seqs %>%
        dplyr::group_by(light_clone) %>%
        dplyr::mutate(
          l_clone_count = n()
        ) %>%  # Count occurrences of each l_clone_id
        dplyr::ungroup() %>%
        # For each cell, keep only the sequence_id corresponding to the most frequent l_clone_id overall
        dplyr::group_by(!!rlang::sym(cell_id)) %>%
        dplyr::filter(light_clone == light_clone[which.max(l_clone_count)]) %>%
        dplyr::ungroup() %>%
        dplyr::select(-l_clone_count)

      # step 3: update original dataframe to only keep filtered shared light chains for these cells,
      # remaining multiplets clones will then be simply resolved as other singlets in the next step
      cells_with_shared_LC <- unique(shared_seqs[[cell_id]])
      exp_db <- exp_db %>%
        dplyr::filter(!(!!rlang::sym(cell_id) %in% cells_with_shared_LC)) %>%
        dplyr::bind_rows(shared_seqs_filtered)

      db_filtered <- dplyr::bind_rows(exp_db, singlets_db)
    } else {db_filtered <- singlets_db}
  }

  # filter on bin, quality (complete_vdj as productive already removed above) then umi_counts/consensus_counts
  db_filtered <- db_filtered %>%
    dplyr::group_by(!!rlang::sym(cell_id)) %>%
    dplyr::arrange(desc(unfiltered_bin), desc(!!rlang::sym(complete_vdj)), desc(!!rlang::sym(umi_count)), desc(!!rlang::sym(consensus_count))) %>%
    dplyr::mutate(
      # save remaining secondary contigs (productive, >20% of total reads and different V/J/CDR3)
      across(all_of(second_columns), ~ lead(.), .names = "second_{.col}")
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cell_id)

  # export final density plots on filtered datasets
  db_filtered <- db_filtered %>%
    dplyr::group_by(!!rlang::sym(split.by)) %>%
    dplyr::mutate(
      filtered_bin = binContigs(log10(umi_count),
                                plot = output,
                                title = paste0(dplyr::first(!!rlang::sym(split.by)), " - ", paste(chain, collapse = "-")," Density plot"),
                                file = paste0(output_folder, analysis_name, "_", dplyr::first(!!rlang::sym(split.by)), "_filtered_", paste(chain, collapse = "-"),"_density_plot.pdf"))[["bins"]]) %>%
    ungroup()

  # do a bit of QC, comparing final results to original SevenBridges filtering
  if(any(db_to_filter[[assay]] == "BD")){
    db_BD_filtered <- db_to_filter %>%
      dplyr::filter(assay == "BD") %>%
      dplyr::group_by(!!rlang::sym(cell_id)) %>%
      dplyr::mutate(
        pct_of_reads = !!rlang::sym(consensus_count)/sum(!!rlang::sym(consensus_count)),
      ) %>%
      dplyr::filter(!!rlang::sym(dominant)) %>%
      dplyr::arrange(desc(!!rlang::sym(productive)), desc(!!rlang::sym(complete_vdj)), desc(!!rlang::sym(umi_count)), desc(!!rlang::sym(consensus_count))) %>%
      dplyr::slice(1) %>%
      ungroup()

    colnames(db_BD_filtered) <- paste0("SB_", colnames(db_BD_filtered))
    db_BD_filtered <- db_BD_filtered %>%
      dplyr::rename(!!rlang::sym(cell_id) := !!rlang::sym(paste0("SB_", cell_id)))

    merged <- dplyr::left_join(db_BD_filtered, db_filtered, by = join_by(!!rlang::sym(cell_id)))
    merged <- merged %>%
      dplyr::mutate(
        same_seq = !!rlang::sym(sequence_id) == !!rlang::sym(paste0("SB_", sequence_id)),
        highest_umi_count = !!rlang::sym(umi_count) > !!rlang::sym(paste0("SB_", umi_count)),
        more_complete = (!!rlang::sym(complete_vdj) & !(!!rlang::sym(paste0("SB_", complete_vdj)))),
        more_productive = (!!rlang::sym(productive) & !(!!rlang::sym(paste0("SB_", productive)))),
        highest_freq = pct_of_reads > SB_pct_of_reads,
        toss_up = (pct_of_reads > 0.2 & SB_pct_of_reads > 0.2 &
                     !!rlang::sym(umi_count) == !!rlang::sym(paste0("SB_", umi_count)) &
                     !!rlang::sym(consensus_count) == !!rlang::sym(paste0("SB_", consensus_count)) &
                     !!rlang::sym(complete_vdj) == !!rlang::sym(paste0("SB_", complete_vdj)) &
                     !!rlang::sym(productive) == !!rlang::sym(paste0("SB_", productive)))
      )
    if(use_clone & all(chain %in% c("IGL", "IGK"))){
      merged <- merged %>%
        dplyr::mutate(
          shared_seq = !!rlang::sym(cell_id) %in% cells_with_shared_LC
        )
    }

    # table(merged$dominant, merged$highest_umi_count, merged$same_seq, useNA= "ifany")
    total_cells <- nrow(merged)
    same_seq <- nrow(dplyr::filter(merged, same_seq))

    cells_with_issues <- merged %>%
      dplyr::filter(!same_seq)
    # highest umi_count cases (likely CDR3 step gain)
    best_all <- nrow(dplyr::filter(cells_with_issues, highest_umi_count, more_complete, more_productive))
    highest_umi_count_more_productive <- nrow(dplyr::filter(cells_with_issues, highest_umi_count, !more_complete, more_productive))
    highest_umi_count_more_complete <- nrow(dplyr::filter(cells_with_issues, highest_umi_count, more_complete, !more_productive))
    highest_umi_count_only <- nrow(dplyr::filter(cells_with_issues, highest_umi_count, !more_complete, !more_productive))
    # lowest umi_count cases (likely bin step gain)
    lowest_umi_count_more_complete_productive <- nrow(dplyr::filter(cells_with_issues, !highest_umi_count, more_complete, more_productive))
    lowest_umi_count_more_productive <- nrow(dplyr::filter(cells_with_issues, !highest_umi_count, !more_complete, more_productive))
    lowest_umi_count_more_complete <- nrow(dplyr::filter(cells_with_issues, !highest_umi_count, more_complete, !more_productive))

    if(use_clone & resolve_chain == "light"){
      # other cases
      highest_freq <- nrow(dplyr::filter(cells_with_issues, !highest_umi_count, !more_complete, !more_productive, highest_freq, !shared_seq))
      toss_up <- nrow(dplyr::filter(cells_with_issues, toss_up, !shared_seq))
      true_issue <- nrow(dplyr::filter(cells_with_issues, !highest_umi_count, !more_complete, !more_productive, !highest_freq, !toss_up, !shared_seq))
      shared_seq <- nrow(dplyr::filter(cells_with_issues, !highest_umi_count, shared_seq))

      log_message <- paste0("out of ", total_cells," cells with ", paste(chain, collapse = "-")," contigs submitted from BD rhapsody pipeline:\n\n",
                            same_seq/total_cells*100, "% (n = ",same_seq,") in full agreement with SevenBridges filtration;\n\n",
                            best_all/total_cells*100, "% (n = ",best_all,") selected based on highest umi_count, productive and complete as compared to SevenBridges filtration (likely pooled CDR3);\n\n",
                            highest_umi_count_more_productive/total_cells*100, "% (n = ",highest_umi_count_more_productive,") selected based on highest umi_count and productive as compared to SevenBridges filtration (likely pooled CDR3);\n\n",
                            highest_umi_count_more_complete/total_cells*100, "% (n = ",highest_umi_count_more_complete,") selected based on highest umi_count and complete as compared to SevenBridges filtration (likely pooled CDR3);\n\n",
                            highest_umi_count_only/total_cells*100, "% (n = ",highest_umi_count_only,") selected based on highest umi_count as compared to SevenBridges filtration (likely pooled CDR3);\n\n",
                            lowest_umi_count_more_complete_productive/total_cells*100, "% (n = ",lowest_umi_count_more_complete_productive,") selected based on lowest umi_count but productive and complete as compared to SevenBridges filtration (likely bin step);\n\n",
                            lowest_umi_count_more_productive/total_cells*100, "% (n = ",lowest_umi_count_more_productive,") selected based on lowest umi_count but productive as compared to SevenBridges filtration (likely bin step);\n\n",
                            lowest_umi_count_more_complete/total_cells*100, "% (n = ",lowest_umi_count_more_complete,") selected based on lowest umi_count but complete as compared to SevenBridges filtration (likely bin step);\n\n",
                            highest_freq/total_cells*100, "% (n = ",highest_freq,") selected based highest reads allowing this contig to pass the 20% of total cell read cut-off (likely pooled CDR3);\n\n",
                            shared_seq/total_cells*100, "% (n = ",shared_seq,") selected based on lowest umi_count but shared between groups (all other shared where included in other selection steps);\n\n")

      if(toss_up>0){
        log_message <- paste0(log_message, toss_up/total_cells*100, "% (n = ",toss_up,") where simply toss-up...;\n\n")
      }
      if(true_issue>0){
        log_message <- paste0(log_message, true_issue/total_cells*100, "% (n = ",true_issue,") represent unexplainable selection (check code?);\n\n")
      }

      cat(log_message)

    } else {
      # other cases
      highest_freq <- nrow(dplyr::filter(cells_with_issues, !highest_umi_count, !more_complete, !more_productive, highest_freq))
      toss_up <- nrow(dplyr::filter(cells_with_issues, toss_up))
      true_issue <- nrow(dplyr::filter(cells_with_issues, !highest_umi_count, !more_complete, !more_productive, !highest_freq, !toss_up))

      log_message <- paste0("out of ", total_cells," cells with ", paste(chain, collapse = "-")," contigs submitted:\n\n",
                            same_seq/total_cells*100, "% (n = ",same_seq,") in full agreement with SevenBridges filtration;\n\n",
                            best_all/total_cells*100, "% (n = ",best_all,") selected based on highest umi_count, productive and complete as compared to SevenBridges filtration (likely pooled CDR3);\n\n",
                            highest_umi_count_more_productive/total_cells*100, "% (n = ",highest_umi_count_more_productive,") selected based on highest umi_count and productive as compared to SevenBridges filtration (likely pooled CDR3);\n\n",
                            highest_umi_count_more_complete/total_cells*100, "% (n = ",highest_umi_count_more_complete,") selected based on highest umi_count and complete as compared to SevenBridges filtration (likely pooled CDR3);\n\n",
                            highest_umi_count_only/total_cells*100, "% (n = ",highest_umi_count_only,") selected based on highest umi_count as compared to SevenBridges filtration (likely pooled CDR3);\n\n",
                            lowest_umi_count_more_complete_productive/total_cells*100, "% (n = ",lowest_umi_count_more_complete_productive,") selected based on lowest umi_count but productive and complete as compared to SevenBridges filtration (likely bin step);\n\n",
                            lowest_umi_count_more_productive/total_cells*100, "% (n = ",lowest_umi_count_more_productive,") selected based on lowest umi_count but productive as compared to SevenBridges filtration; (likely bin step)\n\n",
                            lowest_umi_count_more_complete/total_cells*100, "% (n = ",lowest_umi_count_more_complete,") selected based on lowest umi_count but complete as compared to SevenBridges filtration (likely bin step);\n\n",
                            highest_freq/total_cells*100, "% (n = ",highest_freq,") selected based highest reads allowing this contig to pass the 20% of total cell read cut-off (likely pooled CDR3);\n\n")
      if(toss_up>0){
        log_message <- paste0(log_message, toss_up/total_cells*100, "% (n = ",toss_up,") where simply toss-up...;\n\n")
      }
      if(true_issue>0){
        log_message <- paste0(log_message, true_issue/total_cells*100, "% (n = ",true_issue,") represent unexplainable selection (check code?);\n\n")
      }

      cat(log_message)
    }
  }

  final_seq_nb_post_filtering <- nrow(db_filtered)

  final_log_message <- paste0(final_seq_nb_post_filtering," contigs selected among the ",ini_seq_nb_for_filtering," ", paste(chain, collapse = "-")," contigs submitted.\n")
  cat(final_log_message)  # Print to console

  final_log_message <- paste0(final_log_message,
                              "filtration was made following the following rules: Read count of the contig is at least 20% of the total read count from all contigs for the cell-chain,\n",
                              "to break any ties, the contigs are then sorted in order of: productivity, umi_count bin, full length, highest molecule count, highest read count.\n")
  if(resolve_multi_CDR3){
    final_log_message <- paste0(final_log_message,
                                "contigs with identical v_call, j_call and junction_aa inside a given cell_id where pooled and only the best quality sequence was kept, updating umi counts, read counts and c_call accordingly.\n")
  }
  if(use_clone & all(chain %in% c("IGL", "IGK"))){
    final_log_message <- paste0(final_log_message,
                                "contigs shared inside a given clones were selected in priority.\n")
  }

  cat(final_log_message) 

  resolved_db <- db_filtered %>%
    dplyr::bind_rows(db_not_to_filter) %>%
    dplyr::arrange(!!rlang::sym(cell_id), !!rlang::sym(locus))

  if(use_clone){
    resolved_db <- resolved_db %>%
      dplyr::arrange(!!rlang::sym(clone_id))
  }

  if(is.null(split)){
    resolved_db <- resolved_db %>%
      dplyr::select(-!!rlang::sym(split.by))
  }

  return(resolved_db)
}

#### Duplicate of Dowser::resolveLightChains ####
#' Define subgroups within clones based on light chain rearrangements
#'
#' \code{resolveLightChains} resolve light chain V and J subgroups within a clone
#' @param    data         a tibble containing heavy and light chain sequences with clone_id
#' @param    nproc        number of cores for parallelization
#' @param    minseq       minimum number of sequences per clone
#' @param    locus        name of column containing locus values
#' @param    heavy        value of heavy chains in locus column. All other values will be
#'                        treated as light chains
#' @param    id           name of the column containing sequence identifiers.
#' @param    seq          name of the column containing observed DNA sequences. All
#'                        sequences in this column must be multiple aligned.
#' @param    clone        name of the column containing the identifier for the clone. All
#'                        entries in this column should be identical.
#' @param    cell         name of the column containing identifier for cells.
#' @param    v_call       name of the column containing V-segment allele assignments. All
#'                        entries in this column should be identical to the gene level.
#' @param    j_call       name of the column containing J-segment allele assignments. All
#'                        entries in this column should be identical to the gene level.
#' @param    junc_len     name of the column containing the length of the junction as a
#'                        numeric value. All entries in this column should be identical
#'                        for any given clone.
#' @param    nolight      string to use to indicate a missing light chain
#' @param    pad_ends          pad sequences within a clone to same length?
#'
#' @return   a tibble containing the same data as inputting, but with the column clone_subgroup
#' added. This column contains subgroups within clones that contain distinct light chain
#' V and J genes, with at most one light chain per cell.
#' @details
#' 1. Make temporary array containing light chain clones
#' 2. Enumerate all possible V, J, and junction length combinations
#' 3. Determine which combination is the most frequent
#' 4. Assign sequences with that combination to clone t
#' 5. Copy those sequences to return array
#' 6. Remove all cells with that combination from temp array
#' 7. Repeat 1-6 until temporary array zero.
#' If there is more than rearrangement with the same V/J
#' in the same cell, pick the one with the highest non-ambiguous
#' characters. Cells with missing light chains are grouped with their
#' subgroup with the closest matching heavy chain (Hamming distance)
#' then the largest and lowest index subgroup if ties are present.
#'
#' Outputs of the function are
#' 1. clone_subgroup which identifies the light chain VJ rearrangement that sequence belongs to within it's clone
#' 2. clone_subgroup_id which combines the clone_id variable and the clone_subgroup variable by a "_".
#' 3. vj_cell which combines the vj_gene and vj_alt_cell columns by a ",".
#' [note] pure duplicated function from Dowser:resolveLightChains to circumvent an issue with unlist() in R 4.5.0
#'
#' @export

resolveLightChains2 <- function(data, nproc=1, minseq=1,locus="locus",heavy="IGH",
                               id="sequence_id", seq="sequence_alignment",
                               clone="clone_id", cell="cell_id", v_call="v_call", j_call="j_call",
                               junc_len="junction_length", nolight="missing", pad_ends=TRUE){

  subgroup <- "clone_subgroup"

  light <- data[data[[locus]] != heavy,]
  heavy <- data[data[[locus]] == heavy,]

  if(nrow(heavy) == 0){
    stop("No heavy chains found in data")
  }
  if(nrow(light) == 0){
    warning("No light chains found in data! Assigning all sequences to subgroup 1.")
    data[[subgroup]] = 1
    return(data)
  }

  scount <- table(heavy[[clone]])
  big <- names(scount)[scount >= minseq]
  heavy <- dplyr::filter(heavy,(!!rlang::sym(clone) %in% big))

  if(max(table(heavy[[id]])) > 1){
    stop("Sequence IDs in heavy dataframe must be unique!")
  }
  if(max(table(light[[id]])) > 1){
    stop("Sequence IDs in light dataframe must be unique!")
  }

  heavycount = table(heavy[[cell]])
  if(max(heavycount) > 1){
    stop(paste0(sum(heavycount > 1),
                " cells with multiple heavy chains found. Remove before proceeeding"))
  }

  # filter out light chains with missing cell IDs
  missing_cell <- light[is.na(light[[cell]]),]
  if(nrow(missing_cell) > 0){
    warning(paste("removing",nrow(missing_cell),"light chains with missing cell IDs."))
    light <- light[!is.na(light[[cell]]),]
  }

  heavy$vj_gene <- nolight
  heavy$vj_alt_cell <- NA # set these to NA, since they're pretty rare
  heavy[[subgroup]] <- 1
  light$vj_gene <- nolight
  light$vj_alt_cell <- NA
  light[[subgroup]] <- 1
  light[[clone]] <- -1
  paired <- safe_mclapply(unique(heavy[[clone]]),function(cloneid){
    # Get heavy chains within a clone, and corresponding light chains
    # separate heavy chains with (sc) and without (bulk) paired light chains
    hd <- dplyr::filter(heavy,!!rlang::sym(clone) == cloneid)
    ld <- dplyr::filter(light,!!rlang::sym(cell) %in% hd[[!!cell]])
    ld <- dplyr::filter(ld, !is.na(!!rlang::sym(cell)))
    hd_sc <- hd[hd[[cell]] %in% ld[[cell]] & !is.na(hd[[cell]]),] # added is.na(cell) catch
    hd_bulk <- hd[!hd[[cell]] %in% ld[[cell]] | is.na(hd[[cell]]),]
    if(nrow(ld) == 0){
      hd$clone_subgroup_id <- paste0(hd[[clone]],"_",hd[[subgroup]])
      hd$vj_cell <- sapply(1:nrow(hd), function(x){
        if(!is.na(hd$vj_alt_cell[x])){
          paste(hd$vj_gene[x],hd$vj_alt_cell[x],sep=",")
        }else{
          hd$vj_gene[x]
        }
      })
      return(hd)
    }
    ltemp <- dplyr::filter(ld, !is.na(!!rlang::sym(cell)))
    # CGJ 9/17/24
    # make the gene level partiions be only gene level -- no allele
    ltemp$temp_v <- ltemp[[v_call]]
    ltemp$temp_j <- ltemp[[j_call]]
    ltemp[[v_call]] <- alakazam::getGene(ltemp[[v_call]])
    ltemp[[j_call]] <- alakazam::getGene(ltemp[[j_call]])
    ltemp[[clone]] <- -1
    ld <- dplyr::tibble()
    lclone <- 1
    while(nrow(ltemp) > 0){
      #expand ambiguous V/J calls
      lvs <- strsplit(ltemp[[v_call]],split=",")
      ljs <- strsplit(ltemp[[j_call]],split=",")
      jlens <- ltemp[[junc_len]]
      # get all combinations of V/J calls for each light chain
      combos <-
        lapply(1:length(lvs),function(w)
          unlist(lapply(lvs[[w]],function(x)
            unlist(lapply(ljs[[w]],function(y)
              lapply(jlens[[w]], function(z)paste0(x,":",y,";",z)))))))

      # get unique combinations per cell
      cells <- unique(ltemp[[cell]])
      cellcombos <- lapply(cells,function(x)
        unique(unlist(combos[ltemp[[cell]] == x])))
      #lcounts <- table(unlist(lapply(cellcombos,function(x)x)))
      #lcounts <- table(unlist(cellcombos,function(x)x))
      lcounts <- table(unlist(cellcombos))
      max <- names(lcounts)[which.max(lcounts)]
      cvs <- unlist(lapply(combos,function(x)max %in% x))
      ltemp[cvs,][[subgroup]] <- lclone
      ltemp[cvs,]$vj_gene <- max

      # if a cell has the same combo for two rearrangements, only pick one
      # with the most ACTG characters
      rmseqs <- c()
      cell_counts <- table(ltemp[cvs,][[cell]])
      mcells <- names(cell_counts)[cell_counts > 1]
      for(cellname in mcells){
        # CGJ 1/27/25 -- old way now causes internal dplyr error
        ttemp <- ltemp[cvs & ltemp[[cell]] == cellname, ]
        # ttemp <- dplyr::filter(ltemp,cvs & !!rlang::sym(cell) == cellname)
        ttemp$str_counts <-
          stringr::str_count(ttemp[[seq]],"[A|C|G|T]")
        # keep version with most non-N characters
        keepseq <- ttemp[[id]][which.max(ttemp$str_counts)]
        rmtemp <- ttemp[!ttemp[[id]] == keepseq,]
        rmseqs <- c(rmseqs,rmtemp[[id]])
      }
      # CGJ 1/27/25 -- old way now causes internal dplyr error
      include <- ltemp[cvs & !(ltemp[[id]] %in% rmseqs), ]
      leave <- ltemp[!cvs & (ltemp[[id]] %in% rmseqs), ]
      # include <- dplyr::filter(ltemp, cvs & !(!!rlang::sym(id) %in% rmseqs))
      # leave <- dplyr::filter(ltemp,!cvs | (!!rlang::sym(id) %in% rmseqs))

      # find other cells still in ltemp and add as vj_alt_cell
      mcells <- unique(include[[cell]])
      for(cellname in mcells){
        if(cellname %in% leave[[cell]]){
          include[include[[cell]] == cellname,]$vj_alt_cell <-
            paste(paste0(leave[leave[[cell]] == cellname,][[v_call]],":",
                         leave[leave[[cell]] == cellname,][[j_call]]),
                  collapse=",")
        }
      }
      # CGJ 9/17/24
      # update the include df to have proper v and j call and remove their temp cols
      include[[v_call]] <- include$temp_v
      include[[j_call]] <- include$temp_j
      rm_indx <- which(colnames(include) %in% c("temp_v", "temp_j"))
      include <- include[, -rm_indx]

      ld <- dplyr::bind_rows(ld,include)
      ltemp <- dplyr::filter(ltemp,!(!!rlang::sym(cell) %in% ltemp[cvs,][[!!cell]]))
      lclone <- lclone + 1
    }
    ld[[clone]] <- cloneid
    for(cellname in unique(hd_sc[[cell]])){
      if(cellname %in% ld[[cell]]){
        lclone <- ld[ld[[cell]] == cellname,][[subgroup]]
        hd_sc[hd_sc[[cell]] == cellname,][[subgroup]] <- lclone
        hd_sc[hd_sc[[cell]] == cellname,]$vj_gene <- ld[ld[[cell]] == cellname,]$vj_gene
        hd_sc[hd_sc[[cell]] == cellname,]$vj_alt_cell <-
          ld[ld[[cell]] == cellname,]$vj_alt_cell
      }
    }
    # now get the subgroup_id for the heavy chains lacking paired light chains
    comb <- dplyr::bind_rows(hd_sc,ld)
    comb[[subgroup]] <- as.integer(comb[[subgroup]])
    if(nrow(ld) != 0 & nrow(hd_bulk) != 0){
      for(sequence in 1:nrow(hd_bulk)){
        # CGJ 8/6/24 -- updated to do the padding on the temp sequences
        rating <- unlist(lapply(1:nrow(hd_sc), function(x){
          temp <- rbind(hd_sc[x,], hd_bulk[sequence,])
          if(nchar(temp[[seq]][1] != nchar(temp[[seq]][2]))){
            temp <- alakazam::padSeqEnds(temp[[seq]])
          }
          value <- alakazam::seqDist(temp[1], temp[2])
          return(value)
        }))
        rating <- as.numeric(rating)
        # row number of heavy chain only df with lowest seq dist
        proper_index <- which(rating == min(rating))
        if(length(proper_index) > 1){
          # find the subgroups that belong to the lowest seq dists
          subgroups <- hd_sc[[subgroup]][proper_index]
          if(length(unique(subgroups)) > 1){
            # if there is more than one subgroup find the subgroup sizes of the
            # subgroups being considered
            subgroup_size <- data.frame(clone_subgroup = unique(subgroups))
            subgroup_size$sizes <- unlist(lapply(1:nrow(subgroup_size), function(x){
              nrow(hd_sc[hd_sc[[subgroup]] == subgroup_size$clone_subgroup[x],])
            }))
            # if there is one subgroup that is the largest use it
            if(length(which(subgroup_size$sizes == max(subgroup_size$sizes))) == 1){
              proper_index_value <- subgroup_size$clone_subgroup[
                which(subgroup_size$sizes == max(subgroup_size$sizes))]
            } else {
              # if there are more than one subgroup with the same size use the lower number
              potential_subgroups <- subgroup_size$clone_subgroup[
                which(subgroup_size$sizes == max(subgroup_size$sizes))]
              proper_index_value <- min(potential_subgroups)
            }
          } else{
            proper_index_value <- hd_sc[[subgroup]][proper_index[1]]
          }
        } else{
          proper_index_value <- hd_sc[[subgroup]][proper_index]
        }
        hd_bulk[[subgroup]][sequence] <- proper_index_value
      }
    }
    if(nrow(hd_bulk) != 0){
      comb <- dplyr::bind_rows(comb, hd_bulk)
    }
    comb$clone_subgroup_id <- paste0(comb[[clone]],"_",comb[[subgroup]])
    comb$vj_cell <- sapply(1:nrow(comb), function(x){
      if(!is.na(comb$vj_alt_cell[x])){
        paste(comb$vj_gene[x],comb$vj_alt_cell[x],sep=",")
      }else{
        comb$vj_gene[x]
      }
    })

    size <- c()
    for(subgroups in sort(unique(comb[[subgroup]]))){
      size <- append(size, nrow(comb[comb[[subgroup]] == subgroups,]))
    }
    if(!all(diff(size) <= 0)){
      order_check <- data.frame(table(comb[[subgroup]]))
      colnames(order_check) <- c(subgroup, "size")
      order_check <- order_check[order(-order_check$size),]
      order_check$proper_subgroup <- 1:nrow(order_check)
      comb$new_subgroup <- NA
      for(i in unique(comb[[subgroup]])){
        comb$new_subgroup[comb[[subgroup]] == i] <- order_check$proper_subgroup[order_check[[subgroup]] == i]
      }
      comb <- comb[, -which(names(comb) == subgroup)]
      names(comb)[names(comb) == "new_subgroup"] <- subgroup
    }
    comb
  },mc.cores=nproc)
  paired <- dplyr::bind_rows(paired)
  # remove the junction length from the vj_gene
  paired$vj_gene <- gsub("\\;.*", "", paired$vj_gene)
  return(paired)
}


#### Function to flag B or T cell doublets in single cell data ####
#' check for cases where two IGH or two IGH and two light chains are found in a cell and could represent a B cell doublet.
#'
#' \code{homotypicVDJdoublets} flag potential B cell doublets.
#'
#' @param db            a AIRR formatted dataframe containing heavy and light chain sequences, cell_ids and umi_counts.
#' @param split.by      name of the column in the dataframe to use to split the dataset prior to calculating variable cut-offs
#' @param analysis_name name to use for outputs prefixes.
#' @param use_chains    which chain to use to flag VDJ doublets, can de set to "IGH" or "all" [default = "all"].
#' @param locus         name of column containing locus values.
#' @param seq_type      'Ig' or 'TCR'
#' @param heavy         heavy chains to be used. if set to NULL, will default to IGH for Ig and TRB/TRD for TCR.
#' @param light         light chains to be used. if set to NULL, will default to IGL/IGK for Ig and TRA/TRG for TCR.
#' @param assay         name of column containing assay values.
#' @param scRNAseq.tech vector with names of scRNA-seq technologies used in the assay column, i.e. those from which we have data so far.
#' @param variable_cutoff whether to use variable cutoffs [default = TRUE, cutoffs are automatically calculated based on summary counts in the database to account for potential sequencing bias (1st quartile and max of 10 or 1/10 of median).].
#' @param low_cutoff    cut_off for low probability heavy chain doublets (only used if variable_cutoff = FALSE).
#' @param high_cutoff   cut_off for high probability heavy chain doublets (only used if variable_cutoff = FALSE).
#' @param output        whether to output graphs with umi_counts for dominant versus second IGH VDJ contig and the recap excel workbook. Can be set to none, partial (only pdf + log file with QC parameters) or total (pdf + log file with QC parameters + tsv file with sequences).
#' @param output_folder name of the folder in which graph and recap excel workbooks will be saved [default = "VDJ_QC"].
#' @param cell_id       name of the column containing cell identifier.
#' @param umi_count     name of the column containing the number of unique molecules (UMI) for this contig. Previously called "duplicate_count" in an earlier AIRR standard
#' @param save_plot     format to save the final ggplot2 object. "pdf" or "png". if set to any other value, no plot will be saved.
#' @param verbose       whether to print recap informations to the console
#'
#' @return   the original data frame with for the following additional columns: "is.VDJ_doublet" and "is.IGH_doublet.confidence" (use_chains = "IGH") or "is.IGH_doublet.confidence", "is.IGK_doublet.confidence" and "is.IGL_doublet.confidence" (use_chains = "all")
#' as well as graphs plotting the primary versus "second_" umi_counts by locus for all cell_ids colored by "is."locus"_VDJ_doublet.confidence".
#'
#' @details
#' first group by locus, split.by and calculate variable cutoffs (upper = max(10, 1st quartile of distribution); lower = max(5, median(distribution)/10))
#' then calculate doublet.confidences for all locus: a high confidence doublet is defined as umi_count > 10 and !(umi_count>3*second_umi_count & second_umi_count < upper_cutoff), for low_confidence doublets, we move to second_umi_count < lower_cutoff.
#' then is.VDJ_doublet is defined as that cell must be high for at least one of the chain (and low for the others if multiple chains are used).
#' dependencies: dplyr, tidyr, ggplot2
#'
#' @export
#'
#' @import dplyr

homotypicVDJdoublets <- function(db,
                                 split.by = NULL,
                                 analysis_name = "All_sequences",
                                 use_chains = c("all", "heavy"),
                                 locus = "locus", 
                                 seq_type = c("Ig", "TCR"),
                                 heavy = NULL,
                                 light = NULL, 
                                 assay = "assay",
                                 scRNAseq.tech = c("10X", "BD"),
                                 variable_cutoff = TRUE, 
                                 low_cutoff = 10, 
                                 high_cutoff = 250,
                                 output = TRUE, 
                                 output_folder = "VDJ_QC",
                                 cell_id = "cell_id", 
                                 umi_count = "umi_count",
                                 save_plot = c("pdf", "png"),
                                 verbose = TRUE){
                          

  seq_type <- match.arg(seq_type)
  
  if(is.null(heavy)){
    if(seq_type == "Ig"){
      heavy <- "IGH"
    }
    if(seq_type == "TCR"){
      heavy <- c("TRB","TRD")
    }
  }
  if(is.null(light)){
    if(seq_type == "Ig"){
      light <- c("IGK","IGL")
    }
    if(seq_type == "TCR"){
      light <- c("TRA","TRG")
    }
  }
  
  # define chains to be used:
  use_chains <- match.arg(use_chains)
  
  if(!use_chains %in% c("heavy", "all")){
    stop("use_chains must be one of heavy or all.")
  } else {
    if(use_chains == "all") {
      chains <- c("heavy", "light")
    } else {
      chains <- "heavy"
    }
  }
  
  db <- db %>%
    dplyr::mutate(
      locus_simplified = ifelse(locus %in% heavy, "heavy", ifelse(locus %in% light, "light", "other"))
    ) %>%
    dplyr::select(-dplyr::any_of(
      c("is.VDJ_doublet", 
        "is.VDJ_doublet.confidence", 
        paste0("is.", paste(heavy, collapse = "-"),"_doublet.confidence"), 
        paste0("is.", paste(light, collapse = "-"),"_doublet.confidence"))
    ))
  
  # filter on scRNA-seq contigs and prepare the dataset:
  flagged_db <- db %>%
    dplyr::filter(
      locus_simplified %in% chains & !!rlang::sym(assay) %in% scRNAseq.tech
    ) 
  
  other_db <- db %>%
    dplyr::filter(
      !(locus_simplified %in% chains & !!rlang::sym(assay) %in% scRNAseq.tech)
    ) 
  
  if(is.null(split.by)){
    split <- NULL
    flagged_db <- flagged_db %>%
      dplyr::mutate(
        split.by = "all_sequences"
      )
    split.by = "split.by"
  }
  
  for(chain in chains){
    locus_db <- flagged_db %>%
      dplyr::filter(locus_simplified == chain)

    if(nrow(locus_db)==0){
      warning(paste0("no ",paste(unlist(mget(chain)), collapse = "-")," contig provided"))
      chains <- chains[!chains == chain]
      }
    if(any(duplicated(locus_db$cell_id))){
      stop(paste0(paste(unlist(mget(chain)), collapse = "-")," contig filtration must be done before runing flagVDJdoublets(), run resolveMultiContigs()."))
      }
    if(!paste0("second_", umi_count) %in% colnames(locus_db)){
      stop(paste0(paste(unlist(mget(chain)), collapse = "-")," contig filtration must be done before runing flagVDJdoublets(), run resolveMultiContigs()."))
      }
    if(paste0("second_", umi_count) %in% colnames(locus_db)){
      if(all(is.na(locus_db[[paste0("second_", umi_count)]]))){
        warning("non secondary ",paste(unlist(mget(chain)), collapse = "-")," contig found, either its a perfect dataset or ", paste(unlist(mget(chain)), collapse = "-")," contig filtration has not been run before runing flagVDJdoublets()")
        }
    }
    rm(locus_db)
  }
  
  # update cutoffs if needed:
  if(variable_cutoff){
    # learn dataset (split.by) and locus specific cutoffs (adapted to library depth):
    # 1st quartile of dominant VDJ umi_count (in 75% of doublets the second VDJ should be above this value) and minimum of 10 and 1/10 of Median (below could be considered ambient RNA, i.e; 1/10 of an average cell).
    flagged_db <- flagged_db %>%
      dplyr::group_by(locus_simplified, !!rlang::sym(split.by)) %>%
      dplyr::mutate(
        upper_cutoff = quantile(!!rlang::sym(umi_count), 0.25, na.rm = TRUE),
        lower_cutoff = ceiling(median(!!rlang::sym(umi_count), na.rm = TRUE)/10)
      ) %>%
      dplyr::ungroup()
  } else { # otherwise, use provided cutoffs
    flagged_db <- flagged_db %>%
      dplyr::mutate(
        upper_cutoff = high_cutoff,
        lower_cutoff = low_cutoff
      )
  }

  flagged_db <- flagged_db %>%
    dplyr::select(!!rlang::sym(cell_id), !!rlang::sym(split.by), locus_simplified,
                  !!rlang::sym(umi_count), !!rlang::sym(paste0("second_", umi_count)), upper_cutoff, lower_cutoff) %>%
    tidyr::pivot_wider(names_from = locus_simplified,
                       values_from = c(!!rlang::sym(umi_count), !!rlang::sym(paste0("second_", umi_count)), upper_cutoff, lower_cutoff)) %>%
    dplyr::rename_with(
      .fn = ~ sub(paste0(umi_count, "_"), "", .),
      .cols = dplyr::any_of(
        c(paste0(umi_count, "_heavy"), 
          paste0("second_", umi_count, "_heavy"),
          paste0(umi_count, "_light"), 
          paste0("second_", umi_count, "_light"))
        )
    ) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(chains),
             ~ ifelse(is.na(get(paste0("second_", dplyr::cur_column()))), "not a doublet",
                      ifelse(!(. > 3*get(paste0("second_", dplyr::cur_column())) & get(paste0("second_", dplyr::cur_column())) < get(paste0("upper_cutoff_", dplyr::cur_column()))), "high",
                             ifelse(!(. > 3*get(paste0("second_", dplyr::cur_column())) & get(paste0("second_", dplyr::cur_column())) < get(paste0("lower_cutoff_", dplyr::cur_column()))), "low", "not a doublet"))),
             .names = "is.{.col}_doublet.confidence")
    )

  if(!"light" %in% chains){
    # must be high for the heavy chain
    flagged_db <- flagged_db %>%
      dplyr::mutate(
        is.VDJ_doublet.confidence = is.heavy_doublet.confidence,
        is.VDJ_doublet = (is.VDJ_doublet.confidence == "high")
      ) %>%
      dplyr::select(
        dplyr::all_of(c(cell_id,
                        "is.VDJ_doublet", 
                        "is.VDJ_doublet.confidence", 
                        "is.heavy_doublet.confidence",
                        "upper_cutoff_heavy",
                        "lower_cutoff_heavy"))
        ) %>%
      dplyr::rename(
        !!rlang::sym(paste0("is.",paste(heavy, collapse = "-"),"_doublet.confidence")) := is.heavy_doublet.confidence
      )
  } else {
    # must be high for at least one of the chain and low for the other
    flagged_db <- flagged_db %>%
      dplyr::mutate(
        is.VDJ_doublet.confidence = ifelse(is.heavy_doublet.confidence == "high" & is.light_doublet.confidence %in% c("high", "low"), "high",
                                           ifelse(is.light_doublet.confidence == "high" & is.heavy_doublet.confidence %in% c("high", "low"), "high", 
                                                  ifelse(is.heavy_doublet.confidence %in% c("high", "low") | is.light_doublet.confidence %in% c("high", "low"), "low", "not a doublet"))),
        is.VDJ_doublet = is.VDJ_doublet.confidence == "high"
      ) %>%
      dplyr::select(
        any_of(c(cell_id,
                 "is.VDJ_doublet", 
                 "is.VDJ_doublet.confidence", 
                 "is.heavy_doublet.confidence", 
                 "is.light_doublet.confidence",
                 "upper_cutoff_heavy",
                 "upper_cutoff_light",
                 "lower_cutoff_heavy",
                 "lower_cutoff_light"))
      ) %>%
      dplyr::rename(
        !!rlang::sym(paste0("is.",paste(heavy, collapse = "-"),"_doublet.confidence")) := is.heavy_doublet.confidence,
        !!rlang::sym(paste0("is.",paste(light, collapse = "-"),"_doublet.confidence")) := is.light_doublet.confidence
      ) 
  }
  
  db <- db %>%
    dplyr::left_join(flagged_db, by = join_by(cell_id)) %>%
    dplyr::mutate(
      upper_cutoff = ifelse(locus_simplified == "heavy", upper_cutoff_heavy, 
                              ifelse(locus_simplified == "light", upper_cutoff_light, NA)),
      lower_cutoff = ifelse(locus_simplified == "heavy", lower_cutoff_heavy, 
                                                      ifelse(locus_simplified == "light", lower_cutoff_light, NA))
    )
  
  if(output){
    vdjQCplot(db,
              use_chains = use_chains,
              seq_type = seq_type,
              heavy = heavy,
              light = light,
              split.by = split.by,
              type = "vdj_doublet",
              variable_cutoff = variable_cutoff,
              plot_cutoff = TRUE,
              output_folder = output_folder,
              analysis_name = analysis_name,
              locus = locus,
              cell_id = cell_id,
              umi_count = umi_count,
              second_umi_count = paste0("second_", umi_count),
              save_plot = save_plot,
              return_plot = FALSE)
  }

  groups <- na.omit(levels(as.factor(db[[split.by]])))
  for(group in groups){
    group_db <- db %>%
      dplyr::filter(!!rlang::sym(split.by) == group, !!rlang::sym(locus) %in% heavy)
    cells_in_group <- nrow(group_db)
    doublets <- nrow(dplyr::filter(group_db, is.VDJ_doublet))
    upper_group_cutoff_heavy <- group_db[[paste0("upper_cutoff_", paste(heavy, collapse = "-"))]][1]
    upper_group_cutoff_light <- group_db[[paste0("upper_cutoff_", paste(light, collapse = "-"))]][1]
    lower_group_cutoff_heavy <- group_db[[paste0("lower_cutoff_", paste(heavy, collapse = "-"))]][1]
    lower_group_cutoff_light <- group_db[[paste0("lower_cutoff_", paste(light, collapse = "-"))]][1]
    group_log_message <- paste0(group, ": out of ", cells_in_group," cells:\n",
                                doublets/cells_in_group*100,"% cells (n = ",doublets,") identified as high probability doublets, using ", paste(unlist(mget(chains)), collapse = "-")," contigs;\n",
                                "upper cutoff used:", upper_group_cutoff_heavy, " for ", paste(heavy, collapse = "-"), " and ",upper_group_cutoff_light, " for ", paste(light, collapse = "-"),"\n",
                                "lower cutoff used:", lower_group_cutoff_heavy, " for ", paste(heavy, collapse = "-"), " and ",lower_group_cutoff_light, " for ", paste(light, collapse = "-"),"\n","\n")
    if(verbose){cat(group_log_message)}
  }
  
  db <- db %>%
    dplyr::select(
      -dplyr::any_of(c("locus_simplified",
                       "upper_cutoff_heavy",
                       "upper_cutoff_light",
                       "lower_cutoff_heavy",
                       "lower_cutoff_light"))
      )
  
  if(nrow(other_db)>0){
    #adding back non BCR or TCR contigs depending on the seq_type
    db <- db %>%
      dplyr::bind_rows(other_db)
  }
  
  return(db)
}


#### Function to identify nonB/T-B/T doublets based on azimuth prediction and VDJ detection ####
#' Predicts nonB-B or nonT-T doublets
#'
#' \code{heterotypicVDJdoublets} Predicts nonB/T-B/T doublets based on azimuth predictions and VDJ contigs umi_counts.
#' @param db            a AIRR formatted dataframe containing heavy and light chain sequences, cell_ids and umi_counts.
#' @param split.by      name of the column in the dataframe to use to split the dataset prior to calculating variable cut-offs
#' @param analysis_name name to use for outputs prefixes.
#' @param use_chains     which chain to use to flag VDJ doublets, can de set to "IGH" or "all" [default = "all"].
#' @param locus         name of column containing locus values.
#' @param heavy         heavy chains to be used. if set to NULL, will default to IGH for Ig and TRB/TRD for TCR.
#' @param light         light chains to be used. if set to NULL, will default to IGL/IGK for Ig and TRA/TRG for TCR.
#' @param assay         name of column containing assay values.
#' @param scRNAseq.tech vector with names of scRNA-seq technologies used in the assay column, i.e. those from which we have data so far.
#' @param variable_cutoff whether to use variable cutoffs [default = TRUE, cutoffs are automatically calculated based on summary counts in the database to account for potential sequencing bias (1st quartile and max of 10 or 1/10 of median).].
#' @param low_cutoff    cut_off for low probability heavy chain doublets (only used if variable_cutoff = FALSE).
#' @param high_cutoff   cut_off for high probability heavy chain doublets (only used if variable_cutoff = FALSE).
#' @param output        whether to output graphs with umi_counts for dominant versus second IGH VDJ contig and the recap excel workbook. Can be set to none, partial (only pdf + log file with QC parameters) or total (pdf + log file with QC parameters + tsv file with sequences).
#' @param output_folder name of the folder in which graph and recap excel workbooks will be saved [default = "VDJ_QC"].
#' @param ref           which ref was used for mapping of scRNA-seq dataset: one of "azimuth.pbmcref", "azimuth.tonsilref", "azimuth.bonemarrowref" [default azimuth datasets as of June 2025] [TODO: add gut dataset?]; 
#'                      to add you own clustering info call it 'clusters' and provide the column ref in ref.column = list("clusters"= "column name") and corresponding B or T cell clusters: ref.Bcelltypes = list("clusters" = c("C1", "C3"...))...
#' @param ref.column which output column from azimuth should be used [default values as of June 2025]
#' @param ref.Bcelltypes which clusters in the azimuth output correspond to B/PC cells populations likely to carry a BCR [default values as of June 2025].
#' @param ref.Tcelltypes which clusters in the azimuth output correspond to T cells populations likely to carry a TCR [default values as of June 2025].
#' @param cell_id       name of the column containing cell identifier.
#' @param umi_count     name of the column containing the number of unique molecules (UMI) for this contig. Previously called "duplicate_count" in an earlier AIRR standard
#' @param save_plot     format to save the final ggplot2 object. "pdf" or "png". if set to any other value, no plot will be saved.
#' @param verbose       whether to print recap informations to the console
#'
#' @return   the original data frame with for two new columns: 'is.nonB_VDJ_doublet' and 'is.nonB_VDJ_doublet.confidence' for 'Ig' or nonT for 'TCR' as well as graphs plotting the heavy and light umi_counts for all cell_ids colored by "is.nonB/T_VDJ_doublet.confidence".
#'
#' @details
#' Will define variable cutoffs based on librairie(s) distribution, and assign a probability of being a true nonB-B doublet for each cell not predicted as being a B cell by azimuth.
#' High probably nonB-B doublets are define as having both heavy and light chain in the top three quartiles of IGH/IGL/IGK expression in the dataset yet not being mapped to a B cell cluster by azimuth.
#'
#' https://azimuth.hubmapconsortium.org/
#' for "bonemarrowref", use predicted.celltype.l2 and ! %in% c("Memory B", "Naive B", "Plasma", "pro B, "pre B", "transitional B")
#' for "pbmcref", use predicted.celltype.l2 and ! %in% c("B intermediate", "B memory", "B naive", "Plasmablast")
#' for "tonsilref v2", use predicted.celltype.l1 and ! %in% c("B activated", "B memory", "B naive", "preGCB", "PB", "PC", "PC/doublet", "preMBC/doublet", "prePB", "Cycling DZ GCB","DZ GCB","DZtoLZ GCB transition","FCRL4/5+ B memory","LZ GCB","LZtoDZ GCB transition")
#'
#' dependencies: dplyr, tidyr, ggplot2 and openxlsx
#'
#' @export
#'
#' @import dplyr

heterotypicVDJdoublets <- function(db,
                                   split.by = NULL,
                                   analysis_name = "All_sequences",
                                   locus = "locus", 
                                   seq_type = c("Ig", "TCR"),
                                   use_chains = c("all", "heavy"),
                                   heavy = NULL, 
                                   light = NULL,
                                   assay = "assay",
                                   scRNAseq.tech = c("10X", "BD"),
                                   output = TRUE,
                                   output_folder = "VDJ_QC/",
                                   variable_cutoff =TRUE,
                                   low_cutoff = 10,
                                   high_cutoff = 250,
                                   ref = c("azimuth.pbmcref", "azimuth.tonsilref", "azimuth.bonemarrowref"),
                                   ref.column = list("azimuth.pbmcref" = "predicted.celltype.l2", 
                                                     "azimuth.tonsilref" = "predicted.cell_type.l1", 
                                                     "azimuth.bonemarrowref" = "predicted.celltype.l2"),
                                   ref.Bcelltypes = list("azimuth.pbmcref" = c("B intermediate", "B memory", "B naive", "Plasmablast"),
                                                         "azimuth.tonsilref" = c("B activated", "B memory", "B naive", "preGCB", "PB", "PC", "PC/doublet", "preMBC/doublet", "prePB", "Cycling DZ GCB","DZ GCB","DZtoLZ GCB transition","FCRL4/5+ B memory","LZ GCB","LZtoDZ GCB transition"),
                                                         "azimuth.bonemarrowref" = c("Memory B", "Naive B", "Plasma", "pro B", "pre B", "transitional B")),
                                   ref.Tcelltypes = list("azimuth.pbmcref" = c("CD4 CTL", "CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM", "Treg", "CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM", "dnT", "gdT", "MAIT"),
                                                         "azimuth.tonsilref" = c("CD4 naive", "CD4 Non-TFH", "CD4 TCM", "CD4 TFH", "CD4 TFH Mem", "CD4 TREG", "CD8 naive", "CD8 T", "CD8 TCM", "Cycling T","dnT","MAIT/TRDV2+ gdT","non-TRDV2+ gdT"),
                                                         "azimuth.bonemarrowref" = c("CD4 Effector", "CD4 Memory", "CD4 Naive", "CD8 Effector 1", "CD8 Effector 2", "CD8 Effector 3", "CD8 Memory", "CD8 Naive", "MAIT", "T proliferating")),
                                   cell_id = "cell_id", 
                                   umi_count = "umi_count", 
                                   save_plot = c("pdf", "png"),
                                   verbose = TRUE){
                             
  seq_type <- match.arg(seq_type)
  cell_type <- ifelse(seq_type == "Ig", "B", "T")
  
  if(is.null(heavy)){
    if(seq_type == "Ig"){
      heavy <- "IGH"
    }
    if(seq_type == "TCR"){
      heavy <- c("TRB","TRD")
    }
  }
  if(is.null(light)){
    if(seq_type == "Ig"){
      light <- c("IGK","IGL")
    }
    if(seq_type == "TCR"){
      light <- c("TRA","TRG")
    }
  }
  
  ref <- match.arg(ref)
  
  if (seq_type == "Ig"){
    if(length(ref)!=1 | !(ref %in% names(ref.Bcelltypes)) | !(ref %in% names(ref.column))){
      stop("no clear reference selected, for now you need to select one of azimuth.pbmcref, azimuth.tonsilref or azimuth.bonemarrowref or import new B cell clusters names associated with you reference using the ref, ref.column and ref.Bcelltypes arguments.")
    }
    if(!ref.column[[ref]] %in% colnames(db)){
      stop("no clear reference column selected, missing ", ref.column[[ref]]," in provided data frame")
    }
    VDJ_celltypes <- ref.Bcelltypes[[ref]]
  } 
  if (seq_type == "TCR"){
    if(length(ref)!=1 | !(ref %in% names(ref.Tcelltypes)) | !(ref %in% names(ref.column))){
      stop("no clear reference selected, for now you need to select one of azimuth.pbmcref, azimuth.tonsilref or azimuth.bonemarrowref or import new B cell clusters names associated with you reference using the ref, ref.column and ref.Bcelltypes arguments.")
    }
    if(!ref.column[[ref]] %in% colnames(db)){
      stop("no clear reference column selected, missing ", ref.column[[ref]]," in provided data frame")
    }
    VDJ_celltypes <- ref.Tcelltypes[[ref]]
  }
  
  # define chains to be used:
  use_chains <- match.arg(use_chains)
  
  if(!use_chains %in% c("heavy", "all")){
    stop("use_chains must be one of heavy or all.")
  }
  
  # perform first a few check ups:
  if(use_chains == "all") {
    chains <- c("heavy", "light")
  } else {
    chains <- "heavy"
  }

  db <- db %>%
    dplyr::mutate(
      locus_simplified = ifelse(!!rlang::sym(locus) %in% heavy, "heavy", ifelse(!!rlang::sym(locus) %in% light, "light", "other"))
    ) %>%
    dplyr::select(
      -dplyr::any_of(
        c(paste0("is.non",cell_type,"_VDJ_doublet"), 
          paste0("is.non",cell_type,"_VDJ_doublet.confidence"))
        )
      )
  
  # filter on scRNA-seq contigs and prepare the dataset:
  flagged_db <- db %>%
    dplyr::filter(
      locus_simplified %in% chains & !!rlang::sym(assay) %in% scRNAseq.tech
    ) 
  
  other_db <- db %>%
    dplyr::filter(
      !(locus_simplified %in% chains & !!rlang::sym(assay) %in% scRNAseq.tech)
    ) 
  
  if(is.null(split.by)){
    split <- NULL
    flagged_db <- flagged_db %>%
      dplyr::mutate(
        split.by = "all_sequences"
      )
    split.by = "split.by"
  }

  # check if more than one heavy or one light per cell_id:
  for (chain in chains){
    if(any(duplicated(dplyr::filter(flagged_db, locus_simplified == chain)[[cell_id]]))){
      stop("multiple ", paste(unlist(mget(chain)), collapse = "-")," contigs detected from some cells, run resolveMultiContigs() first.")
    }
  }

  # update cutoffs if needed:
  if(variable_cutoff){
    # learn dataset (split.by) and locus specific cutoffs (adapted to library depth):
    # 1st quartile of dominant VDJ umi_count (in 75% of doublets the second VDJ should be above this value) and 1/10 of Median (below could be considered ambient RNA, i.e; 1/10 of an average cell).
    # if homotypicVDJdoublets was run before, we reuse the calculated cutoffs
    if(!(all(c("upper_cutoff", "lower_cutoff") %in% colnames(flagged_db)))){
      flagged_db <- flagged_db %>%
        dplyr::group_by(locus_simplified, !!rlang::sym(split.by)) %>%
        dplyr::mutate(
          upper_cutoff = quantile(!!rlang::sym(umi_count), 0.25, na.rm = TRUE),
          lower_cutoff = ceiling(median(!!rlang::sym(umi_count), na.rm = TRUE)/10)
        ) %>%
        dplyr::ungroup()
    }
  } else { # otherwise, use provided cutoffs
    flagged_db <- flagged_db %>%
      dplyr::mutate(
        upper_cutoff = high_cutoff,
        lower_cutoff = low_cutoff
      )
  }
  
  cutoffs <- flagged_db %>%
    dplyr::group_by(!!rlang::sym(split.by), locus_simplified) %>%
    dplyr::summarise(dplyr::across(
      dplyr::starts_with("upper_cutoff") | dplyr::starts_with("lower_cutoff"),
      ~ max(.x, na.rm = TRUE),
      .names = "{.col}"
    ), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = locus_simplified, values_from = c(upper_cutoff, lower_cutoff)) %>%
    dplyr::rename_with(.fn = ~ paste0(.x, "_max"), .cols = starts_with("upper_cutoff")) %>%
    dplyr::rename_with(.fn = ~ paste0(.x, "_max"), .cols = starts_with("lower_cutoff"))

  # look for non-B cells identified by azimuth and check whether the level of contig umi_count make them likely to be a true doublet:
  non_VDJcelltype_db <- flagged_db %>%
    dplyr::filter(!(!!rlang::sym(ref.column[[ref]]) %in% VDJ_celltypes)) 
  
  non_VDJcelltype_db <- non_VDJcelltype_db %>%
    dplyr::select(!!rlang::sym(cell_id), locus_simplified, !!rlang::sym(umi_count), upper_cutoff, lower_cutoff, !!rlang::sym(split.by)) %>%
    tidyr::pivot_wider(names_from = locus_simplified, 
                       values_from = c(!!rlang::sym(umi_count), upper_cutoff, lower_cutoff),
                       values_fill = setNames(list(0), umi_count)) %>%
    dplyr::left_join(cutoffs, by = split.by) %>%
    dplyr::mutate(across(
      dplyr::starts_with("upper_cutoff") | dplyr::starts_with("lower_cutoff"),
      ~ ifelse(is.na(.x), get(paste0(cur_column(), "_max")), .x)
    )) %>%
    dplyr::select(-ends_with("_max"))
  
  if(!"light" %in% chains){
    # must be high for the heavy chain
    non_VDJcelltype_db <- non_VDJcelltype_db %>%
      dplyr::mutate(
        !!rlang::sym(paste0("is.non", cell_type,"_VDJ_doublet.confidence")) := ifelse(!!rlang::sym(paste0(umi_count, "_heavy")) > upper_cutoff_heavy, "high",
                                                                                      ifelse(!!rlang::sym(paste0(umi_count, "_heavy")) > lower_cutoff_heavy, "low", "ambient RNA")),
        !!rlang::sym(paste0("is.non", cell_type,"_VDJ_doublet")) := (is.nonB_VDJ_doublet.confidence == "high")
      ) %>%
      dplyr::select(
        dplyr::all_of(c(cell_id,
                        paste0("is.non", cell_type,"_VDJ_doublet"),
                        paste0("is.non", cell_type,"_VDJ_doublet.confidence"), 
                        "upper_cutoff_heavy",
                        "lower_cutoff_heavy"))
      ) 
  } else {
    # must be high for at least one of the chain and low for the other
    non_VDJcelltype_db <- non_VDJcelltype_db %>%
      dplyr::mutate(
        !!rlang::sym(paste0("is.non", cell_type,"_VDJ_doublet.confidence")) := ifelse(!!rlang::sym(paste0(umi_count, "_heavy")) > upper_cutoff_heavy & !!rlang::sym(paste0(umi_count, "_light")) > upper_cutoff_light, "high",
                                                                                      ifelse(!!rlang::sym(paste0(umi_count, "_heavy")) > lower_cutoff_heavy & !!rlang::sym(paste0(umi_count, "_light")) > lower_cutoff_light, "low", "ambient RNA")),
        !!rlang::sym(paste0("is.non", cell_type,"_VDJ_doublet")) := (is.nonB_VDJ_doublet.confidence == "high")
      ) %>%
      dplyr::select(
        any_of(c(cell_id,
                 paste0("is.non", cell_type,"_VDJ_doublet"),
                 paste0("is.non", cell_type,"_VDJ_doublet.confidence"),
                 "upper_cutoff_heavy",
                 "upper_cutoff_light",
                 "lower_cutoff_heavy",
                 "lower_cutoff_light"))
      )
  }
  
  VDJcelltype_db <- flagged_db %>%
    dplyr::filter(!!rlang::sym(ref.column[[ref]]) %in% VDJ_celltypes) 
  
  VDJcelltype_db <- VDJcelltype_db %>%
    dplyr::select(!!rlang::sym(cell_id), locus_simplified, !!rlang::sym(umi_count), upper_cutoff, lower_cutoff, !!rlang::sym(split.by)) %>%
    tidyr::pivot_wider(names_from = locus_simplified, values_from = c(!!rlang::sym(umi_count), upper_cutoff, lower_cutoff)) %>%
    dplyr::mutate(
      !!rlang::sym(paste0("is.non", cell_type,"_VDJ_doublet.confidence")) := paste0(cell_type, " cell"),
      !!rlang::sym(paste0("is.non", cell_type,"_VDJ_doublet")) := FALSE
    ) %>%
    dplyr::select(any_of(c(cell_id,
                           paste0("is.non", cell_type,"_VDJ_doublet"), 
                           paste0("is.non", cell_type,"_VDJ_doublet.confidence"))
    ))
  
  flagged_db <- dplyr::bind_rows(non_VDJcelltype_db, VDJcelltype_db)
  
  db <- db %>%
    dplyr::left_join(flagged_db, by = join_by(!!rlang::sym(cell_id)))
  
  if(output){
    # export QC graph (one by part of the dataframe:
    vdjQCplot(db,
              use_chains = use_chains,
              seq_type = seq_type,
              heavy = heavy,
              light = light,
              split.by = split.by,
              type = paste0("non", cell_type,"_vdj_doublet"),
              locus = locus,
              cell_id = "cell_id",
              umi_count = umi_count,
              output_folder = output_folder,
              analysis_name = analysis_name,
              save_plot = save_plot,
              return_plot = FALSE)
  }

  groups <- na.omit(levels(as.factor(db[[split.by]])))
  for(group in groups){
    group_db <- db %>%
      dplyr::filter(!!rlang::sym(split.by) == group, !!rlang::sym(locus) %in% heavy)
    cells_in_group <- nrow(group_db)
    non_VDJ_group_db <- group_db %>%
      dplyr::filter(!(!!rlang::sym(ref.column[[ref]]) %in% VDJ_celltypes))
    non_VDJ_cell_in_group <- nrow(non_VDJ_group_db)
    high_doublets <- nrow(dplyr::filter(non_VDJ_group_db, !!rlang::sym(paste0("is.non", cell_type,"_VDJ_doublet.confidence")) == "high"))
    low_doublets <- nrow(dplyr::filter(non_VDJ_group_db, !!rlang::sym(paste0("is.non", cell_type,"_VDJ_doublet.confidence")) == "low"))
    ambient_RNA <- nrow(dplyr::filter(non_VDJ_group_db, !!rlang::sym(paste0("is.non", cell_type,"_VDJ_doublet.confidence")) == "ambient RNA"))
    upper_group_cutoff_heavy <- group_db[[paste0("upper_cutoff_", paste(heavy, collapse = "-"))]][1]
    upper_group_cutoff_light <- group_db[[paste0("upper_cutoff_", paste(light, collapse = "-"))]][1]
    lower_group_cutoff_heavy <- group_db[[paste0("lower_cutoff_", paste(heavy, collapse = "-"))]][1]
    lower_group_cutoff_light <- group_db[[paste0("lower_cutoff_", paste(light, collapse = "-"))]][1]
    group_log_message <- paste0(group, ": out of ", cells_in_group," cells:\n",
                                non_VDJ_cell_in_group/cells_in_group*100,"% cells (n = ",non_VDJ_cell_in_group,") identified as non-",cell_type, " cells\n",
                                high_doublets/cells_in_group*100,"% cells (n = ",high_doublets,") identified as high probability non-",cell_type ,"/",cell_type," doublets, using ", paste(unlist(mget(chains)), collapse = "-")," contigs;\n",
                                low_doublets/cells_in_group*100,"% cells (n = ",low_doublets,") identified as low probability non-",cell_type ,"/",cell_type," doublets; \n",
                                ambient_RNA/cells_in_group*100,"% cells (n = ",ambient_RNA,") identified as ambient RNA contamination.\n",
                                "upper cutoff used:", upper_group_cutoff_heavy, " for ", paste(heavy, collapse = "-"), " and ",upper_group_cutoff_light, " for ", paste(light, collapse = "-"),"\n",
                                "lower cutoff used:", lower_group_cutoff_heavy, " for ", paste(heavy, collapse = "-"), " and ",lower_group_cutoff_light, " for ", paste(light, collapse = "-"),"\n","\n")
    if(verbose){cat(group_log_message)}
  }
  
  db <- db %>%
    dplyr::select(
      -dplyr::any_of(c("locus_simplified",
                       "upper_cutoff_heavy",
                       "upper_cutoff_light",
                       "lower_cutoff_heavy",
                       "lower_cutoff_light"))
    )
  
  if(nrow(other_db)>0){
    #adding back non BCR or TCR contigs depending on the seq_type
    db <- db %>%
      dplyr::bind_rows(other_db)
  }

  return(db)
}

#### Wrapper function to flag VDJ doublets in single cell data ####
#' calls successively flagBdoublets and flagNonBdoublets
#'
#' \code{flagVDJdoublets} flag potential B and nonB doublets.
#'
#' @param db            a AIRR formatted dataframe containing heavy and light chain sequences, cell_ids and umi_counts.
#' @param analysis_name name to use for outputs prefixes.
#' @param homotypic     whether to run homotypicVDJdoublets()
#' @param heterotypic   whether to run heterotypicVDJdoublets()
#' @param seq_type        type of VDJ sequence ("Ig" or "TCR" to match igblastb requirements)
#' @param heavy         heavy chains to be used. if set to NULL, will default to IGH for Ig and TRB/TRD for TCR.
#' @param light         light chains to be used. if set to NULL, will default to IGL/IGK for Ig and TRA/TRG for TCR.
#' @param output        whether to output graphs with umi_counts for dominant versus second IGH VDJ contig and the recap excel workbook. Can be set to none, partial (only pdf + log file with QC parameters) or total (pdf + log file with QC parameters + tsv file with sequences).
#' @param output_folder name of the folder in which graph and recap excel workbooks will be saved [default = "VDJ_QC"].
#' @param split.by      which column to use to group cells when learning dataset-specific distributions
#' @param save_plot     format to save the final ggplot2 object. "pdf" or "png". if set to any other value, no plot will be saved.
#' @param verbose       whether to print recap informations to the console
#' @param ...           arguments to pass to homotypicVDJdoublets or heterotypicVDJdoublets.
#'
#' @export

flagVDJdoublets <- function(db,
                            analysis_name = "All_sequences",
                            homotypic = TRUE,
                            heterotypic = TRUE,
                            seq_type = c("Ig", "TCR"),
                            heavy = NULL,
                            light = NULL,
                            output = TRUE,
                            output_folder = "VDJ_QC",
                            split.by = NULL,
                            ref = c("azimuth.pbmcref", "azimuth.tonsilref", "azimuth.bonemarrowref"),
                            ref.column = list("azimuth.pbmcref" = "predicted.celltype.l2", 
                                              "azimuth.tonsilref" = "predicted.cell_type.l1", 
                                              "azimuth.bonemarrowref" = "predicted.celltype.l2"),
                            ref.Bcelltypes = list("azimuth.pbmcref" = c("B intermediate", "B memory", "B naive", "Plasmablast"),
                                                  "azimuth.tonsilref" = c("B activated", "B memory", "B naive", "preGCB", "PB", "PC", "PC/doublet", "preMBC/doublet", "prePB", "Cycling DZ GCB","DZ GCB","DZtoLZ GCB transition","FCRL4/5+ B memory","LZ GCB","LZtoDZ GCB transition"),
                                                  "azimuth.bonemarrowref" = c("Memory B", "Naive B", "Plasma", "pro B", "pre B", "transitional B")),
                            ref.Tcelltypes = list("azimuth.pbmcref" = c("CD4 CTL", "CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM", "Treg", "CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM", "dnT", "gdT", "MAIT"),
                                                  "azimuth.tonsilref" = c("CD4 naive", "CD4 Non-TFH", "CD4 TCM", "CD4 TFH", "CD4 TFH Mem", "CD4 TREG", "CD8 naive", "CD8 T", "CD8 TCM", "Cycling T","dnT","MAIT/TRDV2+ gdT","non-TRDV2+ gdT"),
                                                  "azimuth.bonemarrowref" = c("CD4 Effector", "CD4 Memory", "CD4 Naive", "CD8 Effector 1", "CD8 Effector 2", "CD8 Effector 3", "CD8 Memory", "CD8 Naive", "MAIT", "T proliferating")),
                            save_plot = c("pdf", "png"),
                            save_tsv = TRUE,
                            verbose = TRUE,
                            ...){
  
  save_plot <- match.arg(save_plot)
  
  seq_type <- match.arg(seq_type)
  
  if(is.null(heavy)){
    if(seq_type == "Ig"){
      heavy <- "IGH"
    }
    if(seq_type == "TCR"){
      heavy <- c("TRB","TRD")
    }
  }
  if(is.null(light)){
    if(seq_type == "Ig"){
      light <- c("IGK","IGL")
    }
    if(seq_type == "TCR"){
      light <- c("TRA","TRG")
    }
  }
  
  if(!is.null(output_folder)){
    if(!stringr::str_ends(output_folder, "/")){output_folder = paste0(output_folder, "/")}
    if(!dir.exists(output_folder)){
      dir.create(output_folder)
    }
  }
  
  log_file <- paste0(output_folder, analysis_name, "_flagVDJdoublets_logfile.txt")
  open_mode = "wt"
  
  if(homotypic){
    time_and_log({
      db <- homotypicVDJdoublets(db,
                                 seq_type = seq_type,
                                 heavy = heavy,
                                 light = light,
                                 analysis_name = analysis_name,
                                 split.by = split.by,
                                 output = output,
                                 output_folder = output_folder,
                                 save_plot = save_plot,
                                 verbose = verbose,
                                 ...)
    }, verbose = FALSE, time = TRUE, log_file = log_file, log_title = "homotypic VDJ doublets", open_mode = open_mode)
    open_mode = "a"
  }

  if(heterotypic){
    time_and_log({
      db <- heterotypicVDJdoublets(db,
                                   seq_type = seq_type,
                                   heavy = heavy,
                                   light = light,
                                   analysis_name = analysis_name,
                                   split.by = split.by,
                                   output = output,
                                   output_folder = output_folder,
                                   ref = ref,
                                   ref.column = ref.column,
                                   ref.Bcelltypes = ref.Bcelltypes, 
                                   ref.Tcelltypes = ref.Tcelltypes,
                                   verbose = verbose,
                                   ...)
    }, verbose = FALSE, time = TRUE, log_file = log_file, log_title = "heterotypic VDJ doublets", open_mode = open_mode)
  }
  
  time_and_log({
    cat("\n")
    print(sessionInfo())
  }, verbose = FALSE, time = FALSE, log_file = log_file, log_title = "session info", open_mode = open_mode)
    
  return(db)
}


#### Function to run Changeo igblast pipeline on VDJ sequences from sanger, BD or 10X datasets ####
#' calls Immcantation (ChangeO) AssignGenes.py and MakeDb.py scripts from within R and update some columns names depending on selected technology
#'
#' \code{runAssignGenes} old pipeline to run igblast
#'
#' @param db            a data frame containing at least a sequence and a sequence_id column.
#' @param sequence      name of the column containing the original sequence.
#' @param sequence_id   name of the column containing sequence identifier.
#' @param seq_type      type of VDJ sequence ("Ig" or "TCR" to match igblastb requirements)
#' @param igblast_dir   path to igblast database [default = path suggested on installation: https://changeo.readthedocs.io/en/stable/examples/igblast.html]
#' @param imgt_dir      path to imgt-gapped database [default = path suggested on installation: https://changeo.readthedocs.io/en/stable/examples/igblast.html]
#' @param output        whether to output graphs with umi_counts for dominant versus second IGH VDJ contig and the recap excel workbook. If set to FALSE, only the corrected database is returned.
#' @param output_folder name of the folder in which graph and recap excel workbooks will be saved [default = "igblast_results/"].
#' @param log           whether to create a log file
#' @param log_file      name of the log file
#'
#' @return   a named list of AIRR formatted data frame corresponding to the pass and fail results of Igblast reformatted using Makedb and all additional columns from the submitted data frame.
#'
#' @details
#' 1. will first run igblast
#' 2. import igblast results and reformat data frame to return all metadata from the original dataframe
#' dependencies: standalone igblast; Changeo (AssignGenes.py, MakeDb.py); dplyr; readr; Biostrings
#'
#' @export
#'
#' @import dplyr
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings writeXStringSet

runAssignGenes <- function(db,
                           type = c("nt", "aa"),
                           seq_type = c("Ig", "TCR"),
                           organism = c("human", "mouse", "rabbit", "rat", "rhesus_monkey"),
                           igblast_dir = "~/share/igblast/",
                           imgt_dir = "~/share/germlines/imgt/",
                           sequence = "sequence", sequence_id = "sequence_id",
                           output = FALSE, output_folder = "igblast_results",
                           log = FALSE,
                           log_file= NULL){


  #suppressMessages(library(dplyr))

  type <- match.arg(type)
  seq_type <- match.arg(seq_type)
  organism <- match.arg(organism)
  igblast_dir <- ifelse(stringr::str_ends(igblast_dir,"/"), igblast_dir, paste0(igblast_dir, "/"))
  imgt_dir <- ifelse(stringr::str_ends(imgt_dir,"/"), imgt_dir, paste0(imgt_dir, "/"))

  if(!output){
    output_folder = "temp/"
  } else {
    # check if output folder is properly formatted
    if(!stringr::str_ends(output_folder, "/")){output_folder = paste0(output_folder, "/")}
  }

  if(!dir.exists(output_folder)){
    dir.create(output_folder)
  }

  if(sequence != "sequence"){ # make sure columns are renamed properly to adapt to igblast formatting; old column name is kept
    db$sequence = db[[sequence]]
  }

  if(sequence_id != "sequence_id"){ # make sure columns are renamed properly to adapt to igblast formatting; old column name is kept
    db$sequence_id = db[[sequence_id]]
    if(any(duplicated(db$sequence_id))){stop("the following sequence_id are duplicated: ", paste(db$sequence_id[duplicated(db$sequence_id)], collapse = "; "))}
  }

  # export fasta file:
  seq <- db$sequence
  names(seq) <- db$sequence_id
  dna = Biostrings::DNAStringSet(seq)
  Biostrings::writeXStringSet(dna, paste0(output_folder, "All_seq.fasta"))

  loci <- list("Ig" = "ig", "TCR" = "tr")
  loci <- loci[[seq_type]]

  igblast_version <- list("nt" = "igblast", "aa" = "igblast-aa")
  igblast_version <- igblast_version[[type]]

  # Capture IMGT version
  if(file.exists(paste0(igblast_dir, "/IMGT_version.txt"))){
    database_info <- paste0("(",readLines(paste0(igblast_dir, "/IMGT_version.txt")), ")")
  } else {database_info <- NULL}
  log_message <- c(
    paste0("Database: IMGT ", database_info),
    paste0("Organism: ", organism),
    paste0("Loci: ", seq_type)
  )
  cat(log_message, sep = "\n")

  # run igblast:
  fasta_file <- paste0(output_folder, "All_seq.fasta")
  fasta_file <- gsub(" ", "\\\ ", fasta_file, fixed = TRUE) #to remove any blank in file_path
  #system(paste0("AssignGenes.py igblast -s ", fasta_file, " -b ", igblast_dir ," --organism human --loci ig --format blast"))
  messages_igblast <- system2("AssignGenes.py", args = paste0(igblast_version," -s ", fasta_file, " -b ", igblast_dir ," --organism ", tolower(organism)," --loci ", loci," --format blast"), stdout = TRUE, stderr = TRUE)
  cat(paste(messages_igblast, collapse = "\n"))
  if(log){
    if(!is.null(log_file)){cat(paste(messages_igblast, collapse = "\n"), file = log_file)}
  }
  #system(paste0("AssignGenes.py igblast -s ", fasta_file, " -b ", igblast_dir ," --organism human --loci ig --format blast --nproc ", nproc))

  igblast_output_file <- paste0(output_folder, "All_seq_igblast")
  igblast_output_file <- gsub(" ", "\\\ ", igblast_output_file, fixed = TRUE)
  reference_dir <- paste0(imgt_dir, organism, "/vdj/imgt_", organism, "_*.fasta")
  messages_makedb <- system2("MakeDb.py", args = paste0("igblast -i ", igblast_output_file, ".fmt7 -s ", fasta_file, " -r ", reference_dir, " --extended --failed"), stdout = TRUE, stderr = TRUE)
  cat(paste(messages_makedb, collapse = "\n"))
  if(log){
    if(!is.null(log_file)){cat(paste(messages_makedb, collapse = "\n"), file = log_file)}
  }

  # restore original metadata:
  pass_true <- file.exists(paste0(igblast_output_file, "_db-pass.tsv"))
  fail_true <- file.exists(paste0(igblast_output_file, "_db-fail.tsv"))

  if(pass_true){
    VDJ_db <- readr::read_tsv(file = paste0(igblast_output_file, "_db-pass.tsv"), show_col_types = FALSE)
    diff_cols <- setdiff(colnames(db), colnames(VDJ_db))
    if(length(diff_cols)>0){
      db <- db %>%
        dplyr::select(all_of(c("sequence_id", diff_cols)))
      VDJ_db <- suppressMessages(dplyr::left_join(VDJ_db, db))
    }
    if(fail_true){
      failed_VDJ_db <- readr::read_tsv(file = paste0(igblast_output_file, "_db-fail.tsv"), show_col_types = FALSE)
      if(length(diff_cols)>0){
        failed_VDJ_db <- suppressMessages(dplyr::left_join(failed_VDJ_db, db))
      }
    } else {#100% success!!
      failed_VDJ_db <- VDJ_db[0,]
    }
  } else {#100% failure :-(((
    failed_VDJ_db <- readr::read_tsv(file = paste0(igblast_output_file, "_db-fail.tsv"), show_col_types = FALSE)
    diff_cols <- setdiff(colnames(db), colnames(failed_VDJ_db))
    if(length(diff_cols)>0){
      db <- db %>%
        dplyr::select(all_of(c("sequence_id", diff_cols)))
      failed_VDJ_db <- suppressMessages(dplyr::left_join(failed_VDJ_db, db))
    }
    VDJ_db <- failed_VDJ_db[0,]
  }

  results <- list(pass = as.data.frame(VDJ_db),
                  fail = as.data.frame(failed_VDJ_db))

  if(!output){unlink(output_folder, recursive = TRUE)}

  return(results)
}

#### Function to directly run igblast on VDJ sequences from sanger, BD or 10X datasets ####
#' directly calls igblastn from within R and outputs an airr formatted table
#'
#' \code{runIgblastn} run igblast
#'
#' @param db            a data frame containing at least a sequence and a sequence_id column.
#' @param sequence      name of the column containing the original sequence.
#' @param sequence_id   name of the column containing sequence identifier.
#' @param seq_type      type of VDJ sequence ("Ig" or "TCR" to match igblastb requirements)
#' @param organism      organism (any of "human", "mouse", "rhesus_monkey; for other see https://changeo.readthedocs.io/en/stable/examples/igblast.html)
#' @param igblast_dir   path to igblast database [default = path suggested on installation: https://changeo.readthedocs.io/en/stable/examples/igblast.html]
#' @param imgt_dir      path to imgt-gapped database [default = path suggested on installation: https://changeo.readthedocs.io/en/stable/examples/igblast.html]
#' @param output        whether to output graphs with umi_counts for dominant versus second IGH VDJ contig and the recap excel workbook. If set to FALSE, only the corrected database is returned.
#' @param output_folder name of the folder in which graph and recap excel workbooks will be saved [default = "igblast_results/"].
#' @param log           whether to create a log file
#' @param log_file      name for the log file
#' @param nproc_igblast number of proc to use for igblast
#'
#' @return   a named list of AIRR formatted data frame corresponding to the pass and fail results of Igblast reformatted using Makedb and all additional columns from the submitted data frame.
#'
#' @details
#' 1. will first run igblast and request an airr format output
#' 2. import igblast results and add all metadata from the original dataframe
#' dependencies: standalone igblastn; dplyr; readr; Biostrings; parallel,
#' @export
#'
#' @import dplyr
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings writeXStringSet
#' @importFrom Biostrings readDNAStringSet

runIgblastn <- function(db,
                       sequence = "sequence",
                       sequence_id = "sequence_id",
                       seq_type = c("Ig", "TCR"),
                       organism = c("human", "mouse", "rabbit", "rat", "rhesus_monkey"),
                       igblast_dir = "~/share/igblast/",
                       imgt_dir = "~/share/germlines/imgt/",
                       gapV = TRUE,
                       output = FALSE,
                       output_folder = "igblast_results",
                       log = FALSE,
                       log_file= NULL,
                       nproc_igblast = 10){

  #suppressMessages(library(dplyr))

  seq_type <- match.arg(seq_type)
  organism <- match.arg(organism)
  igblast_dir <- ifelse(stringr::str_ends(igblast_dir,"/"), igblast_dir, paste0(igblast_dir, "/"))
  imgt_dir <- ifelse(stringr::str_ends(imgt_dir,"/"), imgt_dir, paste0(imgt_dir, "/"))

  nproc_igblast <- min(nproc_igblast, parallel::detectCores())

  if(!output){
    output_folder = "temp/"
  } else {
    # check if output folder is properly formatted
    if(!stringr::str_ends(output_folder, "/")){output_folder = paste0(output_folder, "/")}
  }

  if(!dir.exists(output_folder)){
    dir.create(output_folder)
  }

  if(sequence != "sequence"){ # make sure columns are renamed properly to adapt to igblast formatting; old column name is kept
    db$sequence = db[[sequence]]
  }

  if(sequence_id != "sequence_id"){ # make sure columns are renamed properly to adapt to igblast formatting; old column name is kept
    db$sequence_id = db[[sequence_id]]
    if(any(duplicated(db$sequence_id))){stop("the following sequence_id are duplicated: ", paste(db$sequence_id[duplicated(db$sequence_id)], collapse = "; "))}
  }

  # export fasta file:
  seq <- db$sequence
  names(seq) <- db$sequence_id
  dna = Biostrings::DNAStringSet(seq)
  Biostrings::writeXStringSet(dna, paste0(output_folder, "All_seq.fasta"))

  # run igblast:
  fasta_file <- paste0(output_folder, "All_seq.fasta")
  fasta_file <- gsub(" ", "\\\ ", fasta_file, fixed = TRUE) # to remove any blank in file_path

  out_file <- gsub(".fasta", "_igblast.tsv", fasta_file)

  # Define the environment variable
  env_vars <- paste0("IGDATA=", igblast_dir)

  # Define the command
  command <- "igblastn"

  # Capture IMGT version
  version_info <- system2(command, args = "-version", stdout = TRUE)
  if(file.exists(paste0(igblast_dir, "/IMGT_version.txt"))){
    database_info <- paste0("(",readLines(paste0(igblast_dir, "/IMGT_version.txt")), ")")
  } else {database_info <- NULL}

  log_message <- c(
    paste0("Running ", version_info[1]),
    paste0("Organism: ", organism),
    paste0("Loci: ", seq_type),
    paste0("Database: IMGT ", database_info),
    paste0("Output format: airr"),
    paste0("Nproc: ", nproc_igblast)
    )
  cat(log_message, sep = "\n")

  database_type <- list("Ig" = "ig", "TCR" = "tr")

  # Define the arguments
  args <- c(
    "-germline_db_V", paste0(igblast_dir, "database/imgt_",organism,"_",database_type[[seq_type]],"_v"),
    "-germline_db_D", paste0(igblast_dir, "database/imgt_",organism,"_",database_type[[seq_type]],"_d"),
    "-germline_db_J", paste0(igblast_dir, "database/imgt_",organism,"_",database_type[[seq_type]],"_j"),
    "-c_region_db", paste0(igblast_dir, "database/imgt_",organism,"_",database_type[[seq_type]],"_c"),
    "-V_penalty", -1,
    "-D_penalty", -2,
    "-J_penalty", -2,
    "-min_D_match", 5,
    "-auxiliary_data", paste0(igblast_dir, "/optional_file/",organism,"_gl.aux"),
    "-domain_system", "imgt",
    "-ig_seqtype", seq_type,
    "-organism", organism,
    "-outfmt", 19,
    "-query", fasta_file,
    "-out", out_file,
    "-num_threads", nproc_igblast
  )

  # Run igblast
  system2(command, args = args, env = env_vars)

  # restore original metadata:
  if(file.exists(out_file)){
    igblast_output <- readr::read_tsv(file = out_file, show_col_types = FALSE)

    diff_cols <- setdiff(colnames(db), colnames(igblast_output))

    if(length(diff_cols)>0){
        igblast_output <- igblast_output %>%
          dplyr::left_join(
            dplyr::select(db, all_of(c("sequence_id", diff_cols))),
            by = join_by("sequence_id")
          )
    }

    if(gapV){
      fasta_files <- list.files(path = paste0(imgt_dir, organism,"/vdj/"), pattern = paste0("^imgt_",organism,"_",toupper(seq_type),".*V\\.fasta$"), full.names = TRUE)
      fasta_list <- lapply(fasta_files, Biostrings::readDNAStringSet)
      references <- c(fasta_list[[1]], fasta_list[[2]], fasta_list[[3]])
      names(references) <- stringr::str_split(names(references), "\\|", simplify = TRUE)[, 2]

      igblast_output <- igblast_output %>%
        rowwise() %>%
        mutate(tmp = list(gapV(seq = sequence_alignment, v_germ_start = v_germline_start, v_germ_end = v_germline_end, v_call = v_call, references = references))) %>%
        unnest_wider(tmp) %>%
        mutate(
          sequence_alignment = gapped_sequence,
          v_germline_start = gapped_v_germ_start,
          v_germline_end = gapped_v_germ_length
        ) %>%
        select(-gapped_sequence, -gapped_v_germ_start, -gapped_v_germ_length) %>%
        ungroup()
    }

    db_pass <- igblast_output %>%
      dplyr::filter(
        !stop_codon,
        vj_in_frame,
        !v_frameshift,
        !is.na(junction))
    db_fail <- igblast_output %>%
      dplyr::filter(
        stop_codon | !vj_in_frame | v_frameshift | is.na(stop_codon) | is.na(vj_in_frame) | is.na(v_frameshift) | is.na(junction)
        )

    results <- list(pass = as.data.frame(db_pass),
                    fail = as.data.frame(db_fail))

    log <- c(
      paste0("Total sequences : ", nrow(db)),
      paste0("Pass : ", nrow(db_pass)),
      paste0("Fail : ", nrow(db_fail))
    )
    cat(log, sep = "\n")

    if(!output){unlink(output_folder, recursive = TRUE)}
    return(results)

  } else {
    if(!output){unlink(output_folder, recursive = TRUE)}
    stop("issue running igblast, no output file")
  }
  #TODO add germline_aligments, imgt gaps, v_germline_length,	d_germline_length,	j_germline_length,	germline_alignment_d_mask to bypass Dowser::createGermline()
}

#### Function to run introduce IMGT style gap in sequence (similar to MakeDb from Changeo) ####
#' introduce IMGT gap in V sequence
#'
#' \code{gapV} introduce IMGT gap in V sequence
#'
#' @param seq           a sequence
#' @param v_germ_start  v_germline_start output from igblast
#' @param v_germ_end    v_germline_end output from igblast
#' @param v_call        v_call output from igblast
#' @param references    list of gapped reference sequence named according to possible v_calls
#'
#' @keywords internal

gapV <- function(seq,
                 v_germ_start,
                 v_germ_end,
                 v_call,
                 references) {

  v_germ_length = v_germ_end - v_germ_start + 1

  # Initialize gapped sequence with leading dots (IMGTV-style gapping)
  seq_imgt <- paste0(strrep(".", as.integer(v_germ_start) - 1), seq)

  vgene <- strsplit(v_call, ",")[[1]][1]

  # Look up the gapped reference
  if (!vgene %in% names(references)) {
    stop(sprintf("%s was not found in the germline repository.", vgene))
  }
  vgap <- references[[vgene]]

  # Iterate over gaps in germline (".") to insert them into query
  gap_positions <- gregexpr("\\.", vgap)[[1]]
  gap_positions <- gap_positions[gap_positions != -1]

  gapcount <- as.integer(v_germ_start) - 1
  for (i in gap_positions) {
    if (i >= v_germ_end + gapcount) {
      break
    }
    # Insert a dot at position `i` (R is 1-indexed)
    seq_imgt <- paste0(
      substr(seq_imgt, 1, i - 1),
      ".",
      substr(seq_imgt, i, nchar(seq_imgt))
    )
    gapcount <- gapcount + 1
  }

  # Populate return list
  imgt_gapped <- tibble(
    gapped_sequence =  seq_imgt,
    gapped_v_germ_start = 1,
    gapped_v_germ_length = v_germ_end + gapcount
  )
  return(imgt_gapped)
}

#### Function to run Blastn on VDJC sequences to infer constant region ####
#' call C regions from VDJC sequences
#'
#' \code{runBlastnC} run Blastn on a Blastable IMGT C region datbase return a corrected c_call
#'
#' @param db            a AIRR formated data frame containing at least a sequence and a sequence_id column.
#' @param seq_type      type of VDJ sequence ("Ig" or "TCR" to match igblastb requirements)
#' @param organism      organism (any of "human", "mouse", "rhesus_monkey; for other see https://changeo.readthedocs.io/en/stable/examples/igblast.html)
#' @param igblast_dir   path to igblast database [default = path suggested on installation: https://changeo.readthedocs.io/en/stable/examples/igblast.html]
#' @param sequence      name of the column containing the original sequence.
#' @param sequence_id   name of the column containing sequence identifier.
#' @param output        whether to output graphs with umi_counts for dominant versus second IGH VDJ contig and the recap excel workbook. If set to FALSE, only the corrected database is returned.
#' @param output_folder name of the folder in which graph and recap excel workbooks will be saved [default = "VDJ_QC"].
#' @param c_call        name of the column containing c_call.
#'
#' @return   the inputted dataframe with an additional c_call column after the j_call column.
#'
#' @details
#' 1. will first run Blastn
#' 2. import Blastn results and select the alignment(s) with the highest scores
#' 3. in case of multiple gene calls with identical score, will return an agregated call (e.g. (IGHG1|IGHG2))
#'
#' @export
#'
#' @import dplyr
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings writeXStringSet
#' @importFrom readr read_csv

runBlastnC <- function(db,
                       seq_type = c("Ig", "TCR"),
                       organism = c("human", "mouse", "rabbit", "rat", "rhesus_monkey"),
                       igblast_dir = "~/share/igblast",
                       sequence = "sequence",
                       sequence_id = "sequence_id",
                       c_call = "c_call",
                       output = FALSE,
                       output_folder = "C_results",
                       nproc_blast = 10){

  #suppressMessages(library(dplyr))

  seq_type <- match.arg(seq_type)
  organism <- match.arg(organism)
  igblast_dir <- ifelse(stringr::str_ends(igblast_dir,"/"), igblast_dir, paste0(igblast_dir, "/"))

  nproc_blast <- min(nproc_blast, parallel::detectCores())

  if(!output){
    output_folder = "temp/"
  } else {
    # check if output folder is properly formatted
    if(!stringr::str_ends(output_folder, "/")){output_folder = paste0(output_folder, "/")}
  }

  if(!dir.exists(output_folder)){
    dir.create(output_folder)
  }

  # export fasta file:
  seq <- db[[sequence]]
  names(seq) <- db[[sequence_id]]
  dna = Biostrings::DNAStringSet(seq)
  fasta_file <- paste0(output_folder, "All_seq.fasta")
  Biostrings::writeXStringSet(dna, fasta_file)

  # Define the command
  command <- "blastn"

  # Define the arguments
  database_type <- list("Ig" = "ig", "TCR" = "tr")

  out_file <- paste0(output_folder, "All_seq_Ig_C_genes_results.csv")

  args <- c(
    "-db", paste0(igblast_dir, "database/imgt_",organism,"_",database_type[[seq_type]],"_c"),
    "-outfmt", "'10 qseqid sseqid score ppos qlen qstart sstart length'",
    "-max_target_seqs", 10,
    "-query", fasta_file,
    "-out", out_file,
    "-task", "blastn",
    "-num_threads", nproc_blast
  )

  # Run blastn:
  system2(command, args = args)

  # import and summarize results:
  c_results <- readr::read_csv(paste0(output_folder, "All_seq_Ig_C_genes_results.csv"),
                               col_names = c("sequence_id", "c_call_top_match", "c_call_alignment_score", "p_pos", "sequence_length", "sequence_align_start", "IGH_align_start", "c_call_alignment_length"), show_col_types = FALSE)
  c_results_filtered <- c_results %>%
    dplyr::filter(IGH_align_start < 3 | sequence_align_start < 3) # allowing a little bit of flexibility. We have one example an IGH_align_start at 2 for a full IGHA1 clone.

  if(nrow(c_results_filtered)>0){
    c_results_filtered <- c_results_filtered %>%
      dplyr::mutate(c_call_top_match = stringr::str_extract(c_call_top_match, "^[^*]+")) %>%  # Extract everything before the first "*")
      dplyr::group_by(sequence_id) %>%
      dplyr::filter(c_call_alignment_score == max(c_call_alignment_score) & c_call_alignment_length > 25) %>%
      dplyr::summarize(
        !!rlang::sym(c_call) := paste(unique(c_call_top_match), collapse = "|"),  # Aggregate unique values of c_call_top_match
        c_call_alignment_score = max(c_call_alignment_score),  # Retain the max c_call_alignment_score score (should be the same across the group)
        c_call_pct_match = max(p_pos),  # Retain the max p_pos score (should be the same across the group)
        c_call_alignment_length = max(c_call_alignment_length), # Retain the max c_call_alignment_length (should be the same across the group - checked)
        .groups = 'drop'  # Drop the grouping after summarization
      )
    replaced_columns <- intersect(colnames(db), c(c_call, "c_call_alignment_score", "c_call_pct_match", "c_call_alignment_length"))
    if(length(replaced_columns)>0){
      for(replaced_column in replaced_columns){
        db[[replaced_column]] <- NULL
      }
    }
    db <- suppressMessages(dplyr::left_join(db, c_results_filtered))
  } else {
    db[[!!rlang::sym(c_call)]] = NA
    db$c_call_alignment_score = NA
    db$c_call_pct_match = NA
    db$c_call_alignment_length = NA
  }

  if(!output){unlink(output_folder, recursive = TRUE)}

  return(db)
}


#### Function to reconstruct a full VDJ sequence (nt and AA) as well as a Fab (AA) ####
#'
#' \code{reconstructFullVDJ} reconstruct full VDJ sequence for AlphaFold modeling
#'
#' @param db          an AIRR formatted dataframe containing heavy and light chain sequences. Should contain only one heavy chain (IGH) per cell_id, if not run resolveMultiHC() first.
#' @param seq_type      type of VDJ sequence ("Ig" or "TCR" to match igblastb requirements)
#' @param organism      organism (any of "human", "mouse", "rhesus_monkey; for other see https://changeo.readthedocs.io/en/stable/examples/igblast.html)
#' @param igblast_dir   path to igblast database [default = path suggested on installation: https://changeo.readthedocs.io/en/stable/examples/igblast.html]
#' @param CH1_AA      named list of CH1 domain for IGH, IGL and IGK c chains (AA); a default one is provided - imported from Uniprot on Feb 2025.
#'
#' @return    an AIRR formatted dataframe containing the following additional collumns:
#' "missing_v_bp", "missing_j_bp", "full_sequence", "full_sequence_aa", "full_sequence_fab_aa", "comments"
#'
#' @details    Will first try to reconstruct the full VDJ sequence,
#' adding missing V and J parts from germline (if less than 27bp, otherwise return error message in comment),
#' reverting N to germline nucleotide, translating to AA and adding the appropriate CH1 domain based on c_call.
#' Currently limited to human Ig sequences.
#'
#' @export
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings translate
#' @importFrom Biostrings DNAStringSet

reconstructFullVDJ <- function(db,
                               igblast_dir = "~/share/igblast",
                               CH1_AA = list(IGHG1 = "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV",
                                             IGHG2 = "ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSNFGTQTYTCNVDHKPSNTKVDKTV",
                                             IGHG3 = "ASTKGPSVFPLAPCSRSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYTCNVNHKPSNTKVDKRV",
                                             IGHG4 = "ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTKTYTCNVDHKPSNTKVDKRV",
                                             IGHA1 = "ASPTSPKVFPLSLCSTQPDGNVVIACLVQGFFPQEPLSVTWSESGQGVTARNFPPSQDASGDLYTTSSQLTLPATQCLAGKSVTCHVKHYTNPSQDVT",
                                             IGHA2 = "ASPTSPKVFPLSLDSTPQDGNVVVACLVQGFFPQEPLSVTWSESGQNVTARNFPPSQDASGDLYTTSSQLTLPATQCPDGKSVTCHVKHYTNSSQDVT",
                                             IGHM = "GSASAPTLFPLVSCENSPSDTSSVAVGCLAQDFLPDSITFSWKYKNNSDISSTRGFPSVLRGGKYAATSQVLLPSKDVMQGTDEHVVCKVQHPNGNKEKNVPLPV",
                                             IGHD = "APTKAPDVFPIISGCRHPKDNSPVVLACLITGYHPTSVTVTWYMGTQSQPQRTFPEIQRRDSYYMTSSQLSTPLQQWRQGEYKCVVQHTASKSKKEIF",
                                             IGHE = "ASTQSPSVFPLTRCCKNIPSNATSVTLGCLATGYFPEPVMVTWDTGSLNGTTMTLPATTLTLSGHYATISLLTVSGAWAKQMFTCRVAHTPSSTDWVDNKTFS",
                                             IGKC = "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC",
                                             IGLC1 = "GQPKGQPKANPTVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADGSPVKAGVETTKPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS",
                                             IGLC2 = "GQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS",
                                             IGLC3 = "GQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHKSYSCQVTHEGSTVEKTVAPTECS",
                                             IGLC6 = "GQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVKVAWKADGSPVNTGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPAECS",
                                             IGLC7 = "GQPKAAPSVTLFPPSSEELQANKATLVCLVSDFNPGAVTVAWKADGSPVKVGVETTKPSKQSNNKYAASSYLSLTPEQWKSHRSYSCRVTHEGSTVEKTVAPAECS")){

  igblast_dir <- ifelse(stringr::str_ends(igblast_dir,"/"), igblast_dir, paste0(igblast_dir, "/"))

  IMGT_j <- Biostrings::readDNAStringSet(paste0(igblast_dir, "fasta/imgt_human_ig_j.fasta"))
  #TODO add possibility to do that for mouse Ig? maybe getting the nt sequenec from imgt database (issue = multiple version possible)

  db$missing_v_bp <- NA
  db$missing_j_bp <- NA
  db$full_sequence <- NA
  db$full_sequence_aa <- NA
  db$full_sequence_fab_aa <- NA
  db$comments <- NA
  for(i in seq_along(db$sequence_id)){
    ##heavy chain:
    if(!is.na(db$v_call[i])){
      db$missing_v_bp[i] <- gregexpr("G|C|A|T", db$sequence_alignment[i])[[1]][1]-1
      db$missing_j_bp[i] <- length(IMGT_j[strsplit(db$j_call[i], ",")[[1]][1]][[1]]) - db$j_germline_end[i]
      if (db$missing_v_bp[i]==0){
        db$full_sequence[i] <- substring(db$sequence[i], db$v_sequence_start[i], db$j_sequence_end[i])
      } else {
        if(db$missing_v_bp[i] < 27){
          db$full_sequence[i] <- paste0(substring(db$germline_alignment[i], 1, db$missing_v_bp[i]),
                                        substring(db$sequence[i], db$v_sequence_start[i], db$j_sequence_end[i]))
          db$comments[i] <- "missing bp in VH sequence, reverted to germline"
        } else {
          db$comments[i] <- "missing 27 or more bp in VH sequence, check manually"
        }
      }
      # replace N with corresponding germline bp and add comment
      # needed to translate to aa
      Ns_aligned <- gregexpr("N", db$sequence_alignment[i])[[1]]
      Ns_aligned <- Ns_aligned[!Ns_aligned == -1]
      Ns_raw <- gregexpr("N", db$full_sequence[i])[[1]]
      Ns_raw <- Ns_raw[!Ns_raw == -1]
      if(length(Ns_aligned) > 0 & !is.na(db$full_sequence[i]) & length(Ns_aligned) == length(Ns_raw)){
        for (j in seq_along(Ns_aligned)){
          # use the fact that order of Ns is the same between full sequence (trimmed) and aligned sequences even if positions are different
          substr(db$full_sequence[i], Ns_raw[j], Ns_raw[j]) <- substr(db$germline_alignment[i], Ns_aligned[j], Ns_aligned[j])
        }
        db$comments[i] <- paste0(db$comments[i], "; 'N' in VH sequence, reverted to germline")
      }
      if(!length(Ns_aligned) == length(Ns_raw)){
        # rare cases when an N base is wrongly inserted in the read sequence and doesn't align to a known position in V or J genes and is thus discarded in IgBlast sequence alignment.
        # for now, we will just require user to analyze that sequence manually.
        db$comments[i] <- paste0(db$comments[i], "; different numbers of 'N' in raw VH sequence as compared to germline (check sequence manually)")
        db$full_sequence[i] <- NA
      }
      if(is.na(db$c_call_alignment_length[i]) & !is.na(db$missing_j_bp[i]) & db$missing_j_bp[i]>0 & !is.na(db$full_sequence[i])){
        # add missing bp at the end of the J region
        db$full_sequence[i] <- paste0(db$full_sequence[i],
                                      substring(IMGT_j[strsplit(db$j_call[i], ",")[[1]][1]][[1]], db$j_germline_end[i]+1, length(IMGT_j[strsplit(db$j_call[i], ",")[[1]][1]][[1]])))
        db$comments[i] <- paste0(db$comments[i], "; missing bp in JH sequence, reverted to germline")
      }

      if(!is.na(db$full_sequence[i])){
        if(length(stringr::str_extract_all(db$full_sequence[i], "[^ACGT]")[[1]]) == 0){
          options(warn=-1)
          db$full_sequence_aa[i] <- as.character(Biostrings::translate(Biostrings::DNAStringSet(db$full_sequence[i])))
          options(warn=0)
          CH1 <- CH1_AA[[tail(stringr::str_split(db$c_call[i], pattern = "\\|")[[1]], n=1)]]
          db$full_sequence_fab_aa[i] <- paste0(db$full_sequence_aa[i], CH1)
        } else {
          db$comments[i] <- paste0(db$comments[i], "; still non A|C|T|G characters in full sequence after reverting to germline, check germline_alignment")
        }
      }
    }
  }
  return(db)
}


#### Function to import Sanger-based VDJ sequencing data ####
#' Importing Sanger-based VDJ sequencing data post igblast formatting and merging it to an existing AIRR formatted data frame [optional]
#'
#' \code{importSangerVDJ} import a list of AIRR data frames, update colnames and merge with provided AIRR data frame
#'
#' @param sanger_files      a data frame with the following minimal columns : "sample_id", "directory", "filename" 
#' @param db                an AIRR formatted data frame containing heavy and light chain sequences to which the imported sanger VDJ data will be appended. if NULL, returning only the merged sanger VDJ data
#' @param orig.ident        column to use to set orig.ident; should be either "sample_id" (defined in sanger_dir, if importing plate by plate) or "sort_id" (defined inside each sanger file), if importing from one big table for all plates).
#' @param rename.columns    which columns to rename; a list of pairs new_name/old_name is expected
#' @param prefix.columns    which columns to rename by adding teh "sanger_" prefix
#' @param rm.columns        which columns to remove from all imported sanger files
#' @param overwrite.columns whether to overwrite columns in the imported sanger files with data provided in the sanger_dir recap file [default = FALSE].
#' @param cell_id           name of the column containing cell_id informations
#' @param directory         name of the column containing cell_id informations
#' @param filename          name of the column containing cell_id informations
#' @param run_igblast       whether to run igblast
#' @param igblast_dir       path to igblast database
#' @param imgt_dir          path to IMGT-gapped database
#' @param update_c_call     whether to update c_call
#' @param na.rm             whether to remove wells without associated cell_id
#'
#' @export
#'
#' @import dplyr

importSangerVDJ <- function(sanger_files, db = NULL,
                            orig.ident = "sample_id",
                            cell_id = "cell_id",
                            sequence = "sequence",
                            sequence_id = "sequence_id",
                            directory = "directory",
                            filename = "filename",
                            rename.columns = list(
                              list(old_name = "cell_type", new_name = "sorted_cell_type")),
                            prefix.columns = c("low10QC_alternate_calls", "low30QC_alternate_calls", "missing_v_bp", "missing_j_bp", "comments"),
                            rm.columns = c("well_id", "alternate_well_id", "plate_id", "germline_alignment_d_mask", "mu_count", "mu_count_cdr_r", "mu_count_cdr_s", "mu_count_fwr_r", "mu_count_fwr_s", "primary_sequence", "raw_sequence"),
                            overwrite.columns = FALSE,
                            run_igblast = FALSE,
                            update_c_call = FALSE,
                            seq_type = c("Ig", "TCR"),
                            organism = c("human", "mouse", "rabbit", "rat", "rhesus_monkey"),
                            igblast_dir = "~/share/igblast",
                            imgt_dir = "~/share/germlines/imgt/",
                            na.rm = TRUE){

  # the following columns will be recreated later: "germline_alignment_d_mask", "mu_count", "mu_count_cdr_r", "mu_count_cdr_s", "mu_count_fwr_r", "mu_count_fwr_s"
  # the following columns are discarded as non useful or redundant: "primary_sequence", "raw_sequence", "well_id", "sort_id", "alternate_well_id", "plate_id", 
  #
  # the following columns are maintained with name changes: "cell_type", "low10QC_alternate_calls", "low30QC_alternate_calls", "missing_v_bp", "missing_j_bp", "comments"
  # the following columns will be imported as is: "specificity" "full_sequence", "full_sequence_aa", "full_sequence_fab_aa" and any collumn not listed below (added by user)

  #library(dplyr)

  seq_type <- match.arg(seq_type)
  organism <- match.arg(organism)
  igblast_dir <- ifelse(stringr::str_ends(igblast_dir,"/"), igblast_dir, paste0(igblast_dir, "/"))
  imgt_dir <- ifelse(stringr::str_ends(imgt_dir,"/"), imgt_dir, paste0(imgt_dir, "/"))

  filenames <- ifelse(!stringr::str_ends(sanger_files[[directory]],"/"), paste0(sanger_files[[directory]], "/", sanger_files[[filename]]),
                      paste0(sanger_files[[directory]], sanger_files[[filename]]))

  add_columns <- setdiff(colnames(sanger_files), c(orig.ident, directory, filename))

  if(any(duplicated(filenames))){
    message("Duplicated sample_id in the provided scSangerBCR-seq template: ", sanger_files[duplicated(sanger_files[["sample_id"]]), "sample_id"])
    message("Will only import the first to avoid having issues with duplicated cell_ids")
    sanger_files <- sanger_files[!duplicated(sanger_files[[orig.ident]]), ]
    filenames <- ifelse(!stringr::str_ends(sanger_files[[directory]],"/"), paste0(sanger_files[[directory]], "/", sanger_files[[filename]]),
                        paste0(sanger_files[[directory]], sanger_files[[filename]]))
  }

  AIRR.list <- lapply(filenames, FUN=function(file){
    if(file.exists(file)){
      airr <- openxlsx::read.xlsx(file, sheet=1, rowNames = FALSE)
      if(nrow(airr)>0){
        airr$sequencing_plate <- sanger_files[match(file, filenames), orig.ident]
        for(col in add_columns){
          if(overwrite.columns){
            airr[[col]] <- sanger_files[match(file, filenames), col]
          } else {
            if(!col %in% colnames(airr)){airr[[col]] <- sanger_files[match(file, filenames), col]}
          }
        }
        airr <- airr %>%
          dplyr::mutate(
            across(c("low10QC_alternate_calls", "low30QC_alternate_calls"), as.character)
          ) # prevents issues if all values are NA in one of these columns (too good qualities...)
        return(airr)
      } else {return(NULL)}
    } else {
      warning("the following file could not be found: ", file)
      return(NULL)
    }
  })
  sanger_VDJ_db <- dplyr::bind_rows(AIRR.list)

  # update sequence_id column:
  sanger_VDJ_db$sequence_id <- paste0(sanger_VDJ_db$cell_id, "_", sanger_VDJ_db$primers) # required for DefineClones()

  # rename or remove a few columns:
  for(pair in rename.columns) {
    if(pair$old_name %in% colnames(sanger_VDJ_db)){
      sanger_VDJ_db <- sanger_VDJ_db %>%
        dplyr::rename(!!sym(pair$new_name) := !!sym(pair$old_name))
    } else { warning(paste0("missing column: ", pair$old_name, " in the imported Ab1toAIRR or import Template"))}
  }
  for(column in prefix.columns) {
    sanger_VDJ_db <- sanger_VDJ_db %>%
      dplyr::rename(!!sym(paste0("sanger_", column)) := !!sym(column))
  }
  for(column in rm.columns) {
    sanger_VDJ_db[[column]] <- NULL
  }
  
  if(!is.null(db)){
    if(any(colnames(db) %in% c("complete_vdj", "consensus_count", "umi_count"))){
      #creating a few columns to adapt to scRNAseq datasets
      sanger_VDJ_db$complete_vdj <- TRUE
      sanger_VDJ_db$consensus_count <- 10
      sanger_VDJ_db$umi_count <- 10
    }
  }

  if(na.rm){
    missing_cell_id <- sanger_VDJ_db %>%
      dplyr::filter(is.na(!!rlang::sym(cell_id)))
    if(nrow(missing_cell_id)>1){
      warning(nrow(missing_cell_id), " sequences without associated cell_id from the following sequencing plates: ", paste(unique(missing_cell_id$sequencing_plate),collapse =", "), "; will be removed from imported dataset(s)")
    }
    sanger_VDJ_db <- sanger_VDJ_db %>%
      dplyr::filter(!is.na(!!rlang::sym(cell_id)))
  }

  if(run_igblast){
    igblast_results <- runAssignGenes(sanger_VDJ_db,
                                      organism = organism,
                                      seq_type = seq_type,
                                      igblast_dir = igblast_dir,
                                      imgt_dir = imgt_dir,
                                      sequence = sequence,
                                      sequence_id = sequence_id)

    sanger_VDJ_db <- igblast_results[["pass"]]

    sanger_VDJ_db$c_call_igblast <- sanger_VDJ_db$c_call
  }
  if(update_c_call){
    sanger_VDJ_db <- runBlastnC(sanger_VDJ_db,
                                organism = organism,
                                seq_type = seq_type,
                                igblast_dir = igblast_dir)
  }

  if(!is.null(db)){
    if("c_call_SB" %in% colnames(db)){
      # creating a few more columns to adapt to Rhapsody VDJ db
      sanger_VDJ_db$cell_type_experimental <- sanger_VDJ_db$sorted_cell_type
      sanger_VDJ_db$putative_cell <- TRUE
      sanger_VDJ_db$dominant <- TRUE
    }
    full_db <- dplyr::bind_rows(db, sanger_VDJ_db)
    return(full_db)
  } else {return(sanger_VDJ_db)}
}


#### Wrapper function to import and QC VDJ data from single cell RNA-seq datasets (both 10X and BD Rhapsody for now) ####
#' Resolve VDJ heavy chains multiplets
#'
#' \code{scImportVDJ} import, run igblast and QC VDJ output or BD scRNA-seq datasets
#'
#' @param seurat        a merged seurat V5 object on which to filter VDJ contigs or simply the metadata of a seurat object (seurat [[]])[better, more lightweight]
#' @param vdj_files     a dataframe listing the data sets to import; first column should be "SB_analysis_id", "10X_lane_id" or defined by the following sample_id argument, other two expected columns should be "directory" and "filename" and will be used to reconstruct the full_path to the file to import, if a "full path" column is provided it will be used. (see attached excel spreadsheets - one for 10X and one for BD)
#' @param sample_id     name of the column containing the run name in the vdj_dir table [default = "SB_analysis_id" or "10X_lane_id"]; will be used to create unique cell_ids matching the one in the merged seurat object
#' @param analysis_name prefix to use for full fasta and tsv files being exported
#' @param import_from_seurat vector of colnames to import from seurat metadata to the AIRR dataframe (option are Cartridge_Nb, donor_id, sample_tag, tissue, cluster...)
#' @param clean_HC      whether to resolve cases of multiple heavy chains and identify doublets (resolveMultiHC())
#' @param split.by      name of the column in the seurat object to use to split the dataset prior to filtering heavy chain [if NULL, switch to input "SB_analysis_id" or "10X_lane_id"]
#' @param tech          name of the scRNA-seq technology used (one of 10X or BD
#' @param seq_type        type of VDJ sequence ("Ig" or "TCR" to match igblastb requirements)
#' @param organism      organism (any of "human", "mouse", "rhesus_monkey; for other see https://changeo.readthedocs.io/en/stable/examples/igblast.html)
#' @param igblast       whether to run standalone IgBlast, can be set to c("filtered heavy" or "all") if not one of these three values, will be skipped with a warning. [default = "filtered heavy" for both 10X and BD: highly recommended for both to avoid issues at the createGermline() or observedMutation() steps due to different references databases used (10X) or missing imgt gaps in the sequence_alignment collumn (Both))]
#' @param igblast_dir   path to igblast database [default = path suggested on installation: https://changeo.readthedocs.io/en/stable/examples/igblast.html]
#' @param imgt_dir      path to imgt database [default = path suggested on installation: https://changeo.readthedocs.io/en/stable/examples/igblast.html]
#' @param update_c_call whether to run runBlastnC to correct c calls made by igblast (issues with calls with similar scores); can be set to c("heavy" or "all") if not one of these three values, will be skipped with a warning. [default = "heavy" for all cases, will be performed onlight chains after heavy chain clustering and light chain multiplets resolution (see scFindBCRClones())]
#' @param imgt_c_dir    path to Blast-able imgt C gene database [default = path suggested in Create_IMGT_human_IG_C_Blastdb.R script]
#' @param cutoff        choice of cutoff: can be set as "variable" (default) or "fixed";
#'                      if fixed, cutoffs should be provided or preset cutoffs will be used,
#'                      if variable [preferred option], cutoffs are automatically calculated based on summary counts in the database to account for potential sequencing bias (1st quartile and max of 10 or 1/10 of median).
#' @param low_cutoff    cut_off for low probability heavy chain doublets.
#' @param high_cutoff   cut_off for high probability heavy chain doublets.
#' @param na.rm         whether to filter sequences without an identified CDR3.
#' @param output        whether to output graphs with umi_counts for dominant versus second IGH VDJ contig and the recap excel workbook. If set to FALSE, only the corrected database is returned.
#' @param output_folder name of the folder in which graph and recap excel workbooks will be saved [default = "VDJ_QC"].
#' @param cell_id       name of the column containing cell identifier.
#' @param locus         name of column containing locus values.
#' @param heavy         value of heavy chains in locus column. All other values will be
#'                      treated as light chains.
#' @param productive    name of the column containing the info whether a given sequence is productive.
#' @param sequence_id   name of the column containing sequence identifier.
#' @param consensus_count   name of the column containing the number of reads for this contig (usually called consensus_count)
#' @param umi_count     name of the column containing the number of unique molecules (UMI) for this contig. Previously called "duplicate_count" in an earlier AIRR standard
#' @param junction      name of the column containing identified junction in nucleotide format.
#' @param junction_aa   name of the column containing identified junction in amino-acid format.
#' @param sequence      name of the column containing the original sequence.
#' @param v_call        name of the column containing V-segment allele assignments. All
#'                      entries in this column should be identical to the gene level.
#' @param d_call        name of the column containing D-segment allele assignments. All
#'                      entries in this column should be identical to the gene level.
#' @param j_call        name of the column containing J-segment allele assignments. All
#'                      entries in this column should be identical to the gene level.
#' @param c_call        name of the column containing Constant region assignments. All
#'                      entries in this column should be identical to the gene level.
#' @param complete_vdj  name of the column containing complete_vdj assignments.
#' @param verbose       whether to write to the console
#'
#' @return   an AIRR formatted dataframe containing the dominant heavy chain selected for all cells, with the columns: is.VDJ_doublets, is.VDJ_doublet.confidence
#' and with the following columns added for cells with two or more heavy chains detected:
#' second_umi_count, second_concensus_count, second_sequence, second_junction, second_junction_aa, second_v_call, second_d_call, second_j_call, second_c_call.
#' Also outputs one graph and one recap excel workbook with filtered dataframe (sheet 1) and QC parameters (sheet 2) for each group of cells identified by the split.by argument [see resolveMultiHC()].
#'
#' @details
#' 1. will first import VDJ datasets listed in vdj_dir input dataframe and filter on cells in provided seurat V5 object;
#' 2. will run igblast on filtered datasets if requested;
#' 3. will run resolveMultiHC() to resolve HC multiplets
#' @import dplyr
#' @import readr
#' @import purrr
#' @importFrom Biostrings readDNAStringSet
#' 
#' @export

scImportVDJ <- function(vdj_files,
                        seurat = NULL,
                        import_from_seurat = NULL,
                        tech = c("BD", "10X"),
                        remove_columns_on_import = list("10X" = c("high_confidence", "raw_clonotype_id",	"raw_consensus_id",	"exact_subclonotype_id",
                                                                  "fwr1_aa", "fwr2_aa", "fwr3_aa", "fwr4_aa", "cdr1_aa", "cdr2_aa", "cdr3_aa"),
                                                        "BD" = c("cell_type_experimental", "high_quality_cell_tcr_bcr", "sequence_aa", "sequence_aa_length" , "sequence_alignment_length",
                                                                 "cdr3_length", "fwr1_aa", "fwr2_aa", "fwr3_aa", "fwr4_aa", "cdr1_aa", "cdr2_aa", "cdr3_aa",
                                                                 "germline_alignment_aa", "v_germline_alignment", "v_germline_alignment_aa",
                                                                 "d_germline_alignment", "d_germline_alignment_aa", "j_germline_alignment", "j_germline_alignment_aa")),
                        split.by = NULL,
                        sample_id = NULL,
                        analysis_name = "All_sequences",
                        seq_type = c("Ig", "TCR"),
                        organism = c("human", "mouse", "rabbit", "rat", "rhesus_monkey"),
                        igblast = c("all", "filtered heavy", "none"),
                        igblast_dir = "~/share/igblast/",
                        imgt_dir = "~/share/germlines/imgt/",
                        update_c_call = c("none", "filtered heavy", "all"),
                        clean_HC = TRUE,
                        cutoff = "variable",
                        low_cutoff = 10,
                        high_cutoff = 250,
                        na.rm = FALSE,
                        output = TRUE,
                        verbose = TRUE,
                        output_folder = "VDJ_QC",
                        cell_id = "cell_id", locus = "locus", productive = "productive", complete_vdj = "complete_vdj",
                        sequence_id = "sequence_id", umi_count = "umi_count", consensus_count = "consensus_count",
                        junction = "junction", junction_aa = "junction_aa", sequence = "sequence",
                        v_call = "v_call", d_call = "d_call", j_call = "j_call", c_call = "c_call"){

  start <- Sys.time()

  #suppressMessages(library(dplyr))

  seq_type <- match.arg(seq_type)
  if(seq_type == "Ig"){
    heavy <- "IGH"
  }
  if(seq_type == "TCR"){
    heavy <- c("TRB","TRD")
    clean_HC <- FALSE
  }
  
  organism <- match.arg(organism)
  igblast_dir <- ifelse(stringr::str_ends(igblast_dir,"/"), igblast_dir, paste0(igblast_dir, "/"))
  imgt_dir <- ifelse(stringr::str_ends(imgt_dir,"/"), imgt_dir, paste0(imgt_dir, "/"))

  if(!is.null(output_folder)){
    if(!stringr::str_ends(output_folder, "/")){output_folder = paste0(output_folder, "/")}
    if(!dir.exists(output_folder)){
      dir.create(output_folder)
    }
  }

  #initiate log file:
  log_file <- paste0(output_folder, analysis_name, "_scImportVDJ.log")
  time_and_log({
    cat("analysis_name: ", analysis_name, "\n")
    cat("tech: ", tech, "\n")
    cat("imported columns from seurat object: ", import_from_seurat, "\n")
    cat("igblast: ", igblast, " contigs\n")
    cat("organism: ", organism, "\n")
    cat("seq_type: ", seq_type, "\n")
    cat("igblast: ", igblast, "\n")
    cat("update_c_call: ", update_c_call, "\n")
    cat("clean_HC: ", clean_HC, "\n")
  }, verbose = FALSE, time = FALSE, log_file = log_file, log_title = "scImportVDJ", open_mode = "wt")
  
  #log_file <- paste0(output_folder, Sys.time(), "_", analysis_name, "_scImportVDJ_logfile.txt")
  #log_connection <- file(log_file, open = "a")  # a for appending or w for erasing onto previous log

  if(length(tech)!=1 | !(tech %in% c("BD", "10X"))){
    stop("no clear technology selected, you need to select one of BD or 10X")
  }

  if(is.null(sample_id)){
    if(tech == "10X"){sample_id <- "CR_sample_id"}
    if(tech == "BD"){sample_id <- "SB_analysis_id"}
  }

  missing_vdj_files_columns <- setdiff(c(sample_id, "full_path", "directory","filename_vdj"), colnames(vdj_files))
  if(!sample_id %in% colnames(vdj_files) | ((!"full_path" %in% colnames(vdj_files) & (!"directory" %in% colnames(vdj_files) | !"filename_vdj" %in% colnames(vdj_files))))){
    stop(paste0("missing the following collumns in the provided vdj_dir dataframe: ", missing_vdj_files_columns))
    }

  igblast <- match.arg(igblast)
  if(!igblast %in% c("filtered heavy", "all")){warning("It is highly recommended to rerun IgBlast on 10X or BD data to avoid issues later in the pipeline (set igblast = TRUE)")}
  update_c_call <- update_c_call[1]

  ## Part1: import all VDJ data; filter based on seurat cell_ids (if present) and run IgBlast if needed (e.g. 10X):
  if(inherits(seurat, "Seurat")){meta <- seurat[[]]} else {meta <- seurat}
  rm(seurat)
  if(!is.null(meta)){
    if(!cell_id %in% colnames(meta))
    {stop("no cell_id column present in the seurat object provided")}
  }

  if(any(duplicated(meta$cell_id))){
    stop("provided Seurat object contains duplicated cell_id, high risk of confusion in later analysis!!")
  } # cell_ids should be in the format: SB_analysis_id_number (e.g. PC-A3-VDJ_1010)

  # update split.by argument
  if(!is.null(split.by)){
    if(is.null(meta)){
      split.by <- "orig.ident"
    } else {
      if(!split.by %in% colnames(meta)){
        warning(paste0(split.by, " column not found in the provided seurat object, defaulting to split.by = orig.ident"))
        split.by <- "orig.ident"
        }
    }
  } else {
    split.by <- "orig.ident"
  }

  # update import_from_seurat argument to force import of the split.by metadata if present
  if(!is.null(split.by) & split.by %in% colnames(meta)){import_from_seurat <- unique(c(import_from_seurat, split.by))}

  n = 1 + sum(clean_HC, (igblast %in% c("filtered heavy", "all")), (update_c_call %in% c("filtered heavy", "all")))
  step = 1

  message("------------")
  message("Part ", step," of ",n ,": importing unfiltered VDJ outputs")
  message("------------")
  
  safe_pblapply <- function(X, FUN, ..., mc.cores = 1) {
    if (requireNamespace("pbapply", quietly = TRUE)) {
      return(pbapply::pblapply(X, FUN, ...))
    } else {
      message("Package 'pbapply' not found. Falling back to lapply().")
      return(lapply(X, FUN, ...))
    }
  }
  
  VDJ.list <- safe_pblapply(vdj_files[[sample_id]], FUN=function(sample){
    if("full_path" %in% colnames(vdj_files)){
      file <- vdj_files[match(sample,vdj_files[[sample_id]]),"full_path"]
    } else {
      file <-  do.call(paste0, vdj_files[match(sample,vdj_files[[sample_id]]),c("directory","filename_vdj")])
      }

    message("reading file for: ", sample)

    #using readr::read_tsv/csv as it is super fast and works directly on compressed files (tsv.gz)
    if(tech == "BD"){
      #VDJ output from SevenBridges is a .tsv.gz file
      data <- readr::read_tsv(file, show_col_types = FALSE)
    }
    if(tech == "10X"){
      #VDJ output from CellRanger is a .csv file
      #but you also need to go fetch the full sequences stored in a separate .fasta file
      #(if using all_contig_annotations.csv, should be all_contig.fasta and should be in the same folder).
      data <- readr::read_csv(file, show_col_types = FALSE)
      fasta_file <- sub("_annotations.csv", ".fasta", file)
      fasta <- Biostrings::readDNAStringSet(fasta_file)
      data$sequence <- as.character(fasta)
    }
    data <- reformatVDJinput(data, tech = tech,
                             cell_id = cell_id, locus = locus, heavy = heavy, productive = productive, complete_vdj = complete_vdj,
                             sequence_id = sequence_id, umi_count = umi_count, consensus_count = consensus_count,
                             junction = junction, junction_aa = junction_aa, sequence = sequence,
                             v_call = v_call, d_call = d_call, j_call = j_call, c_call = c_call,
                             remove_columns = remove_columns_on_import)
    data[[sequence_id]] <- paste0(sample, "_", data$sequence_id)
    data[[cell_id]] <- paste0(sample, "_", data$cell_id)

    #filter for cell_ids in the provided seurat object
    if(!is.null(meta)){
      cells_to_keep <- meta[[cell_id]]
      data <- dplyr::filter(data, (!!rlang::sym(cell_id) %in% cells_to_keep))
      if(nrow(data)==0){
        warning(paste0(sample, " : no VDJ contigs corresponding to cell_ids present in the provided seurat object. Check cell_id format!"))
        time_and_log({
          warning(paste0(sample, " : no VDJ contigs corresponding to cell_ids present in the provided seurat object. Check cell_id format!"))
        }, verbose = FALSE, time = FALSE, log_file = log_file, log_title = paste0(sample, " contigs import"), open_mode = "a")
        }
    }
    data$orig.ident <- sample
    data$assay <- tech
    readr::write_tsv(data, file = paste0(output_folder, sample, "_VDJ_CellQCfiltered_contigs_AIRR.tsv.gz"))

    return(data)
  })

  VDJ_db <- VDJ.list %>%
    purrr::keep(~ nrow(.) > 0) %>% #removing empty tibble to away issues with Error in `dplyr::bind_rows()`:! Can't combine `..X$consensus_count` <double> (normal) and `..5$consensus_count` <character> (no values, i.e. NA).
    dplyr::bind_rows()
  rm(VDJ.list)

  #import additional metadata from provided seurat object:
  if(!is.null(import_from_seurat)){
    missing_collumns <- setdiff(import_from_seurat, colnames(meta))
    if(length(missing_collumns)>0){warning(paste0("missing the following collumns in the provided seurat object: ", missing_collumns))}
    import_from_seurat <- intersect(import_from_seurat, colnames(meta))
    meta <- meta %>%
      dplyr::select(all_of(c("cell_id", import_from_seurat)))

    already_present <- intersect(import_from_seurat,colnames(VDJ_db))
    if(length(already_present)>0){warning(paste0("the following collumns are alreading present in the AIRR dataframe: ", already_present, ". They will be replaced"))}
    VDJ_db <- VDJ_db %>%
      dplyr::select(-all_of(already_present)) %>%
      dplyr::left_join(meta, by = cell_id)
  }

  #check for duplicated sequence_ids:
  if(any(duplicated(VDJ_db$sequence_id))){
    duplicated_ids <- VDJ_db[duplicated(VDJ_db$sequence_id),"sequence_id"]
    duplicated_db <- VDJ_db %>%
      dplyr::filter(sequence_id %in% duplicated_ids)
    duplicated_filename <- paste0(output_folder, analysis_name, "_VDJ_CellQCfiltered_duplicated_sequence_id")
    readr::write_tsv(duplicated_db, file = paste0(duplicated_filename, ".tsv.gz"))
    warning("the following sequence_id are duplicated and will be removed: ", paste(duplicated_ids, collapse = ", "))
    VDJ_db <- VDJ_db[!duplicated(VDJ_db$sequence_id),]
  }

  filename <- paste0(output_folder, analysis_name, "_VDJ_CellQCfiltered")
  if(!clean_HC & !(igblast %in% c("filtered heavy", "all") & !(update_c_call %in% c("filtered heavy", "all")))){#'we only save this file if no other analysis is performed
    readr::write_tsv(VDJ_db, file = paste0(filename, ".tsv.gz"))
  }

  total_submitted_heavy <- nrow(dplyr::filter(VDJ_db, (!!rlang::sym(locus) == heavy)))
  total_submitted_light <- nrow(dplyr::filter(VDJ_db, (!!rlang::sym(locus) != heavy)))
  total_submitted_cells <- length(unique(VDJ_db$cell_id))

  ##Part2 [option 1]: run igblast on all sequences: can take around 15min for 300k+ sequences
  if(igblast == "all"){# not the preferred option as you will use computer power and time on sequences you will discard later but can help getting rid of contigs that won't pass igblast anyway...
    step <- step + 1
    submitted_seqs <- nrow(VDJ_db)
    submitted_cells <- length(unique(VDJ_db$cell_id))
    expected_time <- ifelse(submitted_seqs < 5000, "less than a minute", 
                            ifelse(submitted_seqs < 200000, "a few minutes", "more than 10 minutes"))
    
    message("------------")
    #message("Part ", step," of ",n ,": running IgBlast on all contigs (can take more than 15min for 300k+ sequences... be patient!)")
    message("Part ", step," of ",n ,": running IgBlast on all contigs (",submitted_seqs," contigs, should take ", expected_time, ")")
    message("------------")
    
    time_and_log({
      igblast_results <- runAssignGenes(VDJ_db,
                                        organism = organism,
                                        seq_type = seq_type,
                                        igblast_dir = igblast_dir,
                                        imgt_dir = imgt_dir,
                                        sequence = sequence,
                                        sequence_id = sequence_id)
    }, verbose = verbose, log_file = log_file, log_title = "running AssignGenes", open_mode = "a")

    VDJ_db <- igblast_results[["pass"]]
    VDJ_db$c_call_igblast <- VDJ_db$c_call

    failed_VDJ_db <- igblast_results[["fail"]]

    rm(igblast_results)

    VDJ_db <- dplyr::relocate(VDJ_db, !!rlang::sym(cell_id), .before = !!rlang::sym(sequence_id))
    VDJ_db <- dplyr::relocate(VDJ_db, c(orig.ident, assay, !!rlang::sym(consensus_count), !!rlang::sym(umi_count)), .after = !!rlang::sym(sequence_id))
    VDJ_db <- dplyr::relocate(VDJ_db, !!rlang::sym(c_call), .after = !!rlang::sym(j_call))

    filename_fail <- paste0(filename, "_igblast_db-fail")
    readr::write_tsv(failed_VDJ_db, file = paste0(filename_fail, ".tsv.gz"))

    filename <- paste0(filename, "_igblast_db-pass")
    if(!(update_c_call %in% c("heavy", "all"))|clean_HC){#we only save this file if no other analysis is performed
      readr::write_tsv(VDJ_db, file = paste0(filename, ".tsv.gz"))
    }

    nb_failed_sequences <- nrow(failed_VDJ_db)
    nb_failed_cells <- length(unique(failed_VDJ_db$cell_id))
    
    if(verbose){
      cat("nb submitted: ", submitted_seqs, " contigs from ", submitted_cells, " cells\n")
      cat("nb pass: ", nrow(VDJ_db), " contigs from ", length(unique(VDJ_db$cell_id)), " cells\n")
      cat("nb failed: ", nb_failed_sequences, " contigs from ", nb_failed_cells, " cells\n")
    }
  }

  ## Part2 [option 2]: split VDJ_db based on sample_id provided and perform HC filtering if needed [skip if using filtered VDJ data from either technologies]:
  if(clean_HC){
    step <- step + 1
    message("------------")
    message("Part ", step," of ", n,": splitting by ",split.by ," and performing heavy chain QC")
    message("------------")
    
    time_and_log({
      VDJ_db <- resolveMultiContigs(VDJ_db, 
                                    split.by = split.by, 
                                    seq_type = seq_type,
                                    resolve_chain = "heavy",
                                    assay = "assay", 
                                    resolve_multi_CDR3 = TRUE, 
                                    use_clone = FALSE,
                                    analysis_name = analysis_name,
                                    output = output, 
                                    output_folder = output_folder,
                                    second_columns = c(sequence_id, locus, umi_count, consensus_count, sequence, v_call, d_call, j_call, c_call, junction, junction_aa, productive, complete_vdj),
                                    cell_id = cell_id, 
                                    sequence_id = sequence_id, 
                                    locus = locus, 
                                    consensus_count = consensus_count, 
                                    umi_count = umi_count, 
                                    v_call = v_call, 
                                    j_call = j_call, 
                                    c_call = c_call, 
                                    junction_aa = junction_aa,
                                    productive = productive, 
                                    complete_vdj = complete_vdj)
    }, verbose = verbose, log_file = log_file, log_title = "running resolveMultiHC()", open_mode = "a")
    
    filename <- paste0(filename, "_HCfilter-pass")
    if(!(igblast == "filtered heavy") & !(update_c_call %in% c("heavy", "all"))){#'we only save this file if no other analysis is performed
      readr::write_tsv(VDJ_db, file = paste0(filename, ".tsv.gz"))
    }
  }

  ## Part3: run IgBlast on filtered heavy chain only
  if(igblast == "filtered heavy"){# default option as much faster
    step <- step + 1
    h_db <- dplyr::filter(VDJ_db, (!!rlang::sym(locus) == heavy))
    l_db <- dplyr::filter(VDJ_db, (!!rlang::sym(locus) != heavy))
    
    submitted_seqs <- nrow(h_db)
    submitted_cells <- length(unique(h_db$cell_id))
    expected_time <- ifelse(submitted_seqs < 5000, "less than a minute", 
                            ifelse(submitted_seqs < 200000, "a few minutes", "more than 10 minutes"))
                            
    message("------------")
    message("Part ", step," of ",n ,": running IgBlast on filtered HC contigs (",submitted_seqs," contigs, should take ", expected_time, ")")
    message("------------")
    
    time_and_log({
      igblast_results <- runAssignGenes(h_db,
                                        organism = organism,
                                        seq_type = seq_type,
                                        igblast_dir = igblast_dir,
                                        imgt_dir = imgt_dir,
                                        sequence = sequence,
                                        sequence_id = sequence_id)
    }, verbose = verbose, log_file = log_file, open_mode = "wt")
    
    h_db <- igblast_results[["pass"]]
    h_db$c_call_igblast <- h_db$c_call

    failed_h_db <- igblast_results[["fail"]]

    rm(igblast_results)

    filename <- paste0(filename, "_HCigblast_db-pass")
    if(!(update_c_call %in% c("heavy", "all"))){# we only save this file if no other analysis is performed
      readr::write_tsv(h_db, file = paste0(filename, ".tsv.gz"))
    }
    filename_fail <- paste0(filename, "_HCigblast_db-fail")
    readr::write_tsv(failed_h_db, file = paste0(filename_fail, ".tsv.gz"))

    nb_failed_sequences <- nrow(failed_h_db)
    nb_failed_cells <- length(unique(failed_h_db$cell_id))
    
    if(verbose){
      cat("nb submitted: ", submitted_seqs, " contigs from ", submitted_cells, " cells\n")
      cat("nb pass: ", nrow(h_db), " contigs from ", length(unique(h_db$cell_id)), " cells\n")
      cat("nb failed: ", nb_failed_sequences, " contigs from ", nb_failed_cells, " cells\n")
    }

    VDJ_db <- dplyr::bind_rows(h_db, l_db)
    VDJ_db <- dplyr::relocate(VDJ_db, !!rlang::sym(cell_id), .before = !!rlang::sym(sequence_id))
    VDJ_db <- dplyr::relocate(VDJ_db, c(orig.ident, assay, !!rlang::sym(consensus_count), !!rlang::sym(umi_count)), .after = !!rlang::sym(sequence_id))
    VDJ_db <- dplyr::relocate(VDJ_db, !!rlang::sym(c_call), .after = !!rlang::sym(j_call))
  }

  ## Part4: update c_calls for heavy chains:
  if(update_c_call == "filtered heavy"){
    step <- step + 1
    message("------------")
    message("Part ", step," of ",n ,": updating c_call for filtered heavy chain contigs")
    message("------------")
    
    h_db <- dplyr::filter(VDJ_db, (!!rlang::sym(locus) == heavy))
    l_db <- dplyr::filter(VDJ_db, (!!rlang::sym(locus) != heavy))

    time_and_log({
      h_db <- runBlastnC(h_db,
                         igblast_dir = igblast_dir,
                         organism = organism,
                         seq_type = seq_type,
                         sequence = sequence,
                         sequence_id = sequence_id)
    }, verbose = verbose, log_file = log_file, open_mode = "a")

    h_db <- dplyr::relocate(h_db, !!rlang::sym(c_call), .after = !!rlang::sym(j_call))

    VDJ_db <- dplyr::bind_rows(h_db, l_db)
    filename <- paste0(filename, "_c_call-pass")
    readr::write_tsv(VDJ_db, file = paste0(filename, ".tsv.gz"))
  }
  if(update_c_call == "all"){# not the preferred option as you will use computer power and time on sequences you will discard later...
    step <- step + 1
    message("------------")
    message("Part ", step," of ", n,": updating c_call for all contigs")
    message("------------")
    
    time_and_log({
      VDJ_db <- runBlastnC(VDJ_db,
                           igblast_dir = igblast_dir,
                           organism = organism,
                           seq_type = seq_type,
                           sequence = sequence,
                           sequence_id = sequence_id)
    }, verbose = verbose, log_file = log_file, open_mode = "a")
    
    VDJ_db <- dplyr::relocate(VDJ_db, !!rlang::sym(c_call), .after = !!rlang::sym(j_call))

    filename <- paste0(filename, "_c_call-pass")
    readr::write_tsv(VDJ_db, file = paste0(filename, ".tsv.gz"))
  }

  message("------------")
  if(verbose){cat("submitted: ", total_submitted_heavy, " IGH contigs and ", total_submitted_light, " IGL/IGK contigs (",total_submitted_cells," unique cells).\n")}
  if(verbose && igblast == "filtered heavy"){
    cat(nb_failed_sequences, " IGH contigs failed igblast.\n")
  }
  if(verbose && igblast == "all"){
    cat(nb_failed_sequences, " contigs from ",nb_failed_cells," cells failed igblast.\n")
  }
  n_final_heavy <- nrow(dplyr::filter(VDJ_db, (!!rlang::sym(locus) == heavy)))
  n_final_light <- nrow(dplyr::filter(VDJ_db, (!!rlang::sym(locus) != heavy)))
  n_final_cells <- length(unique(VDJ_db$cell_id))
  if(verbose){cat("final table: ", n_final_heavy, " IGH contigs and ", n_final_light, " IGL/IGK contigs (",n_final_cells," unique cells).\n")}

  time_and_log({
    print(sessionInfo())
  }, verbose = FALSE, time = FALSE, log_file = log_file, log_title = "session info", open_mode = "a")
  
  end <- Sys.time()
  if(verbose){cat(sprintf("Total running time: %.2f %s", end-start, units(difftime(end, start))),"\n")}

  return(VDJ_db)
}


#### Wrapper function to perform clonal grouping on mix paired (single-cell) and unpaired (bulk or sanger) VDJ data ####
#' Full BCR or TCR analysis pipeline for AIRR single-cell repertoire datasets
#'
#' \code{scFindClones} performs clonal clustering of single cell data
#'
#' @param db              an AIRR formatted dataframe containing BCR (heavy and light chain) or TCR sequences. For BCR, should contain only one heavy chain (IGH) per cell_id, if not will run resolveMultiHC() first. For TCR, expects filtered datasets (1 heavy and 1 light chain per cell).
#' @param method          method to use for scoper::hierarchicalClones(), can be one between: identical, hierarchical or spectral
#' @param spectral_method method to use for scoper::spectralClones(), can be one between: novj or vj
#' @param threshold       method to use for scoper::hierarchicalClones() or scoper::spectralClones(),
#'                        if multiple thresholds are provided, they will be used successively and results stored in the h_clone_"threshold" column. Final clone_id column will reflect the last threshold used.
#' @param clean_LC        whether to resolve cases of multiple light chains
#' @param tech            which tech was used, passed to resolveMultiContigs, only tech = "BD" will results in additional QC being performed.
#' @param split_by_light  whether to clean cases of multiple light chains (resolveMultiLC()) and split clone by ligth chain (dowser::resolveLightChains()), the clone_id collumn will be updated and results will also be stored in the l_clone_id_"last_used_threshold" column;
#' @param seq_type        type of VDJ sequence ("Ig" or "TCR" to match igblastb requirements)
#' @param organism      organism (any of "human", "mouse", "rhesus_monkey; for other see https://changeo.readthedocs.io/en/stable/examples/igblast.html)
#' @param update_germline whether to run dowser::createGermlines() to determine consensus clone sequence and create germline for clone after splitting by light chain;
#' @param SHM             whether to calculate mutational load;
#' @param output          whether to output graphs with umi_counts for dominant versus second IGH VDJ contig and the recap excel workbook. If set to FALSE, only the corrected database is returned.
#' @param output_folder   name of the folder in which graph and recap excel workbooks will be saved [default = "VDJ_Clones"]
#' @param igblast         whether to run standalone IgBlast, can be set to c("filtered_light" or "all") if not one of these three values, will be skipped with a warning. [default = "filtered_light" for both 10X and BD: highly recommended for both to avoid issues at the createGermline() or observedMutation() steps due to different references databases used (10X) or missing imgt gaps in the sequence_alignment collumn (Both))]
#' @param igblast_dir     path to igblast database [default = path suggested on installation: https://changeo.readthedocs.io/en/stable/examples/igblast.html]
#' @param imgt_dir        path to imgt-gapped database [default = path suggested on installation: https://changeo.readthedocs.io/en/stable/examples/igblast.html]
#' @param update_c_call   whether to run runBlastnC to correct c calls made by igblast (issues with calls with similar scores); can be set to c("filtered light" or "all") if not one of these three values, will be skipped with a warning. [default = "all", if using "filtered light" should have already been previously performed on heavy chains (see scImportVDJ())]
#' @param imgt_j_nt       path to imgt_human_ig_j.fasta (passed to reconstructFullVDJ())
#' @param analysis_name   name to use for outputs prefixes [default = "All_sequences"]
#' @param cell_id         name of the column containing cell identifier.
#' @param cell_id_scoper  name of the column containing cell identifier to use scoper in single_cell mode. If only heavy chain data, use NULL to avoid a bug in scoper.
#' @param locus           name of column containing locus values.
#' @param heavy           value of heavy chains in locus column. All other values will be
#' @param na.heavy.rm     whether to remove light chain sequences for cells with no heavy chain from the final database [default = TRUE, !!TODO: work on re-importing these sequences in the final recap table if na.heavy.rm = FALSE]
#' @param nproc           number of processor to use for parallel computing (passed on to all Immcantation functions).
#' @param verbose         whether to show messages and cat() ouputs on console, can be "all" or "partial" (messages only)[TODO], all other values will be equivalent to "none"
#' @param consensus_count     name of the column containing the number of reads for this contig (usually called consensus_count)
#' @param umi_count       name of the column containing the number of unique molecules (UMI) for this contig. Previously called "duplicate_count" in an earlier AIRR standard
#' @param productive      name of the column containing the info whether a given sequence is productive.
#' @param sequence_id     name of the column containing sequence identifier.
#' @param sequence        name of the column containing the original sequence.
#' @param junction        name of the column containing identified junction in nucleotide format.
#' @param junction_aa     name of the column containing identified junction in amino-acid format.
#' @param v_call          name of the column containing V-segment allele assignments. All entries in this column should be identical to the gene level.
#' @param d_call          name of the column containing D-segment allele assignments. All entries in this column should be identical to the gene level.
#' @param j_call          name of the column containing J-segment allele assignments. All entries in this column should be identical to the gene level.
#' @param c_call          name of the column containing Constant region assignments. All entries in this column should be identical to the gene level.
#' @param full_seq_aa     whether to reconstruct full VDJ sequence (nt and AA +/- CH1 domain)
#' @param CH1_AA          named list of CH1 domain for IGH, IGL and IGK c chains (AA); passed to reconstructFullVDJ())
#' @param assay           name of column containing assay values.
#' @param shared.tech     list of grouped technologies to identify shared clones between scRNA-seq and Sanger sequencing (export in recap teable (seperate sheet) + additional shared_clone column)
#' @param fields          Character vector of additional columns to use for grouping in scoper::hierarchicalClones() or dowser::resolveLightChains(). Sequences with disjoint values in the specified fields will be considered as separate clones.
#' @param only_heavy      Whether all submitted chains are heavy chains. Automatically set some options accordingly.
#' @param split.by        name of the column used for grouping at the resolveLigth chain step.
#' @param complete_vdj    name of the column containing complete_vdj info.
#' @param junc_len        name of the column containing junction length info.
#' @param sequence_alignment name of the column containing sequence alignment.
#' @param recap.highlight Which column to use to create extra sheets in the recap table
#' @param top_column      Which columns to relocate to the front of the final table
#'
#' @return    an AIRR formatted dataframe containing selected heavy and light chains (when available) for all cells as well as clonal grouping for all cells (clone_id).
#' If igblast and update_c_call are set to "light" or "all", will run igblast and runBlastnC on requested contigs.
#' If split_by_light is set to TRUE, will also perform clonal partition correction based on observed light chain or heavy chain proximity to cell with known light chain.
#' If update_germline and SHM are set to TRUE, will further update the germline alignment columns based on other contigs in the same clonal group, add a germline_d_mask column,
#' and then analyse the number and frequency of mutation in VH and VL genes, if split_by_light is set to TRUE or only the VH gene if not.
#' Creates a bcr_info or tcr_info columns with three possible values: "full" (if both heavy and light contigs are found for a given cell_id), "heavy_only" or "light_only".
#' Creates an expanded_clone (>1 cell_id) and a shared_clone (between scRNAseq and Sanger), see shared.tech argument.
#' Also outputs a graph for minimum distance between and maximal distance inside heavy clones for each tested threshold, all intermediate results for each step in the analysis and a final recap table which include all chosen parameters for the analysis as sheet 2.
#'
#' @details
#' Works as a wrapper function successively calling in house or Immcantation-based functions to perform full clonal analysis of an AIRR-formatted single B cell repertoire.
#' All key parameters are preset and the simple db <- scFindBCRClones(db) line should be sufficient.
#' Provided dataframe should be an AIRR formatted dataframe and should contain only one heavy chain (IGH) per cell_id, if not run resolveMultiHC() first.
#' Not all cells need to have a detected light chain but cells without a full heavy chain with cdr3 will be removed.
#' 1. it first performs clonal grouping based on the heavy chain using the Scoper package functions;
#' 2. cases of multiple light chains are then resolved using identified groups (resolveMultiLC());
#' 3. clonal groups are further split based on light chain info (step 1) or proximity with other heavy chain for cells with missing light chains (dowser::resolveLightChains);
#' 4. finally, we use dowser::createGermlines() and shazam::observedMutations() to determine consensus clone sequence and create germline for clone and perform V gene somatic mutation analysis.
#' dependencies: stringr; dplyr, scoper, dowser, shazam, alakazam, ggplot2, openxlsx, parallel.
#' cite: Jensen C, Sumner J, Kleinstein S, Hoehn K (2024). “Inferring B Cell Phylogenies from Paired H and L Chain BCR Sequences with Dowser.” The Journal of Immunology. doi:10.4049/jimmunol.2300851 https://doi.org/10.4049/jimmunol.2300851, https://doi.org/10.4049/jimmunol.2300851.
#' @export
#'
#' @import dplyr
#' @importFrom scoper hierarchicalClones
#' @importFrom scoper identicalClones
#' @importFrom scoper spectralClones
#' @importFrom dowser createGermlines
#' @importFrom dowser readIMGT
#' @importFrom shazam observedMutations

scFindClones <- function(db,
                         output = TRUE,
                         output_folder = "VDJ_Clones",
                         analysis_name = "All_sequences",
                         seq_type = c("Ig", "TCR"),
                         organism = c("human", "mouse", "rabbit", "rat", "rhesus_monkey"),
                         only_heavy = FALSE,
                         na.heavy.rm = FALSE,
                         igblast = c("none", "filtered_light", "all"),
                         update_c_call = c("all", "filtered light", "none"),
                         clean_LC = TRUE,
                         split_by_light = TRUE,
                         update_germline = TRUE,
                         SHM = TRUE,
                         full_seq_aa = TRUE,
                         #passed to Scoper:
                         method = c("changeo", "hierarchical", "identical", "spectral"), spectral_method = c("novj", "vj"),
                         threshold = c(0.12, 0.15),
                         cell_id_scoper = "cell_id",
                         #passed to resolveMultiContigs():
                         tech = NULL,
                         split.by = NULL,
                         #passed to runIgblastn()/runAssignGenes()/reconstructFullVDJ/createGermline():
                         igblast_dir = "~/share/igblast",
                         imgt_dir = "~/share/germlines/imgt/",
                         CH1_AA = list(IGHG1 = "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV",
                                       IGHG2 = "ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSNFGTQTYTCNVDHKPSNTKVDKTV",
                                       IGHG3 = "ASTKGPSVFPLAPCSRSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYTCNVNHKPSNTKVDKRV",
                                       IGHG4 = "ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTKTYTCNVDHKPSNTKVDKRV",
                                       IGHA1 = "ASPTSPKVFPLSLCSTQPDGNVVIACLVQGFFPQEPLSVTWSESGQGVTARNFPPSQDASGDLYTTSSQLTLPATQCLAGKSVTCHVKHYTNPSQDVT",
                                       IGHA2 = "ASPTSPKVFPLSLDSTPQDGNVVVACLVQGFFPQEPLSVTWSESGQNVTARNFPPSQDASGDLYTTSSQLTLPATQCPDGKSVTCHVKHYTNSSQDVT",
                                       IGHM = "GSASAPTLFPLVSCENSPSDTSSVAVGCLAQDFLPDSITFSWKYKNNSDISSTRGFPSVLRGGKYAATSQVLLPSKDVMQGTDEHVVCKVQHPNGNKEKNVPLPV",
                                       IGHD = "APTKAPDVFPIISGCRHPKDNSPVVLACLITGYHPTSVTVTWYMGTQSQPQRTFPEIQRRDSYYMTSSQLSTPLQQWRQGEYKCVVQHTASKSKKEIF",
                                       IGHE = "ASTQSPSVFPLTRCCKNIPSNATSVTLGCLATGYFPEPVMVTWDTGSLNGTTMTLPATTLTLSGHYATISLLTVSGAWAKQMFTCRVAHTPSSTDWVDNKTFS",
                                       IGKC = "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC",
                                       IGLC1 = "GQPKGQPKANPTVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADGSPVKAGVETTKPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS",
                                       IGLC2 = "GQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS",
                                       IGLC3 = "GQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHKSYSCQVTHEGSTVEKTVAPTECS",
                                       IGLC6 = "GQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVKVAWKADGSPVNTGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPAECS",
                                       IGLC7 = "GQPKAAPSVTLFPPSSEELQANKATLVCLVSDFNPGAVTVAWKADGSPVKVGVETTKPSKQSNNKYAASSYLSLTPEQWKSHRSYSCRVTHEGSTVEKTVAPAECS"),
                         #used for final recap excel sheet
                         shared.tech = list(scRNAseq = c("10X", "BD"),
                                            scSanger = c("scPCR", "scCulture")),
                         recap.highlight = NULL,
                         top_columns = NULL,
                         #global:
                         cell_id = "cell_id",
                         locus = "locus",
                         assay = "assay",
                         umi_count = "umi_count",
                         consensus_count = "consensus_count",
                         productive = "productive",
                         complete_vdj = "complete_vdj",
                         junction = "junction",
                         junction_aa = "junction_aa",
                         junc_len = "junction_length",
                         sequence = "sequence",
                         sequence_id = "sequence_id",
                         sequence_alignment  = "sequence_alignment",
                         v_call = "v_call",
                         d_call = "d_call",
                         j_call = "j_call",
                         c_call = "c_call",
                         nproc = 1,
                         verbose = TRUE,
                         fields = NULL){
                            
  start <- Sys.time()
  #suppressMessages(library(dplyr))
  #suppressMessages(library(scoper))
  #suppressMessages(library(dowser))

  seq_type <- match.arg(seq_type)
  
  if(seq_type == "Ig"){
    heavy <- "IGH"
    light <- c("IGK","IGL")
    chains <- c(heavy, light)
    fct_type <- "BCR"
  }
  if(seq_type == "TCR"){
    heavy <- c("TRB","TRD")
    light <- c("TRA","TRG")
    chains <- c(heavy, light)
    fct_type <- "TCR"
  }
  #TODO any change made here must be also implemented at the resolveMultiContigs() step
  
  organism <- match.arg(organism)
  igblast_dir <- ifelse(stringr::str_ends(igblast_dir,"/"), igblast_dir, paste0(igblast_dir, "/"))
  imgt_dir <- ifelse(stringr::str_ends(imgt_dir,"/"), imgt_dir, paste0(imgt_dir, "/"))

  if(!is.null(output_folder)){
    if(!stringr::str_ends(output_folder, "/")){output_folder = paste0(output_folder, "/")}
    if(!dir.exists(output_folder)){
      dir.create(output_folder)
    }
  }

  if(is.null(db) | !sequence %in% colnames(db) | !sequence_id %in% colnames(db)){
    stop(paste0("no proper VDJ file provided"))
  }

  if(only_heavy){
    cell_id_scoper = NULL
    clean_LC = FALSE
    split_by_light = FALSE
  }

  if(clean_LC){
    missing_counts <- setdiff(c(umi_count, consensus_count), colnames(db))
    if(length(missing_counts)>0){
      warning("missing the following columns in provided VDJ file: ", paste(missing_counts), collapse = ", ", " all set to 1")
      for(column in missing_counts){
        db[column] <- 1
      }
    }
  }

  #initiate log file:
  log_file <- paste0(output_folder, analysis_name, "_scFind", fct_type,"Clones.log")
  time_and_log({
    cat("analysis_name: ", analysis_name, "\n")
    cat("organism: ", organism, "\n")
    cat("seq_type: ", seq_type, "\n")
    cat("igblast: ", igblast, "\n")
    cat("clonal analysis method: ", method, "\n")
    cat("clonal analysis threshold(s): ", threshold, "\n")
    cat("update_c_call: ", update_c_call, "\n")
    cat("clean_LC: ", clean_LC, "\n")
    cat("split_by_light: ", split_by_light, "\n")
    cat("update_germline: ", update_germline, "\n")
    cat("SHM: ", SHM, "\n")
    cat("full_seq_aa: ", full_seq_aa, "\n")
  }, verbose = FALSE, time = FALSE, log_file = log_file, log_title = paste0("scFind", fct_type,"Clones"), open_mode = "wt")
  
  ## run initial QC on imput and arguments provided:
  required_collumns <- c(cell_id, locus, junction, junction_aa, junc_len, productive)
  missing_collumns <- setdiff(required_collumns, colnames(db))
  if(length(missing_collumns)>0) {
    stop(paste0("missing the following collumns: ", missing_collumns))
  }

  method <- match.arg(method)
  if(!method %in% c("changeo", "identical", "hierarchical", "spectral")){
    stop("need to choose a method for clustering among one of the following: identical, hierarchical or spectral (see https://scoper.readthedocs.io/en/stable/ for detailed informations) or changeo for the older version of hierarchicalClones")
  }

  igblast <- match.arg(igblast)
  update_c_call <- match.arg(update_c_call)
  if(!update_c_call %in% c("filtered light", "all")){warning("It is highly recommended to update light chains c_call using Blastn to avoid artificial selection of one c_call by igblast when multiple calls have similar scores (set update_c_call = light or all)")}

  if(SHM & !update_germline){
    if(!"germline_alignment_d_mask" %in% colnames(db)){
      message("attempting to run SHM analysis without updating germline after groupîng and without providing a proper germline_alignment_d_mask column, setting update_germline to TRUE")
      update_germline = TRUE
    } else {warning("running SHM analysis with provided germline_alignment_d_mask column, not recommanded as it will not be updated upon clonal clustering")}
  }

  if(update_germline & is.null(imgt_dir)){
    stop("no reference database provided for germline reconstruction")
  }

  if("clone_id" %in% colnames(db)){
    clone_id_log <- paste0("clone_id collumn already present in the object, renamed preexisting_clone_id to avoid issue in defineClonesScoper().\n\n")
    if(verbose){cat(clone_id_lo)} # Write to console
    
    db <- dplyr::rename(db, "preexisting_clone_id" = clone_id)
  }
  
  ##initiate step counter and filename:
  n <- 1 + sum(clean_LC, (igblast %in% c("light", "all")), (update_c_call %in% c("light", "all")), split_by_light, update_germline, SHM, full_seq_aa)
  # create filename
  filename <- paste0(output_folder, analysis_name)
  step <- 0
  
  ## 0. [option 1] run igblast on all sequences: can take around 15min for 300k+ sequences
  if(igblast == "all"){# not the preferred option as you will use computer power and time on sequences you will discard later but can help getting rid of contigs that won't pass igblast anyway...
    step <- step + 1
    message(
      "------------\n",
      "Part ", step," of ", n,": Running IgBlast on all sequences (takes around 5min for 100k sequences... be patient!)\n",
      "------------\n"
    )
    
    time_and_log({
      igblast_results <- runAssignGenes(db,
                                        igblast_dir = igblast_dir,
                                        imgt_dir = imgt_dir,
                                        sequence = sequence, 
                                        sequence_id = sequence_id,
                                        log = FALSE, 
                                        log_file = log_connection)
      }, verbose = verbose, log_file = log_file, log_title = paste0("Part ", step," of ", n,": Running IgBlast on all sequences"), open_mode = "a")
    
    VDJ_db <- igblast_results[["pass"]]
    VDJ_db$c_call_igblast <- VDJ_db$c_call
    
    failed_VDJ_db <- igblast_results[["fail"]]
    
    rm(igblast_results)
    
    filename <- paste0(filename, "_igblast_db-pass")
    if(!(update_c_call %in% c("heavy", "all")) & !split_by_light){# we only save this file if no other analysis is performed
      readr::write_tsv(VDJ_db, file = paste0(filename, ".tsv.gz"))
    }
    filename_fail <- paste0(filename, "_igblast_db-fail")
    readr::write_tsv(failed_VDJ_db, file = paste0(filename_fail, ".tsv.gz"))
    
    nb_failed_sequences <- nrow(failed_VDJ_db)
    nb_failed_cells <- length(unique(failed_VDJ_db$cell_id))
    
    if(verbose){
      cat("nb submitted: ", nrow(db), " contigs from ", length(unique(db$cell_id)), " cells\n")
      cat("nb pass: ", nrow(VDJ_db), " contigs from ", length(unique(VDJ_db$cell_id)), " cells\n")
      cat("nb failed: ", nb_failed_sequences, " contigs from ", nb_failed_cells, " cells\n")
    }
    
    h_db <- dplyr::filter(VDJ_db, (!!rlang::sym(locus) %in% heavy))
    l_db <- dplyr::filter(VDJ_db, !(!!rlang::sym(locus) %in% heavy))
  } else {
    h_db <- dplyr::filter(db, (!!rlang::sym(locus) %in% heavy))
    l_db <- dplyr::filter(db, !(!!rlang::sym(locus) %in% heavy))
  }

  # check if remaining heavy chain multiplets:
  if(any(duplicated(h_db[[cell_id]]))){
    if(verbose){cat("remaining cell_id with multiple HC: running resolveMultiHC()")}
    time_and_log({
      h_db <- resolveMultiContigs(h_db, 
                                  split.by = split.by, 
                                  seq_type = seq_type,
                                  resolve_chain = "heavy",
                                  assay = "assay", 
                                  resolve_multi_CDR3 = TRUE, 
                                  use_clone = FALSE,
                                  analysis_name = analysis_name,
                                  output = output, 
                                  output_folder = paste0(output_folder, "/VDJ_QC"),
                                  second_columns = c(sequence_id, locus, umi_count, consensus_count, sequence, v_call, d_call, j_call, c_call, junction, junction_aa, productive, complete_vdj),
                                  cell_id = cell_id, 
                                  sequence_id = sequence_id, 
                                  locus = locus, 
                                  consensus_count = consensus_count, 
                                  umi_count = umi_count, 
                                  v_call = v_call, 
                                  j_call = j_call, 
                                  c_call = c_call, 
                                  junction_aa = junction_aa,
                                  productive = productive, 
                                  complete_vdj = complete_vdj)
    }, verbose = verbose, log_file = log_file, log_title = "running resolveMultiHC()", open_mode = "a")
  }

  # check if remaining NAs in CDR3s:
  if(any(is.na(db[[junction]]))){
    failed_junction <- db[FALSE,]
    filename_missing_junction <- paste0(output_folder, analysis_name, "_VDJ_missing_junction")
    
    if(any(is.na(h_db[[junction]]))){
      log_missing_h_junction <- paste0("removing ", table(is.na(h_db[[junction]]))[2]," heavy chain contigs with no identified Junction (see: ", filename_missing_junction, ".tsv.gz).\n")
      if(verbose){cat(log_missing_h_junction)}
      time_and_log({
        cat(log_missing_h_junction) 
      }, verbose = verbose, log_file = log_file, log_title = "removing HC missing CDR3s", open_mode = "a")
      failed_junction <- dplyr::bind_rows(failed_junction, dplyr::filter(h_db, is.na(!!rlang::sym(junction))))
      h_db <- dplyr::filter(h_db, !is.na(!!rlang::sym(junction)))
    }
    
    if(any(is.na(l_db[[junction]]))){
      log_missing_l_junction <- paste0("removing ", table(is.na(l_db[[junction]]))[2]," light chain contigs with no identified Junction (see: ", filename_missing_junction, ".tsv.gz).\n")
      if(verbose){cat(log_missing_l_junction)}
      time_and_log({
        cat(log_missing_l_junction) 
      }, verbose = verbose, log_file = log_file, log_title = "removing HC missing CDR3s", open_mode = "a")
      failed_junction <- dplyr::bind_rows(failed_junction, dplyr::filter(l_db, is.na(!!rlang::sym(junction))))
      l_db <- dplyr::filter(l_db, !is.na(!!rlang::sym(junction)))
    }
    nb_total_NA_junction <- nrow(failed_junction)
    readr::write_tsv(failed_junction, file = paste0(filename_missing_junction, ".tsv.gz"))
  }
  nb_total_h_sequences <- nrow(h_db)
  nb_total_l_sequences <- nrow(l_db)

  VDJ_db <- dplyr::bind_rows(h_db, l_db)
  nb_total_cells <- length(unique(VDJ_db[[cell_id]]))

  rm(db, h_db, l_db)
  
  other_chains_db <- VDJ_db %>%
    dplyr::filter(!locus %in% chains)
  VDJ_db <- VDJ_db %>%
    dplyr::filter(locus %in% chains)

  # Create the bcr_info or tcr_info column
  VDJ_db <- VDJ_db %>%
    dplyr::group_by(!!rlang::sym(cell_id)) %>%
    dplyr::mutate(
      !!rlang::sym(paste0(tolower(fct_type),"_info")) := dplyr::case_when(
        any(!!rlang::sym(locus) %in% heavy) & any(!(!!rlang::sym(locus) %in% heavy)) ~ "full",
        all(!!rlang::sym(locus) %in% heavy) ~ "heavy_only",
        all(!(!!rlang::sym(locus) %in% heavy)) ~ "light_only",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::ungroup()

  # Check / Remove cells without heavy chain.
  # Technically they can't be used in clustering analysis but they won't interfere in the following steps and will end up with a NA value for clone_id.
  LC_only <- dplyr::filter(VDJ_db, !!rlang::sym(paste0(tolower(fct_type),"_info")) == "light_only")
  if(na.heavy.rm & nrow(LC_only)>0){
    if(nrow(LC_only)>0){
      VDJ_db <- dplyr::filter(VDJ_db, !!rlang::sym(paste0(tolower(fct_type),"_info")) != "light_only")
      filename_LC_only <- paste0(output_folder, analysis_name, "_VDJ_LC_only")
      readr::write_tsv(LC_only_cloned_VDJ_db, file = paste0(filename_LC_only, ".tsv.gz"))
      LC_only_log <- paste0(paste0(length(unique(LC_only[[cell_id]])), " cells out of ", nb_total_cells," are missing a heavy chain and will be removed (see: ",filename_LC_only, ".tsv). \n"))
      if(verbose){cat(LC_only_log)}
      time_and_log({
        cat(LC_only_log) 
      }, verbose = verbose, log_file = log_file, log_title = "Removing cells missing HC", open_mode = "a")
    }
  } else {
    if(nrow(LC_only)>0){
      LC_only_log <- cat(paste0(length(unique(LC_only[[cell_id]])), " cells out of ", nb_total_cells," are missing a heavy chain, will end up with clone_id = NA. \n"))
      if(verbose){cat(LC_only_log)}
      time_and_log({
        cat(LC_only_log) 
      }, verbose = verbose, time = FALSE, log_file = log_file, log_title = "Checking cells missing HC", open_mode = "a")
    }
  }
  
  ## 1. perform clonal partitioning: [note: hierarchical = <15s for a database of 3000 heavy chains with n=4 processors (MacBookPro M1)]
  
  if(method == "changeo"){ # similar to scoper::hierarchicalClones() in theory but using the old python implementation of changeo hierarchical clustering
    step <- step + 1
    message(
      "------------\n",
      "Part ", step,  " of ", n,": Starting clonal partitioning using DefineClones.py from Changeo\n",
      "------------\n"
    )
    
    h_db <- VDJ_db %>%
      dplyr::filter(!!rlang::sym(locus) %in% heavy)
    l_db <- VDJ_db %>%
      dplyr::filter(!(!!rlang::sym(locus) %in% heavy))

    out_temp_folder <- paste0(output_folder, "temp/")
    if(!dir.exists(out_temp_folder)){
      dir.create(out_temp_folder)
    }
    DefineClones_input_file <- paste0(out_temp_folder, "All_seq.tsv")
    DefineClones_input_file <- gsub(" ", "\\\ ", DefineClones_input_file, fixed = TRUE) #'to remove any blank in file_path
    DefineClones_output_file <- paste0(out_temp_folder, "All_seq_clone-pass.tsv")
    DefineClones_output_file <- gsub(" ", "\\\ ", DefineClones_output_file, fixed = TRUE) #'to remove any blank in file_path
    
    for(i in seq_along(threshold)){
      used_threshold <- threshold[i]
      if(verbose){cat("using threshold: ", used_threshold, ". ")}
      
      if("clone_id" %in% colnames(h_db)){
        h_db$clone_id <- NULL
        l_db$clone_id <- NULL
      }
      
      time_and_log({
        cat(paste0("using threshold: ", used_threshold, "\n"))
        readr::write_tsv(h_db, file = DefineClones_input_file)
        messages_DefineClones <- system2("DefineClones.py", args = paste0("-d ", DefineClones_input_file, " --act set --model ham --norm len --dist ", used_threshold), stderr = TRUE, stdout = TRUE)
        cat(paste(messages_DefineClones, collapse = "\n"))
        cloned_h_db <- readr::read_tsv(DefineClones_output_file, show_col_types = FALSE)
      }, verbose = verbose, log_file = log_file, log_title = "DefineClones", open_mode = "a")
      
      cloned_h_db[[paste0("h_clone_id_", used_threshold)]] <- cloned_h_db[["clone_id"]]
      
      h_db <- h_db %>%
        dplyr::left_join(
          dplyr::select(cloned_h_db, all_of(c(cell_id, "clone_id", paste0("h_clone_id_", used_threshold)))), by = join_by(!!rlang::sym(cell_id))
        )
      l_db <- l_db %>%
        dplyr::left_join(
          dplyr::select(cloned_h_db, all_of(c(cell_id, "clone_id", paste0("h_clone_id_", used_threshold)))), by = join_by(!!rlang::sym(cell_id))
        )
    }
    cloned_VDJ_db <- h_db %>%
      dplyr::bind_rows(l_db)

    unlink(out_temp_folder, recursive = TRUE)
  }

  if(method == "hierarchical"){ # preferred method for clustering of BCR data, multiple threshold can be provided
    step <- step + 1
    message(
      "------------\n",
      "Part ", step," of ", n,": Starting clonal partitioning using the ", method, " method from Scoper package\n",
      "------------\n"
    )
    
    for(i in seq_along(threshold)){
      used_threshold <- threshold[i]
      if(verbose){cat("using threshold: ", used_threshold, ". ")}
      if("clone_id" %in% colnames(VDJ_db)){
        VDJ_db$clone_id <- NULL
      }
      time_and_log({
        cat(paste0("using threshold: ", used_threshold))
        clonal_analysis <- scoper::hierarchicalClones(VDJ_db,
                                                      threshold = used_threshold,
                                                      method = "nt",
                                                      normalize = "len",
                                                      linkage = "single",
                                                      cell_id = cell_id_scoper,
                                                      locus = locus,
                                                      only_heavy = TRUE,
                                                      split_light = FALSE,
                                                      junction = junction,
                                                      v_call = v_call,
                                                      j_call = j_call,
                                                      clone = "clone_id",
                                                      fields = fields,
                                                      first = FALSE,
                                                      cdr3 = FALSE,
                                                      mod3 = FALSE,
                                                      max_n = 0,
                                                      nproc = nproc,
                                                      verbose = FALSE,
                                                      log = NULL,
                                                      summarize_clones = TRUE)
        if(output){
          pdf(file = paste0(output_folder, analysis_name, "_IGH_clustering_", used_threshold, ".pdf"))
          scoper::plot(clonal_analysis, binwidth=0.02)
          dev.off()
        }
        VDJ_db <- scoper::as.data.frame(clonal_analysis)
        VDJ_db[[paste0("h_clone_id_", used_threshold)]] <- VDJ_db[["clone_id"]]
      }, verbose = verbose, log_file = log_file, log_title = "hierarchicalClones", open_mode = "a")
    }
    cloned_VDJ_db <- VDJ_db
  }
  
  if(method == "identical"){# for TCR mostly
    step <- step + 1
    message(
      "------------\n",
      "Part ", step," of ", n,": Starting clonal partitioning using the ", method, " method from Scoper package\n",
      "------------\n"
    )
    used_threshold <- 0
    time_and_log({
      cloned_VDJ_db <- scoper::identicalClones(VDJ_db,
                                               method = "nt",
                                               cell_id = cell_id_scoper,
                                               locus = locus,
                                               only_heavy = TRUE,
                                               split_light = FALSE,
                                               junction = junction,
                                               v_call = v_call,
                                               j_call = j_call,
                                               clone = "clone_id",
                                               fields = fields,
                                               first = FALSE,
                                               cdr3 = FALSE,
                                               mod3 = FALSE,
                                               max_n = 0,
                                               nproc = nproc,
                                               verbose = FALSE,
                                               log = NULL,
                                               summarize_clones = FALSE)
    cloned_VDJ_db[[paste0("h_clone_id")]] <- cloned_VDJ_db[["clone_id"]]
    }, verbose = verbose, log_file = log_file, log_title = "identicalClones", open_mode = "a")
  }
  if(method == "spectral"){# !!untested yet; memory intensive for big datasets, will only use the first provided threshold...
    message(
      "------------\n",
      "Part ", step," of ", n,": Starting clonal partitioning using the ", method, " method from Scoper package\n",
      "------------\n"
    )
    used_threshold <- threshold[1]
    time_and_log({
      clonal_analysis <- scoper::spectralClones(VDJ_db,
                                                method = spectral_method,
                                                germline = "germline_alignment",
                                                sequence = "sequence_alignment",
                                                junction = junction,
                                                v_call = v_call,
                                                j_call = j_call,
                                                clone = "clone_id",
                                                fields = fields,
                                                cell_id = cell_id_scoper,
                                                locus = locus,
                                                only_heavy = TRUE,
                                                split_light = TRUE,
                                                targeting_model = NULL,
                                                len_limit = NULL,
                                                first = FALSE,
                                                cdr3 = FALSE,
                                                mod3 = FALSE,
                                                max_n = 0,
                                                threshold = used_threshold,
                                                base_sim = 0.95,
                                                iter_max = 1000,
                                                nstart = 1000,
                                                nproc = 1,
                                                verbose = FALSE,
                                                log = NULL,
                                                summarize_clones = TRUE)
      if(output){
        pdf(file = paste0(output_folder, analysis_name, "_IGH_clustering_", used_threshold, ".pdf"))
        scoper::plot(clonal_analysis, binwidth=0.02)
        dev.off()
      }
      cloned_VDJ_db <- scoper::as.data.frame(clonal_analysis)
      cloned_VDJ_db[[paste0("h_clone_id_", used_threshold)]] <- cloned_VDJ_db[["clone_id"]]
    }, verbose = verbose, log_file = log_file, log_title = "identicalClones", open_mode = "a")
  }

  cloned_VDJ_db <- dplyr::relocate(cloned_VDJ_db, clone_id, .after = !!rlang::sym(cell_id))

  filename <- paste0(filename, "_VDJ_clone-pass")
  if(!clean_LC){
    readr::write_tsv(cloned_VDJ_db, file = paste0(filename, ".tsv.gz"))
  }
  rm(VDJ_db)

  ## 2. resolve cases of light chains multiplets and split clonal groups based on light chains:
  if(clean_LC){
    # technically, dowser::resolveLightChains can deal with multiple light chains but choice of dominant light chain doesn't take into account umi_counts...
    # we use her a similar approach to preprocess the best light chain(s) for a given cell_id and then run the preprocess db through dowser::resolveLightChains for final clustering and dealing with cases of missing light chains.
    step = step + 1
    message(
      "------------\n",
      "Part ", step," of ", n,": Selecting best light chain for each cell\n",
      "------------\n"
    )
    time_and_log({
      cloned_VDJ_db <- resolveMultiContigs(cloned_VDJ_db, 
                                           split.by = split.by, 
                                           seq_type = seq_type,
                                           resolve_chain = "light",
                                           assay = assay, 
                                           resolve_multi_CDR3 = TRUE, 
                                           use_clone = TRUE,
                                           analysis_name = analysis_name,
                                           output = output, 
                                           output_folder = output_folder,
                                           second_columns = c(sequence_id, locus, umi_count, consensus_count, sequence, v_call, d_call, j_call, c_call, junction, junction_aa, productive, complete_vdj),
                                           cell_id = cell_id, 
                                           sequence_id = sequence_id, 
                                           locus = locus, 
                                           consensus_count = consensus_count, 
                                           umi_count = umi_count, 
                                           v_call = v_call, 
                                           j_call = j_call, 
                                           c_call = c_call, 
                                           junction = junction, 
                                           junction_aa = junction_aa,
                                           productive = productive, 
                                           complete_vdj = complete_vdj, 
                                           clone_id = "clone_id", 
                                           nproc = nproc)
    }, verbose = verbose, log_file = log_file, log_title = "resolving LC multiplets", open_mode = "a")
    
    filename <- paste0(filename, "_LCresolved")
    if(!(igblast %in% c("light", "all")) & !(update_c_call %in% c("light", "all"))){
      readr::write_tsv(cloned_VDJ_db, file = paste0(filename, ".tsv.gz"))
    }
  }

  ## 3: [option 2] run IgBlast on light chain only.
  if(igblast == "light"){
    step = step + 1
    message(
      "------------\n",
      "Part ", step, " of ", n,": Running IgBlast on light chain contigs\n",
      "------------\n"
    )
    
    h_db <- dplyr::filter(cloned_VDJ_db, (!!rlang::sym(locus) %in% heavy))
    l_db <- dplyr::filter(cloned_VDJ_db, !(!!rlang::sym(locus) %in% heavy))

    nb_light_submitted <- nrow(l_db)
    nb_cells_ini <- length(unique(l_db$cell_id))
      
    time_and_log({
      igblast_results <- runAssignGenes(l_db,
                                        igblast_dir = igblast_dir,
                                        imgt_dir = imgt_dir,
                                        sequence = sequence, 
                                        sequence_id = sequence_id,
                                        log = FALSE, 
                                        log_file = log_connection)
    }, verbose = verbose, log_file = log_file, log_title = "Running IgBlast on light chain contigs", open_mode = "a")
    
    l_db <- igblast_results[["pass"]]
    l_db$c_call_igblast <- l_db$c_call

    failed_l_db <- igblast_results[["fail"]]

    rm(igblast_results)

    filename <- paste0(filename, "_LCigblast_db-pass")
    if(!(update_c_call %in% c("light", "all")) & !split_by_light){# we only save this file if no other analysis is performed
      readr::write_tsv(l_db, file = paste0(filename, ".tsv.gz"))
    }
    filename_fail <- paste0(filename, "_LCigblast_db-fail")
    readr::write_tsv(failed_l_db, file = paste0(filename_fail, ".tsv.gz"))

    nb_failed_sequences <- nrow(failed_l_db)
    nb_failed_cells <- length(unique(failed_l_db$cell_id))
    
    if(verbose){
      cat("nb submitted: ", nb_light_submitted, " contigs from ", nb_cells_ini, " cells\n")
      cat("nb pass: ", nrow(l_db), " contigs from ", length(unique(l_db$cell_id)), " cells\n")
      cat("nb failed: ", nb_failed_sequences, " contigs from ", nb_failed_cells, " cells\n")
    }

    cloned_VDJ_db <- dplyr::bind_rows(h_db, l_db)
    cloned_VDJ_db <- dplyr::relocate(cloned_VDJ_db, !!rlang::sym(cell_id), .before = !!rlang::sym(sequence_id))
    relocate_columns <- intersect(colnames(cloned_VDJ_db), c(orig.ident, assay, !!rlang::sym(consensus_count), !!rlang::sym(umi_count)))
    if(length(relocate_columns)>0){
      cloned_VDJ_db <- dplyr::relocate(cloned_VDJ_db, !!rlang::sym(relocate_columns), .after = !!rlang::sym(sequence_id))
    }
    cloned_VDJ_db <- dplyr::relocate(cloned_VDJ_db, !!rlang::sym(c_call), .after = !!rlang::sym(j_call))
  }
  
  ## 4. update c_call for selected chains: choice of all or light chain only (heavy chains should have analyzed previouslyin this case)
  if(update_c_call == "light"){
    
    # colname c_call_igblast should already be present as c_call has been updated for heavy chain in prior steps
    step = step + 1
    message(
      "------------\n",
      "Part ", step," of ", n,": Updating c_call for filtered light chains contigs\n",
      "------------\n"
    )
    
    h_db <- dplyr::filter(cloned_VDJ_db, (!!rlang::sym(locus) %in% heavy))
    l_db <- dplyr::filter(cloned_VDJ_db, !(!!rlang::sym(locus) %in% heavy))

    time_and_log({
      l_db <- runBlastnC(l_db,
                         igblast_dir = igblast_dir,
                         organism = organism,
                         seq_type = seq_type,
                         sequence = sequence,
                         sequence_id = sequence_id)
    }, verbose = verbose, log_file = log_file, log_title = "Updating c_call for filtered light chains contigs", open_mode = "a")

    l_db <- dplyr::relocate(l_db, !!rlang::sym(c_call), .after = !!rlang::sym(j_call))

    cloned_VDJ_db <- dplyr::bind_rows(h_db, l_db)
    
  }
  if(update_c_call == "all"){
    step <- step + 1
    message(
      "------------\n",
      "Part ", step," of ", n,": Updating c_call for all contigs\n",
      "------------\n"
    )
    
    time_and_log({
      cloned_VDJ_db <- runBlastnC(cloned_VDJ_db,
                                  igblast_dir = igblast_dir,
                                  organism = organism,
                                  seq_type = seq_type,
                                  sequence = sequence,
                                  sequence_id = sequence_id)
    }, verbose = verbose, log_file = log_file, log_title = "Updating c_call for all contigs", open_mode = "a")
    
    cloned_VDJ_db <- dplyr::relocate(cloned_VDJ_db, !!rlang::sym(c_call), .after = !!rlang::sym(j_call))
  }

  ## 5. update clonal groups based on selected light chains:
  if(split_by_light){
    step = step + 1
    message(
      "------------\n",
      "Part ", step," of ", n,": Spliting clonal groups according to selected light chains\n",
      "------------\n"
    )
    
    heavy_cloned_VDJ_db <- cloned_VDJ_db %>%
      dplyr::filter(!!rlang::sym(paste0(tolower(fct_type),"_info")) %in% c("full", "heavy_only"))
    light_only_db <- cloned_VDJ_db %>%
      dplyr::filter(!!rlang::sym(paste0(tolower(fct_type),"_info")) == "light_only")

    #heavy_cloned_VDJ_db <- dowser::resolveLightChains(heavy_cloned_VDJ_db,
    #                                                  nproc = nproc,
    #                                                  minseq = 1,
    #                                                  locus = locus,
    #                                                  heavy = heavy,
    #                                                  id = sequence_id,
    #                                                  seq = "sequence_alignment",
    #                                                  clone = "clone_id",
    #                                                  cell = cell_id,
    #                                                  v_call = v_call,
    #                                                  j_call = j_call,
    #                                                  junc_len = junc_len,
    #                                                  nolight = "missing",
    #                                                  pad_ends = TRUE)

    if(seq_type == "Ig"){
      time_and_log({
        heavy_cloned_VDJ_db <- resolveLightChains2(heavy_cloned_VDJ_db,
                                                   nproc = nproc,
                                                   minseq = 1,
                                                   locus = locus,
                                                   heavy = "IGH",
                                                   id = sequence_id,
                                                   seq = "sequence_alignment",
                                                   clone = "clone_id",
                                                   cell = cell_id,
                                                   v_call = v_call,
                                                   j_call = j_call,
                                                   junc_len = junc_len,
                                                   nolight = "missing",
                                                   pad_ends = TRUE)
      }, verbose = verbose, log_file = log_file, log_title = "Spliting clonal groups according to selected light chains", open_mode = "a")
    }
    if(seq_type == "TCR"){
      time_and_log({
        # needs to adapt to the fact that resolveLightChains() requires a unique value for heavy chain
        heavy_cloned_VDJ_db <- heavy_cloned_VDJ_db %>% 
          dplyr::group_by(!!rlang::sym(cell_id)) %>% 
          dplyr::mutate(heavy_chain = case_when(
            any(!!rlang::sym(locus) == "TRB") ~ "TRB",
            any(!!rlang::sym(locus) == "TRD") ~ "TRD",
            TRUE ~ NA_character_
          )) %>% 
          dplyr::ungroup() %>%
          dplyr::group_by(heavy_chain) %>% 
          dplyr::group_modify(
            ~ resolveLightChains2(.x,
                                  nproc = nproc,
                                  minseq = 1,
                                  locus = locus,
                                  heavy = unique(.x$heavy_chain),
                                  id = sequence_id,
                                  seq = "sequence_alignment",
                                  clone = "clone_id",
                                  cell = cell_id,
                                  v_call = v_call,
                                  j_call = j_call,
                                  junc_len = junc_len,
                                  nolight = "missing",
                                  pad_ends = TRUE)
          ) %>% 
          dplyr::ungroup()
      }, verbose = verbose, log_file = log_file, log_title = "Spliting clonal groups according to selected light chains", open_mode = "a")
    }
    
    cloned_VDJ_db <- heavy_cloned_VDJ_db %>%
      dplyr::bind_rows(light_only_db)

    # update clone_id
    cloned_VDJ_db[["clone_id"]] <- cloned_VDJ_db[["clone_subgroup_id"]]
    # update clone_id
    cloned_VDJ_db[[paste0("l_subgroup_h_clone_id_", used_threshold)]] <- cloned_VDJ_db[["clone_subgroup_id"]]
    cloned_VDJ_db <- dplyr::relocate(cloned_VDJ_db, clone_id, .after = !!rlang::sym(cell_id))
    cloned_VDJ_db <- dplyr::relocate(cloned_VDJ_db, !!rlang::sym(paste0("l_subgroup_h_clone_id_", used_threshold)), .after = !!rlang::sym(paste0("h_clone_id_", used_threshold)))
    # remove some temporary columns added by dowser::resolveLightChains():
    cloned_VDJ_db <- dplyr::select(cloned_VDJ_db, !(c("vj_gene", "vj_alt_cell",	"clone_subgroup",	"clone_subgroup_id")))

    filename <- paste0(filename, "_light-pass")
    readr::write_tsv(cloned_VDJ_db, file = paste0(filename, ".tsv.gz"))
  }

  ## 6. update germline sequences:
  if(update_germline){
    step = step + 1
    message(
      "------------\n",
      "Part ", step," of ", n,": Updating germline alignments based on clonal groups\n",
      "------------\n"
    )
    
    #TODO would need to incorporate TiGER prediction of germline here

    # to avoid issues (notably missing allele as compared to SevenBridges database) it is recommended to regularly update the IMGT database as follow: https://changeo.readthedocs.io/en/stable/examples/igblast.html
    reference_dir <- paste0(imgt_dir, organism, "/vdj/")
    reference_origin <- paste0("user provided: ", reference_dir)
    
    time_and_log({
      cat(paste0("Loading user provided reference database: ", reference_dir, "\n"))
      reference <- dowser::readIMGT(reference_dir)
      
      nb_submitted_seq <- nrow(cloned_VDJ_db)
      
      heavy_cloned_VDJ_db <- cloned_VDJ_db %>%
        dplyr::filter(!is.na(clone_id))
      light_only_db <- cloned_VDJ_db %>%
        dplyr::filter(is.na(clone_id))
      
      cat("Updating germline alignments for cells with clone_id\n")
      # [note] using here force_trim = TRUE to save some sequences that would otherwise fail germline reconstruction. Likely indicate misalignment of the data or PCR hybrid...
      heavy_cloned_VDJ_db <- dowser::createGermlines(heavy_cloned_VDJ_db,
                                                     references = reference,
                                                     locus = locus,
                                                     trim_lengths = TRUE,
                                                     force_trim = TRUE,
                                                     nproc = nproc,
                                                     v_call = v_call,
                                                     d_call = d_call,
                                                     j_call = j_call,
                                                     amino_acid = FALSE,
                                                     id = sequence_id,
                                                     clone = "clone_id",
                                                     na.rm = FALSE,
                                                     fields = fields)
      
      if(nrow(light_only_db)>0){
        cat("Updating germline alignments for cells with only LC (no clone_id)")
        light_only_db <- dowser::createGermlines(light_only_db,
                                                 references = reference,
                                                 locus = locus,
                                                 trim_lengths = TRUE,
                                                 force_trim = TRUE,
                                                 nproc = nproc,
                                                 v_call = v_call,
                                                 d_call = d_call,
                                                 j_call = j_call,
                                                 amino_acid = FALSE,
                                                 id = sequence_id,
                                                 clone = "cell_id",
                                                 na.rm = FALSE,
                                                 fields = fields)
      }
      cloned_VDJ_db <- heavy_cloned_VDJ_db %>%
        dplyr::bind_rows(light_only_db)
      # remove cells which failed gerline reconstruction:
      germ_failed_cloned_VDJ_db <- dplyr::filter(cloned_VDJ_db, is.na(germline_alignment_d_mask))
      nb_failed_germ <- nrow(germ_failed_cloned_VDJ_db)
      nb_failed_germ_cells <- length(unique(germ_failed_cloned_VDJ_db$cell_id))
      
      if(nrow(germ_failed_cloned_VDJ_db)>0){
        filename_germ_failed <- paste0(output_folder, analysis_name, "_VDJ_germ-failed")
        readr::write_tsv(germ_failed_cloned_VDJ_db, file = paste0(filename_germ_failed, ".tsv.gz"))
        germ_failed_log <- paste0(nb_failed_germ, " sequences out of ", nb_submitted_seq," failed germline reconstruction and will be removed (see: ",filename_germ_failed, ".tsv).\n")
        cat(germ_failed_log)
      } else {
        germ_failed_log <- paste0("All ", nb_submitted_seq," sequences passed germline reconstruction. Congrats!!\n")
        cat(germ_failed_log)
      }
    }, verbose = verbose, log_file = log_file, log_title = "Updating germline alignments based on clonal groups", open_mode = "a")
    if(verbose){cat( "\n", germ_failed_log)}
    
    cloned_VDJ_db <- dplyr::filter(cloned_VDJ_db, !is.na(germline_alignment_d_mask))
    filename <- paste0(filename, "_germ-pass")
    readr::write_tsv(cloned_VDJ_db, file = paste0(filename, ".tsv.gz"))
  } else {reference_origin <- NA}
  
  ## 7. calculate somatic mutations::
  if(SHM){
    step = step + 1
    message(
      "------------\n",
      "Part ", step, " of ", n,": running somatic mutation analysis on V genes (you are almost there...)\n",
      "------------\n"
    )
    
    suppressMessages(library(shazam)) # can't make it work without loading the entire package as shazam calls multiple objects from the package
    if(verbose){cat("Calculating mutation counts. ")}
    time_and_log({
      cloned_VDJ_db <- shazam::observedMutations(cloned_VDJ_db, sequenceColumn= sequence_alignment , germlineColumn="germline_alignment_d_mask", regionDefinition=IMGT_V, frequency=FALSE, nproc=nproc)
      cloned_VDJ_db <- shazam::observedMutations(cloned_VDJ_db, sequenceColumn= sequence_alignment, germlineColumn="germline_alignment_d_mask", regionDefinition=IMGT_V, frequency=FALSE, combine=TRUE, nproc=nproc)
    }, verbose = verbose, log_file = log_file, log_title = "Calculating mutation counts", open_mode = "a")
    if(verbose){cat("Calculating mutation frequencies. ")}
    time_and_log({
      cloned_VDJ_db <- shazam::observedMutations(cloned_VDJ_db, sequenceColumn= sequence_alignment, germlineColumn="germline_alignment_d_mask", regionDefinition=IMGT_V, frequency=TRUE, nproc=nproc)
      cloned_VDJ_db <- shazam::observedMutations(cloned_VDJ_db, sequenceColumn= sequence_alignment, germlineColumn="germline_alignment_d_mask", regionDefinition=IMGT_V, frequency=TRUE, combine=TRUE, nproc=nproc)
    }, verbose = verbose, log_file = log_file, log_title = "Calculating mutation frequencies", open_mode = "a")
    
    filename <- paste0(filename, "_shm-pass")
    readr::write_tsv(cloned_VDJ_db, file = paste0(filename, ".tsv.gz"))
  }
  
  ## 8. reconstruct the full VDJ when possible:
  if(full_seq_aa){
    # useful for AlphaFold analysis (!should be run again after clonal assigments)
    # sequences missing too much pb (early VH or Ns in sequence) will failed reconstruction
    
    step = step + 1
    message(
      "------------\n",
      "Part ", step, " of ", n,": Reconstructing full VDJ\n",
      "------------\n"
    )
    
    time_and_log({
      cloned_VDJ_db <- reconstructFullVDJ(cloned_VDJ_db, igblast_dir = igblast_dir, CH1_AA = CH1_AA)
    }, verbose = verbose, log_file = log_file, log_title = "Reconstructing full VDJ", open_mode = "a")
    
    nb_failed_vdj_seq <- nrow(cloned_VDJ_db[is.na(cloned_VDJ_db$full_sequence_aa),])
    nb_failed_vdj_cells <- nrow(unique(cloned_VDJ_db[is.na(cloned_VDJ_db$full_sequence_aa), cell_id]))

    filename <- paste0(filename, "_full-pass")
    readr::write_tsv(cloned_VDJ_db, file = paste0(filename, ".tsv.gz"))
  }

  cloned_VDJ_db <- dplyr::arrange(cloned_VDJ_db, clone_id, !!rlang::sym(cell_id), !!rlang::sym(v_call))

  # Update bcr_info/tcr_info column in case any sequence as been remove along the way:
  cloned_VDJ_db <- cloned_VDJ_db %>%
    group_by(!!rlang::sym(cell_id)) %>%
    mutate(
      !!rlang::sym(paste0(tolower(fct_type),"_info")) := case_when(
        any(!!rlang::sym(locus) %in% heavy) & any(!(!!rlang::sym(locus) %in% heavy)) ~ "full",
        all(!!rlang::sym(locus) %in% heavy) ~ "heavy_only",
        all(!(!!rlang::sym(locus) %in% heavy)) ~ "light_only",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup()

  # add expanded clone info:
  cloned_VDJ_db <- cloned_VDJ_db %>%
    dplyr::group_by(clone_id) %>%
    dplyr::mutate(expanded_clone = n_distinct(!!rlang::sym(cell_id)) > 1) %>%
    dplyr::ungroup()
  
  # reorder some columns:
  top_columns <- unique(c("clone_id", "expanded_clone", "shared_clone", "mu_count", "donor_id", "donor_group", "timepoint", "tissue", "specificity", "assay", top_columns))
  top_columns <- top_columns[top_columns %in% colnames(cloned_VDJ_db)]
  cloned_VDJ_db <- cloned_VDJ_db %>%
    dplyr::relocate(cell_id, .before = "sequence_id") %>%
    dplyr::relocate(all_of(top_columns), .after = "sequence_id")

  # export recap table:
  message(
    "--------------------------------\n",
    "Exporting recap excel workbook :\n",
    "--------------------------------\n"
  )
  
  n_final_heavy <- as.numeric(table(cloned_VDJ_db[[locus]] %in% heavy)[["TRUE"]])
  n_final_cells <- length(unique(cloned_VDJ_db[[cell_id]]))

  if(!only_heavy){
    n_final_light <- as.numeric(table(cloned_VDJ_db[[locus]] %in% heavy)[["FALSE"]])
    h_db <- cloned_VDJ_db %>%
      dplyr::filter(locus %in% heavy)
    n_final_cells_with_full_BCR <- as.numeric(table(h_db[[paste0(tolower(fct_type),"_info")]] == "full")[["TRUE"]])
    pct_full <- round(n_final_cells_with_full_BCR/n_final_cells*100)
  }
  if(is.null(analysis_name)){analysis_name <- "none provided"}
  if(is.null(output_folder)){output_folder <- "working directory"}
  if(is.null(reference_origin)){reference_origin <- "none provided"}
  if(is.null(fields)){fields <- "none"}
  analysis_parameters <- c("analysis name" = analysis_name,
                           "output folder" = output_folder,
                           "nproc" = nproc,
                           "clustering method" = method,
                           "clustering threshold for final clone_id" = used_threshold,
                           "clustering level" = "nucleotide",
                           "clustering normalization" = "length",
                           "split by light" = split_by_light,
                           "germline updated" = update_germline,
                           "reference origin" = reference_origin,
                           "SHM analysis done" = SHM,
                           "additional fields used for grouping" = fields,
                           "submitted IGH contigs" = nb_total_h_sequences,
                           "submitted IGL/K contigs" = nb_total_l_sequences,
                           "submitted unique cells" = nb_total_cells)

  if(method == "spectral"){
    analysis_parameters["clustering level"] <- spectral_method[1]
    analysis_parameters["clustering normalization"] <- NA
  }

  if(!only_heavy){
    final_log_message <- paste0("submitted: ", nb_total_h_sequences , " IGH contigs and ", nb_total_l_sequences, " IGL/IGK contigs with identified CDR3 (", nb_total_cells,  " unique cells).\n")
  } else {
    final_log_message <- paste0("submitted: ", nb_total_h_sequences , " IGH contigs with identified CDR3 (", nb_total_cells,  " unique cells).\n")
  }

  if(igblast == "light"){
    final_log_message <- paste0(final_log_message, nb_failed_sequences, " IGL/K contigs failed igblast.\n")
    analysis_parameters <- c(analysis_parameters, "IGL/K contigs failing igblast" = nb_failed_sequences)
  }
  if(igblast == "all"){
    final_log_message <- paste0(final_log_message, nb_failed_sequences, " contigs from ", nb_failed_cells,  " cells failed igblast.\n")
    analysis_parameters <- c(analysis_parameters, "contigs failing igblast" = nb_failed_sequences)
    analysis_parameters <- c(analysis_parameters, "cells with contig failing igblast" = nb_failed_cells)
  }
  if(update_germline){
    final_log_message <- paste0(final_log_message, nb_failed_germ, " contigs from ", nb_failed_germ_cells,  " cells failed germline reconstruction.\n")
    analysis_parameters <- c(analysis_parameters, "contigs failing germline reconstruction" = nb_failed_germ)
    analysis_parameters <- c(analysis_parameters, "cells with contig failing germline reconstruction" = nb_failed_germ_cells)
  }
  if(full_seq_aa){
    final_log_message <- paste0(final_log_message, nb_failed_vdj_seq, " contigs from ", nb_failed_vdj_cells,  " cells failed full VDJ reconstruction (too many Ns or missing too much of the fwr1 segment).\n")
    analysis_parameters <- c(analysis_parameters, "contigs failing full VDJ reconstruction" = nb_failed_vdj_seq)
    analysis_parameters <- c(analysis_parameters, "cells with contig failing full VDJ reconstruction" = nb_failed_vdj_cells)
  }

  if(!only_heavy){
    final_log_message <- paste0(final_log_message, "final table: ", n_final_heavy, " IGH contigs and ", n_final_light, " IGL/IGK contigs (", n_final_cells, " unique cells). ", pct_full,"% of cells with full BCR. \n")
  } else {
    final_log_message <- paste0(final_log_message, "final table: ", n_final_heavy, " IGH contigs (", n_final_cells, " unique cells). \n")
  }
  if(verbose){cat(final_log_message)}
  time_and_log({
    cat(final_log_message)
  }, verbose = FALSE, time = FALSE, log_file = log_file, log_title = "Final recap", open_mode = "a")
  
  if(!only_heavy){
    analysis_parameters <- c(analysis_parameters,
                             "final IGH contigs" = n_final_heavy,
                             "final IGL/K contigs" = n_final_light,
                             "final unique cells" = n_final_cells,
                             "% with full VDJ" = pct_full)
  } else {
    analysis_parameters <- c(analysis_parameters,
                             "final IGH contigs" = n_final_heavy,
                             "final unique cells" = n_final_cells)
  }
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    OUT <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(OUT, "Cloned VDJ")
    openxlsx::writeData(OUT, sheet= "Cloned VDJ", x=cloned_VDJ_db, colNames = TRUE, rowNames = FALSE)
    openxlsx::addWorksheet(OUT, "Expanded Cloned VDJ")
    openxlsx::writeData(OUT, sheet= "Expanded Cloned VDJ", x=dplyr::filter(cloned_VDJ_db, expanded_clone), colNames = TRUE, rowNames = FALSE)
    if(assay %in% colnames(cloned_VDJ_db)){
      cloned_VDJ_db <- cloned_VDJ_db %>%
        dplyr::group_by(clone_id) %>%
        dplyr::mutate(
          shared_clone = any(assay %in% shared.tech[["scRNAseq"]]) & any(assay %in% shared.tech[["scSanger"]])
        ) %>%
        dplyr::ungroup()
      if(any(cloned_VDJ_db$shared_clone)){
        openxlsx::addWorksheet(OUT, "Shared Cloned VDJ")
        openxlsx::writeData(OUT, sheet= "Shared Cloned VDJ", x=dplyr::filter(cloned_VDJ_db, shared_clone), colNames = TRUE, rowNames = FALSE)
      }
    }
    
    if(!is.null(recap.highlight)){
      missing_highlights <- setdiff(recap.highlight, colnames(cloned_VDJ_db))
      if(length(missing_highlights)>0){
        warning("missing the following recap.highlight column(s): ", paste(missing_highlights, collapse = ", "))
      }
      recap_highlights <- intersect(recap.highlight, colnames(cloned_VDJ_db))
      if(length(recap_highlights)>0){
        for(highlight in recap_highlights){
          levels <- levels(as.factor(cloned_VDJ_db[[highlight]]))
          for(level in levels){
            openxlsx::addWorksheet(OUT, level)
            openxlsx::writeData(OUT, sheet= level, x=dplyr::filter(cloned_VDJ_db, !!rlang::sym(highlight) == level), colNames = TRUE, rowNames = FALSE)
          }
        }
      }
    }
    
    openxlsx::addWorksheet(OUT, "Analysis parameters")
    openxlsx::writeData(OUT, sheet= "Analysis parameters", x=as.data.frame(analysis_parameters), colNames = FALSE, rowNames = TRUE)
    openxlsx::saveWorkbook(OUT, file= paste0(output_folder, analysis_name,"_VDJ_final_recap.xlsx"), overwrite = TRUE)
    readr::write_tsv(cloned_VDJ_db, file = paste0(output_folder, analysis_name,"_VDJ_final.tsv.gz"))
    
    message("see exported recap file: ", getwd(), "/", output_folder, analysis_name,"_VDJ_final_recap.xlsx. \n")
  } else {
    message("Optional: 'openxlsx' not installed — simply saving as tsv files")
    readr::write_tsv(cloned_VDJ_db, file = paste0(output_folder, analysis_name,"_VDJ_final.tsv.gz"))
  }
  
  if(nrow(other_chains_db) > 0){
    cloned_VDJ_db <- cloned_VDJ_db %>%
      dplyr::bind_rows(other_chains_db)
  }
  
  time_and_log({
    print(sessionInfo())
  }, verbose = FALSE, time = FALSE, log_file = log_file, log_title = "session info", open_mode = "a")
  
  end <- Sys.time()
  if(verbose){cat(paste0("Total running time: ", sprintf("%.2f %s", end-start, units(difftime(end, start))), ".\n"))}

  return(cloned_VDJ_db)
}

#### Wrapper function to call scFindClones with BCR data ####
#' Full BCR analysis pipeline for heavy chain curated AIRR single-cell repertoire datasets
#'
#' \code{scFindBCRlones} performs clonal clustering of single cell data
#'
#' @param db              an AIRR formatted dataframe containing heavy and light chain sequences. Should contain only one heavy chain (IGH) per cell_id, if not run resolveMultiHC() first.
#' @param ...             arguments to be passed to scFindClones
#'
#' @return an AIRR formatted dataframe with a new clone_id column and clone attribution for all BCR related data (see scFindClones).
#'
#' @details
#' wrapper to run scFindClones with BCR specific arguments. All other arguments from scFindClones can also be used.
#'
#' @export

scFindBCRClones <- function(db,
                            ...){
  
  cloned_db <- scFindClones(db,
                            seq_type = "Ig",
                            ...)
                                
  return(cloned_db)
}

#### Wrapper function to call scFindClones with TCR data ####
#' Full TCR analysis pipeline for heavy and light chain curated AIRR single-cell repertoire datasets
#'
#' \code{scFindTCRlones} performs clonal clustering of single cell data
#'
#' @param db              an AIRR formatted dataframe containing heavy and light chain sequences. Should contain only one heavy chain (IGH) per cell_id, if not run resolveMultiHC() first.
#' @param igblast         whether to run standalone IgBlast, can be set to c("filtered_light" or "all") if not one of these three values, will be skipped with a warning. [default = "filtered_light" for both 10X and BD: highly recommended for both to avoid issues at the createGermline() or observedMutation() steps due to different references databases used (10X) or missing imgt gaps in the sequence_alignment collumn (Both))]
#' @param method          method to use for scoper::hierarchicalClones(), can be one between: identical, hierarchical or spectral
#' @param threshold       method to use for scoper::hierarchicalClones() or scoper::spectralClones(),
#'                        if multiple thresholds are provided, they will be used successively and results stored in the h_clone_"threshold" column. Final clone_id column will reflect the last threshold used.
#' @param update_c_call   whether to run runBlastnC to correct c calls made by igblast (issues with calls with similar scores); can be set to c("filtered light" or "all") if not one of these three values, will be skipped with a warning. [default = "all", if using "filtered light" should have already been previously performed on heavy chains (see scImportVDJ())]
#' @param clean_LC        whether to resolve cases of multiple light chains
#' @param shared.tech     list of grouped technologies to identify shared clones between scRNA-seq and Sanger sequencing (export in recap teable (seperate sheet) + additional shared_clone column)
#'
#' @return an AIRR formatted dataframe with a new clone_id column and clone attribution for all BCR related data (see scFindClones).
#'
#' @details
#' wrapper to run scFindClones with TCR specific arguments. All other arguments from scFindClones can also be used.
#'
#' @export

scFindTCRClones <- function(db,
                            igblast = c("none", "all"),
                            update_c_call = c("none"), 
                            method = c("changeo", "identical"),
                            threshold = 0,
                            shared.tech = NULL, #for now only expected to be used for scRNAseq datasets
                            clean_LC = FALSE, # also setting up this last two option as FALSE as they haven't yet been properly tested for TCR datasets:
                            ...){
  
  # no need for mutation analysis and clone-based germline inference
  SHM = FALSE
  update_germline = FALSE
  # no need for full VDJ reconstruction (mostly for Alphafold prediction of Ab/Antigen interaction)
  full_seq_aa = FALSE
  
  cloned_db <- scFindClones(db,
                            seq_type = "TCR",
                            igblast = igblast,
                            method = method,
                            threshold = threshold,
                            update_c_call = update_c_call,
                            SHM = SHM,
                            update_germline = update_germline,
                            full_seq_aa = full_seq_aa,
                            shared.tech = shared.tech,
                            ...)
                                
  return(cloned_db)
}

#### Function to add clonotype information to a Seurat object ####
#' Adds AIRR metadata of a seurat object
#'
#' \code{AddAIRRmetadata} adds the immune receptor information to the metadata of a seurat object
#' @param sc        a seurat object on which to input vdj data, should include a cell_id column matching the cell_id column in the vdj_db dataframe.
#' @param vdj_db    an AIRR formatted dataframe containing bcr (heavy and light chains) or tcr (TCRA, TCRB, TCRG or TCRD) sequences. Should contain only one chain for each type per cell_id, if not run resolveMultiHC() first.
#' @param type      either BCR or TCR, by default BCR.
#' @param import    list of columns to import, if set to "all" will import all columns not already included in the seurat object. If "clone_id" is present will calculate clone size and frequencies.
#' @param locus     name of column containing locus values.
#' @param cell_id   name of the column containing cell identifier.
#' @param clone_id  name of the column containing cell identifier.
#' @param split.by  name of column to use to group sequence when calculating clone size and frequencies.
#' @param commun_columns  name of columns commun between different chains. Will not be rename when all other chain specific columns will be renamed "chain_old name"
#' @param chains    description of chain pairing. list names will be used for renaming.
#'
#' @return    a seurat object with additional metadata imported from the vdj_db dataframe.
#'
#' @details
#' The provided df is expected to follow the AIRR format for column nomenclatures, and to contain a "clone_id" column.
#' Cell_ids with duplicated heavy and light chain should also have been excluded during prior steps of the pipeline. By default will be excluded here.
#' Before using AddAIRRmetadata() ensure the cell_ids of the seurat object match the "cell_id" in the AIRR formatted dataframe.
#' By default, if "clone_id" column exists for tcra and tcrb or for heavy_chain this function also calculates the size and frequencies of the clonotypes (as defined in the clone_id column).
#' (use the "groupBy =" argument, by default NULL will return frequency among all cells in the dataset, if set to NULL, no frequency is returned).
#' Grouping in done at the sc level post incorporation of the VDJ info so all collumns from both sc and import can be used.
#'
#' @import dplyr
#' @import purrr
#'
#' @export

addAIRRmetadata <- function(sc, vdj_db = NULL,
                            seq_type = c("Ig", "TCR"),
                            import = c("bcr_info", "tcr_info", "clone_id", "expanded_clone", "sequence_id", "consensus_count", "umi_count", "locus",
                                       "productive", "v_call", "d_call", "j_call", "c_call", "junction", "junction_aa", "junction_length",
                                       "mu_count", "mu_count_cdr_r", "mu_count_cdr_s", "mu_count_fwr_r", "mu_count_fwr_s",
                                       "is.VDJ_doublet", "is.VDJ_doublet.confidence", "is.nonB_VDJ_doublet", "is.nonB_VDJ_doublet.confidence"),
                            locus = "locus",
                            cell_id = "cell_id",
                            clone_id = "clone_id",
                            split.by = NULL,
                            commun_columns = c("bcr_info", "tcr_info", "clone_id", "expanded_clone", "is.VDJ_doublet", "is.VDJ_doublet.confidence", "is.nonB_VDJ_doublet", "is.nonB_VDJ_doublet.confidence"),
                            chains = list(TCR = list(tcra = "TCRA", tcrb = "TCRB", tcrg = "TCRG", tcrd = "TCRD"),
                                         BCR = list(igheavy = "IGH", iglight = c("IGL", "IGK")))){

  #library(dplyr)

  seq_type <- match.arg(seq_type)
  if(!seq_type %in% c("Ig", "TCR")){
    stop("type should be one of TCR or Ig")
  }
  if(seq_type == "Ig"){
    type = "BCR"
    import <- import[!import == "tcr_info"]
    commun_columns <- commun_columns[!commun_columns == "tcr_info"]
  }
  if(seq_type == "TCR"){
    type = "TCR"
    import <- import[!import == "bcr_info"]
    commun_columns <- commun_columns[!commun_columns == "bcr_info"]
  }

  if(!cell_id %in% colnames(vdj_db)){
    stop("missing cell_ids in vdj data")
  }
  if(is.null(split.by)){
    split <- NULL
    vdj_db <- vdj_db %>%
      dplyr::mutate(
        split.by = "all_sequences"
      )
    split.by = "split.by"
  }

  if(inherits(sc, "Seurat")){
    sc_meta <- sc[[]]
  } else {
    stop(sc, " is NOT a Seurat object")
  }
  
  if(!cell_id %in% rownames(sc[[]])){
    sc_meta[[cell_id]] <- rownames(sc[[]])
  }
  if(!any(vdj_db[[cell_id]] %in% sc_meta[[cell_id]])){
    stop("no matching cell_ids between seurat object and vdj data")
  } else {
    cells_to_keep <- sc_meta[[cell_id]]
    vdj_db <- vdj_db %>%
      dplyr::filter(!!rlang::sym(cell_id) %in% cells_to_keep)
  }

  if(any(!import %in% colnames(vdj_db))){
    missing_columns <- import[!import %in% colnames(vdj_db)]
    import <- import[import %in% colnames(vdj_db)]
    warning("the following collumns are missing in the provided vdj data and will not be imported: ", paste(missing_columns, collapse =  "; "))
  }
  if(any(!commun_columns %in% colnames(vdj_db))){
    missing_columns <- commun_columns[!commun_columns %in% colnames(vdj_db)]
    commun_columns <- commun_columns[commun_columns %in% colnames(vdj_db)]
  }

  if(any(import == "all")){
    import <- colnames(vdj_db)
  }
  import <- import[!import == cell_id]

  if(!split.by %in% import){
    import <- c(import, split.by)
  }
  if(!split.by %in% commun_columns){
    commun_columns <- c(commun_columns, split.by)
  }

  selected_chains <- chains[[type]]

  chains_db.lists <- purrr::imap(selected_chains, function(chain, name){

    chain_db <- vdj_db %>%
      dplyr::filter(!!rlang::sym(locus) %in% chain)

    if(is.null(chain_db)){
      warning("missing ", chain," data")
    }
    # flag duplicated cell_id for a given chain:
    # TODO correct this part to filter automatically using ResolveMultichain if needed (when compatible with tcr data)
    if(any(duplicated(chain_db[[cell_id]]))){
      stop("duplicated cell_ids in ",pair[1]," chains database; run ResolveMultiContigs() first")
    }

    chain_db <- chain_db %>%
      dplyr::select(all_of(c(cell_id, import))) %>%
      dplyr::rename_with(
        ~ paste0(name, "_", .), -all_of(c(cell_id, commun_columns))
        )

    return(chain_db)
  })

 sc_vdj_db <- purrr::reduce(chains_db.lists, full_join, by = c(cell_id, commun_columns))

 na_clone_id <- sc_vdj_db %>%
   dplyr::filter(is.na(!!rlang::sym(clone_id))) %>%
   dplyr::mutate(
     clone_size = NA,
     clone_freq = NA
   )

 sc_vdj_db <- sc_vdj_db %>%
   dplyr::filter(!is.na(!!rlang::sym(clone_id))) %>% #Frequencies calculated among non NAs for clone_ids
   dplyr::group_by(!!rlang::sym(split.by)) %>%
   dplyr::mutate(
     group_size = n()
   ) %>%
   dplyr::group_by(!!rlang::sym(clone_id)) %>%
   dplyr::mutate(
     clone_size = n(),  # Count occurrences of clone_id
     clone_freq = n()/group_size # Frequency among group
   ) %>%
   dplyr::rename_with(
     ~ paste0(tolower(type), "_", .), all_of(c(clone_id, "clone_size", "clone_freq"))
   ) %>%
   dplyr::ungroup() %>%
   dplyr::select(!group_size) %>%
   dplyr::bind_rows(na_clone_id)

 if(is.null(split)){
   sc_vdj_db <- sc_vdj_db %>%
     dplyr::select(-!!rlang::sym(split.by))
 }

 sc_vdj_db <- as.data.frame(sc_vdj_db)
 rownames(sc_vdj_db) <- sc_vdj_db$cell_id

 # remove columns common between sc_vdj_db and seurat object:
 seurat_columns <- intersect(colnames(sc), colnames(sc_vdj_db))
 sc_vdj_db <- sc_vdj_db %>%
   dplyr::select(-all_of(seurat_columns))

 # import in seurat object:
 sc <- AddMetaData(sc, sc_vdj_db)
 return(sc)
}


