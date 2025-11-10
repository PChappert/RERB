
#### Function to summarize cell information ####
#' convert an AIRR formatted data frame to prepare for exportAF3json
#'
#'\code{summarizeBCRClones}
#'@param db AIRR formatted data frame with individual rows for heavy (HC) and light chains (LC)
#'@param cell_id name of column containing cell_id values.
#'@param locus name of column containing locus information.
#'@param import_per_chain names of columns to keep for each chain (will appear as "h_col1" and "l_col1" in the exported data frame) 
#'@param import_global names of columns to summarize at the cell level 
#'@param include_meta metadata to keep in the reformatted data frame (default = selected non AIRR metadata)
#'@param selected_meta names of columns to keep in the final data frame
#'@param filter_incomplete_bcr whether to only keep cells with both heavy and light chains
#' @return
#' takes an AIRR formatted data frame with individual rows for heavy (HC) and light chains (LC) and returns a data frame with one row per cell id.
#' columns listed in the import_per_chain argument are renamed as "h_col1" and "l_col1" in the exported data frame.
#' columns listed in the import_global argument are kept as is.
#' any columns listed in the selected_meta argument will be kept as is unless also listed in the default import_per_chain  argument.
#'
#' @details
#' db should be an AIRR formatted database with defined "cell_id" and "locus" columns.
#'
#' @export

summarizeBCRClones <- function(db,
                               cell_id = "cell_id",
                               locus = "locus",
                               import_per_chain = c("v_call", "d_call", "j_call", "junction_aa", "junction_length", "mu_count", 
                                                    "sequence", "full_sequence", "full_sequence_aa",	"full_sequence_fab_aa"),
                               import_global = c("cell_id", "clone_id"), 
                               include_meta = c("selected", "all"),
                               selected_meta = NULL,
                               filter_incomplete_bcr = TRUE){
  
  #suppressMessages(library(dplyr))
  
  #TODO recode function using pivot_wider, as in AddAIRRMetadata?
  
  chains = list(h = "IGH", 
                l = c("IGL", "IGK"))
  
  include_meta <- match.arg(include_meta)
  if(!is.null(selected_meta)){
    include_meta <- "selected"
  }
  
  airr_columns = c("rev_comp", "productive", "locus", "bcr_info",
                   "v_call", "d_call", "j_call", "c_call", "c_call_igblast", "c_call_bd", "c_call_10x", "junction", "junction_aa", 
                   "sequence_alignment", "germline_alignment",
                   "v_cigar", "d_cigar", "j_cigar", "stop_codon", "vj_in_frame", 
                   "junction_length", "np1_length", "np2_length", "v_sequence_start", "v_sequence_end", "v_germline_start", "v_germline_end", 
                   "d_sequence_start", "d_sequence_end", "d_germline_start", "d_germline_end", "j_sequence_start", "j_sequence_end", "j_germline_start", "j_germline_end", 
                   "v_score", "v_identity", "v_support", "d_score", "d_identity", "d_support", "j_score", "j_identity", "j_support", 
                   "fwr1", "fwr2", "fwr3", "fwr4", "cdr1", "cdr2", "cdr3",
                   "v_germline_length", "d_germline_length", "j_germline_length",
                   "h_clone_id", "l_subgroup_h_clone_id", 
                   "c_call_alignment_score", "c_call_pct_match", "c_call_alignment_length", "vj_cell", "germline_alignment_d_mask",     
                   "mu_count_cdr_r", "mu_count_cdr_s", "mu_count_fwr_r", "mu_count_fwr_s", "mu_freq_cdr_r", "mu_freq_cdr_s", "mu_freq_fwr_r", "mu_freq",
                   "missing_v_bp", "missing_j_bp")
  
  ab1_to_airr_columns = c("sort_id", "sequencing_pool_well_id", "sequencing_plate", "sort_well_id", 
                          "sequence_length", "pct_under_30QC_in_trimmed", 
                          "low10QC_original_seq", "low30QC_original_seq",
                          "sanger_low10QC_original_seq", "sanger_low30QC_original_seq",
                          "low10QC_alternate_calls", "low30QC_alternate_calls",
                          "sanger_low10QC_alternate_calls", "sanger_low30QC_alternate_calls", "QC_passed",
                          "sanger_missing_v_bp", "sanger_missing_j_bp", "sanger_comments")
  
  if(filter_incomplete_bcr){
    db <- db %>%
      dplyr::group_by(cell_id) %>%
      dplyr::mutate(
        bcr_info = case_when(
          any(!!rlang::sym(locus) %in% chains[["h"]]) & any(!(!!rlang::sym(locus) %in% chains[["h"]])) ~ "full",
          all(!!rlang::sym(locus) %in% chains[["h"]]) ~ "heavy_only",
          all(!(!!rlang::sym(locus) %in% chains[["h"]])) ~ "light_only",
          TRUE ~ NA_character_
        )) %>%
      dplyr::ungroup() %>%
      dplyr::filter(bcr_info == "full")
  }
  
  if(include_meta == "all"){
    #by default removing all airr columns and ab1toairr columns
    #to keep them use 'selected_meta' argument and include_meta == "selected"
    import_global <- c(import_global, setdiff(colnames(db), c(import_per_chain, import_global, airr_columns, ab1_to_airr_columns)))
  }
  if(include_meta == "selected"){
    meta_by_chain <- intersect(selected_meta, c(airr_columns, ab1_to_airr_columns))
    import_per_chain <- c(import_per_chain, meta_by_chain)
    if(length(meta_by_chain) > 0){
      warning("The following metadata will be imported at the sequence level: ", paste(meta_by_chain, collapse = "; "))
    }
    import_global <- c(import_global, setdiff(selected_meta, c(import_per_chain, import_global, airr_columns, ab1_to_airr_columns)))
  }
  
  cell_infos <- db %>%
    dplyr::filter(locus == "IGH") %>%
    dplyr::select(any_of(c("cell_id", import_global)))
  
  import_per_chain <- intersect(import_per_chain, colnames(db))
  chains_cols <- unlist(lapply(import_per_chain, function(c) c(paste0("h_", c), paste0("l_", c))))
  
  chains_db.lists <- purrr::imap(chains, function(chain, name){
    chain_db <- db %>%
      dplyr::filter(locus %in% chain) %>%
      dplyr::select(all_of(c("cell_id", import_per_chain))) %>%
      dplyr::rename_with(
        ~ paste0(name, "_", .), -all_of(c("cell_id"))
      )
    return(chain_db)
  })
  
  summarized_db <- cell_infos %>%
    dplyr::left_join(
      purrr::reduce(chains_db.lists, full_join, by = c("cell_id")),
      by = join_by("cell_id")
    )
  
  summarized_db <- summarized_db %>%
    dplyr::select(any_of(c(colnames(cell_infos), chains_cols)))
  
  return(summarized_db)
}

#### Function to export a json with requests for Alphafold3 modeling of Antigen/Ig interactions ####
#' export a Alphafold3 server formatted .json file
#'
#' \code{exportAF3json}
#' @param db a dataframe in the summarizeBCRClones() format with column for heavy and light chains full sequence at the aa level (use fullVDJreconstruction())
#' @param antigen_name name of the antigen to model binding to (only used if no request_name is provided)
#' @param antigen_aa sequence of the antigen to model binding to
#' @param request_name a column in the data frame to use for names of Alphafold3 requests (if NULL, defaulting to "cell_id"_"antigen_name"_"syst.date"
#' @param json_name name for the final .json file
#' @param out_folder path of the folder to export to
#' @param cell_id name of column containing cell_id values.
#' @param full_HC_aa name of column containing full AA sequence for the heavy chain.
#' @param full_LC_aa name of column containing full AA sequence for the light chain.
#' @param useStructureTemplate which chain to use the structure template on ("all", for all chains; "none", for no chains; or any combination of "antigen", "HC" and "LC")
#' @return
#' one .json file containing all the requests ready for upload onto the AlphaFold3 server. HC, LC and Ag are submitted as three independant chains.
#'
#' @export

exportAF3json <- function(db,
                          antigen_name,
                          antigen_aa,
                          request_name = NULL,
                          json_name = "job_request",
                          out_folder = "AF3_json_files/",
                          cell_id = "cell_id",
                          full_HC_aa = "h_full_sequence_aa",
                          full_LC_aa = "l_full_sequence_aa",
                          useStructureTemplate = c("all", "antigen", "HC", "LC", "none"),
                          return_db = FALSE){
  
  useStructureTemplate_on <- match.arg(useStructureTemplate, several.ok = TRUE)
  if("all" %in% useStructureTemplate_on){
    useStructureTemplate_on <- c("antigen", "HC", "LC")
  }
  if("none" %in% useStructureTemplate_on){
    useStructureTemplate_on <- c("none")
  }
  
  chains <- c("antigen", "HC", "LC")
  useStructureTemplate_on <- chains %in% useStructureTemplate_on
  names(useStructureTemplate_on) <- chains
  
  if(!dir.exists(out_folder)){
    dir.create(out_folder)
  }
  
  if(is.null(request_name)){
    warning("no provided request name, defaulting to 'cell_id'_'antigen_name'_'date'")
    db$request_name <- paste0(db[[cell_id]], "_", antigen_name,"_", Sys.Date())
  } else {
    if(!request_name %in% colnames(db)){
      warning("provided request name (", request_name,") not in the colnames of provided db, defaulting to 'cell_id'_'antigen_name'_'date'")
      db$request_name <- paste0(db[[cell_id]], "_", antigen_name,"_", Sys.Date())
    } else {
      db <- db %>%
        dplyr::mutate(
          request_name = !!rlang::sym(request_name)
        )
    }
  }
  
  db <- db %>%
    dplyr::filter(!is.na(!!rlang::sym(cell_id)) & !is.na(!!rlang::sym(full_HC_aa)) & !is.na(!!rlang::sym(full_LC_aa)))
  
  json.files <- lapply(1:nrow(db), function(i) {
    list(
      name = db$request_name[i],
      modelSeeds = list(sample(1:1e6, 1)),   # random seed
      sequences = list(
        list(proteinChain = list(
          sequence = antigen_aa,
          count = 1,
          useStructureTemplate = useStructureTemplate_on["antigen"]
        )),
        list(proteinChain = list(
          sequence = db[[full_HC_aa]][i],
          count = 1,
          useStructureTemplate = useStructureTemplate_on["HC"]
        )),
        list(proteinChain = list(
          sequence = db[[full_LC_aa]][i],
          count = 1,
          useStructureTemplate = useStructureTemplate_on["LC"]
        ))
      ),
      dialect = "alphafoldserver",
      version = 1
    )
  })
  
  json_output <- jsonlite::toJSON(json.files, pretty = TRUE, auto_unbox = TRUE)
  
  write(json_output, paste0(out_folder, "/", json_name, ".json"))
  
  if(return_db){
    db <- db %>%
      dplyr::mutate(
        main_folder = out_folder,
        !!rlang::sym("AF3_", antigen_name,"_run_id") := tolower(request_name)
      ) %>%
      dplyr::select(-request_name)
  }
}

#### Function to extract a list of Ag-mAb results from AlphaFold3 ####
#' extract a matrix of interaction of a list of mAbs (HC+LC) onto a given antigen
#'
#' \code{extractAF3contacts} extract AlphaFold3 predicted interactions from multiple HC/LC/Ag models
#' @param db a dataframe with the following column: "run_id" and "main_folder" and "run_folder"
#' @param main_folder name of the column containing the path to the AlphaFold3 result folder (not including the AF3 folder itself)
#' @param run_folder name of the column containing the name of the AlphaFold3 folder, by default same as run_id.
#' @param QC-cutoff which cutoff to use to filter models based on local ipTM scores (if set to NULL, no filtering is done)
#' @param ... extra parameters to pass on to extractSingleAF3contacts, notably prefix [default = "fold_"]; chains [default = c(A = "Ag", B = "HC", C = "LC")]; Ag_start_num [default = 1] and AF3_model [default = 0]
#' @return
#' If which = "all", a list containing a data frame with QC_scores for each run_id and three matrices of AF3 predicted contacts(HC, LC and all)
#' If which = "contacts", only a list of contact matrices is returned.
#' If which = "QC_scores", only the QC scores dat frame is returned.
#' 
#' @export

extractAF3contacts <- function(db,
                               which = c("all", "contacts", "QC_scores"),
                               run_id = "run_id",
                               main_folder = "main_folder",
                               run_folder = "run_id",
                               type_of_run = "alphafoldserver",
                               QC_cutoff = 0.6,
                               verbose = TRUE,
                               ...){
  
  which <- match.arg(which)
  
  missing_collumns <- setdiff(unique(c(run_id, run_folder, main_folder)), colnames(db))
  if(length(missing_collumns)>0){
    stop("the following columns could not be found in the provided dataframe: ", paste(missing_collumns, collapse = "; "))
  }
  db <- db %>% 
    dplyr::mutate(
      main_folder = .data[[main_folder]],
      run_folder = .data[[run_folder]],
      run_id = .data[[run_id]]
      )
  
  # validate type_of_run
  valid_types <- c("local", "alphafoldserver")
  if (is.character(type_of_run) && length(type_of_run) == 1 && type_of_run %in% valid_types) {
    db$type_of_run <- type_of_run
    
  } else if (is.character(type_of_run) && length(type_of_run) == 1 && type_of_run %in% names(db)) {
    db <- db %>% dplyr::mutate(type_of_run = .data[[type_of_run]])
  } else {
    stop("`type_of_run` must be either 'local', 'server', or a column name in `db`.")
  }
  
  predicted_contacts.list <- purrr::pmap(
    c(
      list(
        main_folder = db$main_folder,
        run_folder = db$run_folder,
        type_of_run = db$type_of_run
      ),
      list(...)
    ),
    extractSingleAF3contacts
  )
  names(predicted_contacts.list) <- db[[run_id]]
  
  #deals with missing file cases (NULL is returned)
  predicted_contacts.list <- predicted_contacts.list %>%
    purrr::discard(is.null)
  
  #concatenate individual results
  QC_scores <- as.data.frame(do.call(rbind, lapply(predicted_contacts.list, FUN = function(prediction){
    return(prediction[["QC_scores"]])
  })))
  QC_scores <- rownames_to_column(QC_scores, var = "run_id") 
    
  predicted_HC_contacts <- t(do.call(cbind, lapply(predicted_contacts.list, FUN = function(prediction){
    return(prediction[["predicted_HC_contacts"]])
  })))
  
  predicted_LC_contacts <- t(do.call(cbind, lapply(predicted_contacts.list, FUN = function(prediction){
    return(prediction[["predicted_LC_contacts"]])
  })))
  
  predicted_all_contacts <- pmax(predicted_HC_contacts, predicted_LC_contacts)
  
  #filter results
  if(!is.null(QC_cutoff)){
    run_to_keep <- dplyr::filter(QC_scores, hc_local_ipTMs >= QC_cutoff & lc_local_ipTMs >= QC_cutoff)$run_id
    run_to_discard <- dplyr::filter(QC_scores, hc_local_ipTMs < QC_cutoff | lc_local_ipTMs < QC_cutoff)$run_id
    if(length(run_to_discard)>0 & verbose){
      warning(paste0("the following predictions were of poor quality: ", paste(run_to_discard, collapse = "; ")))
    }
    predicted_HC_contacts <- predicted_HC_contacts[rownames(predicted_HC_contacts) %in% run_to_keep, ]
    predicted_LC_contacts <- predicted_LC_contacts[rownames(predicted_LC_contacts) %in% run_to_keep, ]
    predicted_all_contacts <- predicted_all_contacts[rownames(predicted_all_contacts) %in% run_to_keep, ]
  }
  
  if(run_id != run_folder){
    QC_scores <- QC_scores %>% 
      dplyr::left_join(dplyr::select(db, all_of(c("run_id", "run_folder"))), by = join_by(run_id))
  }
  
  #return results
  if(which == "all"){
    return(list(
      QC_scores = QC_scores,
      predicted_all_contacts = predicted_all_contacts,
      predicted_HC_contacts = predicted_HC_contacts,
      predicted_LC_contacts = predicted_LC_contacts
    ))
  }
  if(which == "contacts"){
    return(list(
      predicted_all_contacts = predicted_all_contacts,
      predicted_HC_contacts = predicted_HC_contacts,
      predicted_LC_contacts = predicted_LC_contacts
    ))
  }
  if(which == "QC_scores"){
    return(QC_scores)
  }
}

#### Function to extract a single Ag-mAb result from AlphaFold3 ####
#' extract interactions from one mAb onto a given antigen
#'
#' \code{extractSingleAF3contacts} extract AlphaFold3 predicted interactions from a single HC/LC/Ag model
#' @param main_folder path to the AlphaFold3 result folder (not including the AF3 folder itself)
#' @param run_folder name of the AlphaFold3 folder
#' @param prefix whether to add a prefix or not to the run folder name (by default AlphaFold server adds a "fold_" prefix to all files but only to the folder itself if results are individually downloaded)
#' @param chains a named vector giving the order of chains used for AlphaFold3 modeling, expected names are "HC", "LC", "Ag" in any given order and not needing to start at A if a fourth protein was used.
#' @param Ag_start_num option to start numbering of AA in Ag chain somewhere else than 1 if part of the protein was truncated
#' @param Af3_model which of the 5 outputted AlphaFold3 models to use, by default model 0 is the best one.
#' @param type_of_run which implementation of AlphaFold3 was used: one of "local" or "alphafoldserver". Naming and organisation of folders and files changes between both.  
#' @return
#' a list with QC scores, predicted HC contacts and predicted LC contacts. 
#'
#' @keywords internal

extractSingleAF3contacts <- function(main_folder,
                                     run_folder,
                                     prefix = "fold_",
                                     chains = c(A = "Ag", B = "HC", C = "LC"),
                                     Ag_start_num = 1,
                                     AF3_model = NULL,
                                     type_of_run = c("alphafoldserver", "local")
                                     ){
                                      
  #check if folder exists:
  full_run_folder <- paste0(main_folder, "/", run_folder)
  if(!dir.exists(full_run_folder)){
    warning(paste0("missing folder for ", run_folder))
    return(NULL)
  }
  
  #define model to use based on the type of run (alphafold server vs local install:
  type_of_run <- match.arg(type_of_run)
  if(!type_of_run %in% c("alphafoldserver", "local")){
    stop("type_of_run should be one of 'alphafoldserver' or 'local'")
  }
    
  if(type_of_run == "alphafoldserver"){
    if(is.null(AF3_model)){
      #using by default model zero
      AF3_model <- "_0"
      #TODO scan summary confidence json and select the best ranking score?
      } else { AF3_model <- paste0("_", AF3_model) }
  } else {
    if(!is.null(AF3_model)){
      #TODO open "run_folder_ranking_scores.csv", get seed (column 1) and redefine folder to use using seed-"795906"seed nb"_sample-"AF3_model" 
    }
  }
  
  #define file name and check if they exist:
  if(type_of_run == "alphafoldserver") {
    json_confidence_filename <- paste0(full_run_folder, "/", prefix, run_folder, "_summary_confidences", AF3_model,".json")
    data_filename <- paste0(full_run_folder, "/", prefix, run_folder, "_full_data", AF3_model,".json")
    request_filename <- paste0(full_run_folder,"/", prefix, run_folder, "_job_request.json")
  } else if(type_of_run == "local") {
    json_confidence_filename <- paste0(full_run_folder, "/", run_folder, "_summary_confidences", AF3_model,".json")
    data_filename <- paste0(full_run_folder, "/", run_folder, "_confidences", AF3_model,".json")
    request_filename <- paste0(full_run_folder,"/", run_folder, "_data", AF3_model,".json")
  }
  
  json_confidence <- jsonlite::fromJSON(txt = json_confidence_filename, flatten = TRUE)
  json_data <- jsonlite::fromJSON(txt = data_filename, flatten = TRUE)
  
  #extract ipTMs  
  global_ipTM <- json_confidence$iptm
  local_ipTMs <- json_confidence$chain_pair_iptm
  
  if(!ncol(local_ipTMs) == length(chains)){
    stop("issue with : ", full_run_folder, ", ",length(colnames(local_ipTMs)), " chains in request file instead of 3 (Ag, HC and LC)")
  }
  colnames(local_ipTMs) <- paste0(tolower(chains), "_local_ipTMs")
  rownames(local_ipTMs) <- paste0(tolower(chains), "_local_ipTMs")
  
  ipTMs <- c(global_ipTM = global_ipTM, local_ipTMs["ag_local_ipTMs", c("hc_local_ipTMs", "lc_local_ipTMs")])
    
  #extract contacts:
  contact_mat <- as.data.frame(json_data$contact_probs)
    
  colnames(contact_mat) <- paste0(json_data$token_chain_ids, json_data$token_res_ids)
  rownames(contact_mat) <- paste0(json_data$token_chain_ids, json_data$token_res_ids)
    
  chains_aa <- list()
  if(type_of_run == "alphafoldserver") {
    request <- jsonlite::fromJSON(txt = request_filename, flatten = FALSE)
    for (i in seq_along(chains)){
      chains_aa[[i]] <- request$sequences[[1]]$proteinChain[i,1]
    }
    names(chains_aa) <- chains
  } else if(type_of_run == "local") {
    request <- jsonlite::fromJSON(txt = request_filename, flatten = TRUE)
    for (i in seq_along(chains)){
      chains_aa[[i]] <- request$sequences$protein.sequence[i]
    }
    names(chains_aa) <- chains
  }
    
  Ag_aa_seq <- unlist(strsplit(chains_aa[["Ag"]], split = ""))
  Ag_numbered_aa_seq <- paste0(Ag_aa_seq, (seq_along(Ag_aa_seq)+Ag_start_num-1))
  
  HCvsAg_contacts <- contact_mat[grepl(names(chains[chains == "Ag"]), colnames(contact_mat)), grepl(names(chains[chains == "HC"]), rownames(contact_mat))]
  predicted_HC_contacts <- apply(HCvsAg_contacts, 1, max)
  names(predicted_HC_contacts) <- Ag_numbered_aa_seq
  
  LCvsAg_contacts <- contact_mat[grepl(names(chains[chains == "Ag"]), colnames(contact_mat)), grepl(names(chains[chains == "LC"]), rownames(contact_mat))]
  predicted_LC_contacts <- apply(LCvsAg_contacts, 1, max)
  names(predicted_LC_contacts) <- Ag_numbered_aa_seq
  
  return(list(QC_scores = ipTMs, predicted_HC_contacts = predicted_HC_contacts, predicted_LC_contacts = predicted_LC_contacts))
}

#### Function to plot an enriched Heatmap for AlphaFold3 predicted contacts on a given antigen ####
#' plot a heatmap with additional annotations on side top and /or bottom
#'
#' \code{AF3contactsHeatmap} plot AlphaFold3 predicted interactions 
#' @param contact_mat a contact matrix in the format provided by extractsAF3contacts 
#' @param cluster_rows 	If the value is a logical, it controls whether to make cluster on rows. The value can also be a hclust or a dendrogram which already contains clustering. Check https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#clustering .
#' @param show_row_dend whether to plot dendrogram
#' @param row_split A vector or a column in row_annot data frame by which the rows are split if cluster_rows is a logical. But if cluster_rows is a clustering object, split can only be a single number indicating to split the dendrogram by cutree.
#' @param show_column_dend whether to display column dendrogram too (passed to ComplexHeatmap)
#' @param cluster_columns whether to cluster columns too (passed to ComplexHeatmap)
#' @param column_split A vector by which the rows are split if cluster_rows is a logical. But if cluster_rows is a clustering object, split can only be a single number indicating to split the dendrogram by cutree.
#' @param border whether to add border to the heatmap
#' @param heatmap_gap size of gap to be introduced (default = unit(1, "mm")
#' @param hclust_method method for hclust if dendrogram is to be internally calculated (default = "ward.D2")
#' @param col color to use for the heatmap (default = circlize::colorRamp2(c(0, 0.25, 0.5, 1), c("white", "cornflowerblue", "yellow", "red")))
#' @param legend_name name to use for the heatmap color legend (default = "AlphaFold contact probability")
#' @param show_row_names whether to show rownames (of note highlight_row_names will supersede that call)
#' @param row_names_side which side of the heatmap to display rownames (left or right (default))
#' @param show_column_names whether to show colnames
#' @param column_names_side which side of the heatmap to display colnames (top or bottom (default))
#' @param row_names_gp gpar for rownames 
#' @param column_names_gp gpar for colnames
#' @param row_title title for rownames
#' @param title title for the heatmap
#' @param title_gp gpar for heatmap title
#' @param height global height of the heatmap
#' @param row_height individual height of row (setting a value for height will supersede row_height)
#' @param width global width of the heatmap
#' @param col_width individual width of col (setting a value for width will supersede row_height)
#' @param average_contacts whether to plot average contacts for each individual clusters, should point to a column in row_annot. 
#' @param annot_col a list of all colors scheme related to any annotation on the enriched heatmap
#' @param row_annot all row annotations
#' @param row_annot_name_side where to print row annotation names (top (default) or bottom)
#' @param column_annot_name_side where to print column annotation names (right (default) or left)
#' @param row_annot_borders whether to add a border to rowannotations
#' @param highlight_row_names which rownames to highlight (only these will be printed on the right of the heatmap)
#' @param highlight_labels_gp gpar for highlighted row names
#' @param highlight_padding padding for highlighted row names
#' @param top_annot all annotations to be added at the top of the heatmap, order of provide annotation will define the order of plotting from top to bottom.
#' @param bottom_annot all annotations to be added at the bottom of the heatmap, order of provide annotation will define the order of plotting from top to bottom.
#' @param annot_type a list of type of annotation related to all annotation added at the top or bottom of the heatmap, should be one of 'structure' (block will be create) or 'sites' (individual sites will be colored)
#' @param column_annot_gap gap between two top or bottom annotations
#' @param barplot_annot_height height for individual barplot annotations (related to plot_average_contacts)
#' @param site_annot_height height for binding sites annotations
#' @param structure_annot_height a named list of height for structure annotations, if not provided for a given structure, will default internally to site_annot_height for that structure.
#' @param annotation_name_gp gpar for highlighted top or bottom annotation names
#' @param save_pdf whether to save a pdf file
#' @param filename name of saved pdf file (without the ".pdf")
#' @param file_width width of the saved pdf 
#' @param file_height height of the saved pdf 
#' @param plot_row_dendrogram whether to also save the dendrogram
#' @param return_plot whether to return the ComplexHeatmap grop
#' @param ... additional arguments to pass to ComplexHeatmap::Heatmap
#' 
#' @return
#' If return_plot = TRUE, a complexHeatmap grob.
#' If save_pdf = TRUE will automatically save a pdf with the enriched Heatmap (page 1) as well as the dendrogram for all rows (page 2) if plot_row_dendrogram = TRUE.
#' 
#' @export

contactsHeatmap <- function(contact_mat,
                            cluster_rows = FALSE,
                            show_row_dend = FALSE,
                            row_split = NULL,
                            cluster_columns = FALSE,
                            show_column_dend = FALSE,
                            column_split = NULL,
                            border = TRUE,
                            heatmap_gap = grid::unit(1, "mm"),
                            hclust_method = "ward.D2",
                            col = circlize::colorRamp2(c(0, 0.25, 0.5, 1), c("white", "cornflowerblue", "yellow", "red")),
                            legend_name = "AlphaFold contact probability",
                            show_row_names = TRUE,
                            row_names_side = "right",
                            show_column_names = TRUE,
                            column_names_side = "bottom",
                            row_names_gp = grid::gpar(fontsize = 6),
                            column_names_gp = grid::gpar(fontsize = 6),
                            row_title = NULL,
                            title = "VH/VL_predicted_residues",
                            title_gp = grid::gpar(fontsize = 20, fontface = "bold"),
                            height = NULL,
                            row_height = grid::unit(2, "mm"),
                            width = NULL,
                            col_width = grid::unit(2, "mm"),
                            average_contacts = NULL,
                            annot_col = list(),
                            row_annot = list(),
                            row_annot_name_side = "top",
                            column_annot_name_side = "right",
                            row_annot_borders = FALSE,
                            highlight_row_names = NULL,
                            highlight_labels_gp = grid::gpar(fontsize = 9), 
                            highlight_padding = grid::unit(1, "mm"),
                            top_annot = list(),
                            bottom_annot = list(),
                            annot_type = list(),
                            structure_type = list(),
                            print_structure_name = TRUE,
                            column_annot_gap = grid::unit(2, "mm"),
                            barplot_annot_height = grid::unit(4 , "mm"),
                            site_annot_height = grid::unit(4, "mm"),
                            structure_annot_height = list(),
                            annotation_name_gp = grid::gpar(fontsize = 9),
                            save_pdf = TRUE,
                            filename = "all_contacts_heatmap",
                            file_width = 20, 
                            file_height = 20,
                            plot_row_dendrogram = TRUE,
                            plot_column_dendrogram = TRUE,
                            return_plot = TRUE,
                            ...){
                               
  if(is.null(contact_mat)){
    stop("no contact matrix provided")
  }
  
  #define row annotation:
  if(length(row_annot)>0){
    common_rownames <- rownames(contact_mat)[rownames(contact_mat) %in% rownames(row_annot)]
    missing_rownames <- rownames(contact_mat)[!rownames(contact_mat) %in% rownames(row_annot)]
    row_annot <- row_annot[common_rownames, ]
    if(length(missing_rownames)>0){
      warning("Provided row annotation dataframe is missing information regarding the following rownames: ", paste(missing_rownames, collapse = "; "))
      missing_rows <- as.data.frame(matrix(NA, nrow = length(missing_rownames), ncol = ncol(row_annot),
                                           dimnames = list(missing_rownames, colnames(row_annot))))
      row_annot <- row_annot %>%
        dplyr::bind_rows(missing_rows)
    }
    if(row_annot_borders){
      row_anno_gp <- grid::gpar(col = "black")
    } else {row_anno_gp <- NULL}
    row_ha = ComplexHeatmap::rowAnnotation(df = row_annot[rownames(contact_mat),], 
                                           col = annot_col[names(annot_col) %in% colnames(row_annot)], 
                                           gp = row_anno_gp, 
                                           border = border,
                                           annotation_name_side = row_annot_name_side)
  } else {
    row_ha <- NULL
    if(!is.null(average_contacts)){
      warning("The 'average_contacts' argument must point to a column in 'row_annot'. No row annotations dataframe provided, 'average_contacts' will be set to NULL.")
      average_contacts <- NULL
    }
  }
  
  #define row clustering and splits (dendrogram or annotation-based)
  if(isFALSE(cluster_rows)){
    show_row_dend = FALSE
    plot_row_dendrogram = FALSE
    if(row_split %in% colnames(row_annot)){
      row_split = row_annot[[row_split]]
      k_rows = length(levels(as.factor(row_split)))
    } else if(is.vector(row_split)){
      k_rows = length(levels(as.factor(row_split)))
    } else {
      if(!is.null(row_split)){
        warning("row_split should be a vector or a column in row_annot when setting cluster_rows = FALSE")
        row_split = NULL
      }
      k_rows = NULL
    }
  } else if(isTRUE(cluster_rows)){
    #define dendrogram and use it to cluster rows and also rename it as dend to be able to plot it in the recap pdf
    h_rows <- hclust(dist(contact_mat), method = hclust_method)
    cluster_rows <-  as.dendrogram(h_rows)
    if(!is.numeric(row_split)){
      if(!is.null(row_split)){
        warning("cluster_rows should be a numeric value if a dendrogram is used to cluster rows")
        row_split = NULL
      }
    }
    k_rows = row_split
    
  } else if(class(cluster_rows) %in% c("dendrogram", "hclust")){
    if(class(cluster_rows) == "hclust"){
      cluster_rows <- as.dendrogram(cluster_rows)
    } 
    if(!is.numeric(row_split)){
      if(!is.null(row_split)){
        warning("cluster_rows should be a numeric value if a dendrogram is used to cluster rows")
        row_split = NULL
      }
    }
    k_rows = row_split
  } else {
    warning("cluster_rows should be one of logical value or dendrogram/hclust object: set to FALSE")
    cluster_rows = FALSE
    show_row_dend = FALSE
    plot_row_dendrogram = FALSE
    
    if(row_split %in% colnames(row_annot)){
      row_split = row_annot[[row_split]]
      k_rows = length(levels(as.factor(row_split)))
    } else if(is.vector(row_split)){
      k_rows = length(levels(as.factor(row_split)))
    } else {
      if(!is.null(row_split)){
        warning("row_split should be a vector or a column in row_annot when setting cluster_rows = FALSE")
        row_split = NULL
      }
      k_rows = NULL
    }
  }
  
  #define columns clustering and splits (dendrogram or vector-based)
  if(isFALSE(cluster_columns)){
    show_column_dend = FALSE
    plot_column_dendrogram = FALSE
    if(is.vector(column_split)){
      k_columns = length(levels(as.factor(column_split)))
    } else {
      if(!is.null(column_split)){
        warning("row_split should be a vector or a column in row_annot when setting cluster_rows = FALSE")
        column_split = NULL
      }
      k_columns = NULL
    }
  } else if(isTRUE(cluster_columns)){
    #define dendrogram and use it to cluster rows and also rename it as dend to be able to plot it in the recap pdf
    h_columns <- hclust(dist(contact_mat), method = hclust_method)
    cluster_rows <-  as.dendrogram(h_columns)
    if(!is.numeric(column_split)){
      if(!is.null(column_split)){
        warning("cluster_rows should be a numeric value if a dendrogram is used to cluster rows")
        column_split = NULL
      }
    }
    k_columns = column_split
    
  } else if(class(cluster_columns) %in% c("dendrogram", "hclust")){
    if(class(cluster_columns) == "hclust"){
      cluster_columns <- as.dendrogram(cluster_columns)
    } 
    if(!is.numeric(column_split)){
      if(!is.null(column_split)){
        warning("cluster_rows should be a numeric value if a dendrogram is used to cluster rows")
        column_split = NULL
      }
    }
    k_columns = column_split
  } else {
    warning("cluster_rows should be one of logical value or dendrogram/hclust object: set to FALSE")
    cluster_columns = FALSE
    show_column_dend = FALSE
    plot_column_dendrogram = FALSE
    
    if(is.vector(column_split)){
      k_columns = length(levels(as.factor(column_split)))
    } else {
      if(!is.null(column_split)){
        warning("row_split should be a vector or a column in row_annot when setting cluster_rows = FALSE")
        column_split = NULL
      }
      k_columns = NULL
    }
  }
  
  #define bottom annotation:
  if(!is.null(average_contacts)){
    if(average_contacts %in% colnames(row_annot)){
      classes <- levels(droplevels(as.factor(row_annot[[average_contacts]])))
      base_mean.list <- lapply(classes, FUN = function(class){
        rownames_in_class <- rownames(row_annot[!is.na(row_annot[[average_contacts]]) & row_annot[[average_contacts]] == class, ])
        rownames_in_class <- rownames_in_class[rownames_in_class %in% rownames(contact_mat)]
        if(length(rownames_in_class)>1){
          colMeans(as.matrix(contact_mat[rownames_in_class,]))
        } else {
          #cases of only one VH/VL pair in a class
          contact_mat[rownames_in_class,]
        }
      })
      names(base_mean.list) <- classes
      base_col <- annot_col[[average_contacts]]
      missing_base_col <- setdiff(classes, names(base_col))
      if (length(missing_base_col) > 0) {
        warning("Missing colors for the following ", average_contacts, " levels: ", paste(missing_base_col, collapse = "; "),"; assigning automatically.")
        add_base_col <- RColorBrewer::brewer.pal(n = 12, "Paired")[1:length(missing_base_col)]
        names(add_base_col) <- missing_base_col
        base_col <- c(base_col, add_base_col)
        annot_col[[average_contacts]] <- base_col
      }
      base_fill <- base_col
    } else {
      warning("missing ", average_contacts," in row_annot dataframe.")
      if(!is.null(nrow(contact_mat))){
        base_mean.list <- list("Average_binding" = colMeans(contact_mat))
      } else {
        base_mean.list <- list("Average_binding" = contact_mat)
      }
      base_col <- c("Average_binding" = "black")
      base_fill <- c("Average_binding" = "black")
    }
  } else {
    if(!is.null(nrow(contact_mat))){
      base_mean.list <- list("Average_binding" = colMeans(contact_mat))
    } else {
      base_mean.list <- list("Average_binding" = contact_mat)
    }
    base_col <- c("Average_binding" = "black")
    base_fill <- c("Average_binding" = "black")
  }
   
  #define size of the final plot based on row_height, col_width and number of gaps (k_rows, k_cols)
  if(is.null(height)){
    if(!is.null(nrow(contact_mat))){
      if(is.null(k_rows)){
        height <- nrow(contact_mat) * row_height
      } else {
        height <- (nrow(contact_mat) * row_height) + ((k_rows-1)* heatmap_gap)
      }
    } else {
      #case of only one row
      height <- row_height
    }
  }
  if(is.null(width)){
    if(!is.null(ncol(contact_mat))){
      width <- ncol(contact_mat) * col_width
    } else {
      #case of only one row
      width <- length(contact_mat) * col_width
    }
  }
  
  if(!is.null(average_contacts)){
    # Build barplot annotations iteratively
    barplot_annos <- lapply(seq_along(base_mean.list), function(i) {
      name <- names(base_mean.list)[i]
      ComplexHeatmap::anno_barplot(
        base_mean.list[[i]],
        height = barplot_annot_height,
        gp = grid::gpar(col = base_col[[name]], fill = base_fill[[name]]),
        border = FALSE,
        axis_param = list(at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), gp = grid::gpar(fontsize = 6))
      )
    })
    names(barplot_annos) <- names(base_mean.list)
  } else {
    barplot_annos <- NULL
  }
  
  #implement whether or not to print structure names:
  if(isTRUE(print_structure_name)){
    structures_names <- names(annot_type[annot_type == "structure"])
    if(length(structures_names)>0){
      print_structure_name <- rep(TRUE, length(structures_names))
      names(print_structure_name) <- structures_names
    }
  } else {
    structures_names <- names(annot_type[annot_type == "structure"])
    if(isFALSE(print_structure_name)){
      if(length(structures_names)>0){
        print_structure_name <- rep(FALSE, length(structures_names))
        names(print_structure_name) <- structures_names
      }
    } else {
      missing_st_name <- structures_names[!structures_names %in% names(print_structure_name)]
      if(length(missing_st_name)>0){
        add_structure_name <- rep(FALSE, length(missing_st_name))
        names(add_structure_name) <- missing_st_name
        print_structure_name <- c(print_structure_name, add_structure_name)
      }
    }
  }
  
  #create top and bottom annotations:
  bottom_annos <- generate_tb_annotations(annots = bottom_annot,
                                         annot_type = annot_type,
                                         annot_col = annot_col,
                                         structure_type = structure_type,
                                         structure_annot_height = structure_annot_height,
                                         site_annot_height = site_annot_height,
                                         annotation_name_gp = annotation_name_gp,
                                         print_structure_name = print_structure_name,
                                         column_annot_gap = column_annot_gap,
                                         column_annot_name_side = column_annot_name_side,
                                         barplot_annos = barplot_annos)
  
  top_annos <- generate_tb_annotations(annots = top_annot,
                                      annot_type = annot_type,
                                      annot_col = annot_col,
                                      structure_type = structure_type,
                                      structure_annot_height = structure_annot_height,
                                      site_annot_height = site_annot_height,
                                      annotation_name_gp = annotation_name_gp,
                                      print_structure_name = print_structure_name,
                                      column_annot_gap = column_annot_gap,
                                      column_annot_name_side = column_annot_name_side,
                                      barplot_annos = NULL)

  #plot heatmap:
  HM <- ComplexHeatmap::Heatmap(contact_mat,
                                name = legend_name, 
                                cluster_rows = cluster_rows,
                                show_row_dend = show_row_dend,
                                row_split = row_split, 
                                gap = heatmap_gap,
                                border = border,
                                cluster_columns = cluster_columns,
                                show_column_dend = show_column_dend,
                                column_split = column_split,
                                col = col,
                                left_annotation = row_ha, 
                                top_annotation = top_annos$ha,
                                bottom_annotation = bottom_annos$ha,
                                show_row_names = show_row_names,
                                row_names_side = row_names_side,
                                show_column_names = show_column_names,
                                column_names_side = column_names_side,
                                row_names_gp = row_names_gp,
                                column_names_gp = column_names_gp,
                                row_title = row_title,
                                column_title = title,
                                column_title_gp = title_gp,
                                height = height,
                                width = width,
                                ...)
                
  
  if(!is.null(highlight_row_names)){
    HM <- HM +
      ComplexHeatmap::rowAnnotation(link = ComplexHeatmap::anno_mark(at = which(rownames(contact_mat) %in% highlight_row_names), 
                                     labels = rownames(contact_mat)[rownames(contact_mat) %in% highlight_row_names], 
                                     labels_gp = highlight_labels_gp, 
                                     padding = highlight_padding))
  }
  
  if(save_pdf){
    pdf(paste0(filename, ".pdf"), width = grid::unit(file_width, "cm"), height = grid::unit(file_height, "cm"))
    plot(HM)
    if(plot_row_dendrogram){
      plot(dendextend::circlize_dendrogram(dend))
    }
    dev.off()
  }
  
  if(return_plot){
    return(list(
      heatmap = HM,
      annot_col = list(
        row = annot_col[names(annot_col) %in% colnames(row_annot)],
        top = top_annos$col,
        bottom = bottom_annos$col
      ),
      extra_legends = list(top_structures = top_annos$extra_legends,
                           bottom_structures = bottom_annos$extra_legends)))
  }
}

#### Function to plot an enriched Heatmap for AlphaFold3 predicted contacts on a given antigen ####
#' plot a heatmap with additional annotations on side top and /or bottom
#'
#' \code{AF3contactsHeatmap} plot AlphaFold3 predicted interactions 
#' @param annots all annotations to be added, order of provide annotation will define the order of plotting from top to bottom.
#' @param annot_type a list of type of annotation related to all annotation added at the top or bottom of the heatmap, should be one of 'structure' (block will be create) or 'sites' (individual sites will be colored)
#' @param structure_type a list of type of structures related to all annotation added at the top or bottom of the heatmap,
#' @param annot_col a list of all colors scheme related to any annotation on the enriched heatmap
#' @param column_annot_gap gap between two top or bottom annotations
#' @param barplot_annot_height height for individual barplot annotations (related to plot_average_contacts)
#' @param site_annot_height height for binding sites annotations
#' @param structure_annot_height height for structure annotations
#' @param annotation_name_gp gpar for highlighted top or bottom annotation names
#' @return
#' a complete heatmap annotation panel for top or bottom annotation of contacts heatmap
#' 
#' @keywords internal

generate_tb_annotations <- function(annots,
                                    annot_type,
                                    annot_col = list(),
                                    structure_type = list(),
                                    structure_annot_height = list(),
                                    site_annot_height = grid::unit(5, "mm"),
                                    annotation_name_gp = NULL,
                                    print_structure_name = NULL,
                                    column_annot_gap = grid::unit(2, "mm"),
                                    column_annot_name_side = "right",
                                    barplot_annos = list())
                                    {
                                    
  # build annotation + color info
  results <- purrr::imap(annots, function(annot, name) {
    type <- annot_type[[name]]
    if (is.null(type)) {
      warning("No type info for the following annotation: ", name)
      return(NULL)
    }
    
    if (type == "sites") {
      col_map <- annot_col[[name]]
      if (!is.null(col_map)) {
        missing_sites_col <- setdiff(unique(annot), names(col_map))
        if (length(missing_sites_col) > 0) {
          warning("Missing colors for ", name, " sites; assigning automatically.")
          add_col <- RColorBrewer::brewer.pal(n = 12, "Paired")[1:length(missing_sites_col)]
          names(add_col) <- missing_sites_col
          col_map <- c(col_map, add_col)
        }
      }
      
      list(
        anno = annot,
        col = col_map,
        type = "sites"
      )
      
    } else if (type == "structure") {
      annotations_height <- if (!is.null(structure_annot_height) &&
                                is.list(structure_annot_height) &&
                                name %in% names(structure_annot_height)) {
        structure_annot_height[[name]]
      } else {
        site_annot_height
      }
      
      if (!grid::is.unit(annotations_height)) {
        # if numeric, treat as mm; adjust as you prefer
        if (is.numeric(annotations_height) && length(annotations_height) == 1 && !is.na(annotations_height)) {
          annotations_height <- grid::unit(annotations_height, "mm")
        } else {
          annotations_height <- site_annot_height
        }
      }
      
      gp <- if (!is.null(annotation_name_gp)) {
        modifyList(grid::gpar(col = "black"), annotation_name_gp)
      } else {
        grid::gpar(col = "black")
      }
      
      structure_fill_col <- annot_col[[name]]
      
      if (is.character(structure_fill_col) && length(structure_fill_col) == 1 && is.null(names(structure_fill_col))) {
        structure_fill_col <- c("Helix" = structure_fill_col, "Strand" = structure_fill_col)
      }
      
      if (is.null(structure_fill_col)) {
        structure_fill_col <- c("Helix" = "moccasin", "Strand" = "lightblue")
      }
      
      if(name %in% names(structure_type)){
        struct_types <- unique(unlist(structure_type[[name]]))
        missing_struct_col <- setdiff(struct_types, names(structure_fill_col))
        if (length(missing_struct_col) > 0) {
          warning("Missing colors for ", name, " structure; assigning automatically.")
          add_struct_col <- RColorBrewer::brewer.pal(n = 9, "Pastel1")[1:length(missing_struct_col)]
          names(add_struct_col) <- missing_struct_col
          structure_fill_col <- c(structure_fill_col, add_struct_col)
        }
        
        anno_obj <- ComplexHeatmap::anno_block(
          align_to = annot,
          labels = names(annot),
          panel_fun = function(index, level) {
            struct_type <- structure_type[[name]][[level]]
            fill_col <- as.character(structure_fill_col[struct_type])
            if (is.null(fill_col) || length(fill_col) == 0 || is.na(fill_col)) {
              fill_col <- "moccasin"
            }
            grid::grid.rect(gp = grid::gpar(fill = fill_col, col = "black"))
            if (isTRUE(print_structure_name[name])) {
              grid::grid.text(level, 0.5, 0.5, rot = 0, gp = gp)
            }
          },
          height = annotations_height
        )
      } else {
        anno_obj <- ComplexHeatmap::anno_block(
          align_to = annot,
          labels = names(annot),
          panel_fun = function(index, level) {
            fill_col <- as.character(structure_fill_col[1])
            if (is.null(fill_col) || length(fill_col) == 0 || is.na(fill_col)) {
              fill_col <- "moccasin"
            }
            grid.rect(gp = gpar(fill = fill_col, col = "black"))
            if (isTRUE(print_structure_name[name])) {
              grid::grid.text(level, 0.5, 0.5, rot = 0, gp = gp)
            }
          },
          height = annotations_height
        )
      }
      
      list(
        anno = anno_obj,
        col = structure_fill_col,
        type = "structure"
      )
      
    } else {
      warning("Unsupported annotation type for: ", name)
      return(NULL)
    }
  }) %>% purrr::compact()
  
  annos <- purrr::map(results, "anno")
  annos_col   <- purrr::map(results, "col")
  annos_types <- purrr::map_chr(results, "type")
  
  annos <- c(annos, barplot_annos) %>% purrr::compact()
  
  if (is.null(column_annot_gap) || !grid::is.unit(column_annot_gap)) {
    column_annot_gap <- grid::unit(2, "mm")
  }
  if (is.null(site_annot_height) || !grid::is.unit(site_annot_height)) {
    site_annot_height <- grid::unit(5, "mm")
  }
  
  # Combine into HeatmapAnnotation
  if (length(annos) > 0) {
    full_ha <- do.call(
      ComplexHeatmap::HeatmapAnnotation,
      c(
        annos,
        list(
          col = annos_col,
          gap = column_annot_gap,
          annotation_name_side = column_annot_name_side,
          simple_anno_size = site_annot_height,
          annotation_name_gp = annotation_name_gp,
          annotation_label = names(annos)
        )
      )
    )
  } else {
    full_ha <- NULL
  }
  
  structure_legends <- purrr::imap(
    results,
    function(res, name) {
      if (res$type == "structure" && !is.null(res$col)) {
        ComplexHeatmap::Legend(
          labels = names(res$col),
          legend_gp = grid::gpar(fill = unname(res$col)),
          title = paste(name, "Structure Type")
        )
      } else {
        NULL
      }
    }
  ) %>% purrr::compact()
  
  list(
    ha = full_ha,
    col = annos_col,
    type = annos_types,
    extra_legends = structure_legends
  )
}


#### Function to draw an enriched Heatmap from the results of contactsHeatmap ####
#' Draw a heatmap object with optional extra legends (no duplicates)
#' 
#' \code{heatmap_draw} plot a heatmap with additional legends 
#' @param ht_obj a list with $heatmap (a Heatmap or HeatmapList)
#'        and optionally $extra_legends (a list of Legend objects)
#' @return Invisibly returns what `ComplexHeatmap::draw()` returns.

heatmap_draw <- function(ht_obj,
                         heatmap_legend_side = "bottom",
                         annotation_legend_side = "bottom",
                         heatmap_legend_direction = "horizontal",
                         annotation_legend_direction = "horizontal",
                         padding = grid::unit(c(5, 20, 5, 5), "mm"),
                         ...) {
  
  ht <- ht_obj$heatmap
  if (is.list(ht) && length(ht) == 1 && inherits(ht[[1]], "Heatmap")) {
    ht <- ht[[1]]
  }
  
  if (!inherits(ht, c("Heatmap", "HeatmapList"))) {
    stop("`ht_obj$heatmap` must be a ComplexHeatmap object (Heatmap or HeatmapList), not a list.")
  }
  
  legends <- unlist(ht_obj$extra_legends)
  if (is.null(legends)) {
    return(ComplexHeatmap::draw(ht))
  }
  
  legend_titles <- vapply(legends, function(lgd) {
    if (inherits(lgd, "Legend")) lgd@title else NA_character_
  }, character(1))
  
  unique_idx <- !duplicated(legend_titles)
  unique_legends <- legends[unique_idx]
  
  ComplexHeatmap::draw(
    ht,
    annotation_legend_list = unique_legends,
    heatmap_legend_side = heatmap_legend_side,
    annotation_legend_side = annotation_legend_side,
    #heatmap_legend_direction = heatmap_legend_direction,
    #annotation_legend_direction = annotation_legend_direction,
    padding = padding, # extra space on right
    ...)
}
