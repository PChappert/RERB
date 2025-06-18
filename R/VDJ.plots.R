
#### Rapid function to plot histogram (with or without density) ####
#' Plot a basic histogram
#'
#' \code{plotHisto}
#' @param df        a data frame with at least the feature column
#' @param feature   which column to use for plotting
#' @param density   whether to add a density curve
#' @param label_column column used for fill color
#' @param vline     value for added vline
#'
#' @return    a simple histogram
#'
#' @details
#' returns a simple histogram
#'
#' @export

plotHisto <- function(df, density = TRUE, feature, label_column, vline) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Optional: 'ggplot2' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(ggplot2))
  if(isTRUE(density)){
    plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
      geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
      geom_density(alpha=0.7) +
      geom_vline(xintercept = vline, color="black", linetype="dashed", size=1) +
      labs(x=feature, y = "Density")
    plt + guides(fill=guide_legend(title=label_column))
  } else {
    plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
      geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
      geom_vline(xintercept = vline, color="black", linetype="dashed", size=1) +
      labs(x=feature, y = "Density")
    plt + guides(fill=guide_legend(title=label_column))
  }
}

#### Function to plot Donut Plots  ####
#' Plots multiple donut plots and save them as pdf
#'
#' \code{DonutPlotClonotypes3D} Plots multiple donut plots and save them as pdf
#' @param db        an AIRR formatted dataframe containing bcr (heavy and light chains) or tcr (TCRA, TCRB, TCRG or TCRD) sequences. Should contain only one chain for each type per cell_id, if not run resolveMultiHC() first.
#' @param split.by  name of column to use to group sequence when calculating clone size and frequencies.
#' @param ...       arguments to be passed to DonutPlotClonotypes()
#'
#' @return multiple pdf donut plots
#'
#' @details
#' db should be an AIRR formated database with a defined "clone_id" column
#' Defines clones to plot and colors based on clones properties:
#' Expanded clone = clone size >1 in one donor
#' Shared clone = clone found in at least two donors
#' Persisting clone = clone found in at least two time points and/or pop
#' Color choice:
#' Single Greys scale based on clone size in that donor (to account for sampling): "white" = singlets, "black"= biggest clone in all donors
#' Color for shared clones between donors
#' highlight = c("expanded", "top5")
#' Consistent color throughout time-points for persisting-expanded clones
#' Grey for non-persisting expanded clones
#' White for non expanded clones and regroup all sequences found only once in the dataset (non-persisting/non-expanded)
#'
#' @import dplyr
#' @importFrom purrr walk
#'
#' @export

DonutPlotClonotypes3D <- function(db,
                                  split.by = NULL,
                                  prefix = NULL,
                                  ...){

  #library(dplyr, quiet = TRUE)

  if(is.null(split.by)){
    db$origin <- "all"
    DonutPlotClonotypes(db,
                        split.by = "origin",
                        ...)
    return(plots.list)
  }
  if(length(split.by) == 1){
    if(!split.by %in% colnames(db)){
      stop("split.by column does not exist in provided dataframe")
    }
    DonutPlotClonotypes(db,
                        split.by = split.by,
                        ...)
    return(plots.list)
  }
  if(length(split.by) > 1){
    if(any(!split.by %in% colnames(db))){
      stop("not all split.by columns exist in the provided dataframe, missing: ", paste(split.by[!split.by %in% colnames(db)], collapse = ", "))
    }
    origin <- split.by[length(split.by)]
    split.by <- split.by[-length(split.by)]
    
    if(any(is.na(db[[origin]]))){
      warning("missing values for '", origin, "' for some samples, will be set to 'unknown'")
      db <- db %>%
        dplyr::mutate(
          !!rlang::sym(origin) := ifelse(is.na(!!rlang::sym(origin)), "unknown", !!rlang::sym(origin))
        )
    }
    
    plot_group <- function(data) {
      # Extract group name
      group_name <- unique(data[[split.by[1]]])
      if(length(split.by) >1){
        for(i in seq(2, length(split.by)))
        group_name <- paste0(group_name, "_", unique(data[[split.by[i]]]))
      }
      
      # Plot data
      DonutPlotClonotypes(data,
                          split.by = origin,
                          prefix = paste(c(prefix, group_name), collapse = "_"),
                          ...)
    }
    db <- db %>%
      dplyr::group_by(!!!rlang::syms(split.by)) %>%
      dplyr::group_split() %>%
      purrr::walk(plot_group)
  }
}

#### Function to plot Donut Plots  ####
#' Plots individual donut plots and save as pdf
#'
#' \code{DonutPlotClonotypes} Plots individual donut plots and save as pdf
#' @param db        an AIRR formatted dataframe containing bcr (heavy and light chains) or tcr (TCRA, TCRB, TCRG or TCRD) sequences. Should contain only one chain for each type per cell_id, if not run resolveMultiHC() first.
#' @param split.by  name of column to use to group sequence when calculating clone size and frequencies.
#' @param prefix    prefix to use for saved files
#' @param use_chain which chain to use [default: "IGH"], each cell should only have one contig for this chain
#' @param locus     name of column containing locus values.
#' @param cell_id   name of the column containing cell identifier.
#' @param clone_id  name of the column containing cell identifier.
#' @param groups_to_plot which groups to plot
#' @param col.line  color for lines in donut plot
#' @param plots_folder name for export folder [default: "Donut_plots"]
#' @param highlight whether to highlight shared clones or clones based on clone_size or clone_rank. Also possible to pre-classify cells prior and set that column in highlight + a color scheme in highlight_col
#' @param highlight_col color scheme to use for highlighed clones
#' @param external_bar whether to add an external bar for "expanded" or "top5" clones.
#' @param productive  name of column containing productive calls.
#' @param productive_only whether to exclude non productive sequences [default: TRUE]
#' @param height    height of saved plot
#' @param width     width of saved plot
#'
#' @return a seurat object with additional metadata imported from the vdj_db dataframe.
#'
#' @import dplyr
#'
#' @export

DonutPlotClonotypes <- function(db,
                                split.by = NULL,
                                prefix = NULL,
                                use_chain = "IGH",
                                locus = "locus",
                                cell_id = "cell_id",
                                clone_id = "clone_id",
                                groups_to_plot = "all",
                                col.line = "black",
                                plots_folder = "Donut_plots",
                                antigen,
                                highlight = c("shared", "clone_size", "clone_rank"),
                                #also possible to pre-classify cells prior
                                highlight_col = list("shared" = c("Paired", "random"),
                                                     "clone_size" = c("Set1"),
                                                     "clone_rank" = c("Set1")),
                                external_bar = c("expanded", "top5", "none"),
                                productive = "productive",
                                productive_only = TRUE,
                                height = 3,
                                width = 3){

  if (!requireNamespace("circlize", quietly = TRUE)) {
    message("Optional: 'circlize' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(circlize))
  
  if (!requireNamespace("grid", quietly = TRUE)) {
    message("Optional: 'grid' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(grid))
  
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    message("Optional: 'RColorBrewer' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(RColorBrewer))
  
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    message("Optional: 'ComplexHeatmap' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(ComplexHeatmap))

  if(isFALSE(dir.exists(plots_folder))){
    dir.create(plots_folder)
  }

  if(!any(use_chain %in% c("IGH", "IGL", "IGK"))){
    stop("use_chain should be one of IGH, IGL or IGK")
  }
  
  if(!clone_id %in% colnames(db)){
    stop(paste0("missing",  clone_id, "collumn"))
  }

  if(!is.null(prefix)){
    if(isFALSE(dir.exists(paste0(plots_folder, "/", prefix)))){
      dir.create(paste0(plots_folder, "/", prefix))
    }
  }

  highlight <- match.arg(highlight) 
  if(!highlight %in% c("shared", "clone_size", "clone_rank") & !highlight %in% colnames(db)){
    stop("highlight column not properly defined")
  }
  external_bar <- match.arg(external_bar) 
  if(!external_bar %in% c("expanded", "top5", "none")){
    warning("external_bar parameter not properly defined, should be one of expanded, top5 or none, defaulting to none")
    external_bar = "none"
  }

  if(is.null(split.by)){
    Plot_db$origin <- "all"
    origin = "origin"
  }
  if(length(split.by)>1){
    stop("if grouping on more than one parameter, use DonutPlotClonotypes3D")
  }
  if(length(split.by)==1){
    origin = split.by
  }

  Plot_db <- db %>%
    dplyr::filter(!!rlang::sym(locus) %in% use_chain)
  
  if(any(duplicated(Plot_db[[cell_id]]))){
    stop("duplicated cell_id: ", paste(Plot_db[duplicated(Plot_db[[cell_id]]), cell_id], collapse = ", "))
  }

  ##Generate clone repartition table (based on origin and clone_id columns):
  Clones_by_groups_to_plot <- as.data.frame.matrix(table(Plot_db[[clone_id]], Plot_db[[origin]]))
  Clones_by_groups_to_plot <- Clones_by_groups_to_plot %>%
    dplyr::mutate(
      shared = rowSums(dplyr::select(., 1:length(levels(as.factor(Plot_db[[origin]])))) != 0) > 1,
      overall_clone_size = rowSums(dplyr::select(., where(is.numeric))),
      clone_id = rownames(.)
    ) %>%
    dplyr::arrange(desc(overall_clone_size))

  h_col <- highlight_col[[highlight]][1]

  ##Define color for plotting shared clones:
  if(highlight == "shared"){
    #each shared clone gets its own color
    Shared_clones <- Clones_by_groups_to_plot[Clones_by_groups_to_plot$shared, clone_id]
    if(h_col %in% rownames(brewer.pal.info)){
      col_shared <- colorRampPalette(brewer.pal(n = brewer.pal.info[match(h_col, rownames(brewer.pal.info)), "maxcolors"], name = h_col))(length(Shared_clones))
    } else {
      if (!requireNamespace("randomcoloR", quietly = TRUE)) {
        message("To use random color, you need to install: 'randomcoloR' — defaulting to 'Paired' palette")
        h_col <- "Paired"
        col_shared <- colorRampPalette(brewer.pal(n = brewer.pal.info[match(h_col, rownames(brewer.pal.info)), "maxcolors"], name = h_col))(length(Shared_clones))
      } else {
        col_shared <- randomcoloR::randomColor(length(Shared_clones))
      }
    }
    names(col_shared) <- Shared_clones
  }

  ##create list of dataframes for plotting of circosplots for each time points:
  origins <- levels(as.factor(Plot_db[[origin]]))
  Clones_by_groups_to_plot.list <- lapply(origins, FUN = function(y){
    data <- Clones_by_groups_to_plot %>%
      dplyr::select(all_of(c(y, "shared", "overall_clone_size", "clone_id"))) %>%
      dplyr::rename_with(~ "clone_size_in_group", .cols = 1) %>% #rename first column (y)
      dplyr::filter(clone_size_in_group > 0) %>%
      dplyr::arrange(desc(clone_size_in_group)) %>%
      dplyr::mutate(
        expanded = clone_size_in_group > 1,
        clone_rank = row_number()
      ) 
    
    n_unique_seq <- nrow(dplyr::filter(data, !expanded & !shared))
    
    if(n_unique_seq>0){
      data <- data %>%
        dplyr::filter(expanded | shared) %>%
        dplyr::bind_rows(tibble(clone_id = "unique", 
                                clone_size_in_group = n_unique_seq,
                                expanded = FALSE,
                                shared = FALSE)) 
    }
    
    if(highlight == "shared"){
      data <- data %>%
        dplyr::mutate(
          color = ifelse(shared == TRUE, col_shared[clone_id],
                         ifelse(expanded == TRUE, "grey", "white"))
        )
    }
    if(highlight == c("clone_rank")){
      #'each expanded clone gets its own color
      Expanded_clones <- data[data$expanded, clone_id]
      if(length(Expanded_clones)>0){
        if(h_col %in% rownames(brewer.pal.info)){
          col_expanded <- colorRampPalette(brewer.pal(n = brewer.pal.info[match(h_col, rownames(brewer.pal.info)), "maxcolors"], name = h_col))(length(Expanded_clones))
        } else {
          warning("highlight_col should be one of ColorBrewer Palettes, defaulting to Set1")
          col_expanded <- colorRampPalette(brewer.pal(n = 9, name = "Set1"))(length(Expanded_clones))
        }
        names(col_expanded) <- Expanded_clones
        data <- data %>%
          dplyr::mutate(
            color = ifelse(expanded == TRUE, col_expanded[clone_id], "white")
          )
      } else {
        data$color = "white"
      }
    }
    if(highlight == "clone_size"){
      #'color is defined based on clone size
      diff_sizes <- data[!duplicated(data$clone_size_in_group) & data$clone_size_in_group>1,]$clone_size_in_group
      if(length(diff_sizes)>0){
        if(h_col %in% rownames(brewer.pal.info)){
          col_expanded <- colorRampPalette(brewer.pal(n = brewer.pal.info[match(h_col, rownames(brewer.pal.info)), "maxcolors"], name = h_col))(length(diff_sizes))
        } else {
          warning("highlight_col should be one of ColorBrewer Palettes, defaulting to Set1")
          col_expanded <- colorRampPalette(brewer.pal(n = 9, name = "Set1"))(diff_sizes)
        }
        names(col_expanded) <- diff_sizes
        data <- data %>%
          dplyr::mutate(
            color = ifelse(expanded == TRUE, col_expanded[clone_size_in_group], "white")
          )
      } else {
        data$color ="white"
      }
    }  
    return(data)
  })
  names(Clones_by_groups_to_plot.list) <- unlist(origins)
  
  ##Plot Circosplots:
  for (i in seq_along(Clones_by_groups_to_plot.list)){
    Clones <- Clones_by_groups_to_plot.list[[i]]
    df_name <- names(Clones_by_groups_to_plot.list)[i]
    nb_seq <- sum(as.numeric(Clones$clone_size_in_group))
    
    #df for sectors to plot
    df1 <- as.data.frame(Clones$clone_size_in_group)
    names(df1) <- "nb_cells"
    df1$xmin <- 0
    df1$xmax <- df1$nb_cells
    df1$color <- Clones$color
    
    #Plot Circosplot
    if(!is.null(prefix)){
      filename = paste0(plots_folder, "/", prefix, "/", prefix, "_", names(Clones_by_groups_to_plot.list)[i], "_",  highlight,".pdf")
    } else {
      filename = paste0(plots_folder, "/", names(Clones_by_groups_to_plot.list)[i], "_",  highlight,".pdf")
    }
    
    #TODO add possibility to export plot or save as png
    # Function to draw and capture a circlize plot
    #library(gridGraphics)
    #library(gridExtra)
    #circos_plot_as_grob <- function(df) {
      #code for circosplot (what is currently inside pdf())
      #then capture the base R plot as a grid object
    #  grid.echo()
    #  grid.grab()
    #}
    # Example with 2 circlize plots arranged side by side at the end
    #g1 <- circos_plot_as_grob(df1)
    #g2 <- circos_plot_as_grob(df2)
    #grid.arrange(g1, g2, ncol = 2)
    
    
    pdf(file = filename, height = height, width = width)
    par(mar=rep(0,4))
    
    circos.clear()
    
    circos.par(cell.padding=c(0,0,0,0), track.margin=c(0, 0.02), 
               start.degree = 90, gap.degree = 0,
               canvas.xlim = c(-1.5,1.5),
               canvas.ylim = c(-1.5,1.5))
    
    circos.initialize(factors = rownames(df1), x = df1$nb_cells, xlim = cbind(df1$xmin, df1$xmax))
    
    #plot track with clone sizes:
    circos.trackPlotRegion(ylim = c(0, 1), factors = rownames(df1), track.height=0.6,
                           #panel.fun for each sector
                           panel.fun = function(x, y) {
                             name = get.cell.meta.data("sector.index")
                             i = get.cell.meta.data("sector.numeric.index")
                             xlim = get.cell.meta.data("xlim")
                             ylim = get.cell.meta.data("ylim")
                             circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                                         col = df1$color[i], border= "black", lty = 0.5)
                           })
    
    
    #'add nb sequences in the center:
    text(0,0, nb_seq)
    
    #'highlight expanded or top_5 clones:
    if(external_bar == "expanded"){
      #add outside track highlighting expanded clones:
      expanded_clones <- Clones[Clones$expanded == TRUE,]
      nb_expanded_clones <- length(expanded_clones$clone_id)
      size_of_last_expanded_clone <- as.numeric(expanded_clones[nb_expanded_clones,]$clone_size_in_group)
      pct_expanded <- round((sum(as.numeric(expanded_clones$clone_size_in_group))/nb_seq*100), digits =0)
      if (!nb_expanded_clones == 0){
        start_theta <- circlize(0, 0, sector.index = 1, track.index = 1)[1,"theta"]
        end_theta <- circlize(size_of_last_expanded_clone, 0, sector.index = nb_expanded_clones, track.index = 1)[1,"theta"]
        draw.sector(end_theta, max(start_theta), rou1 = 1.01, rou2 = 1.06, clock.wise = FALSE, col = "black")
      }
      #add text:
      text(0.45,1.15, paste0(pct_expanded, "%"))
    }
    
    if(external_bar == "top5"){
      #'add outside track highlighting top5 clones:
      #'accounts for cases with less than five clones (excluding unique sequences of course), no unique sequences or only unique sequences
      if("unique" %in% Clones$expanded){
        nb_top5_clones <- min(5, (length(Clones$clone_id)-1)) 
      } else {nb_top5_clones <- min(5, (length(Clones$clone_id)))}
      
      if(nb_top5_clones > 0){
        top5_clones <- Clones[1:nb_top5_clones,]
        size_of_last_top5_clone <- as.numeric(top5_clones[nb_top5_clones,]$clone_size_in_group)
        pct_top5 <- round((sum(as.numeric(top5_clones$clone_size_in_group))/nb_seq*100), digits =0)
        start_theta <- circlize(0, 0, sector.index = 1, track.index = 1)[1,"theta"]
        end_theta <- circlize(size_of_last_top5_clone, 0, sector.index = nb_top5_clones, track.index = 1)[1,"theta"]
        draw.sector(end_theta, max(start_theta), rou1 = 1.01, rou2 = 1.06, clock.wise = FALSE, col = "#002147", border = "#002147")
        #add text:
        text(0.45,1.15, paste0(pct_top5, "%"), col="#002147")
      } else {text(0.45,1.15, paste0("0%"), col="#002147")}
    }
    
    ##draw legend
    if(highlight == "clone_size"){
      if(nrow(Clones[!duplicated(Clones$clone_size_in_group) & Clones$expanded == TRUE,]) == 0){
        unique <- Clones[!duplicated(Clones$clone_size_in_group),]
        unique$clone_size_in_group <- 1
      } else {
        unique <- Clones[!duplicated(Clones$clone_size_in_group) & Clones$expanded == TRUE,]
        if(isTRUE("unique" %in% Clones$expanded)){
          unique <- rbind(unique, c("1", "FALSE", "1", "unique", "FALSE", "white"))
        }
      }
      lgd_clone = ComplexHeatmap::Legend(labels = unique$clone_size_in_group, 
                                         legend_gp = gpar(fill = unique$color), 
                                         title = "clone size",
                                         border = "black")
      
      ComplexHeatmap::draw(lgd_clone, x = unit(0.2, "in"), y = unit(2.4, "in"), just = "left")
    }
    
    dev.off()
  }
}

#### Function to plot hexmap from repertoire data ####
#' Plots multiple hexmap plots and save as pdf
#'
#' \code{HexmapClonotypes3D} Plots multiple donut plots and save as pdf
#' @param db        an AIRR formatted dataframe containing bcr (heavy and light chains) or tcr (TCRA, TCRB, TCRG or TCRD) sequences. Should contain only one chain for each type per cell_id, if not run resolveMultiHC() first.
#' @param split.by  name of column to use to group sequence when calculating clone size and frequencies.
#' @param ordered   whether to order clones based on size and highlighing parameter
#' @param highlight name of the column use to color each clone [default: "c_call].
#' @param highlight_col color scheme to use for coloring clones
#' @param prefix    prefix to use for saved files
#' @param use_chain which chain to use [default: "IGH"], each cell should only have one contig for this chain
#' @param locus     name of column containing locus values.
#' @param cell_id   name of the column containing cell identifier.
#' @param clone_id  name of the column containing cell identifier.
#' @param groups_to_plot which groups to plot
#' @param col.line  color for lines in donut plot
#' @param plots_folder name for export folder [default: "Donut_plots"]
#' @param productive  name of column containing productive calls.
#' @param productive_only whether to exclude non productive sequences [default: TRUE]
#' @param radius    shape radius, passed to packcircles
#' @param padding   padding, passed to packcircles
#' @param shape     shape to use [default: hex]
#' @param seed      random seed for plotting
#' @param save_plot whether to save the plots [default: TRUE]
#' @param save_as   format to save the plot ("pdf", or "png") [default: pdf]
#' @param height    height of saved plot
#' @param width     width of saved plot
#' @param return_plot whether to return the plots [default: FALSE]
#' @param return_coords whether to return the plots coordinates [default: FALSE]
#'
#' @return a hexmap plot
#'
#' @import dplyr
#'
#' @export

HexmapClonotypes3D <- function(db,
                               split.by = NULL,
                               ordered = TRUE,
                               highlight = "c_call",
                               highlight_col = list("c_call" = c("IGHM"= "green", "IGHD" = "lightgreen", "IGHA1"= "orange", "IGHA2"= "red",
                                                                 "IGHG1" = "blue", "IGHG2" = "lightblue", "IGHG3" = "lightblue", "IGHG4" = "lightblue",
                                                                 "IGHE"= "brown")),
                               prefix = NULL,
                               plots_folder = "Hexbin_plots",
                               use_chain = "IGH",
                               locus = "locus",
                               cell_id = "cell_id",
                               clone_id = "clone_id",
                               productive = "productive",
                               productive_only = TRUE,
                               radius = 1,
                               padding = 0.2,
                               shape = c("hex", "circle"),
                               seed = 42,
                               save_plot = TRUE,
                               save_as = c("pdf", "png"),
                               height = 6,
                               width = 6,
                               return_plot = FALSE,
                               return_coords = FALSE){

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    message("Optional: 'RColorBrewer' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(RColorBrewer))
  
  save.as <- match.arg(save.as)

  if(isFALSE(dir.exists(plots_folder))){
    dir.create(plots_folder)
  }

  if(!any(use_chain %in% c("IGH", "IGL", "IGK"))){
    stop("use_chain should be one of IGH, IGL or IGK")
  }
  Plot_db <- db %>%
    dplyr::filter(!!rlang::sym(locus) %in% use_chain)

  if(!clone_id %in% colnames(db)){
    stop(paste0("missing",  clone_id, "collumn"))
  }

  if(!is.null(prefix)){
    if(isFALSE(dir.exists(paste0(plots_folder, "/", prefix)))){
      dir.create(paste0(plots_folder, "/", prefix))
    }
  }

  if(highlight %in% names(highlight_col)){
    palette <- highlight_col[[highlight]]
    Plot_db$origin <- ifelse(Plot_db[[highlight]] %in% names(palette), Plot_db[[highlight]], NA)
  } else {
    #palette <- NULL
    levels <- levels(as.factor(Plot_db[[highlight]]))
    palette <- colorRampPalette(brewer.pal(n = 12, name = "Paired"))(length(levels))
    names(palette) <- levels
    Plot_db$origin <- Plot_db[[highlight]]
  }

  if(is.null(split.by)){
    if(!is.null(prefix)){
      title = paste0(prefix, "All_sequences_clones_by_", highlight)
      filename = paste0(plots_folder, "/", prefix, "/", prefix, "All_sequences_by_", highlight,".pdf")
    } else {
      title = paste("All_sequences_clones_by_", highlight)
      filename = paste0(plots_folder, "/All_sequences_by_", highlight,".pdf")
    }
    p <- HexmapClonotypes(Plot_db,
                          ordered = ordered,
                          radius = radius,
                          padding = padding,
                          fill_col = "origin",
                          cell_id = cell_id,
                          clone_id = clone_id,
                          title = title,
                          shape = shape,
                          palette = palette,
                          seed = seed,
                          return_coords = FALSE)
    pdf(file = filename, height = 6, width = 6)
    plot(p)
    dev.off()
  } else {
    if(any(!split.by %in% colnames(db))){
      stop("not all split.by columns exist in the provided dataframe, missing: ", paste(split.by[!split.by %in% colnames(db)], collapse = ", "))
    }

    plot_group <- function(data) {
      # Extract group name
      group_name <- unique(data[[split.by[1]]])
      if(length(split.by) >1){
        for(i in seq(2, length(split.by)))
          group_name <- paste0(group_name, "_", unique(data[[split.by[i]]]))
      }
      # Plot data
      p <- HexmapClonotypes(data,
                            ordered = ordered,
                            radius = radius,
                            padding = padding,
                            fill_col = "origin",
                            cell_id = cell_id,
                            clone_id = clone_id,
                            title = title,
                            shape = shape,
                            palette = palette,
                            seed = seed,
                            return_coords = FALSE)
      
      if(save_plot){
        if(!is.null(prefix)){
          title = paste0(prefix, "_", group_name, "_clones_by_", highlight)
          filename = paste0(plots_folder, "/", prefix, "/", prefix, "_", group_name, "_by_", highlight,".pdf")
        } else {
          title = paste(group_name, "_clones_by_", highlight)
          filename = paste0(plots_folder, "/", group_name, "_by_", highlight,".pdf")
        }
        if(save_as == "pdf"){
          pdf(file = filename, height = height, width = width)
          plot(p)
          dev.off()
        }
        if(save_as == "png"){
          png(filename = filename, height = height, width = width)
          plot(p)
          dev.off()
        }
      }
      return(p)
    }
    
    groups <- Plot_db %>% 
      dplyr::group_by(!!!rlang::syms(split.by)) %>% 
      dplyr::group_nest()
    
    plots <- purrr::map(groups$data, plot_group)
    
    if(return_plot){
      group_names <- groups %>% 
        mutate(id = purrr::pmap_chr(., ~ paste(..., sep = "_"))) %>% 
        pull(id)
      
      names(plots) <- group_names
      return(plots)
    }
  }
}

#### Function to plot hexmap from repertoire data ####
#' Plots individual hexmap plots and save as pdf
#'
#' \code{HexmapClonotypes} Plots individual hexmap plots and save as pdf
#' @param db        an AIRR formatted dataframe containing bcr (heavy and light chains) or tcr (TCRA, TCRB, TCRG or TCRD) sequences. Should contain only one chain for each type per cell_id, if not run resolveMultiHC() first.
#' @param split.by  name of column to use to group sequence when calculating clone size and frequencies.
#' @param ordered   whether to order clones based on size and highlighing parameter
#' @param highlight name of the column use to color each clone [default: "c_call].
#' @param highlight_col color scheme to use for coloring clones
#' @param prefix    prefix to use for saved files
#' @param use_chain which chain to use [default: "IGH"], each cell should only have one contig for this chain
#' @param locus     name of column containing locus values.
#' @param cell_id   name of the column containing cell identifier.
#' @param clone_id  name of the column containing cell identifier.
#' @param groups_to_plot which groups to plot
#' @param col.line  color for lines in donut plot
#' @param plots_folder name for export folder [default: "Donut_plots"]
#' @param productive  name of column containing productive calls.
#' @param productive_only whether to exclude non productive sequences [default: TRUE]
#' @param radius    shape radius, passed to packcircles
#' @param padding   padding, passed to packcircles
#' @param shape     shape to use [default: hex]
#' @param seed      random seed for plotting
#' @param save_plot whether to save the plots [default: TRUE]
#' @param save_as   format to save the plot ("pdf", or "png") [default: pdf]
#' @param height    height of saved plot
#' @param width     width of saved plot
#' @param return_plot whether to return the plots [default: FALSE]
#' @param return_coords whether to return the plots coordinates [default: FALSE]
#'
#' @return a hexmap plot
#'
#' @import dplyr
#' @import tibble
#' @import purrr
#'
#' @export

HexmapClonotypes <- function(data,
                             ordered = TRUE,
                             radius = 1,
                             padding = 0.2,
                             fill_col = "origin",
                             cell_id = "cell_id",
                             clone_id = "clone_id",
                             title = "All_sequences",
                             shape = c("hex", "circle"),
                             palette = NULL,
                             seed = 42,
                             max_colors = 10,
                             return_coords = FALSE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Optional: 'ggplot2' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(ggplot2))
  
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    message("Optional: 'RColorBrewer' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(RColorBrewer))
  
  if (!requireNamespace("ggforce", quietly = TRUE)) {
    message("Optional: 'ggforce' not installed — skipping plot.")
    return(invisible(NULL))
  }
  
  if (!requireNamespace("packcircles", quietly = TRUE)) {
    message("Optional: 'packcircles' not installed — skipping plot.")
    return(invisible(NULL))
  }
  
  if (!requireNamespace("scales", quietly = TRUE)) {
    message("Optional: 'scales' not installed — skipping plot.")
    return(invisible(NULL))
  }
  
  shape <- match.arg(shape)

  stopifnot(all(c(cell_id, clone_id, fill_col) %in% colnames(data)))
  set.seed(seed)

  if(any(duplicated(data[[cell_id]]))){
    stop("duplicated ", cell_id, "s, fiilter before running HexmapClonotypes()")
  }

  if(ordered){
    data <- data %>%
      dplyr::arrange(fill_col)
  }

  # Step 1: Compute dominant origin per clone and add it to all sequence rows
  dominant_origin_lookup <- data %>%
    dplyr::count(!!rlang::sym(clone_id), !!rlang::sym(fill_col)) %>%
    dplyr::group_by(!!rlang::sym(clone_id)) %>%
    dplyr::slice_max(n, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(dominant_origin = !!rlang::sym(fill_col)) %>%
    dplyr::select(!!rlang::sym(clone_id), dominant_origin)

  # Add dominant origin to all rows
  data <- data %>%
    left_join(dominant_origin_lookup, by = clone_id)

  # Now count total number of sequences per clone
  clone_sizes <- data %>%
    dplyr::count(!!rlang::sym(clone_id), name = "n") %>%
    dplyr::left_join(dominant_origin_lookup, by = clone_id) %>%
    dplyr::arrange(desc(n), desc(dominant_origin))

  # Pack within each clone
  pack_within_clone <- function(n) {
    # Compute number of rows/cols to fit `n` tightly
    cols <- ceiling(sqrt(n))
    rows <- ceiling(n / cols)

    # Hex spacing: flat-topped hexes
    dx <- radius * 3/2       # horizontal distance between centers
    dy <- sqrt(3) * radius   # vertical distance between centers

    # Generate honeycomb grid positions
    coords <- expand.grid(row = 0:(rows - 1), col = 0:(cols - 1)) %>%
      slice(1:n) %>%
      mutate(
        dx = col * dx,
        dy = row * dy + ifelse(col %% 2 == 1, dy / 2, 0),
        r = radius
      )

    # Center layout around (0, 0)
    coords$dx <- coords$dx - mean(coords$dx)
    coords$dy <- coords$dy - mean(coords$dy)

    # Estimate bounding radius for spacing clones apart
    R <- max(sqrt(coords$dx^2 + coords$dy^2)) + radius + padding

    list(coords = coords, R = R)
  }

  packed_clones <- lapply(clone_sizes$n, pack_within_clone)
  clone_sizes$R <- sapply(packed_clones, function(x) x$R)

  # Arrange clone centroids in circular layout (largest in center)
  theta <- seq(0, 2*pi, length.out = nrow(clone_sizes) + 1)[-1]
  radial_dist <- scales::rescale(rev(clone_sizes$R), to = c(0, 1))
  centroid_layout <- packcircles::circleProgressiveLayout(clone_sizes$R, sizetype = "radius") %>%
    as_tibble() %>%
    rename(x_c = x, y_c = y)

  centroid_layout$clone_id <- clone_sizes$clone_id
  centroid_layout$R <- clone_sizes$R
  centroid_layout$n <- clone_sizes$n


  # Add positions to each sequence
  data <- data %>%
    dplyr::left_join(centroid_layout, by = "clone_id") %>%
    dplyr::arrange(desc(n))

  coords_list <- vector("list", nrow(centroid_layout))

  for (i in seq_len(nrow(centroid_layout))) {
    n <- centroid_layout$n[i]
    clone <- centroid_layout$clone_id[i]
    subset <- data %>%
      dplyr::filter(!!rlang::sym(clone_id) == clone)

    coords <- packed_clones[[i]]$coords
    coords[[cell_id]] <- subset[[cell_id]]
    coords$x <- coords$dx + centroid_layout$x_c[i]
    coords$y <- coords$dy + centroid_layout$y_c[i]
    coords[[fill_col]] <- subset[[fill_col]]

    coords_list[[i]] <- coords
  }

  plot_data <- bind_rows(coords_list)

  # Step 5: Build hexagon coordinates manually
  generate_hex_coords <- function(x_center, y_center, r, id, fill_val) {
    angle <- seq(0, 2 * pi, length.out = 7)
    tibble(
      x = x_center + r * cos(angle),
      y = y_center + r * sin(angle),
      group = id,
      fill = fill_val
    )
  }

  if (shape == "hex") {
    hex_df <- purrr::pmap_dfr(
      list(plot_data$x, plot_data$y, plot_data$r, plot_data[[cell_id]], plot_data[[fill_col]]),
      function(x, y, r, id, fill_val) {
        angle <- seq(0, 2 * pi, length.out = 7)
        tibble(
          x = x + r * cos(angle),
          y = y + r * sin(angle),
          group = id,
          fill = fill_val
        )
      }
    )
  } else {
    hex_df <- plot_data %>%
      mutate(x0 = x, y0 = y, group = cell_id, fill = .data[[fill_col]])
  }

  # Step 6: Assign colors
  hex_df$fill <- factor(hex_df$fill)
  if (is.null(palette)) {
    fill_levels <- levels(hex_df$fill)
    colors <- hue_pal()(min(length(fill_levels), max_colors))
    names(colors) <- fill_levels
  } else {
    colors <- palette
  }

  # Step 7: Plot
  p <- ggplot(hex_df, aes(x = x, y = y, group = group, fill = fill)) +
    {
      if (shape == "hex") geom_polygon(color = NA, alpha = 0.9)
      else ggforce::geom_circle(aes(x0 = x0, y0 = y0, r = r),
                                color = NA, alpha = 0.9)
    } +
    coord_equal() +
    scale_fill_manual(values = colors) +
    theme_void() +
    labs(title = title, fill = fill_col)

  if (return_coords) {
    return(list(plot = p, coords = plot_data))
  } else {
    return(p)
  }
}


#### Function to plot CDR3 logo from repertoire data ####
#' Plots individual CDR3 logo plots and save as pdf
#'
#' \code{plotCDR3logo} Plots individual hexmap plots and save as pdf
#' @param db            an AIRR formatted dataframe containing bcr (heavy and light chains) or tcr (TCRA, TCRB, TCRG or TCRD) sequences. Should contain only one chain for each type per cell_id, if not run resolveMultiHC() first.
#' @param split.by      name of column to use to group sequences.
#' @param locus         name of the column containing locus identifier.
#' @param chain         name of the chain to plot. [default: IGH]
#' @param plots_folder  path to folder for saving plots
#' @param min_size      minimun size of clones to plot
#' @param clone_id      name of the column containing clones identifier.
#' @param junction_type nt or aa
#' @param save_plot     whether to save the plot as a pdf
#' @param return_plot   whether to return the ggplot object
#' @param ...       additional arguments to pass to ggseqlogo
#'
#' @return a ggseqlogo plot
#'
#' @import dplyr
#'
#' @export

plotCDR3logo <- function(db,
                         split.by = NULL,
                         locus = "locus",
                         use_chain = "IGH",
                         plots_folder = "VDJ_Clones/CDR3_logo",
                         min_size = 1,
                         clone_id = "clone_id",
                         junction_type = "aa",
                         save_plot = TRUE,
                         return_plot = FALSE,
                         ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Optional: 'ggplot2' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(ggplot2))
  
  if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
    message("Optional: 'ggseqlogo' not installed — skipping plot.")
    return(invisible(NULL))
  }
  junction <- paste0("junction_", junction_type)
  if(!junction %in% colnames(db)){
    stop("missing '",junction,"' column in provided dataframe")
  }
  
  plot_cdr3 <- function(data){
    junctions <- data[[junction]]
    l = nchar(junctions[1])
    g = ggseqlogo::ggseqlogo(junctions, seq_type = junction_type, method = "prob", ...)
    p <- list(g, l)
    names(p) <- c("cdr3_logo", "cdr3_length")
    return(p)
  }
  
  # make sure clone_id is always at the end
  split.by <- c(split.by[split.by != clone_id], clone_id)
  
  groups <- db %>% 
    dplyr::filter(locus == use_chain) %>%
    dplyr::group_by(!!!rlang::syms(split.by)) %>% 
    dplyr::mutate(
      group_size = n()
    ) %>%
    dplyr::filter(group_size>=min_size) %>%
    dplyr::group_nest()
  
  plots <- purrr::map(groups$data, plot_cdr3)
  
  group_names <- groups %>% 
    dplyr::select(-data) %>%
    dplyr::mutate(id = purrr::pmap_chr(., ~ paste(..., sep = "_"))) %>% 
    dplyr::pull(id)
  
  names(plots) <- group_names
  
  if(save_plot){
    if(isFALSE(dir.exists(plots_folder))){
      dir.create(plots_folder)
    }
    for(i in seq_along(plots)){
      g <- plots[[i]][[1]]
      l <- plots[[i]][[2]]
      name <- names(plots)[i]
      ggsave(g, filename = paste0(plots_folder, "/", name,"_CDR3_logo.pdf"), width=((2+l)/4), height=2)
    }
  }
  if(return_plot){
    return(plots)
  }
}
                             








