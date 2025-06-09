
#### Rapid function to plot histogram (with or without density) ####

plot_multi_histogram <- function(df, density = TRUE, feature, label_column, vline) {
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

#### Function to plot Pie Charts ####
#'db should be an AIRR formated database with a defined "clone_id" column
#'Define clones to plot and colors based on clones properties:
#'Expanded clone = clone size >1 in one donor
#'Shared clone = clone found in at least two donors
#'Persisting clone = clone found in at least two time points and/or pop
#'Color choice:
#'Single Greys scale based on clone size in that donor (to account for sampling): "white" = singlets, "black"= biggest clone in all donors
#'Color for shared clones between donors 
#'highlight = c("expanded", "top5")
#'Consistent color throughout time-points for persisting-expanded clones
#'Grey for non-persisting expanded clones
#'White for non expanded clones and regroup all sequences found only once in the dataset (non-persisting/non-expanded) 
#'requires packages : ComplexHeatmap, circlize, grid and dplyr

DonutPlotClonotypes3D <- function(db, 
                                  split.by = NULL, 
                                  prefix = NULL,
                                  use_chain = "IGH",
                                  locus = "locus",
                                  cell_id = "cell_id",
                                  clone_id = "clone_id",
                                  groups_to_plot = "all", 
                                  col.line = "black", 
                                  plots_folder = "Donut_plots", 
                                  highlight = c("shared", "clone_size", "clone_rank"), 
                                  #'also possible to pre-classify cells prior and set that column in highlight + a color scheme in highlight_col
                                  highlight_col = list("shared" = c("Paired", "random"),
                                                       "clone_size" = c("Set1"),
                                                       "clone_rank" = c("Set1")),
                                  external_bar = c("expanded", "top5", "none"),
                                  productive_only = TRUE){
  
  library(dplyr, quiet = TRUE)
  
  if(is.null(split.by)){
    db$origin <- "all"
    DonutPlotClonotypes(db,
                        split.by = "origin",
                        prefix = prefix,
                        use_chain = use_chain,
                        locus = locus,
                        cell_id = cell_id,
                        clone_id = clone_id,
                        groups_to_plot = groups_to_plot, 
                        col.line = col.line, 
                        plots_folder = plots_folder, 
                        highlight = highlight, 
                        highlight_col = highlight_col,
                        external_bar = external_bar,
                        productive_only = productive_only)
  }
  if(length(split.by) == 1){
    if(!split.by %in% colnames(db)){
      stop("split.by column does not exist in provided dataframe")
    }
    DonutPlotClonotypes(db,
                        split.by = split.by,
                        prefix = prefix,
                        use_chain = use_chain,
                        locus = locus,
                        cell_id = cell_id,
                        clone_id = clone_id,
                        groups_to_plot = groups_to_plot, 
                        col.line = col.line, 
                        plots_folder = plots_folder, 
                        highlight = highlight, 
                        highlight_col = highlight_col,
                        external_bar = external_bar,
                        productive_only = productive_only)
  }
  if(length(split.by) > 1){
    if(any(!split.by %in% colnames(db))){
      stop("not all split.by columns exist in the provided dataframe, missing: ", paste(split.by[!split.by %in% colnames(db)], collapse = ", "))
    }
    origin <- split.by[length(split.by)]
    split.by <- split.by[-length(split.by)]
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
                          use_chain = use_chain,
                          locus = locus,
                          cell_id = cell_id,
                          clone_id = clone_id,
                          groups_to_plot = groups_to_plot, 
                          col.line = col.line, 
                          plots_folder = plots_folder, 
                          highlight = highlight, 
                          highlight_col = highlight_col,
                          external_bar = external_bar,
                          productive_only = productive_only)
    }
    db <- db %>%
      dplyr::group_by(!!!rlang::syms(split.by)) %>%
      dplyr::group_split() %>%
      purrr::walk(plot_group)
  }
}


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
                                #'also possible to pre-classify cells prior 
                                highlight_col = list("shared" = c("Paired", "random"),
                                                     "clone_size" = c("Set1"),
                                                     "clone_rank" = c("Set1")),
                                external_bar = c("expanded", "top5", "none"),
                                productive_only = TRUE){
  
  library(grid, quiet = TRUE)
  library(circlize, quiet = TRUE)
  library(ComplexHeatmap, quiet = TRUE)
  library(RColorBrewer, quiet = TRUE)
  
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
  
  highlight = highlight[1]
  if(!highlight %in% c("shared", "clone_size", "clone_rank") & !highlight %in% colnames(db)){
    stop("highlight column not properly defined")
  }
  external_bar = external_bar[1]
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
    #'each shared clone gets its own color
    Shared_clones <- Clones_by_groups_to_plot[Clones_by_groups_to_plot$shared, clone_id]
    if(h_col %in% rownames(brewer.pal.info)){
      col_shared <- colorRampPalette(brewer.pal(n = brewer.pal.info[match(h_col, rownames(brewer.pal.info)), "maxcolors"], name = h_col))(length(Shared_clones))
    } else {
      col_shared <- randomcoloR::randomColor(length(Shared_clones))
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
    
    pdf(file = filename, height = 3, width = 3)
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

#### working basic hexmapClonotypes ####
#'TODO add documentation

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
                               seed = 42){
  
  library(dplyr, quiet = TRUE)
  library(RColorBrewer, quiet = TRUE)
  
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
      if(!is.null(prefix)){
        title = paste0(prefix, "_", group_name, "_clones_by_", highlight)
        filename = paste0(plots_folder, "/", prefix, "/", prefix, "_", group_name, "_by_", highlight,".pdf")
      } else {
        title = paste(group_name, "_clones_by_", highlight)
        filename = paste0(plots_folder, "/", group_name, "_by_", highlight,".pdf")
      }
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
      pdf(file = filename, height = 6, width = 6)
      plot(p)
      dev.off()
    }
    Plot_db <- Plot_db %>%
      dplyr::group_by(!!!rlang::syms(split.by)) %>%
      dplyr::group_split() %>%
      purrr::walk(plot_group)
  }
}

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
  
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggforce)
  library(packcircles)
  library(scales)
  library(purrr)
  
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
  centroid_layout <- circleProgressiveLayout(clone_sizes$R, sizetype = "radius") %>%
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



  