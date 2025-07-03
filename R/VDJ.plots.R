#### Plotting function for QC graph (flagVDJdoublets and flagNonBdoublets) ####
#' plot QC graph
#'
#' \code{vdjQCplot} plot QC graph
#'
#' @param db            a data table
#' @param seq_type      type of sequences to be plotted: 'Ig' or 'TCR' [default = "Ig"].
#' @param use_chains    which chain(s) to use between "all" or "heavy" only
#' @param heavy         heavy chains to be used. if set to NULL, will default to IGH for Ig and TRB/TRD for TCR.
#' @param light         light chains to be used. if set to NULL, will default to IGL/IGK for Ig and TRA/TRG for TCR.
#' @param split.by      name of the column to use to group cells and learn dataset-specific distributions
#' @param type          whether to plot graphs related to B/B doublets or B/non-B doublets
#' @param output_folder name of the folder in which graph and recap excel workbooks will be saved [default = "VDJ_QC"].
#' @param plot_cutoff   whether to plot cut_off on graph
#' @param variable_cutoff   whether to use dataset-specific cut-offs
#' @param high_cutoff   if fixed cutoffs, value to use to mark high probability heavy chain doublets.
#' @param low_cutfoff   if fixed cutoffs, value to use to mark low-probability doublets
#' @param analysis_name suffix for final plot(s) name
#' @param na_name       name to use for NA values in colour.by
#' @param colour.by     parameter to use for coloring of dots
#' @param colour_code   colors to use for coloring of dots
#' @param locus         name of the column containing locus informations
#' @param heavy         locus value to use for heavy chain identification
#' @param cell_id       name of the column containing cell_ids
#' @param umi_count     name of the column containing umi_count informations
#' @param second_umi_count name of the column containing second_umi_count informations
#' @param save_plot     format to save the final ggplot2 object. "pdf" or "png". if set to any other value, no plot will be saved.
#' @param return_plot   whether to return the ggplot2 object
#'
#' @return   a pdf plot in the output folder if export = "pdf" and/or a ggplot object if export = "graph"
#'
#' @keywords internal
#'
#' @import dplyr

vdjQCplot <- function(db,
                      use_chains = "all",
                      seq_type = c("Ig", "TCR"),
                      heavy = NULL,
                      light = NULL,
                      split.by = NULL,
                      type = c("vdj_doublet", "nonB_vdj_doublet", "nonT_vdj_doublet"),
                      output_folder = NULL,
                      analysis_name = "All_sequences",
                      plot_cutoff = FALSE,
                      variable_cutoff = TRUE,
                      low_cutfoff = 10,
                      high_cutoff = 250,
                      na_name = "unknown",
                      colour.by = list(
                        "vdj_doublet" = "is.VDJ_doublet.confidence",
                        "nonB_vdj_doublet" = "is.nonB_VDJ_doublet.confidence",
                        "nonT_vdj_doublet" = "is.nonT_VDJ_doublet.confidence"
                      ),
                      colour_code = list(
                        "vdj_doublet" = c("high" = "darkred", "low" = "darkorange", "not a doublet" = "darkgreen"),
                        "nonB_vdj_doublet" = c("high" = "darkred", "low" = "darkorange", "ambient RNA" = "darkgreen", "B cell" = "grey"),
                        "nonT_vdj_doublet" = c("high" = "darkred", "low" = "darkorange", "ambient RNA" = "darkgreen", "T cell" = "grey")
                      ),
                      locus = "locus",
                      cell_id = "cell_id",
                      umi_count = "umi_count",
                      second_umi_count = "second_umi_count",
                      save_plot = c("pdf", "png"),
                      return_plot = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Optional: 'ggplot2' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(ggplot2))

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message("Optional: 'patchwork' not installed — skipping plot.")
    return(invisible(NULL))
  }
  save_plot <- match.arg(save_plot)
  seq_type <- match.arg(seq_type)

  if (is.null(heavy)) {
    if (seq_type == "Ig") {
      heavy <- "IGH"
    }
    if (seq_type == "TCR") {
      heavy <- c("TRB", "TRD")
    }
  }

  if (is.null(light)) {
    if (seq_type == "Ig") {
      light <- c("IGK", "IGL")
    }
    if (seq_type == "TCR") {
      light <- c("TRA", "TRG")
    }
  }

  if (!use_chains %in% c("all", "heavy")) {
    stop("use_chains argument must be one of 'all' or 'heavy'")
  } else {
    if (use_chains == "all") {
      chains <- c("heavy", "light")
    } else {
      chains <- "heavy"
    }
  }

  if (!is.null(split.by)) {
    groups <- levels(as.factor(db[[split.by]]))
  } else {
    db <- db %>%
      dplyr::mutate(
        split.by = "all_sequences"
      )
    split.by <- "split.by"
    groups <- "all_sequences"
  }

  colour.group <- colour.by[[type]]
  col <- c(colour_code[[type]], na_name = "blue")

  db <- db %>%
    dplyr::mutate(
      locus_simplified <- ifelse(locus %in% heavy, "heavy", ifelse(locus %in% heavy, "light", "other"))
    ) %>%
    dplyr::filter(locus_simplified %in% chains)

  table(db$locus_simplified)

  # generate recap plots for each sample:
  combined.plots <- lapply(groups, FUN = function(group) {
    # check color groups are well defined
    group_db <- db %>%
      dplyr::filter(!!rlang::sym(split.by) == group) %>%
      dplyr::mutate(
        !!rlang::sym(colour.group) := forcats::fct_na_value_to_level(!!rlang::sym(colour.group), level = na_name),
      )
    group_db[[colour.group]] <- factor(group_db[[colour.group]], levels = names(col))
    group_db <- group_db %>%
      dplyr::arrange(desc(!!rlang::sym(colour.group)))

    # generate dominant versus second umi plots for each chain
    if (type == "vdj_doublet") {
      if (plot_cutoff) {
        if (variable_cutoff) {
          # if not already in the dataframe, learn dataset (split.by) and locus specific cutoffs (adapted to library depth):
          # 1st quartile of dominant VDJ umi_count (in 75% of doublets the second VDJ should be above this value) and minimum of 10 and 1/10 of Median (below could be considered ambient RNA, i.e; 1/10 of an average cell).
          if (!(all(c("upper_cutoff", "lower_cutoff") %in% colnames(group_db)))) {
            group_db <- group_db %>%
              dplyr::group_by(locus_simplified) %>%
              dplyr::mutate(
                upper_cutoff = quantile(!!rlang::sym(umi_count), 0.25, na.rm = TRUE),
                lower_cutoff = ceiling(median(!!rlang::sym(umi_count), na.rm = TRUE) / 10)
              ) %>%
              dplyr::ungroup()
          }
        } else { # otherwise, use provided cutoffs
          group_db <- group_db %>%
            dplyr::mutate(
              upper_cutoff = high_cutoff,
              lower_cutoff = low_cutoff
            )
        }

        plots.list <- lapply(chains, FUN = function(chain) {
          plot_data <- group_db %>%
            dplyr::filter(locus_simplified %in% chain) %>%
            dplyr::mutate(
              !!rlang::sym(second_umi_count) := ifelse(is.na(!!rlang::sym(second_umi_count)), 0.1, !!rlang::sym(second_umi_count))
            )

          max <- max(plot_data[[umi_count]])
          upper_cutoff <- max(na.omit(plot_data[["upper_cutoff"]]))
          lower_cutoff <- max(na.omit(plot_data[["lower_cutoff"]]))

          line1 <- data.frame(x = seq(0.1, 3 * upper_cutoff))
          line1$y <- line1$x / 3
          line2 <- data.frame(x = seq(3 * lower_cutoff, max))
          line2$y <- lower_cutoff
          line3 <- data.frame(x = seq(3 * upper_cutoff, max))
          line3$y <- upper_cutoff

          if (chain == "heavy") {
            full_chain_list <- paste(heavy, collapse = "-")
          } else {
            full_chain_list <- paste(light, collapse = "-")
          }

          p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = !!rlang::sym(umi_count), y = !!rlang::sym(second_umi_count), colour = !!rlang::sym(colour.group))) +
            ggplot2::geom_point() +
            ggplot2::geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
            ggplot2::geom_line(data = line1, aes(x = x, y = y), color = "black", linetype = "dashed") +
            ggplot2::geom_line(data = line2, aes(x = x, y = y), color = "black", linetype = "dashed") +
            ggplot2::geom_line(data = line3, aes(x = x, y = y), color = "black", linetype = "dashed") +
            ggplot2::scale_color_manual(values = col) +
            ggplot2::labs(
              title = paste0(group, " - ", analysis_name),
              x = paste0(full_chain_list, " Dominant contig"),
              y = paste0(full_chain_list, " Secondary contig")
            ) +
            scale_x_log10() +
            scale_y_log10()

          return(p)
        })
        names(plots.list) <- chains
      } else {
        plots.list <- lapply(chains, FUN = function(chain) {
          plot_data <- group_db %>%
            dplyr::filter(locus_simplified %in% chain) %>%
            dplyr::mutate(
              !!rlang::sym(second_umi_count) := ifelse(is.na(!!rlang::sym(second_umi_count)), 0.1, !!rlang::sym(second_umi_count))
            )
          p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = !!rlang::sym(umi_count), y = !!rlang::sym(second_umi_count), colour = !!rlang::sym(colour.group))) +
            ggplot2::geom_point() +
            ggplot2::geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
            ggplot2::scale_color_manual(values = col) +
            ggplot2::labs(
              title = paste0(group, " - ", analysis_name),
              x = paste0(paste(unlist(mget(chain)), collapse = "-"), " Dominant contig"),
              y = paste0(paste(unlist(mget(chain)), collapse = "-"), " Secondary contig")
            ) +
            scale_x_log10() +
            scale_y_log10()
          return(p)
        })
        names(plots.list) <- chains
      }
    }

    # generate heavy versus light plots
    HvsL_plot_db <- group_db %>%
      dplyr::select(!!rlang::sym(cell_id), locus_simplified, !!rlang::sym(umi_count), !!rlang::sym(colour.group)) %>%
      tidyr::pivot_wider(names_from = locus_simplified, values_from = !!rlang::sym(umi_count))

    if ("heavy" %in% colnames(HvsL_plot_db)) {
      HvsL_plot_db <- HvsL_plot_db %>%
        dplyr::mutate(
          heavy = ifelse(is.na(heavy), 0.1, heavy)
        )
    }
    if ("light" %in% colnames(HvsL_plot_db)) {
      HvsL_plot_db <- HvsL_plot_db %>%
        dplyr::mutate(
          light = ifelse(is.na(light), 0.1, light)
        )
    }

    if (type == "vdj_doublet") {
      secondary_plot_db <- group_db %>%
        dplyr::filter(!is.na(!!rlang::sym(second_umi_count))) %>%
        dplyr::select(!!rlang::sym(cell_id), locus_simplified, !!rlang::sym(second_umi_count), !!rlang::sym(colour.group)) %>%
        tidyr::pivot_wider(names_from = locus_simplified, values_from = !!rlang::sym(second_umi_count))

      if ("heavy" %in% colnames(secondary_plot_db)) {
        secondary_plot_db <- secondary_plot_db %>%
          dplyr::mutate(
            heavy = ifelse(is.na(heavy), 0.1, heavy)
          )
      }
      if ("light" %in% colnames(secondary_plot_db)) {
        secondary_plot_db <- secondary_plot_db %>%
          dplyr::mutate(
            light = ifelse(is.na(light), 0.1, light)
          )
      }

      col <- c(col, "primary_umi_counts" = "grey")
      HvsL_plot_db <- HvsL_plot_db %>%
        dplyr::mutate(
          !!rlang::sym(colour.group) := "primary_umi_counts"
        ) %>%
        dplyr::bind_rows(secondary_plot_db)
    }

    p2 <- ggplot2::ggplot(HvsL_plot_db, ggplot2::aes(x = heavy, y = light, colour = !!rlang::sym(colour.group))) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = col) +
      ggplot2::labs(
        title = paste0(group, " - ", analysis_name),
        x = paste0(paste(heavy, collapse = "-"), " - UMI Count"),
        y = paste0(paste(light, collapse = "-"), " - UMI Count")
      ) +
      scale_x_log10() +
      scale_y_log10()

    if (type == "vdj_doublet") {
      plots.list[[length(chains) + 1]] <- p2
      # organize final page
      combined_plot <- plots.list[[1]] + plots.list[[2]]
      if (length(chains) > 1) {
        for (i in 3:(length(chains) + 1)) {
          combined_plot <- combined_plot + plots.list[[i]]
        }
        combined_plot <- combined_plot + patchwork::plot_layout(ncol = 2)
      } else {
        combined_plot <- combined_plot + patchwork::plot_layout(ncol = 1)
      }
    } else {
      combined_plot <- p2 + patchwork::plot_layout(ncol = 1)
    }
    return(combined_plot)
  })
  names(combined.plots) <- groups

  if (save_plot %in% c("pdf", "png")) {
    if (!is.null(output_folder)) {
      if (!stringr::str_ends(output_folder, "/")) {
        output_folder <- paste0(output_folder, "/")
      }
      if (!dir.exists(output_folder)) {
        dir.create(output_folder)
      }
    }
    if (type == "vdj_doublet") {
      filename <- paste0(output_folder, analysis_name, "_VDJ_doublets_QC_plots.pdf")
    }
    if (type == "nonB_vdj_doublet") {
      filename <- paste0(output_folder, analysis_name, "_nonB_VDJ_doublets_QC_plots.pdf")
    }
    if (type == "nonT_vdj_doublet") {
      filename <- paste0(output_folder, analysis_name, "_nonT_VDJ_doublets_QC_plots.pdf")
    }

    if (save_plot == "pdf") {
      pdf(filename, width = 11.69, height = 8.27)
      for (group in groups) {
        print(combined.plots[[group]])
      }
      dev.off()
    }

    if (save_plot == "png") {
      png(filename, width = 11.69, height = 8.27)
      for (group in groups) {
        print(combined.plots[[group]])
      }
      dev.off()
    }
  }
  if (return_plot) {
    return(combined.plots)
  }
}

#### Rapid function to plot histogram (with or without density) ####
#' Plot a basic histogram
#'
#' \code{plotHistogram}
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

plotHistogram <- function(df, density = TRUE, feature, label_column, vline) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Optional: 'ggplot2' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(ggplot2))
  if (isTRUE(density)) {
    plt <- ggplot(df, aes(x = eval(parse(text = feature)), fill = eval(parse(text = label_column)))) +
      geom_histogram(alpha = 0.7, position = "identity", aes(y = ..density..), color = "black") +
      geom_density(alpha = 0.7) +
      geom_vline(xintercept = vline, color = "black", linetype = "dashed", size = 1) +
      labs(x = feature, y = "Density")
    plt + guides(fill = guide_legend(title = label_column))
  } else {
    plt <- ggplot(df, aes(x = eval(parse(text = feature)), fill = eval(parse(text = label_column)))) +
      geom_histogram(alpha = 0.7, position = "identity", aes(y = ..density..), color = "black") +
      geom_vline(xintercept = vline, color = "black", linetype = "dashed", size = 1) +
      labs(x = feature, y = "Density")
    plt + guides(fill = guide_legend(title = label_column))
  }
}

#### Function to plot Donut Plots  ####
#' Plots multiple donut plots and save them as pdf
#'
#' \code{DonutPlotClonotypes} Plots multiple donut plots and save them as pdf
#' @param db        an AIRR formatted dataframe containing bcr (heavy and light chains) or tcr (TCRA, TCRB, TCRG or TCRD) sequences. Should contain only one chain for each type per cell_id, if not run resolveMultiHC() first.
#' @param prefix    prefix to use for saved files
#' @param split.by  name of column to use to group sequence when calculating clone size and frequencies.
#' @param ...       arguments to be passed to SingleDonutPlotClonotypes()
#'
#' @return
#' saved pdf for as many donut plots as there are groups in the dataframe as defined by the split.by argument.
#' if the split.by argument points to more than one column, will group the dataset based on the first element(s) of the split.by argument
#' and keep only the last as split.by argument of the SingleDonutPlotClonotypes function.
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
#' @importFrom purrr map
#'
#' @export

DonutPlotClonotypes <- function(db,
                                split.by = NULL,
                                prefix = NULL,
                                return_plot = FALSE,
                                ...) {
  # library(dplyr, quiet = TRUE)

  if (is.null(split.by)) {
    db$origin <- "all"
    SingleDonutPlotClonotypes(db,
      split.by = "origin",
      ...
    )

    return(plots.list)
  }
  if (length(split.by) == 1) {
    if (!split.by %in% colnames(db)) {
      stop("split.by column does not exist in provided dataframe")
    }
    SingleDonutPlotClonotypes(db,
      split.by = split.by,
      ...
    )

    return(plots.list)
  }
  if (length(split.by) > 1) {
    if (any(!split.by %in% colnames(db))) {
      stop("not all split.by columns exist in the provided dataframe, missing: ", paste(split.by[!split.by %in% colnames(db)], collapse = ", "))
    }
    origin <- split.by[length(split.by)]
    split.by <- split.by[-length(split.by)]

    if (any(is.na(db[[origin]]))) {
      warning("missing values for '", origin, "' for some samples, will be set to 'unknown'")
      db <- db %>%
        dplyr::mutate(
          !!rlang::sym(origin) := ifelse(is.na(!!rlang::sym(origin)), "unknown", !!rlang::sym(origin))
        )
    }

    plot_group <- function(data) {
      # Extract group name
      group_name <- unique(data[[split.by[1]]])
      if (length(split.by) > 1) {
        for (i in seq(2, length(split.by))) {
          group_name <- paste0(group_name, "_", unique(data[[split.by[i]]]))
        }
      }

      # Plot data
      if (!return_plot) {
        SingleDonutPlotClonotypes(data,
          split.by = origin,
          prefix = paste(c(prefix, group_name), collapse = "_"),
          return_plot = FALSE,
          ...
        )
      } else {
        grobs <- SingleDonutPlotClonotypes(data,
          split.by = origin,
          prefix = paste(c(prefix, group_name), collapse = "_"),
          return_plot = TRUE,
          ...
        )
        return(grobs)
      }
    }

    grouped <- db %>%
      dplyr::group_by(!!!syms(split.by))

    split_data <- grouped %>%
      dplyr::group_split(.keep = TRUE)

    group_names <- grouped %>%
      dplyr::group_keys() %>%
      tidyr::unite("group_label", everything(), sep = "_") %>%
      dplyr::pull(group_label)

    grobs.list <- purrr::map(split_data, plot_group)
    names(grobs.list) <- group_names

    # grobs.list <- db %>%
    #  dplyr::group_by(!!!rlang::syms(split.by)) %>%
    #  dplyr::group_split() %>%
    #  purrr::map(plot_group)

    if (return_plot) {
      return(grobs.list)
    }
  }
}

#### Function to plot Donut Plots  ####
#' Plots individual donut plots and save as pdf
#'
#' \code{SingleDonutPlotClonotypes} Plots individual donut plots and save as pdf
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
#' @param highlight whether to highlight shared clones or clones based on clone_size or clone_rank.
#' @param highlight_col color scheme to use for highlighed clones
#' @param external_bar whether to add an external bar for "expanded" or "top5" clones.
#' @param productive  name of column containing productive calls.
#' @param productive_only whether to exclude non productive sequences [default: TRUE]
#' @param height    height of saved plot
#' @param width     width of saved plot
#'
#' @return saved pdf for as many donut plots as there are groups in the dataframe as defined by the split.by argument
#'
#' @import dplyr
#'
#' @export

SingleDonutPlotClonotypes <- function(db,
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
                                      highlight_col = NULL,
                                      external_bar = c("expanded", "top5", "none"),
                                      productive = "productive",
                                      productive_only = TRUE,
                                      height = 3,
                                      width = 3,
                                      save_plot = TRUE,
                                      save_as = c("pdf", "png"),
                                      return_plot = FALSE) {
  save_as <- match.arg(save_as)
  if (save_plot) {
    if (!save_as %in% c("pdf", "png")) {
      warning("The only saving option are 'pdf' and 'png', defaulting to 'pdf'.")
      save_as <- "pdf"
    }
    if (!stringr::str_ends(plots_folder, "/")) {
      plots_folder <- paste0(plots_folder, "/")
    }
    if (isFALSE(dir.exists(plots_folder))) {
      dir.create(plots_folder)
    }
  }

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

  if (!any(use_chain %in% c("IGH", "IGL", "IGK"))) {
    stop("use_chain should be one of IGH, IGL or IGK")
  }

  if (!clone_id %in% colnames(db)) {
    stop(paste0("missing", clone_id, "collumn"))
  }

  if (!is.null(prefix)) {
    if (isFALSE(dir.exists(paste0(plots_folder, prefix)))) {
      dir.create(paste0(plots_folder, prefix))
    }
  }

  highlight <- match.arg(highlight)
  if (!highlight %in% c("shared", "clone_size", "clone_rank")) {
    stop("highlight not properly defined :should be one of 'shared', 'clone_size', 'clone_rank' or a column in the dataframe with an associated named color vector in highlight_col")
  }
  # default color scheme:
  if (is.null(highlight_col)) {
    pallettes <- list(
      "shared" = RColorBrewer::brewer.pal(n = 9, name = "Paired"),
      "clone_size" = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 5, name = "Set1"))(10),
      "clone_rank" = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 5, name = "Set1"))(25)
    )
    if (highlight %in% c("shared", "clone_size", "clone_rank")) {
      highlight_col <- pallettes[[highlight]]
    } else {
      warning("no provided color palette: defaulting to 'Set1'")
      highlight_col <- RColorBrewer::brewer.pal(n = 9, name = "Paired")
    }
  }
  external_bar <- match.arg(external_bar)
  if (!external_bar %in% c("expanded", "top5", "none")) {
    warning("external_bar parameter not properly defined, should be one of expanded, top5 or none, defaulting to none")
    external_bar <- "none"
  }

  if (is.null(split.by)) {
    Plot_db$origin <- "all"
    origin <- "origin"
  }
  if (length(split.by) > 1) {
    stop("if grouping on more than one parameter, use DonutPlotClonotypes3D")
  }
  if (length(split.by) == 1) {
    origin <- split.by
  }

  Plot_db <- db %>%
    dplyr::filter(!!rlang::sym(locus) %in% use_chain)

  if (any(duplicated(Plot_db[[cell_id]]))) {
    stop("duplicated cell_id: ", paste(Plot_db[duplicated(Plot_db[[cell_id]]), cell_id], collapse = ", "))
  }

  ## Generate clone repartition table (based on origin and clone_id columns):
  Clones_by_groups_to_plot <- as.data.frame.matrix(table(Plot_db[[clone_id]], Plot_db[[origin]]))
  Clones_by_groups_to_plot <- Clones_by_groups_to_plot %>%
    dplyr::mutate(
      shared = rowSums(dplyr::select(., 1:length(levels(as.factor(Plot_db[[origin]])))) != 0) > 1,
      overall_clone_size = rowSums(dplyr::select(., where(is.numeric))),
      clone_id = rownames(.)
    ) %>%
    dplyr::arrange(desc(overall_clone_size))

  ## Define color for plotting shared clones:
  if (highlight == "shared") {
    # each shared clone gets its own color
    shared_clones <- Clones_by_groups_to_plot[Clones_by_groups_to_plot$shared, clone_id]
    col_shared <- grDevices::colorRampPalette(highlight_col)(length(shared_clones))
    # if(length(shared_clones) > length(highlight_col)){
    #  col_shared <- grDevices::colorRampPalette(highlight_col)(length(shared_clones))
    # } else {col_shared <- highlight_col[1:(length(shared_clones))]}
    # optional: used random color
    # if (!requireNamespace("randomcoloR", quietly = TRUE)) {col_shared <- randomcoloR::randomColor(length(Shared_clones))}
    names(col_shared) <- shared_clones
  }

  ## create list of dataframes for plotting of circosplots for each time points:
  origins <- levels(as.factor(Plot_db[[origin]]))
  Clones_by_groups_to_plot.list <- lapply(origins, FUN = function(y) {
    data <- Clones_by_groups_to_plot %>%
      dplyr::select(all_of(c(y, "shared", "overall_clone_size", "clone_id"))) %>%
      dplyr::rename_with(~"clone_size_in_group", .cols = 1) %>% # rename first column (y)
      dplyr::filter(clone_size_in_group > 0) %>%
      dplyr::arrange(desc(clone_size_in_group)) %>%
      dplyr::mutate(
        expanded = clone_size_in_group > 1,
        clone_rank = row_number()
      )

    if (highlight == "shared") {
      n_unique_seq <- nrow(dplyr::filter(data, !expanded & !shared))
      if (n_unique_seq > 0) {
        data <- data %>%
          dplyr::filter(expanded | shared) %>%
          dplyr::bind_rows(tibble(
            clone_id = "unique",
            clone_size_in_group = n_unique_seq,
            expanded = FALSE,
            shared = FALSE
          ))
      }
      data <- data %>%
        dplyr::mutate(
          color = ifelse(shared == TRUE, col_shared[clone_id],
            ifelse(expanded == TRUE, "grey", "white")
          )
        )
    }
    if (highlight != "shared") {
      n_unique_seq <- nrow(dplyr::filter(data, !expanded))
      if (n_unique_seq > 0) {
        data <- data %>%
          dplyr::filter(expanded) %>%
          dplyr::bind_rows(tibble(
            clone_id = "unique",
            clone_size_in_group = n_unique_seq,
            expanded = FALSE
          ))
      }
      if (highlight == "clone_rank") {
        # each expanded clone gets its own color
        expanded_clones <- data[data$expanded, clone_id]
        if (length(expanded_clones) > length(highlight_col)) {
          extra_basic_color <- RColorBrewer::brewer.pal(n = ceiling(length(expanded_clones) / 5), name = "Set1")[-(1:5)]
          highlight_col <- c(highlight_col, grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = extra_basic_color, name = "Set1"))(length(extra_basic_color) * 5))
          col_expanded <- highlight_col[1:(length(expanded_clones))]
        } else {
          col_expanded <- highlight_col[1:(length(expanded_clones))]
        }
        names(col_expanded) <- expanded_clones
        data <- data %>%
          dplyr::mutate(
            color = ifelse(expanded == TRUE, col_expanded[clone_id], "white")
          )
      }
      if (highlight == "clone_size") {
        # color is defined based on clone size
        diff_sizes <- unique(data[data$expanded, ]$clone_size_in_group)
        diff_sizes <- rev(diff_sizes[order(diff_sizes)])
        if (length(diff_sizes) > length(highlight_col)) {
          extra_basic_color <- RColorBrewer::brewer.pal(n = ceiling(length(diff_sizes) / 5), name = "Set1")[-(1:5)]
          highlight_col <- c(highlight_col, grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = extra_basic_color, name = "Set1"))(length(extra_basic_color) * 5))
          col_expanded <- highlight_col[1:(length(diff_sizes))]
        } else {
          col_expanded <- highlight_col[1:(length(diff_sizes))]
        }
        names(col_expanded) <- as.character(diff_sizes)
        data <- data %>%
          dplyr::mutate(
            color = ifelse(expanded == TRUE, col_expanded[as.character(clone_size_in_group)], "white")
          )
      }
    }
    return(data)
  })
  names(Clones_by_groups_to_plot.list) <- unlist(origins)

  ## Plot Circosplots:
  if (save_plot) {
    for (i in seq_along(Clones_by_groups_to_plot.list)) {
      Clones <- Clones_by_groups_to_plot.list[[i]]
      df_name <- names(Clones_by_groups_to_plot.list)[i]
      nb_seq <- sum(as.numeric(Clones$clone_size_in_group))

      # df for sectors to plot
      df1 <- as.data.frame(Clones$clone_size_in_group)
      names(df1) <- "nb_cells"
      df1$xmin <- 0
      df1$xmax <- df1$nb_cells
      df1$color <- Clones$color

      # Plot Circosplot
      if (!is.null(prefix)) {
        filename <- paste0(plots_folder, prefix, "/", prefix, "_", names(Clones_by_groups_to_plot.list)[i], "_", highlight, ".pdf")
      } else {
        filename <- paste0(plots_folder, names(Clones_by_groups_to_plot.list)[i], "_", highlight, ".pdf")
      }

      if (save_as == "pdf") {
        pdf(file = filename, height = height, width = width)
      }
      if (save_as == "png") {
        png(file = filename, height = height, width = width)
      }

      par(mar = rep(0, 4))

      circos.clear()

      circos.par(
        cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0.02),
        start.degree = 90, gap.degree = 0,
        canvas.xlim = c(-1.5, 1.5),
        canvas.ylim = c(-1.5, 1.5)
      )

      circos.initialize(factors = rownames(df1), x = df1$nb_cells, xlim = cbind(df1$xmin, df1$xmax))

      # plot track with clone sizes:
      circos.trackPlotRegion(
        ylim = c(0, 1), factors = rownames(df1), track.height = 0.6,
        # panel.fun for each sector
        panel.fun = function(x, y) {
          name <- get.cell.meta.data("sector.index")
          i <- get.cell.meta.data("sector.numeric.index")
          xlim <- get.cell.meta.data("xlim")
          ylim <- get.cell.meta.data("ylim")
          circos.rect(
            xleft = xlim[1], ybottom = ylim[1], xright = xlim[2], ytop = ylim[2],
            col = df1$color[i], border = "black", lty = 0.5
          )
        }
      )

      # add nb sequences in the center:
      text(0, 0, nb_seq)

      # add title:
      title(names(Clones_by_groups_to_plot.list)[i], line = -1, cex.main = 1.5)

      # highlight expanded or top_5 clones:
      if (external_bar == "expanded") {
        # add outside track highlighting expanded clones:
        expanded_clones <- Clones[Clones$expanded == TRUE, ]
        nb_expanded_clones <- length(expanded_clones$clone_id)
        size_of_last_expanded_clone <- as.numeric(expanded_clones[nb_expanded_clones, ]$clone_size_in_group)
        pct_expanded <- round((sum(as.numeric(expanded_clones$clone_size_in_group)) / nb_seq * 100), digits = 0)
        if (!nb_expanded_clones == 0) {
          start_theta <- circlize(0, 0, sector.index = 1, track.index = 1)[1, "theta"]
          end_theta <- circlize(size_of_last_expanded_clone, 0, sector.index = nb_expanded_clones, track.index = 1)[1, "theta"]
          draw.sector(end_theta, max(start_theta), rou1 = 1.01, rou2 = 1.06, clock.wise = FALSE, col = "black")
        }
        # add text:
        text(0.45, 1.15, paste0(pct_expanded, "%"))
      }

      if (external_bar == "top5") {
        # add outside track highlighting top5 clones:
        # accounts for cases with less than five clones (excluding unique sequences of course), no unique sequences or only unique sequences
        if ("unique" %in% Clones$expanded) {
          nb_top5_clones <- min(5, (length(Clones$clone_id) - 1))
        } else {
          nb_top5_clones <- min(5, (length(Clones$clone_id)))
        }

        if (nb_top5_clones > 0) {
          top5_clones <- Clones[1:nb_top5_clones, ]
          size_of_last_top5_clone <- as.numeric(top5_clones[nb_top5_clones, ]$clone_size_in_group)
          pct_top5 <- round((sum(as.numeric(top5_clones$clone_size_in_group)) / nb_seq * 100), digits = 0)
          start_theta <- circlize(0, 0, sector.index = 1, track.index = 1)[1, "theta"]
          end_theta <- circlize(size_of_last_top5_clone, 0, sector.index = nb_top5_clones, track.index = 1)[1, "theta"]
          draw.sector(end_theta, max(start_theta), rou1 = 1.01, rou2 = 1.06, clock.wise = FALSE, col = "#002147", border = "#002147")
          # add text:
          text(0.45, 1.15, paste0(pct_top5, "%"), col = "#002147")
        } else {
          text(0.45, 1.15, paste0("0%"), col = "#002147")
        }
      }

      ## draw legend
      if (highlight == "clone_size") {
        if (nrow(Clones[!duplicated(Clones$clone_size_in_group) & Clones$expanded == TRUE, ]) == 0) {
          unique <- Clones[!duplicated(Clones$clone_size_in_group), ]
          unique$clone_size_in_group <- 1
        } else {
          unique <- Clones[!duplicated(Clones$clone_size_in_group) & Clones$expanded == TRUE, ]
          if (isTRUE("unique" %in% Clones$expanded)) {
            unique <- rbind(unique, c("1", "FALSE", "1", "unique", "FALSE", "white"))
          }
        }
        lgd_clone <- ComplexHeatmap::Legend(
          labels = unique$clone_size_in_group,
          legend_gp = gpar(fill = unique$color),
          title = "clone size",
          border = "black"
        )

        ComplexHeatmap::draw(lgd_clone, x = unit(0.2, "in"), y = unit(2.4, "in"), just = "left")
      }
      dev.off()
    }
  }

  if (return_plot) {
    # for visualisation only as it is currently impossible to directly export a circlize plot
    # TODO check circlizePlus package to plot circos plots directly in ggplot2

    if (!requireNamespace("magick", quietly = TRUE)) {
      message("Optional: 'magick' not installed — skipping plot.")
      return(invisible(NULL))
    }
    if (!requireNamespace("pdftools", quietly = TRUE)) {
      message("Optional: 'pdftools' not installed — skipping plot.")
      return(invisible(NULL))
    }
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
      message("Optional: 'gridExtra' not installed — skipping plot.")
      return(invisible(NULL))
    }
    if (!requireNamespace("grid", quietly = TRUE)) {
      message("Optional: 'grid' not installed — skipping plot.")
      return(invisible(NULL))
    }

    grobs.list <- list()
    for (i in seq_along(Clones_by_groups_to_plot.list)) {
      if (!is.null(prefix)) {
        filename <- paste0(plots_folder, prefix, "/", prefix, "_", names(Clones_by_groups_to_plot.list)[i], "_", highlight, ".pdf")
      } else {
        filename <- paste0(plots_folder, names(Clones_by_groups_to_plot.list)[i], "_", highlight, ".pdf")
      }

      if (save_as == "pdf") {
        img <- magick::image_read_pdf(filename, density = 150)
        tmp_png <- tempfile(fileext = ".png")
        magick::image_write(img, tmp_png)
        img_raster <- png::readPNG(tmp_png)
      }
      if (save_as == "png") {
        img_raster <- png::readPNG(filename)
      }
      grob <- grid::rasterGrob(img_raster, interpolate = TRUE)
      grobs.list[[i]] <- grob
    }
    names(grobs.list) <- names(Clones_by_groups_to_plot.list)
    full_combined_grob <- gridExtra::arrangeGrob(grobs = grobs.list, ncol = length(grobs.list))
    return(full_combined_grob)
  }
}

#### Function to plot hexmap from repertoire data ####
#' Plots multiple hexmap plots and save as pdf
#'
#' \code{HexmapClonotypes} Plots multiple donut plots and save as pdf
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

HexmapClonotypes <- function(db,
                             split.by = NULL,
                             ordered = TRUE,
                             highlight = "c_call",
                             highlight_col = list("c_call" = c(
                               "IGHM" = "green", "IGHD" = "lightgreen", "IGHA1" = "orange", "IGHA2" = "red",
                               "IGHG1" = "blue", "IGHG2" = "lightblue", "IGHG3" = "lightblue", "IGHG4" = "lightblue",
                               "IGHE" = "brown"
                             )),
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
                             return_coords = FALSE) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    message("Optional: 'RColorBrewer' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(RColorBrewer))

  save_as <- match.arg(save_as)

  if (!any(use_chain %in% c("IGH", "IGL", "IGK"))) {
    stop("use_chain should be one of IGH, IGL or IGK")
  }
  Plot_db <- db %>%
    dplyr::filter(!!rlang::sym(locus) %in% use_chain)

  if (!clone_id %in% colnames(db)) {
    stop(paste0("missing", clone_id, "collumn"))
  }

  if (!stringr::str_ends(plots_folder, "/")) {
    plots_folder <- paste0(plots_folder, "/")
  }
  if (save_plot) {
    if (isFALSE(dir.exists(plots_folder))) {
      dir.create(plots_folder)
    }
    if (!is.null(prefix)) {
      if (isFALSE(dir.exists(paste0(plots_folder, prefix)))) {
        dir.create(paste0(plots_folder, prefix))
      }
    }
  }

  if (highlight %in% names(highlight_col)) {
    palette <- highlight_col[[highlight]]
    Plot_db$origin <- ifelse(Plot_db[[highlight]] %in% names(palette), Plot_db[[highlight]], NA)
  } else {
    # palette <- NULL
    levels <- levels(as.factor(Plot_db[[highlight]]))
    palette <- colorRampPalette(brewer.pal(n = 12, name = "Paired"))(length(levels))
    names(palette) <- levels
    Plot_db$origin <- Plot_db[[highlight]]
  }

  if (is.null(split.by)) {
    if (!is.null(prefix)) {
      title <- paste0(prefix, "All_sequences\n (n=", nrow(Plot_db), ")")
      filename <- paste0(plots_folder, prefix, "/", prefix, "All_sequences_by_", highlight, ".pdf")
    } else {
      title <- paste("All_sequences\n (n=", nrow(Plot_db), ")")
      filename <- paste0(plots_folder, "All_sequences_by_", highlight, ".pdf")
    }
    p <- SingleHexmapClonotypes(Plot_db,
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
      return_coords = FALSE
    )

    if (save_plot) {
      if (save_as == "pdf") {
        pdf(file = filename, height = height, width = width)
        plot(p)
        dev.off()
      }
      if (save_as == "png") {
        png(filename = filename, height = height, width = width)
        plot(p)
        dev.off()
      }
    }
  } else {
    if (any(!split.by %in% colnames(db))) {
      stop("not all split.by columns exist in the provided dataframe, missing: ", paste(split.by[!split.by %in% colnames(db)], collapse = ", "))
    }

    groups <- Plot_db %>%
      dplyr::group_by(!!!rlang::syms(split.by)) %>%
      dplyr::group_nest()

    groups$group_name <- groups[[split.by[1]]]
    if (length(split.by) > 1) {
      for (i in seq(2, length(split.by))) {
        groups$group_name <- paste0(groups$group_name, "_", groups[[split.by[i]]])
      }
    }

    plots <- list()

    for (i in 1:length(groups$group_name)) {
      group_name <- groups$group_name[i]
      data <- groups$data[[i]]
      if (!is.null(prefix)) {
        title <- paste0(prefix, "_", group_name, "\n (n=", nrow(data), ")")
        filename <- paste0(plots_folder, prefix, "/", prefix, "_", group_name, "_by_", highlight, ".pdf")
      } else {
        title <- paste(group_name, "\n (n=", nrow(data), ")")
        filename <- paste0(plots_folder, group_name, "_by_", highlight, ".pdf")
      }
      # Plot data
      p <- SingleHexmapClonotypes(data,
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
        return_coords = FALSE
      )
      if (save_plot) {
        if (save_as == "pdf") {
          pdf(file = filename, height = height, width = width)
          plot(p)
          dev.off()
        }
        if (save_as == "png") {
          png(filename = filename, height = height, width = width)
          plot(p)
          dev.off()
        }
      }
      plots[[i]] <- p
    }
    names(plots) <- groups$group_name

    if (return_plot) {
      return(plots)
    }
  }
}

#### Function to plot hexmap from repertoire data ####
#' Plots individual hexmap plots and save as pdf
#'
#' \code{SingleHexmapClonotypes} Plots individual hexmap plots and save as pdf
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

SingleHexmapClonotypes <- function(data,
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

  if (any(duplicated(data[[cell_id]]))) {
    stop("duplicated ", cell_id, "s, fiilter before running HexmapClonotypes()")
  }

  if (ordered) {
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
    dx <- radius * 3 / 2 # horizontal distance between centers
    dy <- sqrt(3) * radius # vertical distance between centers

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
  theta <- seq(0, 2 * pi, length.out = nrow(clone_sizes) + 1)[-1]
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
    dplyr::tibble(
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
        dplyr::tibble(
          x = x + r * cos(angle),
          y = y + r * sin(angle),
          group = id,
          fill = fill_val
        )
      }
    )
  } else {
    hex_df <- plot_data %>%
      dplyr::mutate(x0 = x, y0 = y, group = cell_id, fill = .data[[fill_col]])
  }

  # Step 6: Assign colors
  hex_df$fill <- factor(hex_df$fill)
  if (is.null(palette)) {
    fill_levels <- levels(hex_df$fill)
    colors <- scales::hue_pal()(min(length(fill_levels), max_colors))
    names(colors) <- fill_levels
  } else {
    colors <- palette
  }

  # Step 7: Plot
  p <- ggplot2::ggplot(hex_df, ggplot2::aes(x = x, y = y, group = group, fill = fill)) +
    {
      if (shape == "hex") {
        ggplot2::geom_polygon(color = NA, alpha = 0.9)
      } else {
        ggforce::geom_circle(aes(x0 = x0, y0 = y0, r = r),
          color = NA, alpha = 0.9
        )
      }
    } +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_void() +
    ggplot2::labs(title = title, fill = fill_col)

  if (return_coords) {
    return(list(plot = p, coords = plot_data))
  } else {
    return(p)
  }
}


#### Function to plot CDR3 logos from repertoire data ####
#' Plots multiple individual CDR3 logo plots and save as pdf
#'
#' \code{plotCDR3logos} Batch plot individual CDR3 logo plots and save as pdf
#' @param db            an AIRR formatted dataframe containing bcr (heavy and light chains) or tcr (TCRA, TCRB, TCRG or TCRD) sequences. Should contain only one chain for each type per cell_id, if not run resolveMultiHC() first.
#' @param split.by      name of column to use to group sequences.
#' @param locus         name of the column containing locus identifier.
#' @param use_chains    name of the chain to plot. [default: heavy]
#' @param seq_type      'Ig' or 'TCR'
#' @param heavy         heavy chains to be used. if set to NULL, will default to IGH for Ig and TRB/TRD for TCR.
#' @param light         light chains to be used. if set to NULL, will default to IGL/IGK for Ig and TRA/TRG for TCR.
#' @param trim_junction whether to plot only the cdr3 when a junction is provided
#' @param junction      name of the column were junction sequences are stored. a default list for AIRR junction and junction_aa column is provided.
#' @param junction_type nt or aa
#' @param method        pass to ggseqlogo, method to used for plotting logo [default: "prob"]
#' @param plots_folder  path to folder for saving plots
#' @param min_size      minimun size of clones to plot
#' @param clone_id      name of the column containing clones identifier.
#' @param plot_all      whether to plot all sequence, avoids a warning message if no locus information is provided.
#' @param save_plot     whether to save the plot as a pdf
#' @param save_as       whether to save as 'pdf' or 'png'
#' @param return_plot   whether to return the ggplot object
#'
#' @return a ggseqlogo plot
#'
#' @import dplyr
#'
#' @export

plotCDR3logos <- function(db,
                          split.by = NULL,
                          locus = "locus",
                          seq_type = c("Ig", "TCR"),
                          use_chains = c("heavy", "light", "all"),
                          heavy = NULL,
                          light = NULL,
                          trim_junction = FALSE,
                          junction = list("aa" = "junction_aa", "dna" = "junction"),
                          junction_type = c("aa", "dna"),
                          method = c("prob", "bits"),
                          plots_folder = "VDJ_Clones/CDR3_logo",
                          min_size = 1,
                          clone_id = "clone_id",
                          plot_all = FALSE,
                          save_plot = TRUE,
                          save_as = c("pdf", "png"),
                          return_plot = FALSE) {
  seq_type <- match.arg(seq_type)

  if (is.null(heavy)) {
    if (seq_type == "Ig") {
      heavy <- "IGH"
    }
    if (seq_type == "TCR") {
      heavy <- c("TRB", "TRD")
    }
  }
  if (is.null(light)) {
    if (seq_type == "Ig") {
      light <- c("IGK", "IGL")
    }
    if (seq_type == "TCR") {
      light <- c("TRA", "TRG")
    }
  }

  # define chains to be used:
  use_chains <- match.arg(use_chains)

  if (!use_chains %in% c("heavy", "light")) {
    stop("use_chains must be one of heavy or light.")
  } else {
    if (use_chains == "heavy") {
      chains <- heavy
    } else {
      chains <- light
    }
  }

  if (!locus %in% colnames(db)) {
    if (!plot_all) {
      warning("No locus information provided, plotting all sequences")
    }
    db$locus <- "all"
    chains <- "all"
  }

  junction_type <- match.arg(junction_type)
  if (!junction_type %in% c("aa", "dna")) {
    stop("Choosen 'junction_type' should be one of 'aa' or 'dna', defaulting to 'prob'.")
  }
  if (is.list(junction)) {
    junction_column <- junction[[junction_type]]
  } else {
    junction_column <- junction[1]
  }
  if (!junction_column %in% colnames(db)) {
    stop("missing '", junction_column, "' column in provided dataframe")
  }

  if (is.null(split.by)) {
    split <- FALSE
    db <- db %>%
      dplyr::mutate(
        All_sequences = "all_sequences"
      )
    split.by <- "All_sequences"
  } else {
    split <- TRUE
  }

  # if clone_id is in split.by, make sure clone_id is always at the end and returns a warning if not all cdr3s in a clone are the same length
  if (clone_id %in% split.by) {
    split.by <- c(split.by[split.by != clone_id], clone_id)
    warn_length <- TRUE
  }

  groups <- db %>%
    dplyr::filter(!!rlang::sym(locus) %in% chains) %>%
    dplyr::group_by(!!!rlang::syms(split.by)) %>%
    dplyr::mutate(
      group_size = n()
    ) %>%
    dplyr::filter(group_size >= min_size) %>%
    dplyr::group_nest()

  plots <- purrr::map(groups$data, ~ plotCDR3logo(
    junctions = .x[[junction_column]],
    trim_junction = trim_junction,
    junction_type = junction_type,
    method = method,
    warn_length = warn_length,
    return = "all"
  ))

  if (split) {
    group_names <- groups %>%
      dplyr::select(-data) %>%
      dplyr::mutate(id = purrr::pmap_chr(., ~ paste(..., sep = "_"))) %>%
      dplyr::pull(id)
  } else {
    group_names <- "All_sequences"
  }

  names(plots) <- group_names

  if (save_plot) {
    if (!stringr::str_ends(plots_folder, "/")) {
      plots_folder <- paste0(plots_folder, "/")
    }
    if (isFALSE(dir.exists(plots_folder))) {
      dir.create(plots_folder)
    }
    for (i in seq_along(plots)) {
      g <- plots[[i]][[1]]
      l <- plots[[i]][[2]]
      name <- names(plots)[i]

      save_as <- match.args(save_as)
      if (save_as == "pdf") {
        ggsave(g, filename = paste0(plots_folder, "/", name, "_CDR3_logo.pdf"), width = ((2 + l) / 4), height = 2)
      }
      if (save_as == "png") {
        png(filename = paste0(plots_folder, "/", name, "_CDR3_logo.png"), width = ((2 + l) / 4), height = 2)
        plot(g)
        dev.off()
      }
    }
  }
  if (return_plot) {
    return(plots)
  }
}

#### Function to plot a single CDR3 logo from repertoire data ####
#' returns one CDR3 logo plot
#'
#' \code{plotCDR3logo} returns one CDR3 logo plot
#' @param junctions     a vector of sequences.
#' @param junction_type nt or aa
#' @param trim_junction whether to plot only the cdr3 when a junction is provided
#' @param method        pass to ggseqlogo, method to used for plotting logo [default: "prob"]
#' @param warn_length   whether to send a warning if sequences are of different lengths
#' @param return        whether to simply return the ggplot object or a list containing the ggplot object [default: "plot"] as well as information regarding the sequences lengths ["all"]
#'
#' @return a ggseqlogo plot or a named list with a ggseqlogo plot, a recap table of the length of all sequences provided and a vector of all sequences provided
#'
#' @import dplyr
#'
#' @export

plotCDR3logo <- function(junctions,
                         junction_type = c("aa", "dna"),
                         trim_junction = FALSE,
                         method = c("prob", "bits"),
                         warn_length = FALSE,
                         return = c("plot", "all")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Optional: 'ggplot2' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(ggplot2))

  if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
    message("Optional: 'ggseqlogo' not installed — skipping plot.")
    return(invisible(NULL))
  }

  return <- match.arg(return)

  method <- match.arg(method)
  if (!method %in% c("prob", "bits")) {
    warning("Choosen 'method' for ggseqlogo should be one of 'prob' or 'bits', defaulting to 'prob'.")
    method <- "prob"
  }

  junction_type <- match.arg(junction_type)
  if (!junction_type %in% c("aa", "dna")) {
    stop("Choosen 'junction_type' should be one of 'aa' or 'dna', defaulting to 'prob'.")
  }

  lengths <- as.data.frame(table(nchar(junctions))) %>%
    dplyr::rename(
      nb_seq = Freq,
      cdr3_length = Var1
    )
  if (nrow(lengths) == 1) {
    if (trim_junction) {
      if (junction_type == "aa") {
        junctions <- Biostrings::AAStringSet(junctions)
        # removing first and last AA from junction
        junctions <- subseq(junctions, start = 2, end = -2)
        junctions <- as.character(junctions)
      }
      if (junction_type == "dna") {
        junctions <- Biostrings::DNAStringSet(junctions)
        # removing first three and last three AA from junction
        junctions <- subseq(junctions, start = 4, end = -4)
        junctions <- as.character(junctions)
      }
    }
  }
  if (nrow(lengths) > 1) {
    method <- "bits" # using Prob doesn't take '-' generated by gaps in the alignments into consideration
    if (warn_length) {
      warning("presence of clones with different cdr3 length")
    }
    if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
      message("Alignement needed for cdr3s/junctions of variable length: 'msa' not installed — skipping plot. [try: BiocManager::install('msa')]")
      return(invisible(NULL))
    }
    if (junction_type == "aa") {
      junctions <- Biostrings::AAStringSet(junctions)
      if (trim_junction) {
        # removing first and last AA from junction
        junctions <- subseq(junctions, start = 2, end = -2)
      }
      junctions <- msa::msa(junctions, method = "ClustalW")
      junctions <- as.character(junctions)
    }
    if (junction_type == "dna") {
      junctions <- Biostrings::DNAStringSet(junctions)
      if (trim_junction) {
        # removing first three and last three AA from junction
        junctions <- subseq(junctions, start = 4, end = -4)
      }
      junctions <- msa::msa(junctions, method = "ClustalW")
      junctions <- as.character(junctions)
    }
  }
  g <- ggseqlogo::ggseqlogo(junctions, seq_type = junction_type, method = method)
  p <- list(g, lengths, junctions)
  names(p) <- c("cdr3_logo", "cdr3_length", "cdr3s")
  if (return == "plot") {
    return(g)
  }
  if (return == "all") {
    return(p)
  }
}

#### Function to plot shared clones between groups using a circosplot representation  ####
#' Plots individual donut plots and save as pdf
#'
#' \code{CircosClonotypes} Plots individual donut plots and save as pdf
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
#' @param highlight whether to highlight shared clones or clones based on clone_size or clone_rank.
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

CircosClonotypes <- function(db,
                             split.by = NULL,
                             groups_to_plot = NULL,
                             clone_id = "clone_id",
                             order = TRUE,
                             margin = 0.02,
                             track.height = 0.1,
                             label = TRUE,
                             col = NULL,
                             col.line = "black",
                             highlight = NULL,
                             highlight_col = "red") {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    message("Optional: 'circlize' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(circlize))

  # TODO needs to be updated

  # library(data.table)

  if (any(!split.by %in% colnames(db))) {
    stop("the following columns could not be found in the provided dataframe: ", split.by[!split.by %in% colnames(db)])
  }

  if (is.null(groups_to_plot)) {
    groups_to_plot <- levels(as.factor(db[, split.by]))
  }

  if (is.null(col)) {
    library("RColorBrewer")
    col <- brewer.pal(n = length(groups_to_plot), name = "Set1")
  }
  if (col.line == "match") {
    col.line <- col
  } else {
    col.line <- rep(col.line, length(groups_to_plot))
  }

  # generate clone repartition table:
  Clones_by_groups_to_plot <- as.data.frame.matrix(table(
    sc$clone_id,
    sc@meta.data[, group_by]
  ))
  Clones_by_groups_to_plot <- Clones_by_groups_to_plot[, groups_to_plot]
  Clones_by_groups_to_plot$overall_clone_size <- rowSums(Clones_by_groups_to_plot)
  Clones_by_groups_to_plot$clone_id <- rownames(Clones_by_groups_to_plot)

  if (isTRUE(order)) {
    Clones_by_groups_to_plot <- arrange(Clones_by_groups_to_plot, desc(Clones_by_groups_to_plot$overall_clone_size))
  }

  # dataframe defining the main sectors:
  df1 <- as.data.frame(colSums(Clones_by_groups_to_plot[, groups_to_plot]))
  colnames(df1) <- "nb_cells"
  df1$xmin <- 0
  df1$xmax <- df1$nb_cells + 1

  # list of dataframes defining the clone order and length of each sectors:
  groups <- as.list(groups_to_plot)
  df2.list <- lapply(X = groups, FUN = function(i) {
    df <- Clones_by_groups_to_plot[Clones_by_groups_to_plot[, i] >= 1, c("clone_id", i)]
    colnames(df) <- c("clone_id", "nb_cells")
    df <- arrange(df, desc(df$"nb_cells"))
    df$clone_segments_start <- cumsum(df$nb_cells) - df$nb_cells + 1
    df$clone_segments_end <- cumsum(df$nb_cells)
    return(df)
  })
  names(df2.list) <- groups_to_plot

  # list of dataframes defining the clonal relationships to plot as links:
  df3.list <- lapply(X = groups[1:length(groups) - 1], FUN = function(i, g = groups_to_plot) {
    min <- match(i, g) + 1
    max <- length(g)
    comparison <- g[min:max]
    df.all <- data.frame(
      clone_id = character(),
      size_orig = integer(),
      start_orig = numeric(),
      size_end = integer(),
      start_end = numeric(),
      stringsAsFactors = FALSE
    )
    for (k in 1:length(comparison)) {
      common_clones_k <- intersect(df2.list[[i]]$clone_id, df2.list[[comparison[k]]]$clone_id)
      if (length(common_clones_k) > 0) {
        is.common <- df2.list[[i]]$clone_id %in% common_clones_k
        df <- df2.list[[i]][is.common, c("clone_id", "nb_cells", "clone_segments_start")]
        colnames(df) <- c("clone_id", "size_orig", "start_orig")
        df$orig <- i
        df$end <- comparison[k]
        df$size_end <- df2.list[[comparison[k]]]$nb_cells[match(df$clone_id, df2.list[[comparison[k]]]$clone_id)]
        df$start_end <- df2.list[[comparison[k]]]$clone_segments_start[match(df$clone_id, df2.list[[comparison[k]]]$clone_id)]
        df.all <- rbind(df.all, df)
      }
    }
    return(df.all)
  })
  names(df3.list) <- unlist(groups[1:length(groups) - 1])
  df3 <- do.call(rbind, df3.list)

  ### Basic circos graphic parameters
  par(mar = rep(1, 4))
  plot <- circos.clear()

  plot <- circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0, margin), start.degree = 90, gap.degree = 1)

  ## creating sector for each celltype based on nb of cells in a clone per cell_type:
  plot <- circos.initialize(factors = rownames(df1), x = df1$nb_cells + 1, xlim = cbind(df1$xmin, df1$xmax))

  ## plotting tracks:
  plot <- circos.trackPlotRegion(
    ylim = c(0, 1), factors = rownames(df1), track.height = track.height,
    # panel.fun for each sector
    panel.fun = function(x, y) {
      # select details of current sector
      name <- get.cell.meta.data("sector.index")
      i <- get.cell.meta.data("sector.numeric.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")

      # text direction (dd) and adjusments (aa)
      theta <- circlize(mean(xlim), 1.3)[1, 1] %% 360
      dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
      aa <- c(1, 0.5)
      if (theta < 90 || theta > 270) aa <- c(0, 0.5)

      # plot population labels
      if (isTRUE(label)) {
        circos.text(x = mean(xlim), y = 1.7, labels = name, facing = dd, cex = 0.8, adj = aa)
      }

      # plot main sector
      circos.rect(
        xleft = xlim[1], ybottom = ylim[1], xright = xlim[2], ytop = ylim[2],
        col = col[i], border = "black", lty = 1
      )
    }
  )

  ## add segments based on clone size:
  for (i in seq_along(groups)) {
    pop <- groups[[i]]
    plot <- circos.segments(x0 = df2.list[[pop]]$clone_segments_end, x1 = df2.list[[pop]]$clone_segments_end, y0 = 0.25, y1 = 0.75, sector.index = pop, col = "black", lwd = 0.5)
  }

  ## plot links:
  for (i in 1:nrow(df3)) {
    plot <- circos.link(
      sector.index1 = df3$orig[i], point1 = c(df3$start_orig[i], df3$start_orig[i] + df3$size_orig[i] - 1),
      sector.index2 = df3$end[i], point2 = c(df3$start_end[i], df3$start_end[i] + df3$size_end[i] - 1),
      col = col.line[match(df3$orig[i], groups_to_plot)], lwd = 0.5
    )
  }

  if (!is.null(highlight)) {
    df4.list <- lapply(X = groups[1:length(groups) - 1], FUN = function(i, g = groups_to_plot) {
      min <- match(i, g) + 1
      max <- length(g)
      comparison <- g[min:max]
      df.all <- data.frame(
        clone_id = character(),
        size_orig = integer(),
        start_orig = numeric(),
        size_end = integer(),
        start_end = numeric(),
        stringsAsFactors = FALSE
      )
      for (k in 1:length(comparison)) {
        common_clones_k <- intersect(df2.list[[i]]$clone_id, df2.list[[comparison[k]]]$clone_id)
        highlight_clones_k <- intersect(common_clones_k, highlight)
        if (length(highlight_clones_k) > 0) {
          is.common <- df2.list[[i]]$clone_id %in% highlight_clones_k
          df <- df2.list[[i]][is.common, c("clone_id", "nb_cells", "clone_segments_start")]
          colnames(df) <- c("clone_id", "size_orig", "start_orig")
          df$orig <- i
          df$end <- comparison[k]
          df$size_end <- df2.list[[comparison[k]]]$nb_cells[match(df$clone_id, df2.list[[comparison[k]]]$clone_id)]
          df$start_end <- df2.list[[comparison[k]]]$clone_segments_start[match(df$clone_id, df2.list[[comparison[k]]]$clone_id)]
          df.all <- rbind(df.all, df)
        }
      }
      return(df.all)
    })
    names(df4.list) <- unlist(groups[1:length(groups) - 1])
    df4 <- do.call(rbind, df4.list)

    # plot highlight links
    for (i in 1:nrow(df4)) {
      plot <- circos.link(
        sector.index1 = df4$orig[i], point1 = c(df4$start_orig[i], df4$start_orig[i] + df4$size_orig[i] - 1),
        sector.index2 = df4$end[i], point2 = c(df4$start_end[i], df4$start_end[i] + df4$size_end[i] - 1),
        col = highlight_col, lwd = 0.5
      )
    }
  }
  return(plot)
}


#### Plot frequencies of VH/VL pairing ####
#' Plot frequencies of VH/VL pairing
#'
#' \code{plotVGenePairing} Plot frequencies of VH/VL pairing
#' @param db        an AIRR formatted dataframe containing bcr (heavy and light chains) or tcr (TCRA, TCRB, TCRG or TCRD) sequences. Should contain only one chain for each type per cell_id, if not run resolveMultiHC() first.
#' @param split.by  name of column to use to group sequence.
#' @param groups_to_plot which groups to plot
#' @param prefix    prefix to use for saved files
#' @param plots_folder name for export folder [default: "Vgene_plots"]
#' @param downsample_clones whether to reduce clone to one value
#' @param locus     name of column containing locus values.
#' @param heavy     which chain to use as heavy chain [default: "IGH"], each cell should only have one contig for this chain
#' @param light     which chain to use as hlights chain [default: "IGH"], each cell should only have one contig for this chain
#' @param v_call    name of column containing v calls.
#' @param clone_id  name of column containing clone affiliations.
#' @param plot_freq whether to plot frequenceis (TRUE) or absolute counts (FALSE) [#TODO]
#' @param plot_significance whether to highlight significant differences between group 1 and other groups (red borders for dots)
#' @param ref       which group in groups_to_plot to use as ref. [default: first group in groups_to_plot]
#' @param save_plot whether to save the plot.
#' @param save_as   options for plot saving format (pdf or png)
#' @param return_plot whether to return the plot as a ggplot object
#'
#' @return a ggplot object with associated frequency and binomial test results for inputed data.
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @export

plotVGenePairing <- function(db,
                             split.by = NULL, # max one column
                             groups_to_plot = NULL,
                             prefix = NULL,
                             plots_folder = "Vgene_plots",
                             downsample_clones = TRUE,
                             locus = "locus",
                             seq_type = c("Ig", "TCR"),
                             v_call = "v_call",
                             clone_id = "clone_id",
                             plot_freq = TRUE,
                             plot_significance = FALSE,
                             ref = NULL,
                             save_plot = FALSE,
                             save_as = c("pdf", "png"),
                             height = NULL,
                             width = NULL,
                             return_plot = TRUE,
                             ncol = 2) {
  save_as <- match.arg(save_as)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Optional: 'ggplot2' not installed — skipping plot.")
    return(invisible(NULL))
  }
  suppressMessages(library(ggplot2))

  seq_type <- match.arg(seq_type)

  if (seq_type == "Ig") {
    heavy <- "IGH"
    light <- c("IGK", "IGL")
  }
  if (seq_type == "TCR") {
    heavy <- c("TRB", "TRD")
    light <- c("TRA", "TRG")
  }

  if (!clone_id %in% colnames(db)) {
    stop(paste0("missing", clone_id, "collumn"))
  }

  if (is.null(split.by)) {
    split.by <- "all"
    db$all <- "all sequences"
  }

  if (!is.null(groups_to_plot)) {
    if (any(!groups_to_plot %in% levels(as.factor(db[[split.by]])))) {
      warning("The following group(s) to plot are not found in the ", split.by, " column, will be removed.")
      groups_to_plot <- groups_to_plot[groups_to_plot %in% levels(as.factor(db[[split.by]]))]
    }
    db <- db %>%
      dplyr::filter(!!rlang::sym(split.by) %in% groups_to_plot) %>%
      dplyr::mutate(
        !!rlang::sym(split.by) := factor(!!rlang::sym(split.by), levels = groups_to_plot)
      )
  }

  # Reformat the dataframe
  db <- db %>%
    dplyr::mutate(
      v_gene = alakazam::getGene(v_call, first = TRUE, collapse = TRUE, strip_d = TRUE, omit_nl = FALSE, sep = ","),
      locus_simplified = ifelse(locus %in% heavy, "VH", ifelse(locus %in% light, "VL", NA))
    )

  if (any(is.na(db$locus_simplified))) {
    n_NA <- length(db$locus_simplified[is.na(db$locus_simplified)])
    db <- db %>%
      dplyr::filter(is.na(locus_simplified))
    warning(n_NA, " contigs removed as non ", paste(heavy, collapse = ", "), " or ", paste(light, collapse = ", "))
  }

  # [optional] Select only only representant per clone:
  if (downsample_clones) {
    db <- db %>%
      dplyr::select(cell_id, v_gene, locus_simplified, !!rlang::sym(split.by), clone_id) %>%
      dplyr::group_by(clone_id, locus_simplified) %>%
      dplyr::summarise(
        cell_id = unique(cell_id)[1],
        !!rlang::sym(split.by) := unique(!!rlang::sym(split.by))[1],
        v_gene = unique(v_gene)[1],
        locus_simplified = unique(locus_simplified)[1],
        .groups = "drop"
      )
  } else {
    db <- db %>%
      dplyr::select(cell_id, v_gene, locus_simplified, !!rlang::sym(split.by))
  }

  # Calculate pairing frequencies in cells with paired VH and VL:
  db_wide <- db %>%
    tidyr::pivot_wider(names_from = locus_simplified, values_from = v_gene) %>%
    tidyr::drop_na(VH, VL, !!rlang::sym(split.by))

  # pair_counts <- df_wide %>%
  #  dplyr::count(heavy, light)

  pair_freq <- db_wide %>%
    dplyr::count(!!rlang::sym(split.by), VH, VL) %>%
    dplyr::group_by(!!rlang::sym(split.by)) %>%
    dplyr::mutate(freq = n / sum(n)) %>%
    dplyr::ungroup()

  # Ensure full VH-VL grid for both groups (fill missing with 0 freq)
  all_combinations <- tidyr::expand_grid(
    !!rlang::sym(split.by) := unique(pair_freq[[split.by]]),
    VH = unique(pair_freq$VH),
    VL = unique(pair_freq$VL)
  )

  pair_freq_complete <- all_combinations %>%
    dplyr::left_join(pair_freq, by = c(split.by, "VH", "VL")) %>%
    tidyr::replace_na(list(freq = 0))

  groups <- levels(as.factor(pair_freq_complete[[split.by]]))
  if (any(is.null(ref), !ref %in% groups)) {
    ref <- groups[1]
  }

  # Create expected frequency matrix from the first group in the list:
  expected_db <- pair_freq_complete %>%
    dplyr::filter(!!rlang::sym(split.by) == ref) %>%
    dplyr::group_by(VH, VL) %>%
    dplyr::mutate(
      total = ifelse(is.na(n), 0, n),
      expected_p = freq
    )

  # TODO add possibility to import external Expected_freq_matrix (and then do binomial test on all groups)
  # Convert to a matrix for lookup
  expected_freq_matrix <- expected_db %>%
    dplyr::select(VH, VL, expected_p) %>%
    tidyr::pivot_wider(names_from = VL, values_from = expected_p) %>%
    tibble::column_to_rownames("VH") %>%
    as.matrix()

  # Apply binomial tests to Group 2
  # Total counts for group B
  db_b_list <- lapply(groups[!groups == ref], FUN = function(group) {
    total_in_group <- pair_freq_complete %>%
      dplyr::filter(!!rlang::sym(split.by) == group) %>%
      dplyr::mutate(
        total = ifelse(is.na(n), 0, n)
      ) %>%
      dplyr::pull(total)

    db_group <- pair_freq_complete %>%
      dplyr::filter(!!rlang::sym(split.by) == group) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        n = ifelse(is.na(n), 0, n),
        expected_p = expected_freq_matrix[VH, VL],
        pval = binom.test(x = n, n = sum(total_in_group), p = expected_p, alternative = "two.sided")$p.value
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        pval_adj = p.adjust(pval, method = "fdr"),
        "Significance (p < 0.05)" = pval_adj < 0.05
      )
    return(db_group)
  })

  db_plot <- pair_freq_complete %>%
    dplyr::filter(!!rlang::sym(split.by) == groups[1]) %>%
    dplyr::mutate(pval = NA, pval_adj = NA, "Significance (p < 0.05)" = FALSE) %>%
    dplyr::bind_rows(db_b_list)

  if (plot_significance) {
    p <- ggplot2::ggplot(db_plot %>% dplyr::filter(freq > 0) %>% dplyr::arrange(freq), ggplot2::aes(x = VL, y = VH)) +
      ggplot2::geom_point(ggplot2::aes(size = freq, fill = freq, color = !!rlang::sym("Significance (p < 0.05)")), shape = 21, stroke = 0.5) +
      ggplot2::scale_color_manual(
        values = c("TRUE" = "red", "FALSE" = "black"),
        guide = "legend"
      )
  } else {
    p <- ggplot2::ggplot(db_plot %>% dplyr::filter(freq > 0) %>% dplyr::arrange(freq), ggplot2::aes(x = VL, y = VH)) +
      ggplot2::geom_point(ggplot2::aes(size = freq, fill = freq), shape = 21, stroke = 0.5)
  }
  p <- p +
    # scale_fill_gradientn(
    #  colors = c("white", "cornflowerblue", "yellow", "red"),
    #  values = scales::rescale(c(0, 0.25, 0.5, 1)),
    #  limits = c(0, max(pair_freq_complete$freq)),
    #  name = "Frequency"
    # ) +
    ggplot2::scale_fill_viridis_c(
      option = "D",
      direction = 1,
      limits = c(0, max(pair_freq_complete$freq)),
      name = "Frequency"
    ) +
    ggplot2::scale_size(range = c(2, 6), name = "Frequency") +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::facet_wrap(vars(!!rlang::sym(split.by)), ncol = ncol) +
    ggplot2::labs(
      title = "VH/VL Pairing",
      x = "VL", y = "VH"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.text.y = ggplot2::element_text(size = 6),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 6),
      legend.position = "right"
    )

  if (save_plot) {
    if (!stringr::str_ends(plots_folder, "/")) {
      plots_folder <- paste0(plots_folder, "/")
    }
    if (isFALSE(dir.exists(plots_folder))) {
      dir.create(plots_folder)
    }
    if (!is.null(prefix)) {
      filename <- paste0(plots_folder, prefix, "_Vgene-pairing.pdf")
    } else {
      filename <- paste0(plots_folder, "Vgene-pairing.pdf")
    }
    if (any(is.null(height), is.null(width))) {
      width <- 8 * ncol + 2
      height <- 7 * ceiling(length(groups) / 2)
    }
    if (save_as == "pdf") {
      pdf(file = filename, height = height, width = width)
      plot(p)
      dev.off()
    }
    if (save_as == "png") {
      png(filename = filename, height = height, width = width)
      plot(p)
      dev.off()
    }
  }

  if (return_plot) {
    return(p)
  }
}
