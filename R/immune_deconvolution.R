#' Run Immune Cell Deconvolution Analysis
#'
#' Performs immune cell deconvolution using multiple methods including mMCPcounter for mouse data.
#'
#' @param expr_data Expression data (genes x samples), preferably TPM or FPKM
#' @param metadata Sample metadata
#' @param condition_column Column name for condition comparison
#' @param species Species ("mouse" or "human", default: "mouse")
#' @param method Deconvolution method ("mMCPcounter", "MCPcounter", or "all")
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return List containing immune cell scores and analysis results
#' @export
run_immune_deconvolution <- function(expr_data, metadata, condition_column, species = "mouse",
                                   method = "mMCPcounter", output_dir, experiment_name) {

  create_output_dir(output_dir)
  message("Starting immune cell deconvolution analysis...")

  # Run deconvolution based on method
  if (method == "mMCPcounter" || method == "all") {
    if (species == "mouse") {
      mmc_results <- run_mmcpcounter_analysis(expr_data, metadata, condition_column,
                                            output_dir, experiment_name)
    } else {
      message("mMCPcounter is designed for mouse data. Consider using MCPcounter for human data.")
      mmc_results <- NULL
    }
  } else {
    mmc_results <- NULL
  }

  # Add other methods here in the future (MCPcounter, xCell, etc.)

  # Create comprehensive visualizations
  if (!is.null(mmc_results)) {
    create_immune_comparison_plots(mmc_results$scores, metadata, condition_column,
                                 output_dir, experiment_name)

    # Statistical analysis
    stats_results <- perform_immune_statistical_analysis(mmc_results$scores, metadata,
                                                        condition_column, output_dir, experiment_name)
  } else {
    stats_results <- NULL
  }

  message("Immune deconvolution analysis completed!")

  return(list(
    mMCPcounter = mmc_results,
    statistical_results = stats_results
  ))
}

#' Perform Statistical Analysis of Immune Scores
#'
#' Compares immune cell abundances between conditions and returns a table of p-values.
#'
#' @param immune_scores Immune scores matrix (cell types x samples)
#' @param metadata Sample metadata
#' @param condition_column Column name for condition comparison
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return Data frame with immune cell, test used, p-value, and adjusted p-value
#' @keywords internal
perform_immune_statistical_analysis <- function(immune_scores, metadata, condition_column,
                                                output_dir, experiment_name) {

  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package 'tidyr' is required")

  # Transpose immune_scores to samples x cell types
  immune_df <- as.data.frame(t(immune_scores))
  immune_df$Sample <- gsub("^X", "", rownames(immune_df))

  # Match metadata
  sample_order <- match(immune_df$Sample, metadata$SampleID)
  metadata_ordered <- metadata[sample_order, ]
  metadata_ordered <- metadata_ordered[!is.na(metadata_ordered$SampleID), ]

  immune_df_matched <- immune_df[immune_df$Sample %in% metadata_ordered$SampleID, ]
  immune_df_matched[[condition_column]] <- metadata_ordered[[condition_column]][
    match(immune_df_matched$Sample, metadata_ordered$SampleID)
  ]

  cell_types <- colnames(immune_df_matched)[!colnames(immune_df_matched) %in% c("Sample", condition_column)]

  results <- lapply(cell_types, function(cell) {
    scores <- immune_df_matched[[cell]]
    condition <- immune_df_matched[[condition_column]]
    n_groups <- length(unique(condition))

    if (n_groups == 2) {
      test <- "wilcox.test"
      test_res <- wilcox.test(scores ~ condition)
      p_val <- test_res$p.value
    } else {
      test <- "kruskal.test"
      test_res <- kruskal.test(scores ~ condition)
      p_val <- test_res$p.value
    }

    data.frame(
      Cell_type = cell,
      Test = test,
      P_value = p_val
    )
  })

  results_df <- dplyr::bind_rows(results)
  results_df$Adj_P_value <- p.adjust(results_df$P_value, method = "BH")

  # Save results
  output_file <- file.path(output_dir, paste0(experiment_name, "_immune_stats.tsv"))
  write.table(results_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

  message("Immune statistical analysis saved to: ", output_file)

  return(results_df)
}



#' Run mMCPcounter Analysis
#'
#' Performs mouse immune cell deconvolution using mMCPcounter.
#'
#' @param expr_data Expression data (genes x samples)
#' @param metadata Sample metadata
#' @param condition_column Condition column name
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return List with mMCPcounter scores and analysis
#' @keywords internal
run_mmcpcounter_analysis <- function(expr_data, metadata, condition_column, output_dir, experiment_name) {

  if (!requireNamespace("mMCPcounter", quietly = TRUE)) {
    stop("Package 'mMCPcounter' is required but not installed")
  }

  message("Running mMCPcounter for mouse immune cell deconvolution...")

  # Prepare expression data with gene symbols as rownames
  expr_symbols <- expr_data
  rownames(expr_symbols) <- sub(".*_", "", rownames(expr_symbols))  # Extract gene symbol part

  # Run mMCPcounter - expects log2(TPM) or similar
  mmc_scores <- mMCPcounter::mMCPcounter.estimate(expr_symbols, features = "ENSEMBL.ID")

  # Convert to data frame and clean up
  mmc_df <- as.data.frame(t(mmc_scores))  # Transpose so samples are rows
  mmc_df$Sample <- rownames(mmc_df)
  mmc_df$Sample <- gsub("^X", "", mmc_df$Sample)  # Remove X prefix if present

  # Match with metadata
  metadata_matched <- match_metadata_to_expression_immune(mmc_df, metadata)

  # Save results
  output_file <- file.path(output_dir, paste0(experiment_name, "_mMCPcounter_scores.tsv"))
  write.table(mmc_scores, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)

  message("mMCPcounter analysis completed. Results saved to: ", output_file)
  message("Estimated abundances for ", nrow(mmc_scores), " immune cell types in ",
          ncol(mmc_scores), " samples")

  return(list(
    scores = mmc_scores,
    scores_df = mmc_df,
    metadata = metadata_matched
  ))
}

#' Match Metadata to Immune Scores
#'
#' Helper function to match metadata with immune deconvolution scores.
#'
#' @param immune_df Immune scores data frame
#' @param metadata Sample metadata
#' @param sample_id_col Sample ID column name
#' @return Matched metadata
#' @keywords internal
match_metadata_to_expression_immune <- function(immune_df, metadata, sample_id_col = "SampleID") {

  sample_ids <- unique(immune_df$Sample)
  metadata_matched <- metadata[match(sample_ids, metadata[[sample_id_col]]), ]

  # Remove missing matches
  valid_matches <- !is.na(metadata_matched[[sample_id_col]])
  metadata_matched <- metadata_matched[valid_matches, ]

  return(metadata_matched)
}

#' Create Immune Comparison Plots
#'
#' Creates comprehensive visualization plots for immune cell deconvolution results.
#'
#' @param immune_scores Immune cell scores matrix (cell types x samples)
#' @param metadata Sample metadata
#' @param condition_column Condition column name
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return None (creates plots)
#' @keywords internal
create_immune_comparison_plots <- function(immune_scores, metadata, condition_column,
                                         output_dir, experiment_name) {

  # Create boxplot comparison
  create_immune_boxplots(immune_scores, metadata, condition_column, output_dir, experiment_name)

  # Create heatmap
  create_immune_heatmap(immune_scores, metadata, condition_column, output_dir, experiment_name)

  # Create radar plot
  create_immune_radar_plot(immune_scores, metadata, condition_column, output_dir, experiment_name)

  # Create correlation plot
  #create_immune_correlation_plot(immune_scores, output_dir, experiment_name)
}

#' Create Immune Cell Boxplots
#'
#' Creates boxplots comparing immune cell abundances between conditions.
#'
#' @param immune_scores Immune scores matrix
#' @param metadata Sample metadata
#' @param condition_column Condition column name
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return ggplot object
#' @keywords internal
create_immune_boxplots <- function(immune_scores, metadata, condition_column, output_dir, experiment_name) {

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("ggpubr", quietly = TRUE)) {
    message("Skipping immune boxplots - required packages not available")
    return()
  }

  # Prepare data
  immune_df <- as.data.frame(t(immune_scores))
  immune_df$Sample <- gsub("^X", "", rownames(immune_df))

  # Match metadata
  sample_order <- match(immune_df$Sample, metadata$SampleID)
  metadata_ordered <- metadata[sample_order, ]
  metadata_ordered <- metadata_ordered[!is.na(metadata_ordered$SampleID), ]

  # Add condition to immune data
  immune_df_matched <- immune_df[immune_df$Sample %in% metadata_ordered$SampleID, ]
  immune_df_matched[[condition_column]] <- metadata_ordered[[condition_column]][
    match(immune_df_matched$Sample, metadata_ordered$SampleID)
  ]

  # Convert to long format
  immune_long <- tidyr::pivot_longer(
    immune_df_matched,
    cols = -c("Sample", all_of(condition_column)),
    names_to = "Cell_type",
    values_to = "Score"
  )

  # Create plot
  p <- ggplot2::ggplot(immune_long, ggplot2::aes_string(x = condition_column, y = "Score",
                                                        fill = condition_column)) +
    ggplot2::geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.2, size = 0.7, alpha = 0.6) +
    ggplot2::facet_wrap(~ Cell_type, scales = "free_y", ncol = 3) +
    ggplot2::scale_fill_brewer(type = "qual", palette = "Set1") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::labs(
      title = "Immune Cell Abundances by Condition",
      x = NULL,
      y = "Immune Score",
      fill = condition_column
    )

  # Add statistical comparisons
  if (length(unique(immune_long[[condition_column]])) == 2) {
    p <- p + ggpubr::stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE
    )
  }

  # Save plot
  output_file <- file.path(output_dir, paste0(experiment_name, "_immune_boxplots.pdf"))
  ggplot2::ggsave(output_file, p, width = 12, height = 10)

  message("Immune boxplots saved to: ", output_file)
  return(p)
}

#' Create Immune Cell Heatmap
#'
#' Creates a heatmap showing immune cell abundances across samples.
#'
#' @param immune_scores Immune scores matrix
#' @param metadata Sample metadata
#' @param condition_column Condition column name
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return None (creates plot)
#' @keywords internal
create_immune_heatmap <- function(immune_scores, metadata, condition_column, output_dir, experiment_name) {

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    message("Skipping immune heatmap - ComplexHeatmap not available")
    return(invisible(NULL))
  }

  # ----------------------------
  # 1. Clean & preprocess data
  # ----------------------------

  # Remove zero-variance rows before scaling
  row_sds <- apply(immune_scores, 1, sd, na.rm = TRUE)
  immune_scores <- immune_scores[row_sds > 0, , drop = FALSE]

  if (nrow(immune_scores) == 0) {
    message("No rows with non-zero variance. Skipping heatmap.")
    return(invisible(NULL))
  }

  # Z-score scaling by row
  immune_scaled <- t(scale(t(immune_scores)))

  # Clean sample names (e.g. X prefixes from R)
  colnames(immune_scaled) <- gsub("^X", "", colnames(immune_scaled))

  # Remove rows or columns with any NA/NaN/Inf after scaling
  immune_scaled <- immune_scaled[complete.cases(immune_scaled), , drop = FALSE]
  immune_scaled <- immune_scaled[, apply(immune_scaled, 2, function(x) all(is.finite(x))), drop = FALSE]

  if (nrow(immune_scaled) == 0 || ncol(immune_scaled) == 0) {
    message("Matrix is empty after cleaning. Skipping heatmap.")
    return(invisible(NULL))
  }

  # ----------------------------
  # 2. Match metadata and samples
  # ----------------------------
  common_samples <- intersect(colnames(immune_scaled), metadata$SampleID)

  if (length(common_samples) == 0) {
    message("No matching samples between immune_scores and metadata. Skipping heatmap.")
    return(invisible(NULL))
  }

  immune_scaled <- immune_scaled[, common_samples, drop = FALSE]
  metadata_ordered <- metadata[match(common_samples, metadata$SampleID), ]

  # ----------------------------
  # 3. Column annotation (conditions)
  # ----------------------------
  unique_conditions <- unique(metadata_ordered[[condition_column]])
  n_conditions <- length(unique_conditions)
  palette <- RColorBrewer::brewer.pal(max(3, n_conditions), "Set1")[1:n_conditions]
  names(palette) <- unique_conditions

  ha <- ComplexHeatmap::HeatmapAnnotation(
    Condition = metadata_ordered[[condition_column]],
    col = list(Condition = palette)
  )

  # ----------------------------
  # 4. Create and save heatmap
  # ----------------------------
  output_file <- file.path(output_dir, paste0(experiment_name, "_immune_heatmap.pdf"))

  pdf(output_file, width = 10, height = 8)
  ht <- ComplexHeatmap::Heatmap(
    immune_scaled,
    name = "Immune Score\n(Z-score)",
    col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    top_annotation = ha,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_side = "left",
    column_title = "Immune Cell Abundances"
  )
  draw(ht)
  dev.off()

  message("Immune heatmap saved to: ", output_file)
  invisible(output_file)
}

#' Create Immune Radar Plot
#'
#' Creates radar plots showing immune cell profiles for different conditions.
#'
#' @param immune_scores Immune scores matrix
#' @param metadata Sample metadata
#' @param condition_column Condition column name
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return None (creates plot)
#' @keywords internal
create_immune_radar_plot <- function(immune_scores, metadata, condition_column, output_dir, experiment_name) {

  if (!requireNamespace("fmsb", quietly = TRUE)) {
    message("Skipping radar plot - fmsb package not available")
    return()
  }

  # Prepare data
  immune_df <- as.data.frame(t(immune_scores))
  immune_df$Sample <- gsub("^X", "", rownames(immune_df))

  # Match metadata
  sample_order <- match(immune_df$Sample, metadata$SampleID)
  metadata_ordered <- metadata[sample_order, ]
  metadata_ordered <- metadata_ordered[!is.na(metadata_ordered$SampleID), ]

  immune_df_matched <- immune_df[immune_df$Sample %in% metadata_ordered$SampleID, ]
  immune_df_matched[[condition_column]] <- metadata_ordered[[condition_column]][
    match(immune_df_matched$Sample, metadata_ordered$SampleID)
  ]

  # Calculate means for each condition
  conditions <- unique(immune_df_matched[[condition_column]])
  cell_types <- colnames(immune_df_matched)[!colnames(immune_df_matched) %in% c("Sample", condition_column)]

  radar_data <- data.frame(row.names = cell_types)

  for (condition in conditions) {
    condition_data <- immune_df_matched[immune_df_matched[[condition_column]] == condition, cell_types, drop = FALSE]
    radar_data[[condition]] <- colMeans(condition_data, na.rm = TRUE)
  }

  # fmsb requires max and min rows
  radar_data <- rbind(
    max = apply(radar_data, 1, max, na.rm = TRUE),
    min = apply(radar_data, 1, min, na.rm = TRUE),
    radar_data
  )


  radar_data_t <- t(radar_data[-c(1,2), ])  # remove max/min rows, transpose
  colnames(radar_data_t) <- rownames(radar_data)[-c(1,2)]  # set column names to cell types

  # Add max and min rows as first two rows
  radar_data_ready <- rbind(
    max = apply(radar_data_t, 2, max, na.rm = TRUE),
    min = apply(radar_data_t, 2, min, na.rm = TRUE),
    radar_data_t
  )
  radar_data_ready<-as.data.frame(radar_data_ready)
  # Plot
  output_file <- file.path(output_dir, paste0(experiment_name, "_immune_radar_plot.pdf"))
  pdf(output_file, width = 8, height = 8)


  # Plot
  fmsb::radarchart(
    radar_data_ready,
    axistype = 1,
    pcol = RColorBrewer::brewer.pal(nrow(radar_data_ready) - 2, "Set1"),
    pfcol = scales::alpha(RColorBrewer::brewer.pal(nrow(radar_data_ready) - 2, "Set1"), 0.4),
    plwd = 2,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "grey30",
    caxislabels = seq(0, max(radar_data_ready[1, ]), length.out = 5),
    cglwd = 0.8,
    vlcex = 0.8
  )

  legend("topright", legend = rownames(radar_data_ready)[-c(1,2)],
         col = RColorBrewer::brewer.pal(nrow(radar_data_ready) - 2, "Set1"),
         lty = 1, lwd = 2, bty = "n")
  dev.off()

  message("Immune radar plot saved to: ", output_file)
}
