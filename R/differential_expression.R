#' Run Limma-Voom Differential Expression Analysis
#'
#' Performs differential expression analysis using the limma-voom pipeline.
#'
#' @param expr_data Expression count data (genes x samples)
#' @param sample_groups Character vector of sample group labels
#' @param contrast_matrix Contrast matrix from makeContrasts
#' @param design_matrix Design matrix from model.matrix
#' @param experiment_name Name for the experiment (used in output files)
#' @param output_dir Output directory path
#' @return List containing tfit (treat results) and efit (eBayes results)
#' @export
run_limma_voom <- function(expr_data, sample_groups, contrast_matrix, design_matrix,
                           experiment_name, output_dir,
                           lfc_threshold = 1,
                           fdr_threshold = 0.05) {

  # Check required packages
  required_pkgs <- c("edgeR", "limma")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed")
    }
  }

  # Create output directory
  create_output_dir(output_dir)
  print(head(expr_data))

  # EdgeR filtering and normalization
  d0 <- edgeR::DGEList(expr_data)
  keep.exprs <- edgeR::filterByExpr(d0, group = sample_groups)
  d0 <- d0[keep.exprs, ]
  d0 <- edgeR::calcNormFactors(d0)

  # Voom transformation
  y <- limma::voom(d0, design_matrix)

  # Linear model fitting
  vfit <- limma::lmFit(y, design_matrix)
  vfit <- limma::contrasts.fit(vfit, contrasts = contrast_matrix)
  efit <- limma::eBayes(vfit)

  # Use user-defined LFC for treat
  tfit <- limma::treat(vfit, lfc = lfc_threshold)

  # Extract results (topTreat on tfit, not efit)
  results_efit <- limma::topTreat(tfit, coef = 1, n = Inf)
  results_df <- as.data.frame(results_efit)

  # Save results
  output_file <- file.path(output_dir, paste0(experiment_name, "_DE_results.tsv"))
  write.table(results_df, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)

  # Use user-defined FDR for the message
  n_sig <- sum(results_df$adj.P.Val < fdr_threshold, na.rm = TRUE)
  message("Differential expression results saved to: ", output_file)
  message("Found ", n_sig, " significant genes (adj.P.Val < ", fdr_threshold, ")")

  return(list(
    tfit = tfit,
    efit = efit,
    results = results_df,
    lfc_threshold = lfc_threshold,
    fdr_threshold = fdr_threshold
  ))
}


#' Run Complete Differential Expression Pipeline
#'
#' Runs the complete differential expression analysis pipeline for two groups.
#'
#' @param counts_data Count expression data (genes x samples)
#' @param metadata Sample metadata data frame
#' @param group1_condition Name of condition 1 (treatment group)
#' @param group2_condition Name of condition 2 (control group)
#' @param condition_column Column name in metadata containing condition information (default: "condition")
#' @param sample_id_column Column name in metadata containing sample IDs (default: "SampleID")
#' @param experiment_name Name for the experiment
#' @param output_dir Output directory
#' @param fdr_threshold FDR threshold for significance (default: 0.1)
#' @param lfc_threshold Log fold change threshold (default: 1)
#' @return List containing DE results and summary statistics
#' @export
run_differential_expression <- function(counts_data, metadata, group1_condition, group2_condition,
                                      condition_column = "condition", sample_id_column = "SampleID",
                                      experiment_name, output_dir, fdr_threshold = 0.1, lfc_threshold = 1) {

  message("Starting differential expression analysis: ", experiment_name)

  # Match metadata to expression data
  metadata_ordered <- match_metadata_to_expression(counts_data, metadata, sample_id_column)

  # Get sample names for each group
  group1_samples <- metadata_ordered[[sample_id_column]][metadata_ordered[[condition_column]] == group1_condition]
  group2_samples <- metadata_ordered[[sample_id_column]][metadata_ordered[[condition_column]] == group2_condition]

  # Clean sample names and subset data
  group1_samples_clean <- paste0( group1_samples)
  group2_samples_clean <- paste0( group2_samples)

  # Check if samples exist in data
  all_samples <- c(group1_samples_clean, group2_samples_clean)
  missing_samples <- setdiff(all_samples, colnames(counts_data))
  if (length(missing_samples) > 0) {
    warning("Missing samples in expression data: ", paste(missing_samples, collapse = ", "))
  }

  # Subset expression data
  available_samples <- intersect(all_samples, colnames(counts_data))
  expr_subset <- counts_data[, available_samples, drop = FALSE]

  # Create group labels
  sample_groups <- ifelse(colnames(expr_subset) %in% group1_samples_clean,
                         group1_condition, group2_condition)
  sample_groups <- factor(sample_groups, levels = c(group2_condition, group1_condition))

  # Create design matrix
  design <- model.matrix(~0 + sample_groups)
  colnames(design) <- levels(sample_groups)

  # Create contrast
  contrast_string <- paste0(group1_condition, " - ", group2_condition)
  contrast_matrix <- limma::makeContrasts(contrasts = contrast_string, levels = design)

  # Run limma-voom
  de_results <- run_limma_voom(
    expr_data = expr_subset,
    sample_groups = sample_groups,
    contrast_matrix = contrast_matrix,
    design_matrix = design,
    experiment_name = experiment_name,
    output_dir = output_dir
  )

  # Summarize results
  results_summary <- summarize_de_results(de_results$efit, fdr_threshold, lfc_threshold)

  # Save gene lists
  save_gene_lists(de_results$efit, experiment_name, output_dir, fdr_threshold, lfc_threshold)

  message("Differential expression analysis completed!")

  return(list(
    de_results = de_results,
    summary = results_summary,
    sample_info = list(
      group1 = group1_samples,
      group2 = group2_samples,
      group1_condition = group1_condition,
      group2_condition = group2_condition
    )
  ))
}

#' Summarize Differential Expression Results
#'
#' Creates a summary of differential expression results.
#'
#' @param efit eBayes fit object from limma
#' @param fdr_threshold FDR threshold for significance
#' @param lfc_threshold Log fold change threshold
#' @return List with summary statistics
#' @keywords internal
summarize_de_results <- function(efit, fdr_threshold = 0.1, lfc_threshold = 1) {
  # Get results
  results <- limma::decideTests(efit, p.value = fdr_threshold, lfc = lfc_threshold)

  # Count up/down regulated genes
  n_up <- sum(results == 1)
  n_down <- sum(results == -1)
  n_total <- nrow(efit)

  summary_stats <- list(
    total_genes = n_total,
    significant_genes = n_up + n_down,
    upregulated = n_up,
    downregulated = n_down,
    percent_significant = round((n_up + n_down) / n_total * 100, 2),
    fdr_threshold = fdr_threshold,
    lfc_threshold = lfc_threshold
  )

  message("Summary: ", n_up, " up, ", n_down, " down out of ", n_total, " genes")

  return(summary_stats)
}

#' Save Gene Lists
#'
#' Saves lists of up- and down-regulated genes to files.
#'
#' @param efit eBayes fit object from limma
#' @param experiment_name Name of the experiment
#' @param output_dir Output directory
#' @param fdr_threshold FDR threshold for significance
#' @param lfc_threshold Log fold change threshold
#' @return None (invisible)
#' @keywords internal
save_gene_lists <- function(efit, experiment_name, output_dir, fdr_threshold, lfc_threshold) {
  # Get significant genes
  results <- limma::decideTests(efit, p.value = fdr_threshold, lfc = lfc_threshold)

  up_genes <- rownames(results)[results == 1]
  down_genes <- rownames(results)[results == -1]

  # Save to files
  up_file <- file.path(output_dir, paste0(experiment_name, "_UP_genes.tsv"))
  down_file <- file.path(output_dir, paste0(experiment_name, "_DOWN_genes.tsv"))

  write.table(up_genes, file = up_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(down_genes, file = down_file, quote = FALSE, row.names = FALSE, col.names = FALSE)

  message("Gene lists saved: ", basename(up_file), ", ", basename(down_file))

  invisible()
}
