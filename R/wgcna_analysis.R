#' Run WGCNA Co-expression Network Analysis
#'
#' Performs Weighted Gene Co-expression Network Analysis (WGCNA) to identify gene modules.
#'
#' @param expr_data Expression data matrix (genes x samples), log-transformed
#' @param metadata Sample metadata data frame
#' @param trait_column Column name in metadata for trait of interest
#' @param n_genes Number of most variable genes to use (default: 3000)
#' @param min_module_size Minimum module size (default: 30)
#' @param output_dir Output directory for plots and results
#' @param experiment_name Name for the experiment
#' @param n_threads Number of threads for WGCNA (default: 2)
#' @return List containing module assignments, eigengenes, and trait correlations
#' @importFrom dynamicTreeCut cutreeDynamic
#' @export

run_wgcna_analysis <- function(expr_data, metadata, trait_column, n_genes = 3000,
                              min_module_size = 30, output_dir, experiment_name, n_threads = 2) {

  # Check required packages
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required but not installed")
  }

  create_output_dir(output_dir)

  # Enable WGCNA threads
  WGCNA::enableWGCNAThreads(nThreads = n_threads)
  options(stringsAsFactors = FALSE)

  message("Starting WGCNA analysis with ", n_genes, " most variable genes")

  # Select most variable genes
  gene_vars <- apply(expr_data, 1, var, na.rm = TRUE)
  n_genes_actual <- min(n_genes, sum(!is.na(gene_vars)))
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_genes_actual]

  # Prepare data (samples x genes for WGCNA)
  datExpr <- t(expr_data[top_genes, ])

  # Check data quality
  gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    message("Removing genes/samples with too many missing values")
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  }

  # Choose soft-thresholding power
  message("Determining soft-thresholding power...")
  sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = seq(1, 15, by = 1), verbose = 2)
  softPower <- ifelse(!is.na(sft$powerEstimate), sft$powerEstimate, 6)
  message("Using soft power = ", softPower)

  # Plot soft threshold selection
  plot_soft_threshold_selection(sft, output_dir, experiment_name)

  # Build network and identify modules
  message("Building adjacency matrix and detecting modules...")
  adj <- WGCNA::adjacency(datExpr, power = softPower)
  TOM <- WGCNA::TOMsimilarity(adj)
  dissTOM <- 1 - TOM

  # Hierarchical clustering
  geneTree <- hclust(as.dist(dissTOM), method = "average")

  # Module identification
  dynamicMods <- dynamicTreeCut::cutreeDynamic(
    dendro = geneTree,
    distM = dissTOM,
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    minClusterSize = min_module_size
  )

  dynamicColors <- WGCNA::labels2colors(dynamicMods)
  names(dynamicColors) <- colnames(datExpr)
  # Calculate module eigengenes
  MEList <- WGCNA::moduleEigengenes(datExpr, colors = dynamicColors)
  MEs <- MEList$eigengenes

  # Order modules by clustering
  MEDiss <- 1 - cor(MEs)
  METree <- hclust(as.dist(MEDiss), method = "average")
  MEs <- MEs[, METree$order]

  message("Identified ", length(unique(dynamicColors)), " modules")

  # Create module-trait associations
  trait_data <- create_trait_matrix(metadata, trait_column)
  module_trait_results <- analyze_module_trait_relationships(MEs, trait_data, output_dir, experiment_name)

  # Create visualizations
  plot_dendrogram_and_modules <- function(geneTree, moduleColors, output_dir, experiment_name) {

  output_file <- file.path(output_dir, paste0(experiment_name, "_gene_dendrogram_modules.pdf"))

  pdf(output_file, width = 12, height = 6)
  WGCNA::plotDendroAndColors(
    geneTree, moduleColors, "Module colors",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05,
    main = "Gene dendrogram and module colors"
  )
  dev.off()

  message("Dendrogram plot saved to: ", output_file)
}

#' Plot Module-Trait Heatmap
#'
#' Creates a heatmap showing correlations between modules and traits.
#'
#' @param moduleTraitCor Module-trait correlation matrix
#' @param moduleTraitP Module-trait p-value matrix
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return None (creates plot)
#' @keywords internal
plot_module_trait_heatmap <- function(moduleTraitCor, moduleTraitP, output_dir, experiment_name) {

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    # Fallback to base R heatmap
    output_file <- file.path(output_dir, paste0(experiment_name, "_module_trait_heatmap.pdf"))
    pdf(output_file, width = 8, height = 6)

    # Create text matrix for significance
    textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitP, 1), ")", sep = "")
    dim(textMatrix) <- dim(moduleTraitCor)

    # Plot heatmap
    par(mar = c(6, 8.5, 3, 3))
    WGCNA::labeledHeatmap(
      Matrix = moduleTraitCor,
      xLabels = colnames(moduleTraitCor),
      yLabels = rownames(moduleTraitCor),
      ySymbols = rownames(moduleTraitCor),
      colorLabels = FALSE,
      colors = WGCNA::blueWhiteRed(50),
      textMatrix = textMatrix,
      setStdMargins = FALSE,
      cex.text = 0.8,
      zlim = c(-1, 1),
      main = paste("Module-trait relationships")
    )
    dev.off()
  } else {
    # Use ComplexHeatmap for better visualization
    output_file <- file.path(output_dir, paste0(experiment_name, "_module_trait_heatmap.pdf"))

    # Create significance annotation
    sig_matrix <- ifelse(moduleTraitP < 0.05, "*", "")
    sig_matrix[moduleTraitP < 0.01] <- "**"
    sig_matrix[moduleTraitP < 0.001] <- "***"

    pdf(output_file, width = 8, height = 6)
    ht <- ComplexHeatmap::Heatmap(
      moduleTraitCor,
      name = "Correlation",
      col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
      cell_fun = function(j, i, x, y, w, h, col) {
        grid::grid.text(sig_matrix[i, j], x, y, gp = grid::gpar(fontsize = 12))
      },
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      column_title = "Module-Trait Correlations",
      heatmap_legend_param = list(title = "Correlation")
    )
    print(ht)
    dev.off()
  }

  message("Module-trait heatmap saved to: ", output_file)
}

#' Plot Module Eigengenes
#'
#' Creates plots showing module eigengene expression patterns.
#'
#' @param MEs Module eigengenes
#' @param trait_data Trait data
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return None (creates plot)
#' @keywords internal
plot_module_eigengenes <- function(MEs, trait_data, output_dir, experiment_name) {

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE)) {
    message("Skipping eigengene plots - required packages not available")
    return()
  }

  # Prepare data for plotting
  ME_data <- data.frame(MEs)
  ME_data$sample <- rownames(ME_data)

  # Add trait information
  if (ncol(trait_data) == 1) {
    ME_data$trait <- trait_data[, 1]
    trait_name <- colnames(trait_data)[1]
  } else {
    # Use first trait if multiple
    ME_data$trait <- trait_data[, 1]
    trait_name <- colnames(trait_data)[1]
  }

  # Convert to long format
  ME_long <- tidyr::pivot_longer(ME_data, cols = starts_with("ME"),
                                names_to = "Module", values_to = "Eigengene")
  ME_long$Module <- gsub("ME", "", ME_long$Module)

  # Create boxplot
  p1 <- ggplot2::ggplot(ME_long, ggplot2::aes(x = factor(trait), y = Eigengene, fill = factor(trait))) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.6) +
    ggplot2::facet_wrap(~Module, scales = "free_y", ncol = 4) +
    ggplot2::labs(
      title = "Module Eigengenes by Trait",
      x = trait_name,
      y = "Module Eigengene",
      fill = trait_name
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  # Save plot
  output_file <- file.path(output_dir, paste0(experiment_name, "_module_eigengenes.pdf"))
  ggplot2::ggsave(output_file, p1, width = 12, height = 8)

  message("Module eigengene plot saved to: ", output_file)
}

#' Analyze Gene Significance and Module Membership
#'
#' Calculates gene significance and module membership for trait association.
#'
#' @param datExpr Expression data
#' @param trait Trait vector
#' @param moduleColors Module assignments
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return Data frame with gene significance and module membership
#' @keywords internal
analyze_gene_significance <- function(datExpr, trait, moduleColors, output_dir, experiment_name) {

  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required")
  }

  # Calculate gene significance
  geneTraitSignificance <- as.numeric(cor(datExpr, trait, use = "p"))
  GSPvalue <- as.numeric(WGCNA::corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr)))

  # Calculate module membership
  modNames <- substring(names(table(moduleColors)), 3)
  geneModuleMembership <- as.data.frame(cor(datExpr, WGCNA::moduleEigengenes(datExpr, moduleColors)$eigengenes, use = "p"))
  MMPvalue <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))

  names(geneModuleMembership) <- paste("MM", modNames, sep = "")
  names(MMPvalue) <- paste("p.MM", modNames, sep = "")

  # Create results data frame
  gene_results <- data.frame(
    gene = colnames(datExpr),
    module = moduleColors,
    geneTraitSignificance = geneTraitSignificance,
    GSPvalue = GSPvalue,
    geneModuleMembership,
    MMPvalue
  )

  # Save results
  output_file <- file.path(output_dir, paste0(experiment_name, "_gene_significance_MM.tsv"))
  write.table(gene_results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

  # Create scatter plot for most correlated module
  trait_cor <- cor(WGCNA::moduleEigengenes(datExpr, moduleColors)$eigengenes, trait, use = "p")
  best_module <- names(which.max(abs(trait_cor)))
  best_module_clean <- gsub("ME", "", best_module)

  plot_gene_significance_vs_mm(gene_results, best_module_clean, output_dir, experiment_name)

  message("Gene significance analysis completed")
  return(gene_results)
}

#' Plot Gene Significance vs Module Membership
#'
#' Creates scatter plot of gene significance vs module membership.
#'
#' @param gene_results Gene significance results
#' @param module_name Module name to plot
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return None (creates plot)
#' @keywords internal
plot_gene_significance_vs_mm <- function(gene_results, module_name, output_dir, experiment_name) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Skipping GS vs MM plot - ggplot2 not available")
    return()
  }

  # Filter genes in the module
  module_genes <- gene_results[gene_results$module == module_name, ]
  mm_col <- paste0("MM", module_name)

  if (!mm_col %in% colnames(gene_results)) {
    message("Module membership column not found for module: ", module_name)
    return()
  }

  # Create plot
  p <- ggplot2::ggplot(module_genes, ggplot2::aes_string(x = mm_col, y = "geneTraitSignificance")) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "red") +
    ggplot2::labs(
      title = paste("Gene Significance vs Module Membership -", module_name, "Module"),
      x = paste("Module Membership in", module_name, "module"),
      y = "Gene significance for trait"
    ) +
    ggplot2::theme_minimal()

  # Save plot
  output_file <- file.path(output_dir, paste0(experiment_name, "_GS_vs_MM_", module_name, ".pdf"))
  ggplot2::ggsave(output_file, p, width = 8, height = 6)

  message("GS vs MM plot saved to: ", output_file)
}

#' Save WGCNA Results
#'
#' Saves all WGCNA results to files.
#'
#' @param moduleColors Module color assignments
#' @param MEs Module eigengenes
#' @param module_trait_results Module-trait correlation results
#' @param gene_module_analysis Gene significance analysis results
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return None (invisible)
#' @keywords internal
save_wgcna_results <- function(moduleColors, MEs, module_trait_results, gene_module_analysis,
                              output_dir, experiment_name) {

  # Save module assignments
  module_file <- file.path(output_dir, paste0(experiment_name, "_module_assignments.tsv"))

  module_df <- data.frame(
    gene = colnames(datExpr),
    module = moduleColors
  )

  # Convert moduleTraitCor and moduleTraitP to data frames
  moduleCor <- as.data.frame(module_trait_results$moduleTraitCor)
  moduleP   <- as.data.frame(module_trait_results$moduleTraitP)

  # Add module names as a column (remove "ME" prefix to match colors if needed)
  moduleCor$module <- gsub("^ME", "", rownames(moduleCor))
  moduleP$module   <- gsub("^ME", "", rownames(moduleP))


  # Merge module correlation info with each gene
  geneModuleTraitDF <- merge(
    module_df,
    moduleCor,
    by.x = "module",
    by.y = "module",
    all.x = TRUE,
    sort = FALSE
  )

  # Optionally add p-values
  geneModuleTraitDF <- merge(
    geneModuleTraitDF,
    moduleP,
    by.x = "module",
    by.y = "module",
    suffixes = c("_cor", "_p"),
    all.x = TRUE,
    sort = FALSE
  )

  # Save
  write.table(
    geneModuleTraitDF,
    file = file.path(output_dir, paste0(experiment_name, "_genes_with_module_trait.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  # Save module eigengenes
  me_file <- file.path(output_dir, paste0(experiment_name, "_module_eigengenes.tsv"))
  write.table(MEs, file = me_file, sep = "\t", quote = FALSE, row.names = TRUE)

  # Module summary
  module_summary <- data.frame(
    module = names(table(moduleColors)),
    size = as.numeric(table(moduleColors))
  )
  summary_file <- file.path(output_dir, paste0(experiment_name, "_module_summary.tsv"))
  write.table(module_summary, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

  message("WGCNA results saved to multiple files in: ", output_dir)

  invisible()
}

#' Get Hub Genes from WGCNA Modules
#'
#' Identifies hub genes (highly connected genes) within each module.
#'
#' @param datExpr Expression data used in WGCNA
#' @param moduleColors Module assignments
#' @param n_hub_genes Number of hub genes to return per module (default: 10)
#' @return Data frame with hub genes for each module
#' @export
get_hub_genes <- function(datExpr, moduleColors, n_hub_genes = 10) {

  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required")
  }

  # Calculate module membership for each gene
  MEs <- WGCNA::moduleEigengenes(datExpr, moduleColors)$eigengenes
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))

  # Get hub genes for each module
  modules <- unique(moduleColors)
  hub_genes_list <- list()

  for (module in modules) {
    if (module == "grey") next  # Skip unassigned genes

    # Get genes in this module
    module_genes <- names(moduleColors)[moduleColors == module]

    # Get module membership values
    me_col <- paste0("ME", module)
    if (me_col %in% colnames(geneModuleMembership)) {
      mm_values <- geneModuleMembership[module_genes, me_col]

      # Get top hub genes
      top_genes <- names(sort(abs(mm_values), decreasing = TRUE))[1:min(n_hub_genes, length(mm_values))]

      hub_genes_list[[module]] <- data.frame(
        gene = top_genes,
        module = module,
        module_membership = mm_values[top_genes],
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine results
  hub_genes_df <- do.call(rbind, hub_genes_list)
  rownames(hub_genes_df) <- NULL

  return(hub_genes_df)
}
  plot_dendrogram_and_modules(geneTree, dynamicColors, output_dir, experiment_name)
  plot_module_trait_heatmap(module_trait_results$moduleTraitCor, module_trait_results$moduleTraitP,
                           output_dir, experiment_name)
  plot_module_eigengenes(MEs, trait_data, output_dir, experiment_name)

  # Gene significance and module membership analysis
  if (ncol(trait_data) == 1) {
    gene_module_analysis <- analyze_gene_significance(datExpr, trait_data[,1], dynamicColors,
                                                    output_dir, experiment_name)
  } else {
    gene_module_analysis <- NULL
    message("Skipping gene significance analysis (multiple traits detected)")
  }

  # Save results
  save_wgcna_results(dynamicColors, MEs, module_trait_results, gene_module_analysis,
                    output_dir, experiment_name)

  # Clean up
  rm(adj, TOM, dissTOM)
  gc()

  results <- list(
    moduleColors = dynamicColors,
    moduleEigengenes = MEs,
    moduleTraitCor = module_trait_results$moduleTraitCor,
    moduleTraitP = module_trait_results$moduleTraitP,
    geneTree = geneTree,
    softPower = softPower,
    datExpr = datExpr,
    trait_data = trait_data
  )

  if (!is.null(gene_module_analysis)) {
    results$geneSignificance <- gene_module_analysis
  }

  message("WGCNA analysis completed!")
  return(results)
}

#' Plot Soft Threshold Selection
#'
#' Creates plots to help select the soft-thresholding power for WGCNA.
#'
#' @param sft Soft threshold analysis results from pickSoftThreshold
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return None (creates plots)
#' @keywords internal
plot_soft_threshold_selection <- function(sft, output_dir, experiment_name) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  # Create data frame for plotting
  plot_data <- data.frame(
    Power = sft$fitIndices$Power,
    SFT_R_sq = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
    mean_k = sft$fitIndices$mean.k.
  )

  # Scale-free topology fit
  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Power, y = SFT_R_sq)) +
    ggplot2::geom_point() +
    ggplot2::geom_text(ggplot2::aes(label = Power), vjust = -0.5, size = 3) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = "Scale-free Topology Model Fit",
      x = "Soft Threshold (power)",
      y = "Scale Free Topology Model Fit, signed R-squared"
    ) +
    ggplot2::theme_minimal()

  # Mean connectivity
  p2 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Power, y = mean_k)) +
    ggplot2::geom_point() +
    ggplot2::geom_text(ggplot2::aes(label = Power), vjust = -0.5, size = 3) +
    ggplot2::labs(
      title = "Mean Connectivity",
      x = "Soft Threshold (power)",
      y = "Mean Connectivity"
    ) +
    ggplot2::theme_minimal()

  # Save plots
  output_file <- file.path(output_dir, paste0(experiment_name, "_soft_threshold_selection.pdf"))
  pdf(output_file, width = 12, height = 5)
  print(gridExtra::grid.arrange(p1, p2, ncol = 2))
  dev.off()

  message("Soft threshold selection plots saved to: ", output_file)
}

#' Create Trait Matrix
#'
#' Converts trait information to numeric matrix for WGCNA analysis.
#'
#' @param metadata Sample metadata
#' @param trait_column Column name containing trait information
#' @return Numeric matrix of traits
#' @keywords internal
create_trait_matrix <- function(metadata, trait_column) {

  trait_values <- metadata[[trait_column]]

  if (is.factor(trait_values) || is.character(trait_values)) {
    # Convert categorical to numeric (0/1 for each level)
    trait_levels <- unique(trait_values)
    trait_matrix <- matrix(0, nrow = length(trait_values), ncol = length(trait_levels))
    colnames(trait_matrix) <- paste0(trait_column, "_", trait_levels)

    for (i in seq_along(trait_levels)) {
      trait_matrix[trait_values == trait_levels[i], i] <- 1
    }
  } else {
    # Numeric trait
    trait_matrix <- matrix(as.numeric(trait_values), ncol = 1)
    colnames(trait_matrix) <- trait_column
  }

  rownames(trait_matrix) <- rownames(metadata)
  return(trait_matrix)
}

#' Analyze Module-Trait Relationships
#'
#' Calculates correlations between module eigengenes and traits.
#'
#' @param MEs Module eigengenes
#' @param trait_data Trait data matrix
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return List with correlation matrix and p-values
#' @keywords internal
analyze_module_trait_relationships <- function(MEs, trait_data, output_dir, experiment_name) {

  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required")
  }

  # Calculate correlations
  moduleTraitCor <- cor(MEs, trait_data, use = "pairwise.complete.obs")
  moduleTraitP <- WGCNA::corPvalueStudent(moduleTraitCor, nrow(MEs))

  # Save correlation results
  cor_file <- file.path(output_dir, paste0(experiment_name, "_module_trait_correlations.tsv"))
  cor_results <- data.frame(
    module = rownames(moduleTraitCor),
    moduleTraitCor,
    moduleTraitP,
    check.names = FALSE
  )
  write.table(cor_results, file = cor_file, sep = "\t", quote = FALSE, row.names = FALSE)

  return(list(moduleTraitCor = moduleTraitCor, moduleTraitP = moduleTraitP))
}

#' Perform Functional Enrichment Analysis for Significant WGCNA Modules
#'
#' Identifies significant modules and performs GO/KEGG/Reactome enrichment analysis
#' and creates STRING network visualizations.
#'
#' @param wgcna_results Results from run_wgcna_analysis()
#' @param significance_threshold P-value threshold for significant modules (default: 0.05)
#' @param min_correlation_threshold Minimum absolute correlation for significance (default: 0.3)
#' @param organism Organism for enrichment analysis ('human', 'mouse', etc.) (default: 'human')
#' @param output_dir Output directory for results
#' @param experiment_name Experiment name prefix
#' @param species_id STRING species ID (default: 9606 for human)
#' @return List containing enrichment results and network data
#' @export
analyze_significant_modules <- function(wgcna_results,
                                       significance_threshold = 0.05,
                                       min_correlation_threshold = 0.3,
                                       organism = "human",
                                       output_dir,
                                       experiment_name,
                                       species_id = 9606) {

  # Check required packages
  required_packages <- c("clusterProfiler", "org.Hs.eg.db", "ReactomePA",
                        "STRINGdb", "ggplot2", "dplyr", "enrichplot")

  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(paste("Required packages not installed:", paste(missing_packages, collapse = ", ")))
  }

  # Create enrichment output directory
  enrichment_dir <- file.path(output_dir, "enrichment_analysis")
  create_output_dir(enrichment_dir)

  message("Starting functional enrichment analysis for significant modules...")

  # Identify significant modules
  significant_modules <- identify_significant_modules(
    wgcna_results,
    significance_threshold,
    min_correlation_threshold
  )

  if (length(significant_modules) == 0) {
    message("No significant modules found with current thresholds")
    return(NULL)
  }

  message("Found ", length(significant_modules), " significant modules: ",
          paste(significant_modules, collapse = ", "))

  # Perform enrichment analysis for each significant module
  enrichment_results <- list()

  for (module in significant_modules) {
    message("Analyzing module: ", module)

    # Get genes in the module
    module_genes <- get_module_genes(wgcna_results, module)

    if (length(module_genes) < 5) {
      message("Skipping module ", module, " - too few genes (", length(module_genes), ")")
      next
    }

    # Convert gene symbols to Entrez IDs if needed
    entrez_genes <- convert_gene_symbols_to_entrez(module_genes, organism)

    if (length(entrez_genes) < 3) {
      message("Skipping module ", module, " - too few genes after ID conversion")
      next
    }

    # Perform enrichment analyses
    module_enrichment <- perform_module_enrichment(
      entrez_genes, module_genes, module, organism, enrichment_dir, experiment_name
    )

    # Create STRING network
    string_network <- create_string_network(
      module_genes, module, species_id, enrichment_dir, experiment_name
    )

    # Store results
    enrichment_results[[module]] <- list(
      genes = module_genes,
      entrez_ids = entrez_genes,
      enrichment = module_enrichment,
      network = string_network
    )
  }

  # Create summary report
  create_enrichment_summary_report(enrichment_results, significant_modules,
                                  enrichment_dir, experiment_name)

  message("Functional enrichment analysis completed!")
  return(enrichment_results)
}

#' Identify Significant Modules
#'
#' Identifies modules that are significantly correlated with traits.
#'
#' @param wgcna_results WGCNA analysis results
#' @param significance_threshold P-value threshold
#' @param min_correlation_threshold Minimum correlation threshold
#' @return Vector of significant module names
#' @keywords internal
identify_significant_modules <- function(wgcna_results, significance_threshold, min_correlation_threshold) {

  moduleTraitCor <- wgcna_results$moduleTraitCor
  moduleTraitP <- wgcna_results$moduleTraitP

  # Find modules with significant correlations
  significant_mask <- (moduleTraitP < significance_threshold) &
                     (abs(moduleTraitCor) > min_correlation_threshold)

  # Get module names (remove ME prefix)
  significant_modules <- c()
  for (i in 1:nrow(significant_mask)) {
    if (any(significant_mask[i, ])) {
      module_name <- gsub("^ME", "", rownames(moduleTraitCor)[i])
      if (module_name != "grey") {  # Exclude unassigned genes
        significant_modules <- c(significant_modules, module_name)
      }
    }
  }

  return(unique(significant_modules))
}

#' Get Genes in a Module
#'
#' Extracts gene names for a specific module.
#'
#' @param wgcna_results WGCNA analysis results
#' @param module Module name
#' @return Vector of gene names
#' @keywords internal
get_module_genes <- function(wgcna_results, module) {
  moduleColors <- wgcna_results$moduleColors
  module_genes <- names(moduleColors)[moduleColors == module]
  module_genes <- sub("_.*", "", module_genes)
  return(module_genes)
}

#' Convert Gene Symbols to Entrez IDs
#'
#' Converts gene symbols to Entrez IDs for enrichment analysis.
#'
#' @param gene_symbols Vector of gene symbols
#' @param organism Organism ('human' or 'mouse')
#' @return Vector of Entrez IDs
#' @keywords internal
convert_gene_symbols_to_entrez <- function(gene_symbols, organism = "human") {

  if (organism == "human") {
    org_db <- org.Hs.eg.db::org.Hs.eg.db
  } else if (organism == "mouse") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("org.Mm.eg.db package required for mouse analysis")
    }
    org_db <- org.Mm.eg.db::org.Mm.eg.db
  } else {
    stop("Organism not supported. Use 'human' or 'mouse'")
  }

  # Convert symbols to Entrez IDs
  entrez_ids <- clusterProfiler::bitr(gene_symbols,
                                     fromType = "SYMBOL",
                                     toType = "ENTREZID",
                                     OrgDb = org_db)

  return(entrez_ids$ENTREZID)
}

#' Perform Module Enrichment Analysis
#'
#' Performs GO, KEGG, and Reactome enrichment analysis for a module.
#'
#' @param entrez_genes Entrez IDs
#' @param gene_symbols Gene symbols
#' @param module Module name
#' @param organism Organism
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return List of enrichment results
#' @keywords internal
perform_module_enrichment <- function(entrez_genes, gene_symbols, module, organism,
                                     output_dir, experiment_name) {

  if (organism == "human") {
    org_db <- org.Hs.eg.db::org.Hs.eg.db
    kegg_organism <- "hsa"
  } else if (organism == "mouse") {
    org_db <- org.Mm.eg.db::org.Mm.eg.db
    kegg_organism <- "mmu"
  }

  enrichment_results <- list()

  # GO Biological Process
  message("  - GO Biological Process enrichment...")
  tryCatch({
    go_bp <- clusterProfiler::enrichGO(gene = entrez_genes,
                                      OrgDb = org_db,
                                      ont = "BP",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.05,
                                      readable = TRUE)
    enrichment_results$GO_BP <- go_bp

    if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
      # Save results
      write.csv(go_bp@result,
               file.path(output_dir, paste0(experiment_name, "_", module, "_GO_BP.csv")),
               row.names = FALSE)

      # Create plot
      if (nrow(go_bp@result) >= 5) {
        p <- enrichplot::dotplot(go_bp, showCategory = 15) +
          ggplot2::ggtitle(paste("GO BP Enrichment -", module, "Module"))
        ggplot2::ggsave(file.path(output_dir, paste0(experiment_name, "_", module, "_GO_BP.pdf")),
                       p, width = 10, height = 8)
      }
    }
  }, error = function(e) {
    message("    GO BP analysis failed: ", e$message)
  })

  # GO Molecular Function
  message("  - GO Molecular Function enrichment...")
  tryCatch({
    go_mf <- clusterProfiler::enrichGO(gene = entrez_genes,
                                      OrgDb = org_db,
                                      ont = "MF",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.05,
                                      readable = TRUE)
    enrichment_results$GO_MF <- go_mf

    if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
      write.csv(go_mf@result,
               file.path(output_dir, paste0(experiment_name, "_", module, "_GO_MF.csv")),
               row.names = FALSE)

      if (nrow(go_mf@result) >= 5) {
        p <- enrichplot::dotplot(go_mf, showCategory = 15) +
          ggplot2::ggtitle(paste("GO MF Enrichment -", module, "Module"))
        ggplot2::ggsave(file.path(output_dir, paste0(experiment_name, "_", module, "_GO_MF.pdf")),
                       p, width = 10, height = 8)
      }
    }
  }, error = function(e) {
    message("    GO MF analysis failed: ", e$message)
  })

  # KEGG Pathway
  message("  - KEGG pathway enrichment...")
  tryCatch({
    kegg <- clusterProfiler::enrichKEGG(gene = entrez_genes,
                                       organism = kegg_organism,
                                       pvalueCutoff = 0.05,
                                       pAdjustMethod = "BH")

    if (!is.null(kegg) && nrow(kegg@result) > 0) {
      # Convert to readable format
      kegg <- clusterProfiler::setReadable(kegg, org_db, keyType = "ENTREZID")
      enrichment_results$KEGG <- kegg

      write.csv(kegg@result,
               file.path(output_dir, paste0(experiment_name, "_", module, "_KEGG.csv")),
               row.names = FALSE)

      if (nrow(kegg@result) >= 5) {
        p <- enrichplot::dotplot(kegg, showCategory = 15) +
          ggplot2::ggtitle(paste("KEGG Enrichment -", module, "Module"))
        ggplot2::ggsave(file.path(output_dir, paste0(experiment_name, "_", module, "_KEGG.pdf")),
                       p, width = 10, height = 8)
      }
    }
  }, error = function(e) {
    message("    KEGG analysis failed: ", e$message)
  })

  # Reactome Pathway (human only)
  if (organism == "human") {
    message("  - Reactome pathway enrichment...")
    tryCatch({
      reactome <- ReactomePA::enrichPathway(gene = entrez_genes,
                                          pvalueCutoff = 0.05,
                                          pAdjustMethod = "BH",
                                          readable = TRUE)

      if (!is.null(reactome) && nrow(reactome@result) > 0) {
        enrichment_results$Reactome <- reactome

        write.csv(reactome@result,
                 file.path(output_dir, paste0(experiment_name, "_", module, "_Reactome.csv")),
                 row.names = FALSE)

        if (nrow(reactome@result) >= 5) {
          p <- enrichplot::dotplot(reactome, showCategory = 15) +
            ggplot2::ggtitle(paste("Reactome Enrichment -", module, "Module"))
          ggplot2::ggsave(file.path(output_dir, paste0(experiment_name, "_", module, "_Reactome.pdf")),
                         p, width = 10, height = 8)
        }
      }
    }, error = function(e) {
      message("    Reactome analysis failed: ", e$message)
    })
  }

  return(enrichment_results)
}

#' Create STRING Network
#'
#' Creates protein-protein interaction network using STRINGdb.
#'
#' @param gene_symbols Gene symbols
#' @param module Module name
#' @param species_id STRING species ID
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return STRING network data
#' @keywords internal
create_string_network <- function(gene_symbols, module, species_id, output_dir, experiment_name) {

  message("  - Creating STRING network...")

  tryCatch({
    # Initialize STRING database
    string_db <- STRINGdb::STRINGdb$new(version = "11.5", species = species_id,
                                       score_threshold = 400, input_directory = "")

    # Prepare gene data
    gene_df <- data.frame(gene = gene_symbols, stringsAsFactors = FALSE)

    # Map to STRING IDs
    mapped_genes <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)

    if (nrow(mapped_genes) < 3) {
      message("    Too few genes mapped to STRING database")
      return(NULL)
    }

    message("    Mapped ", nrow(mapped_genes), " genes to STRING database")

    # Get interactions
    interactions <- string_db$get_interactions(mapped_genes$STRING_id)

    if (nrow(interactions) == 0) {
      message("    No interactions found in STRING database")
      return(NULL)
    }

    # Create network plot
    output_file <- file.path(output_dir, paste0(experiment_name, "_", module, "_STRING_network.pdf"))

    pdf(output_file, width = 12, height = 10)
    string_db$plot_network(mapped_genes$STRING_id,
                          payload_id = mapped_genes$gene,
                          add_link = FALSE,
                          add_summary = TRUE)
    title(main = paste("STRING Network -", module, "Module"),
          sub = paste("Interactions:", nrow(interactions), "genes"))
    dev.off()

    # Save network data
    network_file <- file.path(output_dir, paste0(experiment_name, "_", module, "_STRING_interactions.csv"))

    # Add gene symbols to interactions
    interactions_with_symbols <- merge(interactions,
                                     mapped_genes[, c("STRING_id", "gene")],
                                     by.x = "from", by.y = "STRING_id")
    names(interactions_with_symbols)[names(interactions_with_symbols) == "gene"] <- "from_gene"

    interactions_with_symbols <- merge(interactions_with_symbols,
                                     mapped_genes[, c("STRING_id", "gene")],
                                     by.x = "to", by.y = "STRING_id")
    names(interactions_with_symbols)[names(interactions_with_symbols) == "gene"] <- "to_gene"

    write.csv(interactions_with_symbols, network_file, row.names = FALSE)

    message("    STRING network saved to: ", output_file)

    return(list(
      mapped_genes = mapped_genes,
      interactions = interactions_with_symbols,
      network_file = network_file
    ))

  }, error = function(e) {
    message("    STRING network analysis failed: ", e$message)
    return(NULL)
  })
}

#' Create Enrichment Summary Report
#'
#' Creates a summary report of all enrichment analyses.
#'
#' @param enrichment_results Enrichment results for all modules
#' @param significant_modules Vector of significant module names
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @return None (creates files)
#' @keywords internal
create_enrichment_summary_report <- function(enrichment_results, significant_modules,
                                            output_dir, experiment_name) {

  message("Creating enrichment summary report...")

  # Create summary data frame
  summary_data <- data.frame()

  for (module in names(enrichment_results)) {
    module_data <- enrichment_results[[module]]

    # Count enriched terms for each analysis type
    go_bp_count <- if (!is.null(module_data$enrichment$GO_BP)) nrow(module_data$enrichment$GO_BP@result) else 0
    go_mf_count <- if (!is.null(module_data$enrichment$GO_MF)) nrow(module_data$enrichment$GO_MF@result) else 0
    kegg_count <- if (!is.null(module_data$enrichment$KEGG)) nrow(module_data$enrichment$KEGG@result) else 0
    reactome_count <- if (!is.null(module_data$enrichment$Reactome)) nrow(module_data$enrichment$Reactome@result) else 0

    # Network stats
    network_interactions <- if (!is.null(module_data$network)) nrow(module_data$network$interactions) else 0
    network_genes <- if (!is.null(module_data$network)) nrow(module_data$network$mapped_genes) else 0

    summary_row <- data.frame(
      Module = module,
      Total_Genes = length(module_data$genes),
      Mapped_Genes = length(module_data$entrez_ids),
      GO_BP_Terms = go_bp_count,
      GO_MF_Terms = go_mf_count,
      KEGG_Pathways = kegg_count,
      Reactome_Pathways = reactome_count,
      STRING_Interactions = network_interactions,
      STRING_Mapped_Genes = network_genes
    )

    summary_data <- rbind(summary_data, summary_row)
  }

  # Save summary
  summary_file <- file.path(output_dir, paste0(experiment_name, "_enrichment_summary.csv"))
  write.csv(summary_data, summary_file, row.names = FALSE)

  # Create detailed text report
  report_file <- file.path(output_dir, paste0(experiment_name, "_enrichment_report.txt"))

  cat("WGCNA Functional Enrichment Analysis Report\n",
      "==========================================\n\n",
      "Experiment: ", experiment_name, "\n",
      "Analysis Date: ", Sys.Date(), "\n",
      "Significant Modules: ", length(significant_modules), "\n\n",
      file = report_file)

  for (module in names(enrichment_results)) {
    module_data <- enrichment_results[[module]]

    cat("Module: ", module, "\n",
        "Number of genes: ", length(module_data$genes), "\n",
        "Genes mapped to databases: ", length(module_data$entrez_ids), "\n\n",
        file = report_file, append = TRUE)

    # Top enriched terms
    if (!is.null(module_data$enrichment$GO_BP) && nrow(module_data$enrichment$GO_BP@result) > 0) {
      top_go_bp <- head(module_data$enrichment$GO_BP@result$Description, 3)
      cat("Top GO BP terms:\n", paste("  -", top_go_bp, collapse = "\n"), "\n\n",
          file = report_file, append = TRUE)
    }

    if (!is.null(module_data$enrichment$KEGG) && nrow(module_data$enrichment$KEGG@result) > 0) {
      top_kegg <- head(module_data$enrichment$KEGG@result$Description, 3)
      cat("Top KEGG pathways:\n", paste("  -", top_kegg, collapse = "\n"), "\n\n",
          file = report_file, append = TRUE)
    }

    cat("----------------------------------------\n\n", file = report_file, append = TRUE)
  }

  message("Summary report saved to: ", summary_file)
  message("Detailed report saved to: ", report_file)
}

#' Enhanced WGCNA Analysis with Functional Enrichment
#'
#' Runs complete WGCNA analysis followed by functional enrichment analysis.
#'
#' @param expr_data Expression data matrix (genes x samples), log-transformed
#' @param metadata Sample metadata data frame
#' @param trait_column Column name in metadata for trait of interest
#' @param n_genes Number of most variable genes to use (default: 3000)
#' @param min_module_size Minimum module size (default: 30)
#' @param output_dir Output directory for plots and results
#' @param experiment_name Name for the experiment
#' @param n_threads Number of threads for WGCNA (default: 2)
#' @param run_enrichment Whether to run enrichment analysis (default: TRUE)
#' @param significance_threshold P-value threshold for significant modules (default: 0.05)
#' @param min_correlation_threshold Minimum absolute correlation for significance (default: 0.3)
#' @param organism Organism for enrichment analysis (default: 'human')
#' @param species_id STRING species ID (default: 9606 for human)
#' @return List containing WGCNA results and enrichment results
#' @export
run_enhanced_wgcna_analysis <- function(expr_data, metadata, trait_column, n_genes = 3000,
                                       min_module_size = 30, output_dir, experiment_name,
                                       n_threads = 2, run_enrichment = TRUE,
                                       significance_threshold = 0.05,
                                       min_correlation_threshold = 0.3,
                                       organism = "human", species_id = 9606) {

  # Run standard WGCNA analysis
  message("=== Running WGCNA Analysis ===")
  wgcna_results <- run_wgcna_analysis(
    expr_data = expr_data,
    metadata = metadata,
    trait_column = trait_column,
    n_genes = n_genes,
    min_module_size = min_module_size,
    output_dir = output_dir,
    experiment_name = experiment_name,
    n_threads = n_threads
  )

  enrichment_results <- NULL

  if (run_enrichment) {
    message("\n=== Running Functional Enrichment Analysis ===")
    enrichment_results <- analyze_significant_modules(
      wgcna_results = wgcna_results,
      significance_threshold = significance_threshold,
      min_correlation_threshold = min_correlation_threshold,
      organism = organism,
      output_dir = output_dir,
      experiment_name = experiment_name,
      species_id = species_id
    )
  }

  # Return combined results
  return(list(
    wgcna = wgcna_results,
    enrichment = enrichment_results
  ))
}


