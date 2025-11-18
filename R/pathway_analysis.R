#' Run GSEA Analysis
#'
#' Performs Gene Set Enrichment Analysis using fgsea on differential expression results.
#'
#' @param de_results Differential expression results from limma (efit object or data frame)
#' @param species Species for gene sets ("Mus musculus" or "Homo sapiens")
#' @param category MSigDB category (default: "H" for Hallmark)
#' @param subcategory MSigDB subcategory (default: NULL)
#' @param min_size Minimum gene set size (default: 15)
#' @param max_size Maximum gene set size (default: 500)
#' @param n_perm Number of permutations (default: 100000)
#' @return fgsea results data frame
#' @export
run_gsea_analysis <- function(de_results, species = "Mus musculus", category = "H",
                             subcategory = NULL, min_size = 15, max_size = 500,
                             n_perm = 100000) {

  # Check required packages
  required_pkgs <- c("fgsea", "msigdbr", "limma")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed")
    }
  }

  # choose correct species string for msigdbr
  if (species == "MM") {
    species1 <- "Mus musculus"
    category <- "MH"
  } else {
    species1 <- "Homo sapiens"
    category <- "H"
  }


  # Get gene sets
  message("Loading gene sets for ", species, ", category: ", category)
  if (is.null(subcategory)) {
    gene_sets_df <- msigdbr::msigdbr(db_species = species, species = species1, collection = category)
  } else {
    gene_sets_df <- msigdbr::msigdbr(db_species = species, species = species1, collection =  category, subcollection = subcategory)
  }

  gene_sets <- split(gene_sets_df$ensembl_gene, gene_sets_df$gs_name)
  message("Loaded ", length(gene_sets), " gene sets")

  # Prepare rankings
  if (is.data.frame(de_results)) {
    # If already a data frame
    ranks <- de_results$logFC
    names(ranks) <- sub(".*_", "", rownames(de_results))  # Extract ENSEMBL ID
  } else {
    # If limma efit object
    efit_results <- limma::topTreat(de_results, coef = 1, n = Inf)
    ranks <- efit_results$logFC
    names(ranks) <- sub(".*_", "", rownames(efit_results))  # Extract ENSEMBL ID
  }

  # Remove NA values
  ranks <- ranks[!is.na(ranks)]
  message("Created rankings for ", length(ranks), " genes")

  # Run GSEA
  message("Running GSEA with ", n_perm, " permutations...")
  set.seed(123)
  gsea_results <- fgsea::fgsea(
    pathways = gene_sets,
    stats = ranks,
    minSize = min_size,
    maxSize = max_size,
    nPermSimple = n_perm
  )

  # Sort by absolute NES
  gsea_results <- gsea_results[order(-abs(gsea_results$NES)), ]

  message("GSEA completed. Found ", sum(gsea_results$padj < 0.05, na.rm = TRUE),
          " significant pathways (padj < 0.05)")

  return(gsea_results)
}

#' Create GSEA Barplot
#'
#' Creates a barplot visualization of top GSEA results.
#'
#' @param gsea_results Results from fgsea
#' @param n_pathways Number of top pathways to show (default: 20)
#' @param output_file Output PDF file path (optional)
#' @param width Plot width in inches (default: 8)
#' @param height Plot height in inches (default: 6)
#' @return ggplot object
#' @export
plot_gsea_barplot <- function(gsea_results, n_pathways = 20, output_file = NULL,
                             width = 8, height = 6) {

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Required packages (ggplot2, RColorBrewer) not installed")
  }

  # Select top pathways
  top_up <- gsea_results[gsea_results$NES > 0, ][1:min(n_pathways/2, sum(gsea_results$NES > 0)), ]
  top_down <- gsea_results[gsea_results$NES < 0, ][1:min(n_pathways/2, sum(gsea_results$NES < 0)), ]
  plot_data <- rbind(top_up, top_down)
  plot_data <- plot_data[!is.na(plot_data$pathway), ]

  # Clean pathway names
  plot_data$pathway_clean <- gsub("HALLMARK_", "", plot_data$pathway)
  plot_data$pathway_clean <- gsub("_", " ", plot_data$pathway_clean)

  # Create significance indicator
  plot_data$significant <- ifelse(plot_data$padj < 0.05, "Significant", "Not Significant")

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = reorder(pathway_clean, NES), y = NES,
                                               fill = significant)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = c("Significant" = "#E31A1C", "Not Significant" = "#BDBDBD")) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = "GSEA Results - Top Enriched Pathways",
      x = "Pathway",
      y = "Normalized Enrichment Score (NES)",
      fill = "Significance"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10),
      legend.position = "bottom"
    )

  # Save if requested
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = width, height = height)
    message("GSEA barplot saved to: ", output_file)
  }

  return(p)
}

#' Create GSEA Dotplot
#'
#' Creates a dotplot visualization of GSEA results showing NES, significance, and gene set size.
#'
#' @param gsea_results Results from fgsea
#' @param n_pathways Number of top pathways to show (default: 30)
#' @param output_file Output PDF file path (optional)
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @return ggplot object
#' @export
plot_gsea_dotplot <- function(gsea_results, n_pathways = 30, output_file = NULL,
                             width = 10, height = 8) {

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Required packages (ggplot2, RColorBrewer) not installed")
  }

  # Select top pathways by absolute NES
  plot_data <- gsea_results[order(-abs(gsea_results$NES)), ][1:min(n_pathways, nrow(gsea_results)), ]
  plot_data <- plot_data[!is.na(plot_data$pathway), ]

  # Clean pathway names
  plot_data$pathway_clean <- gsub("HALLMARK_", "", plot_data$pathway)
  plot_data$pathway_clean <- gsub("_", " ", plot_data$pathway_clean)

  # Order pathways by NES
  plot_data <- plot_data[order(plot_data$NES), ]
  plot_data$pathway_clean <- factor(plot_data$pathway_clean, levels = plot_data$pathway_clean)

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = NES, y = pathway_clean)) +
    ggplot2::geom_point(ggplot2::aes(size = size, fill = -log10(padj)),
                       shape = 21, stroke = 0.6, color = "black") +
    ggplot2::scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = -log10(0.05),
      name = "-log10(padj)",
      na.value = "grey80"
    ) +
    ggplot2::scale_size_continuous(name = "Gene set size", range = c(2, 8)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = "GSEA Results - Pathway Enrichment",
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      axis.text.y = ggplot2::element_text(size = 9)
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")

  # Save if requested
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = width, height = height)
    message("GSEA dotplot saved to: ", output_file)
  }

  return(p)
}

#' Create GSEA Table Plot
#'
#' Creates a table-style visualization of GSEA results using fgsea's plotGseaTable.
#'
#' @param gsea_results Results from fgsea
#' @param gene_sets Gene sets list used for GSEA
#' @param gene_ranks Named vector of gene rankings
#' @param n_pathways Number of pathways to show (default: 30)
#' @param output_file Output PDF file path (optional)
#' @param width Plot width in inches (default: 11)
#' @param height Plot height in inches (default: 15)
#' @return None (creates plot)
#' @export
create_gsea_table_plot <- function(gsea_results, gene_sets, gene_ranks, n_pathways = 30,
                                  output_file = NULL, width = 11, height = 15) {

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required but not installed")
  }

  # Select top pathways
  top_up <- gsea_results[gsea_results$NES > 0, ][order(-gsea_results$NES[gsea_results$NES > 0]), ]
  top_down <- gsea_results[gsea_results$NES < 0, ][order(gsea_results$NES[gsea_results$NES < 0]), ]

  n_up <- min(n_pathways/2, nrow(top_up))
  n_down <- min(n_pathways/2, nrow(top_down))

  top_pathways_up <- head(top_up$pathway, n_up)
  top_pathways_down <- head(top_down$pathway, n_down)
  top_pathways <- c(top_pathways_up, rev(top_pathways_down))

  # Create plot
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
  }

  p <- fgsea::plotGseaTable(
    pathways = gene_sets[top_pathways],
    stats = gene_ranks,
    fgseaRes = gsea_results,
    gseaParam = 0.5
  )

  if (!is.null(output_file)) {
    dev.off()
    message("GSEA table plot saved to: ", output_file)
  }

  return(invisible(p))
}

#' Run Camera Gene Set Testing
#'
#' Performs competitive gene set testing using limma's camera method.
#'
#' @param voom_object Voom-transformed expression object
#' @param design_matrix Design matrix
#' @param contrast_vector Contrast vector
#' @param gene_sets List of gene sets (indices)
#' @param species Species for gene sets ("Mus musculus" or "Homo sapiens")
#' @param min_size Minimum gene set size (default: 15)
#' @return Camera results data frame
#' @export
run_camera_analysis <- function(voom_object, design_matrix, contrast_vector, gene_sets = NULL,
                               species = "Mus musculus", min_size = 15) {

  if (!requireNamespace("limma", quietly = TRUE) ||
      !requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Required packages (limma, msigdbr) not installed")
  }

  if (species == "MM") {
    species1 <- "Mus musculus"
    category <- "MH"
  } else {
    species1 <- "Homo sapiens"
    category <- "H"
  }

  # If gene_sets not provided, use Hallmark sets
  if (is.null(gene_sets)) {
    message("Loading Hallmark gene sets for ", species)
    h_gene_sets <- msigdbr::msigdbr(db_species = species, species = species1, collection  = "H")
    gene_sets_list <- split(h_gene_sets$ensembl_gene, h_gene_sets$gs_name)

    # Extract ENSEMBL IDs from rownames
    rownames_voom <- rownames(voom_object$E)
    ensembl_ids <- sub(".*_(ENSMUSG[0-9]+).*", "\\1", rownames_voom)

    # Convert to indices
    gene_sets_indices <- lapply(gene_sets_list, function(gs) {
      which(!is.na(ensembl_ids) & ensembl_ids %in% gs)
    })

    # Filter by minimum size
    gene_sets_indices <- gene_sets_indices[sapply(gene_sets_indices, length) >= min_size]
  } else {
    gene_sets_indices <- gene_sets
  }

  message("Running camera analysis with ", length(gene_sets_indices), " gene sets")

  # Run camera
  camera_results <- limma::camera(
    y = voom_object$E,
    index = gene_sets_indices,
    design = design_matrix,
    contrast = contrast_vector
  )

  camera_results_df <- as.data.frame(camera_results)
  camera_results_df$pathway <- rownames(camera_results_df)

  message("Camera analysis completed. Found ", sum(camera_results_df$PValue < 0.05),
          " significant pathways (p < 0.05)")

  return(camera_results_df)
}

#' Save Pathway Analysis Results
#'
#' Saves pathway analysis results to files with proper formatting.
#'
#' @param gsea_results GSEA results data frame
#' @param camera_results Camera results data frame (optional)
#' @param experiment_name Name of the experiment
#' @param output_dir Output directory
#' @return None (invisible)
#' @export
save_pathway_results <- function(gsea_results, camera_results = NULL, experiment_name, output_dir) {

  create_output_dir(output_dir)

  # Save GSEA results
  gsea_file <- file.path(output_dir, paste0(experiment_name, "_GSEA_results.tsv"))
  gsea_results$leadingEdge <- sapply(gsea_results$leadingEdge, function(x) paste(x, collapse = ";"))
  write.table(gsea_results, file = gsea_file, sep = "\t", quote = FALSE, row.names = FALSE)

  message("GSEA results saved to: ", gsea_file)

  # Save Camera results if provided
  if (!is.null(camera_results)) {
    camera_file <- file.path(output_dir, paste0(experiment_name, "_Camera_results.tsv"))
    write.table(camera_results, file = camera_file, sep = "\t", quote = FALSE, row.names = TRUE)
    message("Camera results saved to: ", camera_file)
  }

  invisible()
}
