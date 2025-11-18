# Required packages
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(msigdbr)

# -------------------------------
# Your ssgsea function
# -------------------------------
ssgsea <- function(X, gene_sets, alpha = 0.25, scale = TRUE, norm = FALSE, single = TRUE) {
  row_names <- rownames(X)
  num_genes <- nrow(X)
  gene_sets <- lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  R <- matrixStats::colRanks(X, preserveShape = TRUE, ties.method = 'average')
  
  es <- apply(R, 2, function(R_col) {
    gene_ranks <- order(R_col, decreasing = TRUE)
    
    es_sample <- sapply(gene_sets, function(gene_set_idx) {
      indicator_pos <- gene_ranks %in% gene_set_idx
      indicator_neg <- !indicator_pos
      rank_alpha <- (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos <- cumsum(rank_alpha) / sum(rank_alpha)
      step_cdf_neg <- cumsum(indicator_neg) / sum(indicator_neg)
      step_cdf_diff <- step_cdf_pos - step_cdf_neg
      
      if (scale) step_cdf_diff <- step_cdf_diff / num_genes
      if (single) sum(step_cdf_diff) else step_cdf_diff[which.max(abs(step_cdf_diff))]
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es <- matrix(es, nrow = 1)
  if (norm) es <- es / diff(range(es))
  
  rownames(es) <- names(gene_sets)
  colnames(es) <- colnames(X)
  return(es)
}

# -------------------------------
# Function to collapse duplicate gene symbols
# -------------------------------
collapse_by_symbol <- function(expr_df) {
  # Split gene_name_ENSMUSG -> gene_name
  expr_df$symbol <- sub("_.*", "", rownames(expr_df))
  message("1 this ok")
  expr_df_collapsed <- expr_df %>%
    group_by(symbol) %>%
    summarize(across(where(is.numeric), sum))
  message("2 this ok")
  mat <- as.matrix(expr_df_collapsed[,-1])
  message("3 this ok")
  rownames(mat) <- expr_df_collapsed$symbol
  message("4 this ok")
  return(mat)
}

# -------------------------------
# Function to run ssgsea and plot heatmap
# -------------------------------
run_ssgsea_analysis <- function(
    tpm_file,
    metadata_file,
    gene_sets = NULL,   # list of vectors for custom gene sets
    gtf_file,
    remove_samples = NULL,
    species = "MM",     # "MM" for mouse, "HS" for human
    output_prefix = "ssgsea",
    plot_gene_level = TRUE
) {
  # Load data
  tpm <- read.table(tpm_file, header = TRUE, sep = "\t")
  colnames(tpm) <- c(colnames(tpm)[1:2], clean_sample_names(colnames(tpm)[3:ncol(tpm)]))
  
  metadata <- read.delim(metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Filter protein coding genes
  gtf <- rtracklayer::import(gtf_file)
  pc_tpm <- filter_protein_coding_genes(tpm, gtf)
  
  # Remove unwanted samples
  if (!is.null(remove_samples)) {
    pc_tpm <- pc_tpm[, !colnames(pc_tpm) %in% remove_samples]
    metadata <- metadata[!metadata$SampleID %in% remove_samples, ]
  }
  print(head(pc_tpm))
  # Collapse duplicates to gene symbols
  expr_mat <- collapse_by_symbol(pc_tpm)
  print(head(expr_mat))
  # Log-transform
  expr_mat_log <- log2(expr_mat + 1)
  
  # MSigDB Hallmark loading if gene_sets not provided
  if (is.null(gene_sets)) {
    species1 <- ifelse(species == "MM", "Mus musculus", "Homo sapiens")
    h_gene_sets <- msigdbr::msigdbr(db_species = species, species = species1, collection = "H")
    gene_sets_list <- split(h_gene_sets$gene_symbol, h_gene_sets$gs_name)
  } else {
    gene_sets_list <- gene_sets
  }
  
  print(gene_sets_list)
  # Run ssgsea
  es <- ssgsea(expr_mat_log, gene_sets_list)
  which(!is.finite(es), arr.ind = TRUE)
  es <- es[rowSums(is.finite(es)) > 0, ]
  print(es)
  # -------------------------------
  # Plot heatmap
  # -------------------------------
  # Create annotation column: condition
  
  meta_col <- metadata$condition
  names(meta_col) <- metadata$SampleID
  ha <- HeatmapAnnotation(Condition = meta_col[colnames(es)],
                          col = list(Condition = c("Vehicle" = "lightblue", "Mirdametinib" = "salmon")))
  
  heatmap_plot<- Heatmap(es,
          name = "ssGSEA",
          top_annotation = ha,
          show_row_names = TRUE,
          show_column_names = TRUE,
          cluster_rows = TRUE,
          cluster_columns = TRUE)
  # Optionally save
  if(!is.null(output_prefix)){
    pdf(paste0(output_prefix, "_ssgsea_heatmap.pdf"))
    print(heatmap_plot)
    dev.off()
  }
  
  
  # 2️⃣ Gene-level heatmap
  if(plot_gene_level & !is.null(gene_sets)) {
    genes_to_plot <- unique(unlist(gene_sets_list))
    genes_in_matrix <- genes_to_plot[genes_to_plot %in% rownames(expr_mat_log)]
    expr_subset <- expr_mat_log[genes_in_matrix, , drop=FALSE]
    
    heatmap_genes <- Heatmap(expr_subset,
                             name = "Expression",
                             top_annotation = ha,
                             show_row_names = TRUE,
                             show_column_names = TRUE,
                             cluster_rows = TRUE,
                             cluster_columns = TRUE)
    
    if(!is.null(output_prefix)){
      pdf(paste0(output_prefix, "_genes_heatmap.pdf"))
      print(heatmap_genes)
      dev.off()
    }
  }
  return(list(ssgsea_scores = es))
}

# -------------------------------
# Example: single custom gene set
# -------------------------------
# gene_set1 <- c("Clec4e", "Il1b", "Tlr4")
# run_ssgsea_analysis(tpm_file = "Data/salmon.merged.gene_tpm.tsv",
#                     metadata_file = "Data/testCondition.txt",
#                     gene_sets = list(CustomSet = gene_set1),
#                     gtf_file = "Data/Mus_musculus.GRCm39.111.gtf",
#                     remove_samples = c("4183"),
#                     species = "MM",
#                     output_prefix = "CustomSet")

# -------------------------------
# Example: multiple custom gene sets
# -------------------------------
# gene_sets_multi <- list(
#   Inflammation = c("Clec4e", "Il1b", "Tlr4"),
#   CellCycle = c("Ccna2", "Ccnb1", "Cdk1")
# )
# run_ssgsea_analysis(tpm_file = "Data/salmon.merged.gene_tpm.tsv",
#                     metadata_file = "Data/testCondition.txt",
#                     gene_sets = gene_sets_multi,
#                     gtf_file = "Data/Mus_musculus.GRCm39.111.gtf",
#                     remove_samples = c("4183"),
#                     species = "MM",
#                     output_prefix = "MultiGeneSets")

# -------------------------------
# Example: Hallmark gene sets
# -------------------------------
# run_ssgsea_analysis(tpm_file = "Data/salmon.merged.gene_tpm.tsv",
#                     metadata_file = "Data/testCondition.txt",
#                     gene_sets = NULL,
#                     gtf_file = "Data/Mus_musculus.GRCm39.111.gtf",
#                     remove_samples = c("4183"),
#                     species = "MM",
#                     output_prefix = "Hallmarks")
