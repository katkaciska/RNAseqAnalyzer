#' Filter for Protein-Coding Genes
#'
#' Filters expression data to keep only protein-coding genes using GTF annotation.
#'
#' @param expr_df Data frame with gene expression (must have a 'gene_id' column)
#' @param gtf GTF object loaded with rtracklayer::import()
#' @return Data frame with only protein-coding genes, rownames as geneName_geneID
#' @export
filter_protein_coding_genes <- function(expr_df, gtf) {
  # Check required packages
  if (!requireNamespace("rtracklayer", quietly = TRUE) ||
      !requireNamespace("stringr", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Required packages (rtracklayer, stringr, dplyr, S4Vectors) not installed")
  }

  mcols_gtf <- S4Vectors::mcols(gtf)

  # Decide which biotype column to use
  if ("gene_biotype" %in% colnames(mcols_gtf)) {
    biotype_col <- "gene_biotype"
  } else if ("gene_type" %in% colnames(mcols_gtf)) {
    biotype_col <- "gene_type"
  } else {
    stop("Neither 'gene_biotype' nor 'gene_type' found in GTF metadata.")
  }

  # Filter GTF for protein-coding gene entries
  gene_entries <- gtf[
    gtf$type == "gene" &
      !is.na(mcols_gtf[[biotype_col]]) &
      mcols_gtf[[biotype_col]] == "protein_coding"
  ]

  # Extract minimal annotation
  gene_annotations <- S4Vectors::mcols(gene_entries)[, c("gene_id", biotype_col, "gene_name")]
  gene_annotations_df <- as.data.frame(gene_annotations)

  # Clean gene_id versions in GTF annotations
  gene_annotations_df <- gene_annotations_df %>%
    dplyr::mutate(
      gene_id_clean = stringr::str_replace(gene_id, "\\..*", ""),
      gene_name_gtf = gene_name
    ) %>%
    dplyr::select(gene_id_clean, gene_name_gtf, dplyr::all_of(biotype_col))

  # Clean gene_id versions in expression table
  # If expr_df already has a gene_id column, use it; otherwise assume rownames
  if (!"gene_id" %in% colnames(expr_df)) {
    expr_df <- expr_df %>%
      tibble::rownames_to_column("gene_id")
  }

  expr_df <- expr_df %>%
    dplyr::mutate(gene_id_clean = stringr::str_replace(gene_id, "\\..*", ""))

  # Join expression with GTF annotation on cleaned gene_id
  annotated_expr <- dplyr::inner_join(
    expr_df,
    gene_annotations_df,
    by = "gene_id_clean"
  )

  # Choose gene symbol: fall back to gene_id if missing
  gene_symbol <- annotated_expr$gene_name_gtf
  missing_symbol <- is.na(gene_symbol) | gene_symbol == ""
  gene_symbol[missing_symbol] <- annotated_expr$gene_id_clean[missing_symbol]

  # Create new rownames: GeneSymbol_geneID
  annotated_expr$gene_name_id <- paste0(gene_symbol, "_", annotated_expr$gene_id_clean)
  rownames(annotated_expr) <- annotated_expr$gene_name_id

  # Drop metadata columns but KEEP gene_name_id
  cleaned_expr <- annotated_expr %>%
    dplyr::select(
      -gene_id,
      -gene_id_clean,
      -gene_name_gtf,
      -dplyr::all_of(biotype_col)
    )

  message("Filtered to ", nrow(cleaned_expr), " protein-coding genes")
  print(cleaned_expr)
  return(cleaned_expr)

}
#' Prepare Expression Data for Analysis
#'
#' Prepares expression data by log-transforming and filtering samples.
#'
#' @param expr_data Expression data matrix (genes x samples)
#' @param log_transform Logical, whether to log2 transform (default: TRUE)
#' @param pseudocount Pseudocount to add before log transformation (default: 1)
#' @param remove_samples Character vector of sample names to remove (default: NULL)
#' @return Processed expression matrix
#' @export
prepare_expression_data <- function(expr_data,
                                    log_transform = TRUE,
                                    pseudocount = 1,
                                    remove_samples = NULL) {

  # Separate numeric and non-numeric columns
  numeric_cols <- vapply(expr_data, is.numeric, logical(1))
  if (!all(numeric_cols)) {
    message("Detected non-numeric columns: ",
            paste(colnames(expr_data)[!numeric_cols], collapse = ", "),
            ". These will be kept aside and not log-transformed.")
  }

  expr_num <- as.matrix(expr_data[, numeric_cols, drop = FALSE])
  anno_df  <- NULL
  if (!all(numeric_cols)) {
    anno_df <- expr_data[, !numeric_cols, drop = FALSE]
  }

  # Remove specified samples (columns) in the numeric part
  if (!is.null(remove_samples)) {
    samples_to_keep <- !colnames(expr_num) %in% remove_samples
    expr_num <- expr_num[, samples_to_keep, drop = FALSE]
    message("Removed ", sum(!samples_to_keep), " samples")
  }

  # Log transform numeric expression
  if (log_transform) {
    expr_num <- log2(expr_num + pseudocount)
    message("Applied log2(x + ", pseudocount, ") transformation")
  }

  message("Final data dimensions: ", nrow(expr_num), " genes x ", ncol(expr_num), " samples")

  # Option A: return only numeric matrix for downstream modelling
  return(expr_num)

  # Option B (if you prefer): return both expression and annotation
  # return(list(expr = expr_num, anno = anno_df))
}

#' Clean Sample Names
#'
#' Removes leading 'X' from sample names that R adds to numeric column names.
#'
#' @param sample_names Character vector of sample names
#' @return Character vector of cleaned sample names
#' @export
clean_sample_names <- function(sample_names) {
  gsub("^X", "", sample_names)
}

#' Match Metadata to Expression Data
#'
#' Matches and orders metadata to match the sample order in expression data.
#'
#' @param expr_data Expression data matrix with samples as columns
#' @param metadata Data frame with metadata, must have sample ID column
#' @param sample_id_col Name of the sample ID column in metadata (default: "SampleID")
#' @return Ordered metadata data frame matching expression data column order
#' @export
match_metadata_to_expression <- function(expr_data, metadata, sample_id_col = "SampleID") {
  # Clean sample names from expression data
  sample_ids <- clean_sample_names(colnames(expr_data))

  # Match metadata order to expression data columns
  metadata_ordered <- metadata[match(sample_ids, metadata[[sample_id_col]]), ]

  # Check for missing matches
  missing_samples <- is.na(metadata_ordered[[sample_id_col]])
  if (any(missing_samples)) {
    warning("Could not find metadata for samples: ",
            paste(sample_ids[missing_samples], collapse = ", "))
  }

  # Remove rows with missing metadata
  metadata_ordered <- metadata_ordered[!missing_samples, , drop = FALSE]

  return(metadata_ordered)
}

#' Filter Low Expression Genes
#'
#' Removes genes with low expression across samples.
#'
#' @param expr_data Expression data matrix (genes x samples)
#' @param min_count Minimum count threshold (default: 10)
#' @param min_samples Minimum number of samples that must meet threshold (default: 2)
#' @return Filtered expression matrix
#' @export
filter_low_expression_genes <- function(expr_data, min_count = 10, min_samples = 2) {
  # Count samples meeting threshold for each gene
  keep_genes <- rowSums(expr_data >= min_count) >= min_samples

  # Filter
  filtered_data <- expr_data[keep_genes, , drop = FALSE]

  message("Kept ", sum(keep_genes), " genes out of ", nrow(expr_data),
          " (", round(sum(keep_genes)/nrow(expr_data)*100, 1), "%)")

  return(filtered_data)
}
