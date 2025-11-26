#' Human Cell Type Reference Data
#'
#' A reference dataset containing gene expression profiles for various human cell types.
#' Used as a reference for the predict_cell_types function to predict cell type proportions
#' in spatial transcriptomics data.
#'
#' @format A data frame with rows corresponding to cell types and columns corresponding to genes.
#'   The data frame has a 'cell_type' column that identifies the cell type, and numerous gene columns
#'   with expression values.
#'
#' @source Generated from reference single-cell RNA sequencing datasets
"human_ref"

#' Mouse Cell Type Reference Data
#'
#' A reference dataset containing gene expression profiles for various mouse cell types.
#' Used as a reference for the predict_cell_types function to predict cell type proportions
#' in spatial transcriptomics data.
#'
#' @format A data frame with rows corresponding to cell types and columns corresponding to genes.
#'   The data frame has a 'cell_type' column that identifies the cell type, and numerous gene columns
#'   with expression values.
#'
#' @source Generated from reference single-cell RNA sequencing datasets
"mouse_ref"