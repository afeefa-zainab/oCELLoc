#' @title oCELLoc: Spatial Transcriptomics Cell Type Prediction
#' @description Predicts average cell type proportions for a spatial transcriptomics sample
#' using Lasso regression on average spot expression. Applies specific lambda selection
#' rules and normalizes output proportions.
#' @name oCELLoc
#'
#' @importFrom rlang .data
NULL # Necessary for roxygen2 documentation generation

#' Predict Average Cell Type Proportions for a Sample
#'
#' This function takes spatial transcriptomics data for a single sample (potentially
#' across multiple spots), calculates the average expression, and predicts
#' average cell type proportions using Lasso regression against a reference dataset.
#' It applies an exponential transformation to input data, uses a specific rule for
#' lambda selection (seeking 3-14 non-zero coefficients), filters coefficients,
#' and normalizes the final proportions to sum to 1.
#'
#' @param spatial_data A data.frame or matrix containing spatial gene expression data.
#'        Genes should be in row names, and columns should represent spots/barcodes.
#'        Assumes expression values are log-transformed (e.g., log(CPM+1) or log(TPM+1)).
#' @param reference A data.frame or matrix containing reference expression data.
#'        Genes should be in row names, cell types should be in column names.
#'        Alternatively, a character string specifying a built-in reference ("human" or "mouse").
#' @param sample_name Optional name for the sample (used in plot titles). If NULL, uses "Sample".
#' @param nfolds Number of folds for cross-validation in `cv.glmnet`. (Default: 5)
#' @param transform_input Logical, whether to apply `exp(data) - 1` transformation to the
#'        input spatial data. Set to `FALSE` if data is already in linear scale (e.g., counts, CPM). (Default: TRUE)
#' @param normalize_reference Logical, whether to normalize each cell type in the reference
#'        to have the same total expression. (Default: TRUE)
#' @param lambda_selection_rule Character, method for lambda selection. Options are:
#'        "auto" (use glmnet's default lambda sequence) or "custom" (use custom lambda range). (Default: "auto")
#' @param alpha The elasticnet mixing parameter, where alpha=1 is the lasso (default) and `alpha=0` is ridge.
#' @param lambda_min Minimum lambda value for custom lambda sequence (only used when lambda_selection_rule="custom"). (Default: 0.001)
#' @param lambda_max Maximum lambda value for custom lambda sequence (only used when lambda_selection_rule="custom"). (Default: 1.0)
#' @param lambda_n Number of lambda values in custom sequence (only used when lambda_selection_rule="custom"). (Default: 100)
#' @param min_nonzero Minimum number of desired non-zero coefficients for lambda selection. (Default: 3)
#' @param max_nonzero Maximum number of desired non-zero coefficients for lambda selection. (Default: 14)
#' @param keep_top_n Maximum number of positive coefficients to retain after filtering. If more
#'        coefficients are positive, only the top `keep_top_n` are kept. Set to `Inf` to disable. (Default: 14)
#' @param nonzero_threshold Threshold below which coefficients are considered zero during lambda
#'        selection and final filtering. (Default: 1e-3)
#' @param generate_plots Logical, whether to generate CV and coefficient path plots. (Default: TRUE)
#' @param verbose Logical, whether to print progress messages. (Default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{proportions}: Data frame with columns 'Cell_Type' and 'Proportion'
#'   \item \code{nonzero_celltypes}: Vector of cell type names with non-zero proportions
#'   \item \code{selected_lambda}: The lambda value selected by the algorithm
#'   \item \code{selection_rule}: Whether lambda was selected by "3-14_rule_glmnet", "3-14_rule_custom", or "fallback"
#'   \item \code{common_genes}: Vector of genes used in the analysis
#'   \item \code{cv_plot}: Function to generate cross-validation ggplot (if generate_plots=TRUE)
#'   \item \code{coef_plot}: Function to generate coefficient path ggplot (if generate_plots=TRUE)
#' }
#' Returns `NULL` if processing fails.
#'
#' @examples
#' # Example 1: Using built-in human reference with glmnet lambda sequence
#' # Load example human average expression data
#' load(system.file("extdata", "human_avg_expression.rda", package = "oCELLoc"))
#' 
#' # Run with built-in human reference and glmnet lambda sequence
#' results_human <- predict_cell_types(
#'   spatial_data = human_avg_expression,
#'   reference = "human",
#'   sample_name = "Human_Example",
#'   lambda_selection_rule = "auto"
#' )
#' 
#' # View top results
#' print(head(results_human$proportions, 10))
#' print(results_human$nonzero_celltypes)
#' 
#' 
#' # Example 2: Using built-in mouse reference with custom lambda sequence
#' # Load example mouse average expression data  
#' load(system.file("extdata", "mouse_avg_expression.rda", package = "oCELLoc"))
#' 
#' # Run with built-in mouse reference and custom lambda sequence
#' results_mouse <- predict_cell_types(
#'   spatial_data = mouse_avg_expression,
#'   reference = "mouse",
#'   sample_name = "Mouse_Example",
#'   lambda_selection_rule = "custom",
#'   lambda_min = 0.001,
#'   lambda_max = 0.5,
#'   lambda_n = 50
#' )
#' 
#' # View top results
#' print(head(results_mouse$proportions, 10))
#' print(results_mouse$nonzero_celltypes)
#' 
#' @export
predict_cell_types <- function(spatial_data,
                                reference,
                                sample_name = NULL,
                                nfolds = 5,
                                transform_input = TRUE,
                                normalize_reference = TRUE,
                                lambda_selection_rule = "auto",
                                alpha = 1,
                                lambda_min = 0.001,
                                lambda_max = 1.0,
                                lambda_n = 100,
                                min_nonzero = 3,
                                max_nonzero = 14,
                                keep_top_n = 14,
                                nonzero_threshold = 1e-3,
                                generate_plots = TRUE,
                                verbose = TRUE) {

  # --- 1. Input Validation and Setup ---
  if (is.null(sample_name)) sample_name <- "Sample"
  if (verbose) message("Starting oCELLoc prediction for sample: ", sample_name)
  
  # Validate lambda_selection_rule
  if (!lambda_selection_rule %in% c("auto", "custom")) {
    stop("lambda_selection_rule must be either 'auto' or 'custom'")
  }
  
  # Validate custom lambda parameters
  if (lambda_selection_rule == "custom") {
    if (lambda_min <= 0 || lambda_max <= 0) {
      stop("lambda_min and lambda_max must be positive when using custom lambda selection")
    }
    if (lambda_min >= lambda_max) {
      stop("lambda_min must be less than lambda_max")
    }
    if (lambda_n < 2) {
      stop("lambda_n must be at least 2")
    }
  }

  if (alpha < 0 || alpha > 1) {
    stop("alpha must be between 0 and 1")
  }

  # Wrap main processing in tryCatch for better error handling
  result_data <- tryCatch({
    # --- 2. Process Spatial Data (V) ---
    if (verbose) message("--> Processing spatial data...")
    
    # Convert to matrix if needed
    if (is.data.frame(spatial_data)) {
      V <- as.matrix(spatial_data)
    } else if (is.matrix(spatial_data)) {
      V <- spatial_data
    } else {
      stop("spatial_data must be a data.frame or matrix")
    }
    
    # Check for gene names in row names
    if (is.null(rownames(V))) {
      stop("spatial_data must have gene names as row names")
    }
    
    # Check for numeric data
    if (!is.numeric(V)) {
      stop("spatial_data must contain numeric expression values")
    }
    
    # Replace NA with 0
    V[is.na(V)] <- 0

    # Apply input transformation if requested
    if(transform_input) {
        if (verbose) message("--> Applying exp(data) - 1 transformation...")
        if (any(V > 50, na.rm = TRUE)) {
            warning("Input data contains large values (>50). Are you sure it's log-transformed? Applying exp(data)-1.")
        }
        if (any(V < 0, na.rm = TRUE)) {
            warning("Input data contains negative values before transformation. Setting negative results to 0.")
            V[V < 0] <- 0
        }
        V <- exp(V) - 1
        V[V < 0] <- 0
    } else {
        if (any(V < 0, na.rm = TRUE)) {
            warning("Input data contains negative values and transform_input=FALSE. Setting negative values to 0.")
            V[V < 0] <- 0
        }
    }

    # --- 3. Process Reference Data (W) ---
    if (verbose) message("--> Loading reference data...")
    
    # Handle different reference input types
    if (is.character(reference) && length(reference) == 1) {
      if (reference %in% c("human", "mouse")) {
        # Load built-in reference (already in new format: genes as rows, cell types as columns)
        ref_name <- paste0(reference, "_ref")
        utils::data(list = ref_name, package = "oCELLoc", envir = environment())
        W <- get(ref_name)

      } else {
        stop("Invalid built-in reference. Must be 'human' or 'mouse'")
      }
    } else if (is.data.frame(reference) || is.matrix(reference)) {
      # Handle user-provided reference (expected to be in format: genes as rows, cell types as columns)
      W <- as.matrix(reference)
    } else {
      stop("reference must be a data.frame, matrix, or valid built-in reference name ('human' or 'mouse')")
    }

    if (verbose) message("--> Processing reference data...")
    
    # Check format (genes should be rows, cell types should be columns)
    if (is.null(rownames(W))) {
      stop("Reference data must have gene names as row names")
    }
    if (is.null(colnames(W))) {
      stop("Reference data must have cell type names as column names")
    }
    
    # Replace NA with 0
    W[is.na(W)] <- 0
    
    # Ensure non-negativity
    if (any(W < 0, na.rm = TRUE)) {
       warning("Reference data contains negative values. Setting them to 0.")
       W[W < 0] <- 0
    }
    
    # Normalize reference if requested
    if (normalize_reference) {
      if (verbose) message("--> Normalizing reference cell types...")
      col_sums <- colSums(W, na.rm = TRUE)
      if (any(col_sums <= 0)) {
        warning("Some cell types have zero or negative total expression. Normalization may fail.")
      }
      # Normalize each cell type to have total expression = 1 (arbitrary but consistent)
      target_sum <- 1
      W <- sweep(W, 2, col_sums / target_sum, "/")
      W[is.na(W)] <- 0
    }

    # --- 4. Find Common Genes ---
    if (verbose) message("--> Finding common genes...")
    common_genes <- intersect(rownames(V), rownames(W))
    n_common <- length(common_genes)
    if (verbose) message("    Found ", n_common, " common genes.")

    if (n_common < 50) {
      stop("Insufficient common genes (", n_common, " < 50). Cannot proceed.")
    }

    # Subset matrices to common genes
    V_sub <- V[common_genes, , drop = FALSE]
    W_sub <- W[common_genes, , drop = FALSE] # Genes x CellTypes

    # --- 5. Prepare for Lasso: Calculate Average Expression ---
    if (verbose) message("--> Calculating average expression across spots/columns...")
    y <- rowMeans(V_sub, na.rm = TRUE)
    x <- W_sub # Reference matrix: Genes x CellTypes

    # Check if y has variance
    if (stats::sd(y, na.rm = TRUE) < 1e-9) {
        stop("Average expression vector 'y' has zero variance across common genes. Cannot run Lasso.")
    }

    # --- 6. Run Lasso Regression ---
    if (verbose) message("--> Running Lasso cross-validation with constraints...")
    
    # Determine lambda sequence based on selection rule
    if (lambda_selection_rule == "custom") {
      # Create custom lambda sequence
      lambda_seq <- exp(seq(log(lambda_max), log(lambda_min), length.out = lambda_n))
      if (verbose) message("    Using custom lambda sequence: min=", lambda_min, ", max=", lambda_max, ", n=", lambda_n)
    } else {
      # Use glmnet's default lambda sequence (set lambda to NULL)
      lambda_seq <- NULL
      if (verbose) message("    Using glmnet's default lambda sequence")
    }
    
    # Run cross-validation WITH constraint
    cv_fit <- glmnet::cv.glmnet(x, y, alpha = alpha, lower.limits = 0, lambda = lambda_seq, nfolds = nfolds)
    all_lambda <- cv_fit$lambda
    fit <- cv_fit$glmnet.fit  # Use the fit from cv.glmnet
    
    if(verbose) message("--> Selecting optimal lambda...")
    
    optimal_lambda <- NA
    best_mse <- Inf
    best_H <- NULL
    selection_rule <- paste0(min_nonzero, "-", max_nonzero, "_rule_", lambda_selection_rule)
    
    # Use cv.glmnet's built-in nzero (number of non-zero coefficients)
    nonzero_counts <- cv_fit$nzero
    
    # Find indices where coefficient count is within range
    valid_indices <- which(nonzero_counts >= min_nonzero & nonzero_counts <= max_nonzero)
    
    if (length(valid_indices) > 0) {
        # Among valid lambdas, find the one with minimum MSE
        valid_mse <- cv_fit$cvm[valid_indices]
        best_idx <- valid_indices[which.min(valid_mse)]
        
        optimal_lambda <- all_lambda[best_idx]
        best_mse <- cv_fit$cvm[best_idx]
        
        # Extract coefficients at optimal lambda
        best_H <- stats::coef(fit, s = optimal_lambda)[-1]
        best_H[best_H < 0] <- 0
        
        if(verbose) {
            message("    Found optimal lambda: ", format(optimal_lambda, digits=4), 
                   " with ", nonzero_counts[best_idx], " coefficients, MSE=", format(best_mse, digits=4))
        }
    } else {
        # Fallback if no lambda with desired number of coefficients
        if(verbose) message("    No lambda found within coefficient range [", min_nonzero, ", ", max_nonzero, "]")
        
        # Find lambdas with any non-zero coefficients (up to max_nonzero)
        valid_indices <- which(nonzero_counts > 0 & nonzero_counts <= max_nonzero)
        
        if (length(valid_indices) > 0) {
            # Among valid lambdas, find the one with the most non-zero coefficients
            # If tie, use the one with lowest MSE
            max_nonzero_count <- max(nonzero_counts[valid_indices])
            best_indices <- valid_indices[nonzero_counts[valid_indices] == max_nonzero_count]
            
            if (length(best_indices) == 1) {
                best_idx <- best_indices
            } else {
                # Multiple lambdas with same count, choose one with lowest MSE
                best_idx <- best_indices[which.min(cv_fit$cvm[best_indices])]
            }
            
            optimal_lambda <- all_lambda[best_idx]
            best_H <- stats::coef(fit, s = optimal_lambda)[-1]
            best_H[best_H < 0] <- 0
            selection_rule <- paste0("fallback_", nonzero_counts[best_idx], "_CTs")
            
            warning(paste("No lambda with", min_nonzero, "-", max_nonzero, "CTs for", sample_name, 
                         "- using best available with", nonzero_counts[best_idx], "CTs"))
            
            if(verbose) {
                message("    Using fallback lambda: ", format(optimal_lambda, digits=4))
                message("    Number of non-zero cell types: ", nonzero_counts[best_idx])
            }
        } else {
            # No lambda gives non-zero coefficients, use lambda.min and force some coefficients
            warning(paste("No lambda with any non-zero coefficients for", sample_name, "- using lambda.min"))
            selection_rule <- "fallback"
            optimal_lambda <- cv_fit$lambda.min
            
            # Extract coefficients at lambda.min
            best_H <- stats::coef(fit, s = optimal_lambda)[-1]
            best_H[best_H < 0] <- 0
            
            # Keep top N cell types even if fewer are non-zero
            top_n <- min(keep_top_n, length(best_H))
            top_indices <- order(best_H, decreasing = TRUE)[1:top_n]
            best_H[-top_indices] <- 0
            
            if(verbose) {
                message("    Using fallback lambda.min: ", format(optimal_lambda, digits=4))
                message("    Keeping top ", sum(best_H > 0), " cell types")
            }
        }
    }
    
    # --- 7. Process Final Coefficients ---
    if(verbose) message("--> Normalizing coefficients...")
    
    # Assign names to coefficients
    names(best_H) <- colnames(W_sub)
    H <- best_H
    
    # Ensure non-negativity (redundant but safe)
    H[H < 0] <- 0
    
    # Additional filtering for keep_top_n if not in fallback mode
    if (selection_rule != "fallback" && !is.infinite(keep_top_n)) {
      non_zero_count <- sum(H > nonzero_threshold)
      if (non_zero_count > keep_top_n) {
        top_indices <- order(H, decreasing = TRUE)[1:keep_top_n]
        H[-top_indices] <- 0
        if(verbose) message("    Reduced from ", non_zero_count, " to ", keep_top_n, " cell types")
      }
    }
    
    # Normalize to sum to 1
    total_coef_sum <- sum(H, na.rm = TRUE)
    if (total_coef_sum > 1e-9) {
      H_norm <- H / total_coef_sum
      if(verbose) message("    Coefficients normalized. Sum before: ", format(total_coef_sum, digits=4))
    } else {
      if (verbose) warning("    Sum of positive coefficients is near zero. Setting all proportions to 0.")
      H_norm <- H
      H_norm[] <- 0
    }

    # Generate CV plot if requested
    cv_plot <- NULL
    if (generate_plots) {
        if(verbose) message("--> Preparing CV plot function...")
        tryCatch({
            # Create a function that generates the CV ggplot
            cv_plot <- function() {
              cv_data <- data.frame(
                lambda = cv_fit$lambda,
                mean_mse = cv_fit$cvm,
                se = cv_fit$cvsd
              )
              
              ggplot2::ggplot(cv_data, ggplot2::aes(x = log(.data$lambda), y = .data$mean_mse)) +
                ggplot2::geom_line(color = "blue") +
                ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$mean_mse - .data$se, 
                                                ymax = .data$mean_mse + .data$se), 
                                    alpha = 0.2) +
                ggplot2::geom_vline(xintercept = log(optimal_lambda), 
                                  linetype = "dashed", color = "red") +
                ggplot2::ggtitle(paste("CV Curve -", sample_name)) +
                ggplot2::xlab("log(Lambda)") + 
                ggplot2::ylab("MSE") +
                ggplot2::theme_minimal()
            }
        }, error = function(e) {
            warning("Failed to prepare CV plot function: ", e$message)
            cv_plot <<- NULL
        })
    }

    # Generate coefficient path plot if requested
    coef_plot <- NULL
    if (generate_plots) {
       if(verbose) message("--> Preparing coefficient path plot function...")
       tryCatch({
            # Create a function that generates the coefficient path ggplot
            coef_plot <- function() {
              beta_matrix <- as.matrix(fit$beta)
              
              # Create melted data frame for plotting
              df_beta <- reshape2::melt(beta_matrix)
              colnames(df_beta) <- c("Gene", "Lambda_Index", "Coefficient")
              df_beta$Lambda <- rep(log(fit$lambda), each = nrow(beta_matrix))
              
              # Filter out zero coefficients for plotting
              df_beta_nonzero <- df_beta[abs(df_beta$Coefficient) > 1e-6, ]
              
              if (nrow(df_beta_nonzero) > 0) {
                ggplot2::ggplot(df_beta_nonzero, ggplot2::aes(x = .data$Lambda, y = .data$Coefficient, color = .data$Gene)) +
                  ggplot2::geom_line(show.legend = FALSE, linewidth = 0.8) +
                  ggplot2::geom_vline(xintercept = log(optimal_lambda), 
                                    linetype = "dashed", color = "red") +
                  ggplot2::ggtitle(paste("Lasso Path -", sample_name)) +
                  ggplot2::xlab("log(Lambda)") + 
                  ggplot2::ylab("Coefficient") +
                  ggplot2::theme_minimal()
              } else {
                ggplot2::ggplot() + 
                  ggplot2::ggtitle(paste("Lasso Path - No non-zero coefficients -", sample_name)) +
                  ggplot2::theme_minimal()
              }
            }
       }, error = function(e) {
           warning("Failed to prepare coefficient path plot function: ", e$message)
           coef_plot <<- NULL
       })
    }

    # --- 8. Format Results ---
    if(verbose) message("--> Formatting results...")
    
    # Create results data frame (without redundant row names)
    results_df <- data.frame(
      Cell_Type = names(H_norm), 
      Proportion = as.numeric(H_norm),
      stringsAsFactors = FALSE
    )
    results_df <- results_df[order(results_df$Proportion, decreasing = TRUE), ]
    rownames(results_df) <- NULL  # Remove redundant row names
    
    # Get non-zero cell types
    nonzero_celltypes <- results_df$Cell_Type[results_df$Proportion > 0]
    
    # Create return list
    result_list <- list(
      proportions = results_df,
      nonzero_celltypes = nonzero_celltypes,
      selected_lambda = optimal_lambda,
      selection_rule = selection_rule,
      common_genes = common_genes,
      cv_plot = cv_plot,
      coef_plot = coef_plot
    )
    
    if (verbose) message("--> Prediction for sample ", sample_name, " finished successfully.")
    
    return(result_list)

  }, error = function(e) {
      error_message <- paste("Error processing sample:", sample_name, "-", e$message)
      warning(error_message)
      return(NULL)
  })

  return(result_data)
}