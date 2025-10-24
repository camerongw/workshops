library(Matrix)
library(dbscan)
library(dplyr)
library(ggplot2)
library(reshape2)

# Compute spatial neighbourhood means once per replicate for *all* variables.
# This replaces repeated averaging inside every pairwise comparison.
compute_spatial_means_once <- function(
    input_data,
    reps = 1:length(input_data),
    x = "x", y = "y",
    radius = 50,
    nbins = 2000,
    vars = NULL,
    seed = NULL,
    verbose = FALSE
){
  if (!is.null(seed)) set.seed(seed)
  
  # Infer variables if not provided
  if (is.null(vars)) {
    vars <- setdiff(colnames(input_data[[1]]), c(x, y))
  }
  
  means_list <- vector("list", length(reps))
  meta_list  <- vector("list", length(reps))
  
  for (k in seq_along(reps)) {
    i <- reps[k]
    df <- input_data[[i]]
    
    # Neighbourhoods at given radius (computed once per replicate)
    pos <- dbscan::frNN(df[, c(x, y)], eps = radius)
    
    n <- nrow(df)
    # Sample bins (spots) where we'll compute neighbourhood means
    bins <- sample.int(n, size = min(nbins, n))  # safety if n < nbins
    
    if (verbose) message(sprintf("Rep %d: computing means on %d bins", i, length(bins)))
    
    # For each variable, compute mean over each bin's neighbours
    # Result: nbins x p matrix (p = number of variables)
    var_means <- sapply(
      vars,
      function(v) vapply(bins, function(idx) {
        mean(df[[v]][pos$id[[idx]]])
      }, numeric(1))
    )
    
    # Store as data.frame with variables as columns
    means_df <- as.data.frame(var_means)
    colnames(means_df) <- vars
    
    # Optional meta info (replicate id, bin index)
    meta_df <- data.frame(rep = i, bin = bins)
    
    means_list[[k]] <- means_df
    meta_list[[k]]  <- meta_df
  }
  
  pooled_means <- do.call(rbind, means_list)
  pooled_meta  <- do.call(rbind, meta_list)
  
  list(
    pooled_means = pooled_means,  # rows = all bins across all reps; cols = variables
    pooled_meta  = pooled_meta,
    vars         = vars
  )
}

# Vectorised neighbourhood averaging via sparse matrices:
# builds W (row-normalised adjacency), then M = W %*% X for all variables at once.
compute_spatial_means_sparse <- function(
    input_data,
    reps = 1:length(input_data),
    x = "x", y = "y",
    radius = 50,
    nbins = 10000,
    vars = NULL,
    seed = NULL,
    verbose = FALSE
){
  if (!is.null(seed)) set.seed(seed)
  if (is.null(vars)) vars <- setdiff(colnames(input_data[[1]]), c(x, y))
  
  means_list <- vector("list", length(reps))
  meta_list  <- vector("list", length(reps))
  
  for (k in seq_along(reps)) {
    i <- reps[k]
    df <- input_data[[i]]
    n  <- nrow(df)
    
    pos <- dbscan::frNN(df[, c(x, y)], eps = radius)
    
    # Build sparse adjacency: one row per spot, entries for its neighbours
    rows <- rep(seq_len(n), lengths(pos$id))
    cols <- unlist(pos$id)
    vals <- rep(1, length(cols))
    A    <- sparseMatrix(i = rows, j = cols, x = vals, dims = c(n, n))
    
    # Ensure self-neighbourhood (avoid zero row sums)
    A <- A + Diagonal(n)
    
    # Row-normalise to get mean operator W
    rs <- Matrix::rowSums(A)
    W  <- Diagonal(n, 1/rs) %*% A
    
    # All variables at once
    X <- as.matrix(df[, vars, drop = FALSE])
    M <- as.matrix(W %*% X)  # n x p (means at every spot)
    
    bins <- sample.int(n, size = min(nbins, n))
    means_df <- as.data.frame(M[bins, , drop = FALSE])
    colnames(means_df) <- vars
    
    meta_df <- data.frame(rep = i, bin = bins)
    means_list[[k]] <- means_df
    meta_list[[k]]  <- meta_df
    
    if (verbose) message(sprintf("Rep %d: sparse averaging done (%d bins)", i, length(bins)))
  }
  
  pooled_means <- do.call(rbind, means_list)
  pooled_meta  <- do.call(rbind, meta_list)
  
  list(
    pooled_means = pooled_means,
    pooled_meta  = pooled_meta,
    vars         = vars
  )
}

## downstream pairwise correlations
coloc_analysis_replicates_fast <- function(
    input_data,
    reps = 1:length(input_data),
    x = "x", y = "y",
    nbins = 2000,
    radius = 50,
    verbose = FALSE,
    use_sparse = TRUE,
    compute_p = FALSE # set TRUE if you really need p-values per pair
){
  # 1) Compute neighbourhood means once (choose sparse or simple method)
  pre <- if (use_sparse) {
    compute_spatial_means_sparse(
      input_data = input_data,
      reps = reps, x = x, y = y,
      radius = radius, nbins = nbins,
      verbose = verbose
    )
  } else {
    compute_spatial_means_once(
      input_data = input_data,
      reps = reps, x = x, y = y,
      radius = radius, nbins = nbins,
      verbose = verbose
    )
  }
  
  vars <- pre$vars
  P    <- pre$pooled_means  # rows: all bins across reps; cols: variables
  
  # 2) Full Spearman correlation matrix across variables (computed once)
  rho_mat <- suppressWarnings(cor(P[, vars, drop = FALSE], method = "spearman", use = "pairwise.complete.obs"))
  
  # Optional p-values (loop only over upper triangle)
  p_mat <- matrix(NA_real_, nrow = length(vars), ncol = length(vars),
                  dimnames = list(vars, vars))
  if (compute_p) {
    for (a in seq_along(vars)) {
      for (b in a:length(vars)) {
        if (a == b) { p_mat[a, b] <- NA_real_; next }
        ct <- suppressWarnings(cor.test(P[, vars[a]], P[, vars[b]], method = "spearman"))
        p_mat[a, b] <- ct$p.value
        p_mat[b, a] <- ct$p.value
      }
    }
  }
  
  # 3) Long-format for plotting
  heatmap_data <- expand.grid(celltype1 = vars, celltype2 = vars, stringsAsFactors = FALSE)
  heatmap_data$ρ     <- as.vector(rho_mat)
  heatmap_data$pval  <- if (compute_p) as.vector(p_mat) else NA_real_
  
  # Same-type comparisons set to NA (as you did)
  same <- heatmap_data$celltype1 == heatmap_data$celltype2
  heatmap_data$ρ[same]    <- NA_real_
  heatmap_data$pval[same] <- NA_real_
  
  # 4) Order per hierarchical clustering
  # Use 1 - rho as a distance; replace NAs with 0 on the diagonal just for clustering
  rho_for_clust <- rho_mat
  diag(rho_for_clust) <- 1  # perfect self-correlation
  rho_for_clust[is.na(rho_for_clust)] <- 0
  dist_mat <- as.dist(1 - rho_for_clust)
  hc <- hclust(dist_mat, method = "average")
  
  # Reorder factors
  ord <- hc$order
  heatmap_data$celltype1 <- factor(heatmap_data$celltype1, levels = vars[ord], ordered = TRUE)
  heatmap_data$celltype2 <- factor(heatmap_data$celltype2, levels = vars[ord], ordered = TRUE)
  
  # 5) Plot
  cor_plot <- ggplot(heatmap_data, aes(x = celltype1, y = celltype2, fill = ρ)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "royalblue3", mid = "white", high = "red3",
      midpoint = 0, na.value = "grey80", limits = c(-1, 1)
    ) +
    scale_y_discrete(limits = rev(levels(heatmap_data$celltype2))) +
    xlab(NULL) + ylab(NULL) +
    theme_minimal() +
    scale_x_discrete(position = "top") +
    theme(
      axis.text.x = element_text(hjust = 0, angle = 45),
      legend.position = "right",
      legend.direction = "vertical",
      axis.line = element_line(colour = "black")
    ) +
    coord_fixed()
  
  # Output
  list(
    plot = cor_plot,
    heatmap_data = heatmap_data,
    rho_mat = rho_mat,
    p_mat = if (compute_p) p_mat else NULL,
    pooled_means = P,
    meta = pre$pooled_meta
  )
}

## test
outs <- coloc_analysis_replicates_fast(res)
