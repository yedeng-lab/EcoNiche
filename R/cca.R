#' CCA-Based Niche Analysis Workflow Controller
#'
#' A controller function that orchestrates the CCA-based niche analysis workflow.
#' It dispatches the execution to either the gradient workflow (Steps 1-4)
#' or the group workflow (Steps 1-2, 5-6) based on the specified \code{mode}.
#'
#' @param mode A character string specifying the analysis mode. 
#'   Must be one of \code{"gradient"} or \code{"group"}.
#' @param ... Additional arguments passed to the specific workflow functions 
#'   (\code{\link{cca_workflow_gradient}} or \code{\link{cca_workflow_group}}).
#'
#' @return A list containing the results of the executed workflow steps. 
#'   See \code{\link{cca_workflow_gradient}} or \code{\link{cca_workflow_group}} for structure details.
#' 
#' @examples
#' \dontrun{
#' # Run the gradient workflow
#' res <- cca_workflow(mode = "gradient", otu = my_otu, env = my_env, sel = "Temp")
#' }
#' @export
cca_workflow <- function(mode = c("gradient", "group"), ...) {
  mode <- match.arg(mode)
  if (mode == "gradient") {
    cca_workflow_gradient(...)
  } else {
    cca_workflow_group(...)
  }
}

# -------------------------------------------------------------------
# Step 1: PCA on environmental variables (screening / checking)
# -------------------------------------------------------------------

#' Prepare Environmental Data for CCA (Step 1)
#'
#' Performs standardization and Principal Component Analysis (PCA) on environmental variables.
#' This step assists in selecting constrained variables and covariates for downstream 
#' CCA or partial-CCA.
#'
#' @param env A data.frame of environmental variables. 
#'   Rows must be sample IDs and columns must be environmental variables.
#' @param sel A vector of variable names or column indices to be used as constrained variables.
#' @param constrain Optional. A vector of variable names or column indices to be used as covariates.
#'   (Note: in the legacy workflow this parameter corresponds to 'covariates').
#' @param standardize Logical; if \code{TRUE}, environmental variables are standardized 
#'   (centered and scaled) column-wise. Default is \code{TRUE}.
#' @param galaxy_colnum Logical; if \code{TRUE} and numeric indices are provided, 
#'   they are treated as 1-based indices where column 1 is assumed to be SampleID 
#'   (and thus indices are decremented by 1). Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{env_sd}{The standardized environmental data frame.}
#'   \item{sel_idx}{Indices of selected variables.}
#'   \item{con_idx}{Indices of constrained covariates.}
#'   \item{pca}{The PCA result object (from \code{\link[stats]{prcomp}}).}
#'   \item{combined}{A data frame combining SampleID, PCA axes, and constrained variables.}
#'   \item{cor_res}{A matrix/table of correlations between PCA axes and original variables, 
#'     including significance levels.}
#' }
#' @export
cca_prep_env <- function(env,
                         sel,
                         constrain = NULL,
                         standardize = TRUE,
                         galaxy_colnum = TRUE) {
  env <- as.data.frame(env, check.names = FALSE)
  
  if (is.null(rownames(env))) {
    stop("Step 1 Error: `env` must have row names as sample IDs (e.g., read.table(..., row.names=1)).")
  }
  
  # 1) standardize
  if (isTRUE(standardize)) {
    env_sd <- vegan::decostand(env, method = "standardize", MARGIN = 2)
  } else {
    env_sd <- env
  }
  
  # Helper: Index resolution
  to_idx_kai <- function(x, cn, galaxy_colnum = TRUE) {
    if (is.null(x) || length(x) == 0) return(integer(0))
    
    # Handle character input
    if (is.character(x)) {
      # Allow comma-separated string "2,3,4"
      if (length(x) == 1 && grepl(",", x, fixed = TRUE)) {
        x <- strsplit(x, ",", fixed = TRUE)[[1]]
        x <- trimws(x)
      }
      
      # If all strings are digits, treat as numeric
      if (all(grepl("^\\d+$", x))) {
        x <- as.numeric(x)
      } else {
        idx <- match(x, cn)
        if (anyNA(idx)) {
          bad <- x[is.na(idx)]
          stop(sprintf("Step 1: unknown variable name(s) in selection: %s", paste(bad, collapse = ", ")))
        }
        return(idx)
      }
    }
    
    # Handle numeric input
    if (!is.numeric(x)) stop("Step 1 Error: selection must be variable names or numeric indices.")
    
    # Galaxy compatibility: Column 1 is often SampleID in raw files, so it is forbidden here
    if (isTRUE(galaxy_colnum) && any(x == 1)) {
      stop("Step 1 Error: Column 1 is SampleID; please choose another column.")
    }
    
    if (isTRUE(galaxy_colnum)) {
      x <- x - 1
    }
    
    x <- as.integer(x)
    if (any(x < 1) || any(x > length(cn))) {
      stop("Step 1 Error: selection index out of range after column adjustment.")
    }
    unique(x)
  }
  
  cn <- colnames(env_sd)
  sel_idx <- to_idx_kai(sel, cn, galaxy_colnum = galaxy_colnum)
  con_idx <- to_idx_kai(constrain, cn, galaxy_colnum = galaxy_colnum)
  
  # Requirement: at least 2 variables for PCA in this specific workflow
  if (length(sel_idx) < 2) {
    stop("Step 1 Error: Please select more than two variables in `sel` for PCA.")
  }
  
  # avoid overlap (keeps logic clean; Kai doesn't explicitly check, but overlap can confuse interpretation)
  if (length(intersect(sel_idx, con_idx)) > 0) {
    stop("Step 1: the same variable cannot be in both `sel` and `constrain`.")
  }
  
  # 2) PCA
  env_sd2 <- env_sd[, sel_idx, drop = FALSE]
  pca <- stats::prcomp(t(env_sd2), scale. = TRUE)
  
  # 3) Build combined table: SampleID + PCA axes + constrain vars
  pca_axis <- as.data.frame(pca$rotation, check.names = FALSE)
  combined <- data.frame(
    SampleID = rownames(env_sd),
    pca_axis,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (length(con_idx) > 0) {
    combined <- cbind(combined, env_sd[, con_idx, drop = FALSE])
    colnames(combined) <- c("SampleID", colnames(pca_axis), colnames(env)[con_idx])
  } else {
    colnames(combined) <- c("SampleID", colnames(pca_axis))
  }
  
  # 4) Correlate PCA axis to selected variables
  tmp <- data.frame(matrix(NA_character_, nrow = ncol(env_sd2), ncol = ncol(pca_axis)),
                    stringsAsFactors = FALSE)
  rownames(tmp) <- colnames(env_sd2)
  colnames(tmp) <- colnames(pca_axis)
  
  # Use total comparisons for adjustment conceptualization, 
  # though method="none" is used below to preserve raw p-values in display string.
  ncom <- nrow(t(pca_axis)) * ncol(env_sd2)
  
  for (i in seq_len(ncol(pca_axis))) {
    for (j in seq_len(ncol(env_sd2))) {
      ct <- stats::cor.test(as.numeric(pca_axis[, i]),
                            as.numeric(env_sd2[, j]),
                            method = "pearson")
      
      r_txt <- format(signif(unname(ct$estimate), digits = 3), scientific = FALSE)
      p_txt <- format(signif(stats::p.adjust(ct$p.value, method = "none", n = ncom), digits = 3),
                      scientific = FALSE)
      
      tmp[j, i] <- paste0(r_txt, " (", p_txt, ")")
    }
  }
  
  # Combine importance summary and correlations
  cor_res <- rbind(summary(pca)$importance, tmp)
  list(
    env_sd    = env_sd,
    sel_idx   = sel_idx,
    con_idx   = con_idx,
    pca       = pca,
    combined  = combined,
    cor_res   = cor_res
  )
}


# -------------------------------------------------------------------
# Step 2: CCA / partial-CCA; site scores on position and width axes
# -------------------------------------------------------------------

#' Fit CCA and Partial-CCA Ordination (Step 2)
#'
#' Performs Constrained Correspondence Analysis (CCA) and Partial CCA to obtain
#' site scores for niche position and niche width calculations.
#'
#' @param otu An OTU/species matrix or data frame (rows = taxa, columns = samples).
#' @param env A data.frame of environmental variables (rows = samples).
#' @param sel Variables (names or indices) defining the width axis (used in partial CCA).
#' @param covariates Variables (names or indices) defining the position axis (used in CCA and as covariates in partial CCA).
#' @param standardize Logical; if \code{TRUE}, environmental variables are standardized. Default is \code{TRUE}.
#' @param galaxy_colnum Logical; for numeric indices, whether to treat them as Galaxy-style 1-based indices (skipping column 1).
#'
#' @return A list containing:
#' \describe{
#'   \item{site_width}{Matrix of site scores from the partial CCA (width axes).}
#'   \item{site_pos}{Named vector of site scores from the first axis of the CCA (position axis).}
#'   \item{site_pos_mat}{Matrix of all site scores from the CCA.}
#'   \item{cca}{The CCA model object.}
#'   \item{pcca}{The Partial CCA model object.}
#' }
#' @export
cca_fit_ordination <- function(otu,
                               env,
                               sel,
                               covariates,
                               standardize = TRUE,
                               galaxy_colnum = TRUE) {
  
  otu <- as.data.frame(otu, check.names = FALSE)
  env <- as.data.frame(env, check.names = FALSE)
  
  if (is.null(rownames(otu)) || is.null(colnames(otu))) {
    stop("Step 2 Error: `otu` must have row names as taxa IDs and column names as sample IDs.")
  }
  if (is.null(rownames(env))) {
    stop("Step 2 Error: `env` must have row names as sample IDs.")
  }
  
  # Align samples
  samp <- intersect(colnames(otu), rownames(env))
  if (length(samp) == 0) {
    stop("Step 2 Error: no shared samples between OTU columns and ENV row names")
  }
  otu <- otu[, samp, drop = FALSE]
  env <- env[samp, , drop = FALSE]
  
  # Standardize
  env_sd <- if (isTRUE(standardize)) {
    vegan::decostand(env, method = "standardize", MARGIN = 2)
  } else {
    env
  }
  
  # Helper: Index resolution
  to_idx_kai <- function(x, cn, galaxy_colnum = TRUE) {
    if (is.null(x) || length(x) == 0) return(integer(0))
    
    if (is.character(x)) {
      if (length(x) == 1 && grepl(",", x, fixed = TRUE)) {
        x <- strsplit(x, ",", fixed = TRUE)[[1]]
        x <- trimws(x)
      }
      if (all(grepl("^\\d+$", x))) {
        x <- as.numeric(x)
      } else {
        idx <- match(x, cn)
        if (anyNA(idx)) {
          bad <- x[is.na(idx)]
          stop(sprintf("Step 2: unknown variable name(s): %s", paste(bad, collapse = ", ")))
        }
        return(idx)
      }
    }
    
    if (!is.numeric(x)) stop("Step 2: selection must be variable names or numeric indices.")
    
    if (isTRUE(galaxy_colnum) && any(x == 1)) {
      stop("Step 2: wrong selection. Column 1 is SampleID; please choose another column.")
    }
    if (isTRUE(galaxy_colnum)) {
      x <- x - 1
    }
    
    x <- as.integer(x)
    if (any(x < 1) || any(x > length(cn))) {
      stop("Step 2: selection index out of range after Galaxy adjustment (-1).")
    }
    unique(x)
  }
  
  cn <- colnames(env_sd)
  sel_idx <- to_idx_kai(sel, cn, galaxy_colnum = galaxy_colnum)
  cov_idx <- to_idx_kai(covariates, cn, galaxy_colnum = galaxy_colnum)
  
  if (length(sel_idx) == 0) {
    stop("Step 2: `sel` is empty. Provide at least one constrained variable for partial-CCA width.")
  }
  if (length(cov_idx) == 0) {
    stop("Step 2: `covariates` is empty. Kai Step2 requires covariates for position CCA and partial-CCA.")
  }
  if (length(intersect(sel_idx, cov_idx)) > 0) {
    stop("Step 2: the same variable cannot be in both `sel` and `covariates`.")
  }
  
  # Community matrix: samples as rows for vegan
  comm_t <- t(otu)
  
  # 1. Width axis: partial-CCA
  # Formula conceptually: comm ~ sel | covariates
  C_width <- vegan::cca(
    comm_t,
    env_sd[, sel_idx, drop = FALSE],
    env_sd[, cov_idx, drop = FALSE]
  )
  
  # 2. Position axis: CCA
  # Formula conceptually: comm ~ covariates
  C_pos <- vegan::cca(
    comm_t,
    env_sd[, cov_idx, drop = FALSE]
  )
  
  # Extract site scores
  site_width_mat <- C_width$CCA$wa
  site_pos_mat   <- C_pos$CCA$wa
  
  # Align order to `samp` (though likely already aligned, this is safe)
  site_width_mat <- site_width_mat[samp, , drop = FALSE]
  site_pos_mat   <- site_pos_mat[samp, , drop = FALSE]
  
  # Position axis 1 as named vector
  site_pos <- site_pos_mat[, 1]
  names(site_pos) <- rownames(site_pos_mat)
  
  list(
    site_width   = site_width_mat,  # Matrix: all axes (for width calculation)
    site_pos     = site_pos,        # Vector: axis 1 (for position calculation)
    site_pos_mat = site_pos_mat,    # Matrix: all axes
    cca          = C_pos,
    pcca         = C_width
  )
}


# -------------------------------------------------------------------
# Step 3: species niche width & optimum + position-width plot
# -------------------------------------------------------------------

#' Calculate Species Niche Width and Position (Step 3)
#'
#' Computes the niche width and niche position for each species using the site scores
#' from Step 2. Generates a plot of Niche Width vs. Niche Position.
#'
#' @param otu OTU/species matrix (rows = taxa, columns = samples).
#' @param site_width Site scores for width axis (usually matrix from partial CCA).
#' @param site_pos Site scores for position axis (usually vector from CCA axis 1).
#' @param method Smoothing method for the plot ("lm" or "loess").
#' @param make_plot Logical; if \code{TRUE}, returns a ggplot object.
#' @param top_node Integer; maximum number of most abundant species to use for 
#'   loess smoothing calculation (default 10000). Points are plotted for all valid species.
#'
#' @return A list containing:
#' \describe{
#'   \item{species}{Data frame of species traits (NichePosition, NicheWidth).}
#'   \item{plot}{The ggplot object (if \code{make_plot = TRUE}).}
#' }
#' @export
cca_calc_species <- function(otu,
                             site_width,
                             site_pos,
                             method = c("lm", "loess"),
                             make_plot = TRUE,
                             top_node = 10000) {
  method <- match.arg(method)
  
  otu <- as.data.frame(otu, check.names = FALSE)
  if (is.null(rownames(otu)) || is.null(colnames(otu))) {
    stop("Step 3: `otu` must have row names as taxa IDs and column names as sample IDs.")
  }
  
  # Normalize site_width
  if (is.matrix(site_width) || is.data.frame(site_width)) {
    sw_df <- as.data.frame(site_width, check.names = FALSE)
    if (is.null(rownames(sw_df))) {
      stop("Step 3: `site_width` matrix/data.frame must have row names as sample IDs.")
    }
  } else {
    if (is.null(names(site_width))) {
      stop("Step 3: `site_width` must be a named vector OR a matrix/data.frame with rownames.")
    }
    sw_df <- data.frame(WA = as.numeric(site_width),
                        row.names = names(site_width),
                        check.names = FALSE)
  }
  
  # Normalize site_pos
  if (is.matrix(site_pos) || is.data.frame(site_pos)) {
    sp_df <- as.data.frame(site_pos, check.names = FALSE)
    if (is.null(rownames(sp_df))) {
      stop("Step 3: `site_pos` matrix/data.frame must have row names as sample IDs.")
    }
    sp_vec <- sp_df[[1]]
    names(sp_vec) <- rownames(sp_df)
  } else {
    if (is.null(names(site_pos))) {
      stop("Step 3: `site_pos` must be a named vector OR a matrix/data.frame with rownames.")
    }
    sp_vec <- as.numeric(site_pos)
    names(sp_vec) <- names(site_pos)
  }
  
  # Align samples
  samp <- Reduce(intersect, list(colnames(otu), rownames(sw_df), names(sp_vec)))
  if (length(samp) == 0) {
    stop("Step 3: sample IDs do not match across inputs.")
  }
  
  otu <- otu[, samp, drop = FALSE]
  
  # Prepare site score matrices/vectors aligned to samples
  sc_width <- sw_df[samp, , drop = FALSE]
  sc_width <- as.matrix(sc_width)
  
  # Helper functions for weighted calculations
  sc_pos <- sp_vec[samp]
  
  wa <- function(b, a) {
    sumb <- sum(b)
    if (sumb == 0) return(0)
    sum(a * b) / sumb
  }
  
  wa.sd3 <- function(b, a) {
    sumb <- sum(b)
    if (sumb == 0) return(0)
    
    bprop <- b / sumb
    denom <- (sumb * (1 - sum(bprop^2)))
    if (denom == 0) return(0)
    
    bm <- b / denom
    
    a <- as.matrix(a)  
    wam <- sum(a * bm) / sum(bm)  
    
    niche <- sqrt(sum((a - wam)^2 * bm))
    round(niche, 4)
  }
  
  # Compute traits for all species
  # Columns of t(otu) are species; as.list iterates over species
  niche_width <- sapply(as.list(as.data.frame(t(otu))), wa.sd3, a = sc_width)
  niche_pos   <- sapply(as.list(as.data.frame(t(otu))), wa,     a = sc_pos)
  
  niche_width <- as.numeric(niche_width)
  names(niche_width) <- rownames(otu)
  niche_pos <- as.numeric(niche_pos)
  names(niche_pos) <- rownames(otu)
  
  # Filter invalid results
  n <- which(niche_width > 0)
  if (length(n) == 0) {
    stop("No results for your data. (All NicheWidth <= 0)")
  }
  sp_keep <- names(n)
  
  temp.data1 <- data.frame(
    NichePosition = niche_pos[sp_keep],
    NicheWidth    = niche_width[sp_keep],
    check.names = FALSE
  )
  rownames(temp.data1) <- sp_keep
  
  # Subset for smoothing if method is loess
  top_n <- min(length(sp_keep), top_node)
  top.grp <- names(sort(rowSums(otu[sp_keep, , drop = FALSE]), decreasing = TRUE))[1:top_n]
  
  if (method == "loess") {
    temp.data <- temp.data1[top.grp, , drop = FALSE]
  } else {
    temp.data <- temp.data1
  }
  
  # Plot generation
  mm <- max(temp.data1$NicheWidth, na.rm = TRUE)
  x_anno <- as.numeric(stats::quantile(temp.data1$NichePosition, 0.05, na.rm = TRUE))
  y_anno <- mm - 0.02 * diff(range(temp.data1$NicheWidth, na.rm = TRUE))
  
  if (method == "lm") {
    reg <- summary(stats::lm(temp.data[, 2] ~ temp.data[, 1]))
    rp.r <- reg$adj.r.squared
    rp.p <- reg$coefficients[1, 4]  # Kai script uses intercept P
  } else {
    fit <- stats::loess(temp.data1[, 2] ~ temp.data1[, 1])  # Kai fits on all points
    reg.cor <- stats::cor.test(temp.data1[, 2], stats::fitted(fit))
    rp.r <- (reg.cor$estimate)^2
    rp.p <- reg.cor$p.value
  }
  
  p <- NULL
  if (isTRUE(make_plot)) {
    col_point_fill <- "#0EA5E9"
    col_point_edge <- "#FFFFFF"
    col_smooth     <- "#1E40AF"
    col_ci_fill    <- "#93C5FD"
    col_text       <- "#111827"
    grid_major     <- "#E5E7EB"
    
    label <- paste0(
      "R\u00B2 = ", round(rp.r, 4),
      ",  P ", ifelse(rp.p < 0.001, "< 0.001", paste0("= ", round(rp.p, 3)))
    )
    
    p <- ggplot2::ggplot(data = temp.data, ggplot2::aes(x = NichePosition, y = NicheWidth)) +
      ggplot2::geom_point(
        data = temp.data1,
        ggplot2::aes(x = NichePosition, y = NicheWidth),
        shape = 21, size = 3.8, stroke = 0.9,
        fill = col_point_fill, colour = col_point_edge, alpha = 0.95
      ) +
      ggplot2::geom_smooth(
        method = method, se = TRUE,
        colour = col_smooth, fill = col_ci_fill,
        size = 1.6, alpha = 0.35
      ) +
      ggplot2::xlab("Niche Position") +
      ggplot2::ylab("Niche Width") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(colour = col_text),
        legend.title = ggplot2::element_text(size = 18, face = "bold"),
        legend.text  = ggplot2::element_text(size = 15, face = "bold"),
        axis.title.x = ggplot2::element_text(size = 20, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 20, face = "bold"),
        axis.text.x  = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
        axis.text.y  = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
        panel.grid.major = ggplot2::element_line(colour = grid_major, linewidth = 0.4),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        axis.line = ggplot2::element_line(linewidth = 1.0, colour = "#0F172A"),
        axis.line.x = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        plot.margin = grid::unit(rep(0.25, 4), "line")
      ) +
      ggplot2::annotate(
        "text",
        x = x_anno, y = y_anno,
        label = label,
        size = 5, fontface = "bold",
        colour = col_smooth, hjust = 0
      )
  }
  
  list(species = temp.data1, plot = p)
}

# -------------------------------------------------------------------
# Step 4: environmental gradient vs niche position (sites)
# -------------------------------------------------------------------

#' Plot Site Niche Position vs Environmental Gradient (Step 4)
#'
#' Evaluates the relationship between the environmental gradient and the 
#' aggregated niche position of sites.
#'
#' @param env A sample-by-environment data.frame (rows = samples).
#' @param site_pos Named numeric vector of site scores on the position axis.
#' @param var Environmental variable to plot against (name or index).
#' @param make_plot Logical; if \code{TRUE}, returns a ggplot object.
#' @param galaxy_colnum Logical; for numeric indices, whether to treat them as 1-based indices.
#'
#' @return A list containing:
#' \describe{
#'   \item{data}{Data frame used for plotting (ENV, NichePosition).}
#'   \item{plot}{The ggplot object.}
#' }
#' @export
cca_calc_gradient <- function(env,
                              site_pos,
                              var,
                              make_plot = TRUE,
                              galaxy_colnum = TRUE) {
  env <- as.data.frame(env, check.names = FALSE)
  
  samp <- intersect(rownames(env), names(site_pos))
  if (length(samp) == 0) stop("Step 4 Error: sample IDs do not match across inputs.")
  
  env <- env[samp, , drop = FALSE]
  pos <- site_pos[samp]
  
  # Resolve variable
  cn <- colnames(env)
  
  resolve_var_idx <- function(v, cn, galaxy_colnum = TRUE) {
    if (is.null(v) || length(v) == 0) stop("Step 4 Error: `var` is required.")
    
    if (is.character(v)) {
      v <- trimws(v)
      if (grepl("^\\d+$", v)) {
        v <- as.numeric(v)
      } else {
        idx <- match(v, cn)
        if (is.na(idx)) stop("Step 4 Error: specified environmental variable not found.")
        return(idx)
      }
    }
    
    if (!is.numeric(v)) stop("Step 4 Error: `var` must be a name or numeric index.")
    
    if (isTRUE(galaxy_colnum) && any(v == 1)) {
      stop("Step 4 Error: Wrong selection. Column 1 is SampleID.")
    }
    
    idx <- as.integer(v)
    if (isTRUE(galaxy_colnum)) idx <- idx - 1
    
    if (any(idx < 1) || any(idx > length(cn))) {
      stop("Step 4 Error: variable index out of range.")
    }
    idx[1] 
  }
  
  idx <- resolve_var_idx(var, cn, galaxy_colnum = galaxy_colnum)
  var_name <- cn[idx]
  
  df <- data.frame(
    ENV = env[, idx],
    NichePosition = as.numeric(pos),
    row.names = samp,
    check.names = FALSE
  )
  
  p <- NULL
  if (isTRUE(make_plot)) {
    col_point_fill <- "#0EA5E9"
    col_point_edge <- "#FFFFFF"
    col_text       <- "#111827"
    grid_major     <- "#E5E7EB"
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = NichePosition, y = ENV)) +
      ggplot2::geom_point(
        shape = 21, size = 4.0, stroke = 1.0,
        fill = col_point_fill, colour = col_point_edge, alpha = 0.95
      ) +
      ggplot2::ylab(paste(var_name, "gradient")) +
      ggplot2::xlab(paste(var_name, "niche position")) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(colour = col_text),
        legend.title = ggplot2::element_text(size = 18, face = "bold"),
        legend.text  = ggplot2::element_text(size = 15, face = "bold"),
        axis.title.x = ggplot2::element_text(size = 20, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 20, face = "bold"),
        axis.text.x  = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
        axis.text.y  = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
        panel.grid.major = ggplot2::element_line(colour = grid_major, linewidth = 0.4),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        axis.line = ggplot2::element_line(linewidth = 1.0, colour = "#0F172A"),
        axis.line.x = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        plot.margin = grid::unit(rep(0.2, 4), "line")
      )
  }
  
  list(data = df, plot = p)
}

# -------------------------------------------------------------------
# Step 5: group/sample niche width summary
# -------------------------------------------------------------------

#' Calculate Niche Width Summary by Group (Step 5)
#'
#' Computes aggregated niche width statistics at the sample level and group level.
#' Supports two calculation modes via \code{choice}: "all" (global context) 
#' or group-specific context.
#'
#' @param otu OTU/species matrix (rows = taxa, columns = samples).
#' @param site_width Site scores for width axis (usually matrix).
#' @param group A named vector or factor defining groups. Names must match sample IDs.
#' @param choice Calculation mode. If \code{"all"}, computes niche width using all samples. 
#'   Otherwise, computes within-group niche widths.
#'
#' @return A list containing:
#' \describe{
#'   \item{sample_summary}{Data frame of sample-level statistics (ID, AverageValue, StandardDeviation/Error).}
#'   \item{group_summary}{Data frame of group-level statistics.}
#' }
#' @export
cca_calc_group <- function(otu,
                           site_width,
                           group,
                           choice = "all") {
  
  otu <- as.data.frame(otu, check.names = FALSE)
  if (is.null(rownames(otu)) || is.null(colnames(otu))) {
    stop("Step 5 Error: `otu` must have row names (taxa) and column names (samples).")
  }
  
  # ---- normalize site_width to a data.frame (rows = samples, cols = axes) ----
  if (is.matrix(site_width) || is.data.frame(site_width)) {
    sw_df <- as.data.frame(site_width, check.names = FALSE)
    if (is.null(rownames(sw_df))) {
      stop("Step 5 Error: `site_width` must have row names.")
    }
  } else {
    # fallback: named vector
    if (is.null(names(site_width))) {
      stop("Step 5 Error: `site_width` must be named.")
    }
    sw_df <- data.frame(WA = as.numeric(site_width),
                        row.names = names(site_width),
                        check.names = FALSE)
  }
  
  if (is.null(names(group))) {
    stop("Step 5 Error: `group` must be a named vector (names = sample IDs).")
  }
  
  samp_comm <- colnames(otu)
  
  # Consistency checks
  if (length(group) != length(samp_comm)) {
    stop("Please make group file and community file consistent.")
  }
  if (!all(samp_comm %in% names(group))) {
    stop("Please make group file and community file consistent.")
  }
  if (!all(samp_comm %in% rownames(sw_df))) {
    stop("Step 5: site scores (site_width) and community samples are inconsistent.")
  }
  
  # Align
  group_comm <- group[samp_comm]
  sc <- as.matrix(sw_df[samp_comm, , drop = FALSE])
  
  # Assign samples to groups
  sites <- unique(as.character(group_comm))
  grp <- setNames(vector("list", length(sites)), sites)
  for (i in seq_along(sites)) {
    grp[[i]] <- samp_comm[which(as.character(group_comm) == sites[i])]
  }
  
  wa.sd3 <- function(b, a) {
    sb <- sum(b)
    if (sb == 0) return(0)
    
    bm <- b / (sb * (1 - sum((b / sb)^2)))
    
    a <- as.matrix(a)                 
    wam <- sum(a * bm) / sum(bm)     
    niche <- sqrt(sum((a - wam)^2 * bm))
    
    round(niche, 4)
  }
  
  outSample <- NULL
  outSite <- NULL
  
  if (identical(choice, "all")) {
    
    # 1) Global niche width
    res <- sapply(as.list(as.data.frame(t(otu))), wa.sd3, a = sc)
    niche_all <- as.numeric(res)
    names(niche_all) <- rownames(otu)
    
    # 2) Sample summary
    for (sp in seq_along(samp_comm)) {
      sub_sp <- rownames(otu)[which(otu[, sp] > 0)]
      nic <- niche_all[sub_sp]
      n <- which(nic > 0)
      
      if (length(n) != 0) {
        mu <- mean(nic[n], na.rm = TRUE)
        sd_val <- stats::sd(nic[n], na.rm = TRUE)  
        outSample <- rbind(outSample, cbind(mu, sd_val))
      } else {
        stop("Please check with your data.")
      }
    }
    rownames(outSample) <- samp_comm
    
    # 3) Group summary
    for (i in seq_along(grp)) {
      gname <- names(grp)[i]
      outSite <- rbind(outSite, cbind(mean(outSample[grp[[i]], 1]),
                                      stats::sd(outSample[grp[[i]], 1])))
      rownames(outSite)[nrow(outSite)] <- gname
    }
    
    outSample <- data.frame(
      ID = rownames(outSample),
      AverageValue = outSample[, 1],
      StandardDeviation = outSample[, 2],  
      row.names = NULL,
      check.names = FALSE
    )
    outSite <- data.frame(
      ID = rownames(outSite),
      AverageValue = outSite[, 1],
      StandardDeviation = outSite[, 2],
      row.names = NULL,
      check.names = FALSE
    )
    
  } else {
    
    # Group-wise calculation
    res <- matrix(0, nrow = nrow(otu), ncol = length(grp))
    for (z in seq_along(grp)) {
      sub_samples <- grp[[z]]
      sub_data <- otu[, sub_samples, drop = FALSE]
      res[, z] <- sapply(
        as.list(as.data.frame(t(sub_data))),
        wa.sd3,
        a = sc[colnames(sub_data), , drop = FALSE]   #关键：保留多列轴
      )
    }
    colnames(res) <- names(grp)
    rownames(res) <- rownames(otu)
    niche_sites <- res
    
    for (i in seq_along(grp)) {
      gname <- names(grp)[i]
      nic <- niche_sites[, gname]
      n <- which(nic > 0)
      
      if (length(n) != 0) {
        tmp_mean <- rep(mean(nic[n]), length(grp[[i]]))
        tmp_sd   <- rep(stats::sd(nic[n]), length(grp[[i]]))
        
        outSample <- rbind(outSample, cbind(tmp_mean, tmp_sd))
        outSite   <- rbind(outSite,   cbind(mean(nic[n]), stats::sd(nic[n])))
        rownames(outSite)[nrow(outSite)] <- gname
      }
    }
    
    rownames(outSample) <- as.vector(unlist(grp))
    
    outSample <- data.frame(
      ID = rownames(outSample),
      AverageValue = outSample[, 1],
      StandardError = outSample[, 2],  
      row.names = NULL,
      check.names = FALSE
    )
    outSite <- data.frame(
      ID = rownames(outSite),
      AverageValue = outSite[, 1],
      StandardError = outSite[, 2],
      row.names = NULL,
      check.names = FALSE
    )
  }
  
  list(sample_summary = outSample, group_summary = outSite)
}

# -------------------------------------------------------------------
# Step 6: sample-level niche width vs environmental gradient
# -------------------------------------------------------------------

#' Plot Sample/Group Niche Width vs Environmental Gradient (Step 6)
#'
#' Visualizes the relationship between environmental gradients and niche width
#' at the sample or group level.
#'
#' @param env A sample-by-environment data.frame.
#' @param sample_summary Output from \code{\link{cca_calc_group}} (step 5).
#' @param var Environmental variable to plot against (name or index).
#' @param group Named vector/factor of group labels (required for group plots).
#' @param plot_type Type of plot: "sample", "group", or "both".
#' @param method Smoothing method ("lm" or "loess").
#' @param make_plot Logical; if \code{TRUE}, returns ggplot object(s).
#' @param galaxy_colnum Logical; for numeric indices, whether to treat them as 1-based indices.
#' @param show_ci Logical; whether to show confidence intervals on the smooth line.
#' @param annotate_stats Logical; whether to annotate R-squared and P-value on the plot.
#'
#' @return A list containing data frames and ggplot objects for samples and/or groups.
#' @export
cca_plot_group_env <- function(env,
                               sample_summary,
                               var,
                               group = NULL,
                               plot_type = c("sample", "group", "both"),
                               method = c("lm", "loess"),
                               make_plot = TRUE,
                               galaxy_colnum = TRUE,
                               show_ci = TRUE,
                               annotate_stats = TRUE) {
  plot_type <- match.arg(plot_type)
  method <- match.arg(method)
  
  env <- as.data.frame(env, check.names = FALSE)
  sample_summary <- as.data.frame(sample_summary, check.names = FALSE)
  
  # Normalize ID column to rownames
  if ("ID" %in% colnames(sample_summary)) {
    rownames(sample_summary) <- as.character(sample_summary$ID)
  } else if (is.null(rownames(sample_summary))) {
    stop("Step 6 Error: `sample_summary` must have an `ID` column or row names.")
  }
  
  # Identify columns
  if (!("AverageValue" %in% colnames(sample_summary))) {
    stop("Step 6 Error: `sample_summary` must contain `AverageValue`.")
  }
  mean_col <- "AverageValue"
  
  if ("StandardDeviation" %in% colnames(sample_summary)) {
    sd_col <- "StandardDeviation"
  } else if ("StandardError" %in% colnames(sample_summary)) {
    sd_col <- "StandardError"
  } else {
    stop("Step 6 Error: `sample_summary` must contain SD/SE column.")
  }
  
  # Align samples
  samp <- intersect(rownames(env), rownames(sample_summary))
  if (length(samp) == 0) stop("Step 6 Error: sample IDs do not match.")
  
  env <- env[samp, , drop = FALSE]
  nicheSample <- sample_summary[samp, c(mean_col, sd_col), drop = FALSE]
  
  # Resolve variable
  cn <- colnames(env)
  
  resolve_var_idx <- function(v, cn, galaxy_colnum = TRUE) {
    if (is.null(v) || length(v) == 0) stop("Step 6 Error: `var` is required.")
    
    if (is.character(v)) {
      v <- trimws(v)
      if (grepl("^\\d+$", v)) {
        v <- as.numeric(v)
      } else {
        idx <- match(v, cn)
        if (is.na(idx)) stop("Step 6 Error: variable not found.")
        return(idx)
      }
    }
    
    if (!is.numeric(v)) stop("Step 6 Error: `var` must be a name or numeric index.")
    
    if (isTRUE(galaxy_colnum) && any(v == 1)) {
      stop("Step 6 Error: Column 1 is SampleID.")
    }
    
    idx <- as.integer(v)
    if (isTRUE(galaxy_colnum)) idx <- idx - 1
    
    if (any(idx < 1) || any(idx > length(cn))) {
      stop("Step 6 Error: variable index out of range.")
    }
    idx[1]
  }
  
  idx <- resolve_var_idx(var, cn, galaxy_colnum = galaxy_colnum)
  var_name <- cn[idx]
  
  combined <- cbind(env[, idx, drop = FALSE], nicheSample)
  colnames(combined) <- c("ENV", "NicheWidth", "NicheWdithSD") 
  rownames(combined) <- samp
  df_sample <- as.data.frame(combined, check.names = FALSE)
  
  col_point_fill <- "#F59E0B"  
  col_point_edge <- "#FFFFFF"  
  col_smooth    <- "#6D28D9"  
  col_ci_fill    <- "#A78BFA"  
  
  add_stats_label <- function(p, dfit, xcol, ycol, method, col_text) {
    dfit <- dfit[stats::complete.cases(dfit[, c(xcol, ycol)]), , drop = FALSE]
    if (!isTRUE(annotate_stats) || nrow(dfit) < 3) return(p)
    
    if (method == "lm") {
      fit <- stats::lm(dfit[[ycol]] ~ dfit[[xcol]])
      sm  <- summary(fit)
      r2  <- sm$adj.r.squared
      pvl <- sm$coefficients[2, 4]
    } else {
      fit <- stats::loess(dfit[[ycol]] ~ dfit[[xcol]])
      ct  <- stats::cor.test(dfit[[ycol]], stats::fitted(fit))
      r2  <- (ct$estimate)^2
      pvl <- ct$p.value
    }
    
    label <- paste0(
      "R^2 = ", round(r2, 4),
      ", P ", ifelse(pvl < 0.001, "< 0.001", paste0("= ", round(pvl, 3)))
    )
    
    x_anno <- stats::quantile(dfit[[xcol]], 0.05, na.rm = TRUE)
    y_anno <- max(dfit[[ycol]], na.rm = TRUE) * 0.98
    
    p + ggplot2::annotate(
      "text", x = x_anno, y = y_anno,
      label = label, size = 5, fontface = "bold",
      colour = col_text, hjust = 0
    )
  }
  
  plot_sample <- NULL
  if (isTRUE(make_plot) && plot_type %in% c("sample", "both")) {
    plot_sample <- ggplot2::ggplot(df_sample, ggplot2::aes(x = ENV, y = NicheWidth)) +
      ggplot2::geom_point(
        shape = 21, size = 3.8, stroke = 0.9,
        fill = col_point_fill, colour = col_point_edge, alpha = 0.95
      ) +
      ggplot2::geom_smooth(
        method = method, se = isTRUE(show_ci),
        colour = col_smooth, fill = col_ci_fill,
        linewidth = 1.6, alpha = 0.35
      ) +
      ggplot2::xlab(paste(var_name, "gradient")) +
      ggplot2::ylab("Niche width within each sites") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 20, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 20, face = "bold"),
        axis.text.x  = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
        axis.text.y  = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
        axis.line    = ggplot2::element_line(linewidth = 1.2),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin = grid::unit(rep(0.2, 4), "line")
      )
    
    plot_sample <- add_stats_label(plot_sample, df_sample, "ENV", "NicheWidth", method, col_smooth)
  }
  
  df_group <- NULL
  plot_group <- NULL
  
  if (plot_type %in% c("group", "both")) {
    if (is.null(group)) {
      stop("Step 6: `group` is required when plot_type is 'group' or 'both'.")
    }
    if (is.null(names(group))) {
      stop("Step 6: `group` must be a named vector/factor (names are sample IDs).")
    }
    
    group2 <- as.factor(group[samp])
    if (anyNA(group2)) stop("Step 6: some samples in env/sample_summary are missing in `group` names.")
    
    df_tmp <- data.frame(
      Sample = samp,
      Group = as.character(group2),
      ENV = df_sample$ENV,
      NicheWidth = df_sample$NicheWidth,
      stringsAsFactors = FALSE
    )
    
    env_group <- stats::aggregate(ENV ~ Group, data = df_tmp, FUN = function(x) mean(x, na.rm = TRUE))
    width_group <- stats::aggregate(NicheWidth ~ Group, data = df_tmp, FUN = function(x) mean(x, na.rm = TRUE))
    
    df_group <- merge(env_group, width_group, by = "Group")
    df_group <- df_group[order(df_group$ENV), , drop = FALSE]
    
    if (isTRUE(make_plot)) {
      plot_group <- ggplot2::ggplot(df_group, ggplot2::aes(x = ENV, y = NicheWidth)) +
        ggplot2::geom_point(
          shape = 21, size = 4.2, stroke = 1.0,
          fill = col_point_fill, colour = col_point_edge, alpha = 0.95
        ) +
        ggplot2::geom_smooth(
          method = method, se = isTRUE(show_ci),
          colour = col_smooth, fill = col_ci_fill,
          linewidth = 1.6, alpha = 0.35
        ) +
        ggplot2::xlab(paste0("Environmental gradient (group mean): ", var_name)) +
        ggplot2::ylab("Group mean niche width") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 20, face = "bold"),
          axis.title.y = ggplot2::element_text(size = 20, face = "bold"),
          axis.text.x  = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
          axis.text.y  = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
          axis.line    = ggplot2::element_line(linewidth = 1.2),
          panel.grid.minor = ggplot2::element_blank(),
          plot.margin = grid::unit(rep(0.2, 4), "line")
        )
      
      plot_group <- add_stats_label(plot_group, df_group, "ENV", "NicheWidth", method, col_smooth)
    }
  }
  
  if (plot_type == "sample") {
    return(list(data = df_sample, plot = plot_sample))
  }
  if (plot_type == "group") {
    return(list(data = df_group, plot = plot_group))
  }
  # both
  list(
    data_sample = df_sample,
    plot_sample = plot_sample,
    data_group  = df_group,
    plot_group  = plot_group
  )
}

# -------------------------------------------------------------------
# Two end-to-end workflows: 1-2-3-4 and 1-2-5-6
# -------------------------------------------------------------------

#' Execute Gradient-Based Workflow (Steps 1-4)
#'
#' Runs the complete environmental gradient analysis pipeline:
#' 1. Prepare env data (PCA)
#' 2. Fit CCA/Partial-CCA
#' 3. Calculate species niche traits
#' 4. Analyze gradient vs niche position
#'
#' @inheritParams cca_prep_env
#' @inheritParams cca_fit_ordination
#' @inheritParams cca_calc_species
#' @inheritParams cca_calc_gradient
#' @param var Environmental variable for gradient plotting (passed to Step 4).
#'
#' @return A named list containing outputs of steps 1 through 4.
#' @export
cca_workflow_gradient <- function(otu,
                                  env,
                                  sel,
                                  covariates,
                                  standardize = TRUE,
                                  method = c("lm", "loess"),
                                  var = sel[1],
                                  galaxy_colnum = TRUE,
                                  make_plot = TRUE,
                                  top_node = 10000) {
  
  method <- match.arg(method)
  
  step1 <- cca_prep_env(
    env,
    sel = sel,
    constrain = covariates,
    standardize = standardize,
    galaxy_colnum = galaxy_colnum
  )
  
  step2 <- cca_fit_ordination(
    otu, env,
    sel = sel,
    covariates = covariates,
    standardize = standardize,
    galaxy_colnum = galaxy_colnum
  )
  
  step3 <- cca_calc_species(
    otu,
    site_width = step2$site_width,  
    site_pos   = step2$site_pos,    
    method     = method,
    make_plot  = make_plot,
    top_node   = top_node
  )
  
  step4 <- cca_calc_gradient(
    env,
    site_pos = step2$site_pos,
    var = var,
    make_plot = make_plot,
    galaxy_colnum = galaxy_colnum
  )
  
  list(step1 = step1, step2 = step2, step3 = step3, step4 = step4)
}

#' Execute Group-Based Workflow (Steps 1-2, 5-6)
#'
#' Runs the complete group-based analysis pipeline:
#' 1. Prepare env data (PCA)
#' 2. Fit CCA/Partial-CCA
#' 5. Calculate group niche widths
#' 6. Plot group/sample niche width vs gradient
#'
#' @inheritParams cca_prep_env
#' @inheritParams cca_fit_ordination
#' @inheritParams cca_calc_group
#' @inheritParams cca_plot_group_env
#' @param var Environmental variable for gradient plotting (passed to Step 6).
#'
#' @return A named list containing outputs of steps 1, 2, 5, and 6.
#' @export
cca_workflow_group <- function(otu,
                               env,
                               sel,
                               group,
                               covariates,
                               standardize = TRUE,
                               method = c("lm", "loess"),
                               var = sel[1],
                               choice = "all",
                               galaxy_colnum = TRUE,
                               make_plot = TRUE,
                               plot_type = c("sample", "group", "both"),
                               show_ci = TRUE,
                               annotate_stats = TRUE) {
  
  method <- match.arg(method)
  plot_type <- match.arg(plot_type)
  
  step1 <- cca_prep_env(
    env,
    sel = sel,
    constrain = covariates,
    standardize = standardize,
    galaxy_colnum = galaxy_colnum
  )
  
  step2 <- cca_fit_ordination(
    otu, env,
    sel = sel,
    covariates = covariates,
    standardize = standardize,
    galaxy_colnum = galaxy_colnum
  )
  
  step5 <- cca_calc_group(
    otu = otu,
    site_width = step2$site_width,  
    group = group,
    choice = choice
  )
  
  step6 <- cca_plot_group_env(
    env = env,
    sample_summary = step5$sample_summary,
    var = var,
    group = group,                  
    plot_type = plot_type,          
    method = method,                
    make_plot = make_plot,
    galaxy_colnum = galaxy_colnum,
    show_ci = show_ci,
    annotate_stats = annotate_stats
  )
  
  list(step1 = step1, step2 = step2, step5 = step5, step6 = step6)
}

