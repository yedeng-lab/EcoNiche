# ===========================
# Niche width (Levins / Shannon): samples as states
# ===========================

#' Calculate Niche Width (Samples as States)
#'
#' Computes OTU-level niche breadth by treating each sample as a discrete state.
#' For OTU \eqn{i}, let \eqn{p_{is}} be the proportional abundance of OTU \eqn{i} in sample \eqn{s}
#' (i.e., abundance of OTU \eqn{i} normalized to sum to 1 across all samples).
#'
#' Two niche-width indices are supported via \code{method}:
#' \itemize{
#'   \item \strong{Levins breadth (Levins)}:
#'     \deqn{B_i = 1 / \sum_s p_{is}^2}
#'     Optionally, a standardized breadth ranging from 0 to 1 is returned:
#'     \deqn{B_i' = (B_i - 1)/(K - 1)}
#'     where \eqn{K} is the total number of samples (states).
#'   \item \strong{Shannon breadth (Shannon)}:
#'     \deqn{H_i = - \sum_s p_{is}\log(p_{is})}
#'     Terms with \eqn{p_{is}=0} are excluded from the sum to avoid \eqn{\log(0)}.
#' }
#'
#' @param otu An OTU-by-sample matrix or data frame (rows = OTUs, columns = samples).
#'   Values should be non-negative (counts or relative abundance).
#' @param env Optional sample metadata data frame with row names as sample IDs.
#'   If provided, samples are aligned by the intersection of \code{colnames(otu)}
#'   and \code{rownames(env)}.
#' @param min_occ Minimum number of samples with abundance > 0 required to keep an OTU.
#'   Default is 3.
#' @param min_abund Minimum total abundance required to keep an OTU (sum across samples).
#'   Default is 5.
#' @param standardize Logical; if \code{TRUE} and \code{method} includes \code{"levins"},
#'   returns the standardized breadth \code{levins_Bstd} in addition to the raw Levins breadth.
#'   Default is \code{TRUE}.
#' @param method Character; which niche width index to compute. One of
#'   \code{"levins"}, \code{"shannon"}, or \code{"both"}. Default is \code{"levins"}.
#'
#' @return A data frame with one row per OTU. Columns include:
#' \describe{
#'   \item{OTU}{OTU identifier.}
#'   \item{n_states}{Number of states (samples) used in calculation.}
#'   \item{n_samples}{Number of samples where the OTU is present.}
#'   \item{total_abund}{Total abundance of the OTU across samples.}
#'   \item{levins_B}{Raw Levins niche breadth. Present if \code{method} is \code{"levins"} or \code{"both"}.}
#'   \item{levins_Bstd}{Standardized Levins niche breadth (if \code{standardize=TRUE}).
#'     Present if \code{method} is \code{"levins"} or \code{"both"}.}
#'   \item{shannon_H}{Shannon niche breadth (entropy). Present if \code{method} is \code{"shannon"} or \code{"both"}.}
#' }
#'
#' @export
niche_width_calc <- function(otu,
                        env = NULL,
                        min_occ   = 3L,
                        min_abund = 5,
                        standardize = TRUE,
                        method = c("levins", "shannon", "both")) {
  
  method <- match.arg(method)
  
  otu <- as.data.frame(otu, check.names = FALSE)
  
  if (is.null(rownames(otu)) || is.null(colnames(otu))) {
    stop("`otu` must have row names as OTU IDs and column names as sample IDs.")
  }
  
  # Optional: Align samples if `env` is provided; otherwise use `otu` columns
  if (!is.null(env)) {
    env <- as.data.frame(env, check.names = FALSE)
    if (is.null(rownames(env))) stop("`env` must have row names as sample IDs.")
    samp <- intersect(colnames(otu), rownames(env))
    if (length(samp) < 5L) stop("Too few overlapping samples between `otu` and `env`.")
    otu <- otu[, samp, drop = FALSE]
  } else {
    samp <- colnames(otu)
  }
  
  K <- length(samp)
  if (K < 2L && standardize) {
    stop("Need at least 2 samples to compute standardized Levins width.")
  }
  
  # OTU filtering
  occ  <- rowSums(otu > 0, na.rm = TRUE)
  tot  <- rowSums(otu,     na.rm = TRUE)
  keep <- (occ >= min_occ) & (tot >= min_abund)
  if (!any(keep)) stop("No OTUs remain after filtering. Try lowering `min_occ` / `min_abund`.")
  comm2 <- otu[keep, , drop = FALSE]
  
  levins_list <- lapply(seq_len(nrow(comm2)), function(i) {
    y <- as.numeric(comm2[i, ])
    y[is.na(y)] <- 0
    
    s <- sum(y)
    if (!is.finite(s) || s <= 0) return(NULL)
    
    p <- y / s
    p[!is.finite(p)] <- 0
    
    out <- list(
      OTU         = rownames(comm2)[i],
      n_states    = K,
      n_samples   = sum(y > 0),
      total_abund = s
    )
    
    # Levins (same as spaa: 1/sum(p^2)) :contentReference[oaicite:1]{index=1}
    if (method %in% c("levins", "both")) {
      B  <- 1 / sum(p^2)
      Bn <- if (standardize && K > 1) (B - 1) / (K - 1) else NA_real_
      out$levins_B    <- B
      out$levins_Bstd <- Bn
    }
    
    # Shannon (same as spaa: -sum(p*log(p)), but must use p>0) :contentReference[oaicite:2]{index=2}
    if (method %in% c("shannon", "both")) {
      p_pos <- p[p > 0]
      H <- -sum(p_pos * log(p_pos))
      out$shannon_H <- H
    }
    
    as.data.frame(out, stringsAsFactors = FALSE)
  })
  
  levins_list <- levins_list[!vapply(levins_list, is.null, logical(1))]
  if (!length(levins_list)) stop("All OTUs were excluded during calculations.")
  
  res <- do.call(rbind, levins_list)
  rownames(res) <- NULL
  res
}

# ===========================
# Levins niche width along a composite axis (e.g., CCA1): binned states
# ===========================

#' Calculate Levins Niche Width Along a Gradient (Binned States)
#'
#' Computes OTU-level Levins niche breadth along a continuous composite axis (e.g., CCA1) 
#' by discretizing the axis into bins (states). Abundance is first aggregated within bins 
#' (mean or sum), converted to within-bin proportions \eqn{p_{ij}}, and then used to compute:
#' \deqn{B_i = 1 / \sum_j p_{ij}^2}
#' A standardized breadth is also reported:
#' \deqn{B_i' = (B_i - 1)/(K - 1)}
#' where \eqn{K} is the number of bins.
#'
#' @param otu An OTU-by-sample matrix or data frame.
#' @param env A sample-by-environment data frame with row names as sample IDs.
#' @param axis_var Character string; column name in \code{env} containing the axis values.
#' @param nbin Integer; number of bins used to discretize the axis. Default is 8.
#' @param bin_method Character string; binning strategy. \code{"equal_freq"} (quantile-based) 
#'   or \code{"equal_width"}.
#' @param agg_fun Character string; aggregation method within bins. \code{"mean"} (recommended) 
#'   or \code{"sum"}.
#' @param otu_mode Character string; input data type. \code{"auto"} (detect), \code{"count"}, 
#'   or \code{"relative"}. If \code{"count"}, columns are normalized to relative abundance 
#'   before aggregation to reduce sequencing depth bias.
#' @param min_occ Minimum number of samples with abundance > 0 required to keep an OTU.
#' @param min_abund Minimum total abundance required to keep an OTU.
#'
#' @return A data frame with one row per OTU and the following columns:
#' \describe{
#'   \item{OTU}{OTU identifier.}
#'   \item{axis}{Name of the gradient axis used.}
#'   \item{n_states}{Number of bins (K).}
#'   \item{levins_B}{Raw Levins niche breadth.}
#'   \item{levins_Bstd}{Standardized Levins niche breadth.}
#'   \item{n_samples}{Number of samples where the OTU is present.}
#'   \item{total_abund}{Total abundance of the OTU.}
#' }
#' @export
levins_calc_binned <- function(otu,
                               env,
                               axis_var,                 # 例如 "CCA1" 或你加到 env 的综合轴列名
                               nbin = 8L,
                               bin_method = c("equal_freq", "equal_width"),
                               agg_fun = c("mean", "sum"),
                               otu_mode = c("auto", "count", "relative"),
                               min_occ   = 3L,
                               min_abund = 5) {
  bin_method <- match.arg(bin_method)
  agg_fun    <- match.arg(agg_fun)
  otu_mode   <- match.arg(otu_mode)
  
  otu <- as.data.frame(otu, check.names = FALSE)
  env <- as.data.frame(env, check.names = FALSE)
  
  if (is.null(rownames(otu)) || is.null(colnames(otu))) {
    stop("`otu` must have row names as OTU IDs and column names as sample IDs.")
  }
  if (is.null(rownames(env))) stop("`env` must have row names as sample IDs.")
  if (!(axis_var %in% colnames(env))) stop("`axis_var` was not found in `env`: ", axis_var)
  
  # Align samples
  samp <- intersect(colnames(otu), rownames(env))
  if (length(samp) < 5L) stop("Too few overlapping samples between `otu` and `env`.")
  otu <- otu[, samp, drop = FALSE]
  env <- env[samp, , drop = FALSE]
  
  x <- as.numeric(env[[axis_var]])
  ok <- is.finite(x)
  x <- x[ok]
  otu <- otu[, ok, drop = FALSE]
  samp <- samp[ok]
  
  if (length(x) < 5L) stop("Too few valid samples after removing NA in axis.")
  
  # Automatically detect otu_mode
  if (otu_mode == "auto") {
    vals <- as.numeric(as.matrix(otu))
    vals <- vals[is.finite(vals)]
    frac_integer <- mean(abs(vals - round(vals)) < 1e-10)
    max_val <- max(vals, na.rm = TRUE)
    # If mostly integers and values > 1.5, assume counts; otherwise relative
    otu_mode <- if (frac_integer > 0.98 && max_val > 1.5) "count" else "relative"
  }
  
  # For counts: normalize columns to relative abundance first 
  if (otu_mode == "count") {
    col_tot <- colSums(otu, na.rm = TRUE)
    otu <- sweep(otu, 2, col_tot, "/")
    otu[!is.finite(as.matrix(otu))] <- 0
  }
  
  # Binning: Equal frequency is preferred 
  if (bin_method == "equal_freq") {
    brks <- stats::quantile(x, probs = seq(0, 1, length.out = nbin + 1L), na.rm = TRUE, type = 7)
    brks <- unique(as.numeric(brks))
    if (length(brks) < 3L) stop("Not enough unique axis values for binning. Reduce `nbin`.")
    bin <- cut(x, breaks = brks, include.lowest = TRUE, right = TRUE)
  } else {
    brks <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = nbin + 1L)
    bin <- cut(x, breaks = brks, include.lowest = TRUE, right = TRUE)
  }
  levs <- levels(bin)
  K <- length(levs)
  
  # OTU filtering: Filter by occurrence and total abundance
  occ  <- rowSums(otu > 0, na.rm = TRUE)
  tot  <- rowSums(otu,     na.rm = TRUE)
  keep <- (occ >= min_occ) & (tot >= min_abund)
  if (!any(keep)) stop("No OTUs remain after filtering. Try lowering thresholds.")
  comm2 <- otu[keep, , drop = FALSE]
  
  levins_list <- lapply(seq_len(nrow(comm2)), function(i) {
    y <- as.numeric(comm2[i, ])
    s <- sum(y, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(NULL)
    
    df_tmp <- data.frame(state = bin, abund = y)
    
    agg <- if (agg_fun == "mean") {
      tapply(df_tmp$abund, df_tmp$state, mean, na.rm = TRUE)
    } else {
      tapply(df_tmp$abund, df_tmp$state, sum, na.rm = TRUE)
    }
    # Fill all bins to ensure fixed K for comparability
    agg_full <- setNames(rep(0, K), levs)
    agg_full[names(agg)] <- agg
    agg_full[is.na(agg_full)] <- 0
    
    A_i <- sum(agg_full)
    if (!is.finite(A_i) || A_i <= 0) return(NULL)
    
    p <- agg_full / A_i
    
    B  <- 1 / sum(p^2)
    Bn <- if (K > 1) (B - 1) / (K - 1) else NA_real_
    
    data.frame(
      OTU          = rownames(comm2)[i],
      axis         = axis_var,
      n_states     = K,
      levins_B     = B,
      levins_Bstd  = Bn,
      n_samples    = sum(y > 0, na.rm = TRUE),
      total_abund  = s,
      stringsAsFactors = FALSE
    )
  })
  
  levins_list <- levins_list[!vapply(levins_list, is.null, logical(1))]
  if (!length(levins_list)) stop("All OTUs were excluded during Levins calculations.")
  
  res <- do.call(rbind, levins_list)
  rownames(res) <- NULL
  res
}


#' Calculate Community-Mean Levins Niche Width Along a Gradient
#'
#' Computes the abundance-weighted mean Levins width for each sample, using OTU-level 
#' Levins widths and sample-level OTU abundances. Optionally plots the relationship 
#' between this community metric and an environmental gradient.
#'
#' The community-mean width for sample \eqn{s} is:
#' \deqn{B_s = \sum_i Q_{is} B_i'}
#' where \eqn{Q_{is}} is the relative abundance of OTU \eqn{i} in sample \eqn{s} 
#' and \eqn{B_i'} is the standardized Levins width (\code{levins_Bstd}).
#'
#' @param otu An OTU-by-sample matrix or data frame.
#' @param env A sample-by-environment data frame.
#' @param levins_df A data frame of OTU-level Levins results containing at least
#'   columns \code{OTU} and \code{width_col}.
#' @param grad The environmental gradient to plot against (column name or index in \code{env}).
#' @param width_col Character string; column name in \code{levins_df} to use as the width metric.
#'   Default is \code{"levins_Bstd"}.
#' @param method Character string; smoothing method for the trend line. \code{"lm"} or \code{"loess"}.
#' @param make_plot Logical; whether to return a ggplot object.
#'
#' @return A list containing:
#' \describe{
#'   \item{data}{Data frame with columns Sample, ENV, and CommLevinsWidth.}
#'   \item{plot}{The ggplot object (if \code{make_plot=TRUE}).}
#' }
#' @export
levins_calc_group <- function(otu,
                                  env,
                                  levins_df,
                                  grad,
                                  width_col  = "levins_Bstd",
                                  method     = c("lm", "loess"),
                                  make_plot  = TRUE) {
  method <- match.arg(method)
  
  otu       <- as.data.frame(otu,       check.names = FALSE)
  env       <- as.data.frame(env,       check.names = FALSE)
  levins_df <- as.data.frame(levins_df, check.names = FALSE)
  
  if (is.null(rownames(otu)) || is.null(colnames(otu))) {
    stop("`otu` must have row names as OTU IDs and column names as sample IDs.")
  }
  if (is.null(rownames(env))) {
    stop("`env` must have row names as sample IDs.")
  }
  if (!("OTU" %in% colnames(levins_df))) {
    stop("`levins_df` must contain a column named 'OTU'.")
  }
  if (!(width_col %in% colnames(levins_df))) {
    stop("Column '", width_col, "' was not found in `levins_df`.")
  }
  
  # Align samples
  samp <- intersect(colnames(otu), rownames(env))
  if (length(samp) < 5L) {
    stop("Too few overlapping samples between `otu` and `env`. Check sample IDs.")
  }
  otu <- otu[, samp, drop = FALSE]
  env <- env[samp, , drop = FALSE]
  
  # Parse gradient
  if (is.character(grad)) {
    idx <- match(grad, colnames(env))
    if (is.na(idx)) stop("Gradient variable '", grad, "' was not found in `env`.")
    x <- env[, idx]
    grad_name <- colnames(env)[idx]
  } else {
    idx <- as.integer(grad)
    if (is.na(idx) || idx < 1L || idx > ncol(env)) {
      stop("`grad` is neither a valid column name nor a valid column index.")
    }
    x <- env[, idx]
    grad_name <- colnames(env)[idx]
  }
  x <- as.numeric(x)
  
  # Align OTUs
  common_otus <- intersect(rownames(otu), levins_df$OTU)
  if (!length(common_otus)) {
    stop("No overlapping OTU IDs between `otu` and `levins_df`.")
  }
  
  otu_sub <- otu[common_otus, , drop = FALSE]
  lev2    <- levins_df[match(common_otus, levins_df$OTU), , drop = FALSE]
  
  width_vec <- as.numeric(lev2[[width_col]])
  names(width_vec) <- lev2$OTU
  width_vec <- width_vec[rownames(otu_sub)]
  
  # --------- KEY FIX ----------
  # Remove OTUs with invalid widths from BOTH numerator and denominator.
  keep_w <- is.finite(width_vec) & !is.na(width_vec) & (width_vec >= 0)
  
  # If standardized width, optionally keep within [0,1] (tolerate tiny numeric overshoot)
  if (grepl("Bstd", width_col, fixed = TRUE)) {
    keep_w <- keep_w & (width_vec <= 1 + 1e-8)
  }
  
  if (!any(keep_w)) {
    stop("All overlapping OTUs have invalid/non-finite widths in `", width_col, "`.")
  }
  
  otu_sub  <- otu_sub[keep_w, , drop = FALSE]
  width_vec <- width_vec[keep_w]
  # ----------------------------
  
  num <- colSums(otu_sub * width_vec, na.rm = TRUE)
  den <- colSums(otu_sub,            na.rm = TRUE)
  
  B_comm <- ifelse(den > 0, num / den, NA_real_)
  
  df <- data.frame(
    Sample          = samp,
    ENV             = x,
    CommLevinsWidth = B_comm,
    stringsAsFactors = FALSE
  )
  
  p1 <- NULL
  if (make_plot) {
    df_plot <- df[is.finite(df$ENV) & is.finite(df$CommLevinsWidth), ]
    if (!nrow(df_plot)) stop("No valid samples are available for plotting.")
    
    # Calculate R^2 and P-value
    if (method == "lm") {
      reg <- summary(stats::lm(CommLevinsWidth ~ ENV, data = df_plot))
      rp.r <- reg$adj.r.squared
      rp.p <- reg$coefficients[2, 4]
    } else {
      fit <- stats::loess(CommLevinsWidth ~ ENV, data = df_plot)
      reg.cor <- stats::cor.test(df_plot$CommLevinsWidth, fitted(fit))
      rp.r <- (reg.cor$estimate)^2
      rp.p <- reg.cor$p.value
    }
    
    mm     <- max(df_plot$CommLevinsWidth, na.rm = TRUE)
    x_anno <- as.numeric(stats::quantile(df_plot$ENV, 0.05, na.rm = TRUE))
    y_anno <- mm - 0.02 * diff(range(df_plot$CommLevinsWidth, na.rm = TRUE))
    
    # --- Plot Aesthetics ---
    col_point_fill <- "#0EA5E9"
    col_point_edge <- "#FFFFFF"
    col_smooth     <- "#1E40AF"
    col_ci_fill    <- "#93C5FD"
    col_text       <- "#111827"
    grid_major     <- "#E5E7EB"
    # ------------------------------------
    
    methodReg <- if (method == "lm") "lm" else "loess"
    
    p1 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = ENV, y = CommLevinsWidth)) +
      ggplot2::geom_point(
        shape = 21, size = 4.0, stroke = 1.0,
        fill = col_point_fill, colour = col_point_edge, alpha = 0.95
      ) +
      ggplot2::geom_smooth(
        method = methodReg, se = TRUE,
        colour = col_smooth, fill = col_ci_fill,
        linewidth = 1.6, alpha = 0.35
      ) +
      ggplot2::xlab(paste0(grad_name, " gradient")) +
      ggplot2::ylab(paste0("Community-mean Levins niche width (", grad_name, ")")) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text         = ggplot2::element_text(colour = col_text),
        axis.title.x = ggplot2::element_text(size = 20, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 20, face = "bold"),
        axis.text.x  = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
        axis.text.y  = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
        panel.grid.major = ggplot2::element_line(colour = grid_major, linewidth = 0.4),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        axis.line.x = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank()
      ) +
      ggplot2::annotate(
        "text", x = x_anno, y = y_anno,
        label = paste0(
          "R^2 = ", round(rp.r, 4),
          ", P ", ifelse(rp.p < 0.001, "< 0.001", paste0("= ", round(rp.p, 3)))
        ),
        size = 5, fontface = "bold", colour = col_smooth, hjust = 0
      )
  }
  
  list(data = df, plot = p1)
}

#' Plot Species Niche Position vs Levins Width
#'
#' Merges species niche optima (from any method such as CCA, GAM, or weighted average)
#' with standardized Levins niche width, and plots the position-width relationship
#' showing the correlation (R-squared) and P-value.
#'
#' @param pos_df A data frame containing species niche positions. Must include \code{id_col} and \code{pos_col}.
#' @param levins_df A data frame containing Levins widths. Must include \code{id_col} and \code{width_col}.
#' @param id_col Character string; column name of the species/OTU ID shared by both tables (e.g., \code{"OTU"}).
#' @param pos_col Character string; column name of the niche position in \code{pos_df}. Default \code{"NichePosition"}.
#' @param width_col Character string; column name of the Levins width in \code{levins_df}. Default \code{"levins_Bstd"}.
#' @param method Character string; fitting method for the trend line. \code{"lm"} or \code{"loess"}.
#' @param make_plot Logical; whether to return a ggplot object.
#'
#' @return A list containing:
#' \describe{
#'   \item{data}{Merged data frame used for plotting.}
#'   \item{plot}{The ggplot object (if \code{make_plot=TRUE}).}
#' }
#' @export
levins_plot_pos_width <- function(pos_df,
                                  levins_df,
                                  id_col,
                                  pos_col   = "NichePosition",
                                  width_col = "levins_Bstd",
                                  method    = c("lm", "loess"),
                                  make_plot = TRUE) {
  method    <- match.arg(method)
  pos_df    <- as.data.frame(pos_df,    check.names = FALSE)
  levins_df <- as.data.frame(levins_df, check.names = FALSE)
  
  for (nm in c(id_col, pos_col)) {
    if (!(nm %in% colnames(pos_df))) {
      stop("pos_df is missing required column '", nm, "'.")
    }
  }
  for (nm in c(id_col, width_col)) {
    if (!(nm %in% colnames(levins_df))) {
      stop("levins_df is missing required column '", nm, "'.")
    }
  }
  
  df <- merge(
    pos_df[, c(id_col, pos_col), drop = FALSE],
    levins_df[, c(id_col, width_col,
                  intersect(c("total_abund", "n_samples"), colnames(levins_df))), drop = FALSE],
    by = id_col
  )
  
  if (!nrow(df)) stop("No overlap between tables on '", id_col, "'.")
  
  colnames(df)[colnames(df) == pos_col]   <- "NichePosition"
  colnames(df)[colnames(df) == width_col] <- "LevinsWidth"
  
  rownames(df) <- df[[id_col]]
  
  # --------- KEY FIX ----------
  # Allow LevinsWidth == 0 (valid for standardized Levins width).
  df <- df[
    is.finite(df$NichePosition) &
      is.finite(df$LevinsWidth) &
      df$LevinsWidth >= 0,
    , drop = FALSE
  ]
  
  # If standardized width, keep within [0,1] (tolerate tiny numeric overshoot)
  if (grepl("Bstd", width_col, fixed = TRUE)) {
    df <- df[df$LevinsWidth <= 1 + 1e-8, , drop = FALSE]
  }
  # ----------------------------
  
  if (!nrow(df)) stop("No valid NichePosition / LevinsWidth values remain after filtering.")
  
  # Optionally keep top 10000 by total abundance
  df_top <- df
  if (nrow(df) > 10000 && "total_abund" %in% colnames(df)) {
    ord <- order(df$total_abund, decreasing = TRUE)
    df_top <- df[ord[seq_len(10000)], , drop = FALSE]
  }
  
  mm     <- max(df$LevinsWidth, na.rm = TRUE)
  x_anno <- as.numeric(stats::quantile(df$NichePosition, 0.05, na.rm = TRUE))
  y_anno <- mm - 0.02 * diff(range(df$LevinsWidth, na.rm = TRUE))
  
  if (method == "lm") {
    reg <- summary(stats::lm(LevinsWidth ~ NichePosition, data = df_top))
    rp.r <- reg$adj.r.squared
    rp.p <- reg$coefficients[2, 4]
    methodReg <- "lm"
  } else {
    fit <- stats::loess(LevinsWidth ~ NichePosition, data = df_top)
    reg.cor <- stats::cor.test(df_top$LevinsWidth, fitted(fit))
    rp.r <- (reg.cor$estimate)^2
    rp.p <- reg.cor$p.value
    methodReg <- "loess"
  }
  
  p1 <- NULL
  if (make_plot) {
    # --- Plot Aesthetics ---
    col_point_fill <- "#0EA5E9"
    col_point_edge <- "#FFFFFF"
    col_smooth     <- "#1E40AF"
    col_ci_fill    <- "#93C5FD"
    col_text       <- "#111827"
    grid_major     <- "#E5E7EB"
    # ------------------------------------
    
    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = NichePosition, y = LevinsWidth)) +
      ggplot2::geom_point(
        shape = 21, size = 3.8, stroke = 0.9,
        fill = col_point_fill, colour = col_point_edge, alpha = 0.95
      ) +
      ggplot2::geom_smooth(
        method = methodReg, se = TRUE,
        colour = col_smooth, fill = col_ci_fill,
        linewidth = 1.6, alpha = 0.35
      ) +
      ggplot2::xlab("Niche position") +
      ggplot2::ylab("Levins niche width") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text         = ggplot2::element_text(colour = col_text),
        axis.title.x = ggplot2::element_text(size = 20, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 20, face = "bold"),
        axis.text.x  = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
        axis.text.y  = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
        panel.grid.major = ggplot2::element_line(colour = grid_major, linewidth = 0.4),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        axis.line.x = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank()
      ) +
      ggplot2::annotate(
        "text", x = x_anno, y = y_anno,
        label = paste0(
          "R^2 = ", round(rp.r, 4),
          ", P ", ifelse(rp.p < 0.001, "< 0.001", paste0("= ", round(rp.p, 3)))
        ),
        size = 5, fontface = "bold", colour = col_smooth, hjust = 0
      )
  }
  
  list(data = df, plot = p1)
}
