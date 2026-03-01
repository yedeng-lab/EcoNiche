#' Fit GAM-Based Niche Curves for All OTUs
#'
#' Fits a Generalized Additive Model (GAM) response curve for each OTU along a single 
#' environmental gradient. This function handles both count data (using Negative Binomial 
#' or Poisson families) and relative abundance data (using logit-transformed Gaussian models).
#' It estimates the niche optimum and the 50\% niche breadth (breadth50) for each taxon.
#'
#' @param otu An OTU-by-sample matrix or data frame. Rows should be OTUs/taxa and 
#'   columns should be samples. Values can be counts or relative abundances.
#' @param env A sample-by-environment data frame. Rows must be samples and columns 
#'   must be environmental variables.
#' @param env_var A character string specifying the column name in \code{env} to use 
#'   as the environmental gradient.
#' @param data_type Character string specifying the data type. One of \code{"auto"} 
#'   (heuristic detection), \code{"count"}, or \code{"relative"}.
#' @param count_family Character string specifying the error distribution family for 
#'   count data. One of \code{"nb"} (Negative Binomial, default) or \code{"poisson"}.
#' @param use_offset Logical; whether to use an offset term (log library size) in the 
#'   GAM for count data. Default is \code{TRUE}.
#' @param lib_size Optional named numeric vector of library sizes (sequencing depth) 
#'   for each sample. If \code{NULL} and \code{data_type="count"}, column sums of 
#'   \code{otu} are used.
#' @param min_prev Minimum prevalence threshold (proportion of samples with abundance > 0). 
#'   OTUs below this threshold are skipped. Default is 0.10.
#' @param min_total Minimum total abundance threshold. Relevant for count data. 
#'   Default is 100.
#' @param min_mean Minimum mean abundance threshold. Relevant for relative abundance data. 
#'   Default is 1e-5.
#' @param k_spline Integer; the upper limit on the spline basis dimension for the 
#'   smooth term \code{s(env, k = ...)}. Default is 5.
#' @param n_grid Integer; number of grid points along the gradient used to estimate 
#'   optimum and niche breadth. Default is 200.
#' @param verbose Logical; whether to print progress messages. Default is \code{TRUE}.
#'
#' @return A data frame with one row per OTU and the following columns:
#' \describe{
#'   \item{OTU}{OTU identifier.}
#'   \item{n_nonzero}{Number of samples with non-zero abundance.}
#'   \item{prevalence}{Proportion of samples with non-zero abundance.}
#'   \item{total_counts}{Total counts (for count data) or NA.}
#'   \item{r2}{Adjusted R-squared (for relative abundance models) or NA.}
#'   \item{dev_expl}{Proportion of deviance explained.}
#'   \item{edf}{Effective degrees of freedom of the smooth term.}
#'   \item{F_stat}{Test statistic (F or Chi-sq) for the smooth term.}
#'   \item{p_value}{Approximate p-value for the smooth term.}
#'   \item{optimum_env}{Environmental value at the peak of the fitted curve.}
#'   \item{env50_min}{Lower environmental bound where fitted abundance >= 50\% of peak.}
#'   \item{env50_max}{Upper environmental bound where fitted abundance >= 50\% of peak.}
#'   \item{breadth50}{Niche breadth (env50_max - env50_min).}
#' }
#' @importFrom mgcv gam predict.gam nb
#' @importFrom stats median poisson predict
#' @export
gam_fit_model <- function(otu,
                          env,
                          env_var,
                          data_type  = c("auto", "count", "relative"),
                          # count 模式参数
                          count_family = c("nb", "poisson"),
                          use_offset   = TRUE,
                          lib_size     = NULL,      # For count mode: named vector; if NULL use colSums(otu)
                          # 通用过滤参数
                          min_prev   = 0.10,
                          min_total  = 100,         # Relevant for count mode
                          min_mean   = 1e-5,        # Relevant for relative mode
                          # GAM参数
                          k_spline = 5,
                          n_grid   = 200,
                          verbose  = TRUE) {
  data_type <- match.arg(data_type)
  count_family <- match.arg(count_family)
  
  otu <- as.data.frame(otu, check.names = FALSE)
  env <- as.data.frame(env, check.names = FALSE)
  
  if (is.null(rownames(otu)) || is.null(colnames(otu))) {
    stop("`otu` must have row names as OTU IDs and column names as sample IDs.")
  }
  if (is.null(rownames(env))) {
    stop("`env` must have row names as sample IDs.")
  }
  if (!(env_var %in% colnames(env))) {
    stop("Environmental variable '", env_var, "' was not found in `env`.")
  }
  
  # Align samples
  common_samples <- intersect(colnames(otu), rownames(env))
  if (length(common_samples) < 10) {
    stop("Fewer than 10 overlapping samples between `otu` and `env`. Check sample IDs.")
  }
  otu <- otu[, common_samples, drop = FALSE]
  env <- env[common_samples, , drop = FALSE]
  env_vec_full <- env[common_samples, env_var, drop = TRUE]
  
  # Automatically detect data_type (conservative and interpretable)
  if (data_type == "auto") {
    vals <- as.numeric(as.matrix(otu))
    vals <- vals[is.finite(vals)]
    if (!length(vals)) stop("No finite values found in `otu`.")
    
    frac_integer <- mean(abs(vals - round(vals)) < 1e-10)
    max_val <- max(vals, na.rm = TRUE)
    
    # Heuristic: if mostly integers and values > 1.5, assume counts.
    data_type <- if (frac_integer > 0.98 && max_val > 1.5) "count" else "relative"
  }
  
  if (verbose) message("Detected/selected data_type = ", data_type)
  
  # Count mode: prepare lib_size
  if (data_type == "count") {
    if (is.null(lib_size)) {
      lib_size <- colSums(otu, na.rm = TRUE)
    } else {
      if (is.null(names(lib_size))) stop("`lib_size` must be a named numeric vector.")
      lib_size <- lib_size[common_samples]
    }
    if (any(!is.finite(lib_size)) || any(lib_size <= 0)) {
      stop("`lib_size` contains non-finite or non-positive values.")
    }
  }
  
  # Helper: Safely extract statistics from s.table to avoid index errors
  .extract_smooth_stats <- function(sm) {
    s_table <- sm$s.table
    
    if (is.null(s_table)) {
      return(list(edf = NA_real_, stat = NA_real_, p = NA_real_))
    }
    
    # Ensure it is an object with dimensions
    if (is.null(dim(s_table))) {
      # Rarely degenerates to a vector; return NA to avoid errors
      return(list(edf = NA_real_, stat = NA_real_, p = NA_real_))
    }
    if (nrow(s_table) < 1) {
      return(list(edf = NA_real_, stat = NA_real_, p = NA_real_))
    }
    
    cn <- colnames(s_table)
    
    edf_val <- if (!is.null(cn) && ("edf" %in% cn)) as.numeric(s_table[1, "edf"]) else NA_real_
    
    # Statistic column: Prioritize F, then Chi.sq (common for NB/Poisson), then others
    stat_val <- NA_real_
    if (!is.null(cn)) {
      if ("F" %in% cn) {
        stat_val <- as.numeric(s_table[1, "F"])
      } else if ("Chi.sq" %in% cn) {
        stat_val <- as.numeric(s_table[1, "Chi.sq"])
      } else if ("Chisq" %in% cn) {
        stat_val <- as.numeric(s_table[1, "Chisq"])
      } else if ("chisq" %in% cn) {
        stat_val <- as.numeric(s_table[1, "chisq"])
      }
    }
    
    p_val <- NA_real_
    if (!is.null(cn)) {
      if ("p-value" %in% cn) {
        p_val <- as.numeric(s_table[1, "p-value"])
      } else if ("p.value" %in% cn) {
        p_val <- as.numeric(s_table[1, "p.value"])
      } else if ("Pr(>F)" %in% cn) {
        p_val <- as.numeric(s_table[1, "Pr(>F)"])
      }
    }
    
    list(edf = edf_val, stat = stat_val, p = p_val)
  }
  
  otu_ids  <- rownames(otu)
  res_list <- vector("list", length(otu_ids))
  
  if (verbose) {
    if (data_type == "count") {
      message("Fitting count GAM: y ~ s(env) ",
              if (use_offset) "+ offset(log(lib_size)) " else "",
              "with ", if (count_family == "nb") "Negative Binomial" else "Poisson",
              " family; including zeros.")
    } else {
      message("Fitting relative-abundance GAM: logit(p*) ~ s(env) with Gaussian; including zeros via smoothing.")
    }
  }
  
  fam_obj <- NULL
  if (data_type == "count") {
    fam_obj <- if (count_family == "nb") mgcv::nb() else poisson()
  }
  
  for (i in seq_along(otu_ids)) {
    otu_id <- otu_ids[i]
    y_raw  <- as.numeric(otu[otu_id, ])
    
    # Prevalence filtering
    prev <- mean(y_raw > 0, na.rm = TRUE)
    if (!is.finite(prev) || prev < min_prev) next
    
    keep <- is.finite(y_raw) & is.finite(env_vec_full)
    if (sum(keep) < 10) next
    
    x <- env_vec_full[keep]
    
    if (data_type == "count") {
      y <- y_raw[keep]
      
      tot <- sum(y, na.rm = TRUE)
      if (!is.finite(tot) || tot < min_total) next
      if (sum(y > 0, na.rm = TRUE) < 3) next
      
      dat <- data.frame(y = y, env = x)
      if (use_offset) dat$off <- log(lib_size[keep])
      
      fit <- try({
        if (use_offset) {
          mgcv::gam(y ~ s(env, k = k_spline) + offset(off),
                    data = dat, family = fam_obj, method = "REML")
        } else {
          mgcv::gam(y ~ s(env, k = k_spline),
                    data = dat, family = fam_obj, method = "REML")
        }
      }, silent = TRUE)
      if (inherits(fit, "try-error")) next
      
      sm <- summary(fit)
      st <- .extract_smooth_stats(sm)
      
      env_grid <- seq(min(dat$env), max(dat$env), length.out = n_grid)
      
      # Prediction: Fix a "typical depth" offset for comparability
      if (use_offset) {
        off0 <- log(stats::median(exp(dat$off), na.rm = TRUE))
        pred <- as.numeric(predict(fit,
                                   newdata = data.frame(env = env_grid, off = off0),
                                   type = "response"))
      } else {
        pred <- as.numeric(predict(fit,
                                   newdata = data.frame(env = env_grid),
                                   type = "response"))
      }
      
      if (!all(is.finite(pred)) || max(pred) <= 0) {
        optimum_env <- env50_min <- env50_max <- breadth50 <- NA_real_
      } else {
        optimum_env <- env_grid[which.max(pred)]
        thr <- 0.5 * max(pred)
        idx50 <- which(pred >= thr)
        if (!length(idx50)) {
          env50_min <- env50_max <- breadth50 <- NA_real_
        } else {
          env50_min <- min(env_grid[idx50])
          env50_max <- max(env_grid[idx50])
          breadth50 <- env50_max - env50_min
        }
      }
      
      res_list[[i]] <- data.frame(
        OTU          = otu_id,
        n_nonzero    = sum(dat$y > 0),
        prevalence   = prev,
        total_counts = tot,
        r2           = NA_real_,     # r.sq is not recommended for count mode
        dev_expl     = sm$dev.expl,
        edf          = st$edf,
        F_stat       = st$stat,      # Could be F or Chi.sq; NA if missing
        p_value      = st$p,
        optimum_env  = optimum_env,
        env50_min    = env50_min,
        env50_max    = env50_max,
        breadth50    = breadth50,
        stringsAsFactors = FALSE
      )
      
    } else {
      # Relative abundance mode
      p <- y_raw[keep]
      
      mval <- mean(p, na.rm = TRUE)
      if (!is.finite(mval) || mval < min_mean) next
      
      n <- length(p)
      p_star <- (p * (n - 1) + 0.5) / n
      p_star <- pmin(pmax(p_star, 1e-12), 1 - 1e-12)
      
      y <- log(p_star / (1 - p_star))
      dat <- data.frame(y = y, env = x)
      
      fit <- try(mgcv::gam(y ~ s(env, k = k_spline), data = dat, method = "REML"),
                 silent = TRUE)
      if (inherits(fit, "try-error")) next
      
      sm <- summary(fit)
      st <- .extract_smooth_stats(sm)
      
      env_grid <- seq(min(dat$env), max(dat$env), length.out = n_grid)
      eta <- as.numeric(predict(fit, newdata = data.frame(env = env_grid), type = "link"))
      pred <- 1 / (1 + exp(-eta))
      
      if (!all(is.finite(pred)) || max(pred) <= 0) {
        optimum_env <- env50_min <- env50_max <- breadth50 <- NA_real_
      } else {
        optimum_env <- env_grid[which.max(pred)]
        thr <- 0.5 * max(pred)
        idx50 <- which(pred >= thr)
        if (!length(idx50)) {
          env50_min <- env50_max <- breadth50 <- NA_real_
        } else {
          env50_min <- min(env_grid[idx50])
          env50_max <- max(env_grid[idx50])
          breadth50 <- env50_max - env50_min
        }
      }
      
      res_list[[i]] <- data.frame(
        OTU          = otu_id,
        n_nonzero    = sum(p > 0, na.rm = TRUE),
        prevalence   = prev,
        total_counts = NA_real_,
        r2           = sm$r.sq,
        dev_expl     = sm$dev.expl,
        edf          = st$edf,
        F_stat       = st$stat,   # Usually F here; NA if missing
        p_value      = st$p,
        optimum_env  = optimum_env,
        env50_min    = env50_min,
        env50_max    = env50_max,
        breadth50    = breadth50,
        stringsAsFactors = FALSE
      )
    }
    
    if (verbose && i %% 100 == 0) {
      message("Processed OTU ", i, "/", length(otu_ids))
    }
  }
  
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]
  if (!length(res_list)) {
    warning("No OTUs passed filters or model fitting.")
    return(
      data.frame(
        OTU = character(0),
        n_nonzero = integer(0),
        prevalence = numeric(0),
        total_counts = numeric(0),
        r2 = numeric(0),
        dev_expl = numeric(0),
        edf = numeric(0),
        F_stat = numeric(0),
        p_value = numeric(0),
        optimum_env = numeric(0),
        env50_min = numeric(0),
        env50_max = numeric(0),
        breadth50 = numeric(0)
      )
    )
  }
  
  res <- do.call(rbind, res_list)
  res[order(res$p_value), ]
}

#' Calculate Sample-Level Niche Width Weighted by Abundance
#'
#' Computes the abundance-weighted mean niche breadth for each sample, given an 
#' OTU-by-sample table and OTU-level niche breadth estimates (e.g., from \code{gam_fit_model}).
#' The formula is:
#' \deqn{B_j = \sum_i (abundance_{ij} \times breadth50_i) / \sum_i abundance_{ij}}
#'
#' @param otu An OTU-by-sample matrix or data frame.
#' @param niche_df A data frame containing at least columns for OTU IDs and niche width values.
#' @param otu_col Character string; the column name in \code{niche_df} containing OTU IDs. 
#'   Default is \code{"OTU"}.
#' @param width_col Character string; the column name in \code{niche_df} containing niche width values. 
#'   Default is \code{"breadth50"}.
#' @param weight_mode Character string; how to handle abundance weights. One of \code{"auto"}, 
#'   \code{"counts"}, or \code{"relative"}.
#'
#' @return A data frame with columns \code{sample} and \code{Bw50_abundance_weighted}.
#' @export
gam_calc_sitewidth <- function(otu,
                               niche_df,
                               otu_col   = "OTU",
                               width_col = "breadth50",
                               weight_mode = c("auto", "counts", "relative")) {
  weight_mode <- match.arg(weight_mode)
  
  otu      <- as.data.frame(otu, check.names = FALSE)
  niche_df <- as.data.frame(niche_df, check.names = FALSE)
  
  if (is.null(rownames(otu))) stop("`otu` must have row names as OTU IDs.")
  if (!(otu_col %in% colnames(niche_df))) stop("Column '", otu_col, "' not found in `niche_df`.")
  if (!(width_col %in% colnames(niche_df))) stop("Column '", width_col, "' not found in `niche_df`.")
  
  niche2 <- niche_df[!is.na(niche_df[[width_col]]), , drop = FALSE]
  common_otus <- intersect(rownames(otu), niche2[[otu_col]])
  if (!length(common_otus)) stop("No overlapping OTU IDs between `otu` and `niche_df`.")
  
  otu_sub <- otu[common_otus, , drop = FALSE]
  niche2  <- niche2[match(common_otus, niche2[[otu_col]]), , drop = FALSE]
  
  breadth50 <- niche2[[width_col]]
  names(breadth50) <- niche2[[otu_col]]
  
  # Select weighting mode
  col_tot <- colSums(otu_sub, na.rm = TRUE)
  
  if (weight_mode == "auto") {
    # If column sums are close to 1 (or 100), treat as relative abundance; 
    # otherwise treat as counts and convert to relative weights
    if (median(col_tot, na.rm = TRUE) > 0.9 && median(col_tot, na.rm = TRUE) < 1.1) {
      w <- otu_sub
    } else if (median(col_tot, na.rm = TRUE) > 90 && median(col_tot, na.rm = TRUE) < 110) {
      w <- otu_sub / 100
    } else {
      # counts -> relative weights
      w <- sweep(otu_sub, 2, col_tot, "/")
      w[!is.finite(as.matrix(w))] <- 0
    }
  } else if (weight_mode == "counts") {
    w <- otu_sub
  } else { # "relative"
    w <- sweep(otu_sub, 2, col_tot, "/")
    w[!is.finite(as.matrix(w))] <- 0
  }
  
  num <- colSums(w * breadth50, na.rm = TRUE)
  den <- colSums(w, na.rm = TRUE)
  
  B_sample <- ifelse(den > 0, num / den, NA_real_)
  
  data.frame(
    sample = colnames(otu_sub),
    Bw50_abundance_weighted = B_sample,
    stringsAsFactors = FALSE
  )
}

#' Plot GAM-Based Niche Curve for a Single OTU
#'
#' Fits a GAM response curve for a specific OTU and plots the results.
#' Includes options for confidence intervals and different color palettes.
#'
#' @param otu An OTU-by-sample matrix or data frame.
#' @param env A sample-by-environment data frame.
#' @param env_var Character string; the environmental variable to use as the gradient.
#' @param otu_id Character string; the ID of the OTU to plot (must exist in \code{otu}).
#' @param data_type Character string; \code{"auto"}, \code{"count"}, or \code{"relative"}.
#' @param count_family Character string; \code{"nb"} or \code{"poisson"}.
#' @param use_offset Logical; whether to use library size offset for count data.
#' @param lib_size Optional named vector of library sizes.
#' @param min_mean Minimum mean abundance filter (for relative mode).
#' @param min_prev Minimum prevalence filter.
#' @param k_spline Spline basis dimension.
#' @param n_grid Number of grid points for prediction.
#' @param add_ci Logical; whether to plot the 95\% confidence interval.
#' @param palette Character string; color palette name (\code{"blue"}, \code{"orange"}, 
#'   \code{"green"}, \code{"purple"}, or \code{"viridis"}).
#' @param point_alpha Numeric; transparency alpha for observed data points.
#'
#' @return A list containing:
#' \describe{
#'   \item{data}{Data frame of observed values (env, y).}
#'   \item{grid}{Data frame of fitted values along the gradient.}
#'   \item{fit}{The GAM model object.}
#'   \item{optimum_env}{Estimated niche optimum.}
#'   \item{env50_min}{Lower bound of 50\% niche breadth.}
#'   \item{env50_max}{Upper bound of 50\% niche breadth.}
#'   \item{breadth50}{Niche breadth.}
#'   \item{plot}{The ggplot object.}
#' }
#' @export
gam_plot_species <- gam_plot_species_unified <- function(otu,
                                                         env,
                                                         env_var,
                                                         otu_id,
                                                         data_type = c("auto", "count", "relative"),
                                                         # count 模式
                                                         count_family = c("nb", "poisson"),
                                                         use_offset   = TRUE,
                                                         lib_size     = NULL,     # NULL means colSums(otu)
                                                         # relative 模式
                                                         min_mean     = 1e-5,
                                                         # 通用
                                                         min_prev     = 0.10,
                                                         k_spline     = 5,
                                                         n_grid       = 200,
                                                         add_ci       = TRUE,
                                                         # 配色
                                                         palette = c("blue", "orange", "green", "purple", "viridis"),
                                                         point_alpha = 0.85) {
  
  data_type    <- match.arg(data_type)
  count_family <- match.arg(count_family)
  palette      <- match.arg(palette)
  
  otu <- as.data.frame(otu, check.names = FALSE)
  env <- as.data.frame(env, check.names = FALSE)
  
  if (is.null(rownames(otu)) || is.null(colnames(otu))) {
    stop("`otu` must have row names as OTU IDs and column names as sample IDs.")
  }
  if (is.null(rownames(env))) stop("`env` must have row names as sample IDs.")
  if (!(otu_id %in% rownames(otu))) stop("OTU '", otu_id, "' not found in `otu`.")
  if (!(env_var %in% colnames(env))) stop("Environmental variable '", env_var, "' not found in `env`.")
  
  common_samples <- intersect(colnames(otu), rownames(env))
  if (length(common_samples) < 10) stop("Too few overlapping samples to fit the model.")
  
  otu_sub <- suppressWarnings(as.numeric(as.character(otu[otu_id, common_samples, drop = TRUE])))
  env_vec <- env[common_samples, env_var, drop = TRUE]
  
  keep <- is.finite(otu_sub) & is.finite(env_vec)
  otu_sub <- otu_sub[keep]
  env_vec <- env_vec[keep]
  
  if (length(env_vec) < 10) stop("Too few valid samples after NA removal.")
  
  # Automatically detect type
  if (data_type == "auto") {
    frac_integer <- mean(abs(otu_sub - round(otu_sub)) < 1e-10)
    max_val <- max(otu_sub, na.rm = TRUE)
    data_type <- if (frac_integer > 0.98 && max_val > 1.5) "count" else "relative"
  }
  
  # Color scheme
  cols <- switch(
    palette,
    blue    = list(point="#0EA5E9", line="#1E40AF", ci="#93C5FD", vline="#EF4444", edge="#FFFFFF"),
    orange  = list(point="#FB923C", line="#C2410C", ci="#FDBA74", vline="#2563EB", edge="#FFFFFF"),
    green   = list(point="#34D399", line="#047857", ci="#A7F3D0", vline="#DC2626", edge="#FFFFFF"),
    purple  = list(point="#A78BFA", line="#6D28D9", ci="#DDD6FE", vline="#DC2626", edge="#FFFFFF"),
    viridis = list(point=NULL,      line=NULL,      ci=NULL,      vline="#DC2626", edge="#FFFFFF")
  )
  
  # Build grid
  env_grid <- seq(min(env_vec), max(env_vec), length.out = n_grid)
  
  if (data_type == "count") {
    # Counts: filtering thresholds (optional)
    prev <- mean(otu_sub > 0, na.rm = TRUE)
    if (prev < min_prev) stop("This OTU prevalence is below min_prev; plot not fitted.")
    
    tot <- sum(otu_sub, na.rm = TRUE)
    if (!is.finite(tot) || tot <= 0) stop("This OTU has zero total counts; plot not fitted.")
    
    fam_obj <- if (count_family == "nb") mgcv::nb() else poisson()
    
    if (use_offset) {
      if (is.null(lib_size)) {
        # Use total column sums as sequencing depth (consistent with samples)
        lib_size <- colSums(otu[, common_samples, drop=FALSE], na.rm = TRUE)
      } else {
        if (is.null(names(lib_size))) stop("`lib_size` must be a named vector.")
        lib_size <- lib_size[common_samples]
      }
      lib_size <- lib_size[keep]
      if (any(!is.finite(lib_size)) || any(lib_size <= 0)) stop("Invalid lib_size values.")
      off <- log(lib_size)
      dat <- data.frame(env = env_vec, y = otu_sub, off = off)
      
      fit <- mgcv::gam(y ~ s(env, k = k_spline) + offset(off),
                       data = dat, family = fam_obj, method = "REML")
      
      # Prediction: Fix at "typical sequencing depth"
      off0 <- log(stats::median(exp(off), na.rm = TRUE))
      pred <- predict(fit,
                      newdata = data.frame(env = env_grid, off = off0),
                      se.fit = add_ci,
                      type = "link")
      # link scale -> response scale
      eta <- as.numeric(pred$fit)
      se  <- if (add_ci) as.numeric(pred$se.fit) else NULL
      
      mu  <- as.numeric(predict(fit, newdata=data.frame(env=env_grid, off=off0), type="response"))
      mu_u <- mu_l <- NULL
      if (add_ci) {
        # Approx: +/- 1.96*se on link scale, then back to response
        eta_u <- eta + 1.96*se
        eta_l <- eta - 1.96*se
        # For NB/Poisson, link is log by default
        mu_u <- exp(eta_u)
        mu_l <- exp(eta_l)
      }
      
      y_label <- "Expected count (typical depth)"
      grid_df <- data.frame(env=env_grid, fit=mu,
                            fit_upper=if (add_ci) mu_u else NA_real_,
                            fit_lower=if (add_ci) mu_l else NA_real_)
    } else {
      dat <- data.frame(env = env_vec, y = otu_sub)
      fit <- mgcv::gam(y ~ s(env, k = k_spline),
                       data = dat, family = fam_obj, method = "REML")
      
      pred <- predict(fit,
                      newdata = data.frame(env = env_grid),
                      se.fit = add_ci,
                      type = "link")
      
      eta <- as.numeric(pred$fit)
      se  <- if (add_ci) as.numeric(pred$se.fit) else NULL
      mu  <- as.numeric(predict(fit, newdata=data.frame(env=env_grid), type="response"))
      mu_u <- mu_l <- NULL
      if (add_ci) {
        mu_u <- exp(eta + 1.96*se)
        mu_l <- exp(eta - 1.96*se)
      }
      
      y_label <- "Expected count"
      grid_df <- data.frame(env=env_grid, fit=mu,
                            fit_upper=if (add_ci) mu_u else NA_real_,
                            fit_lower=if (add_ci) mu_l else NA_real_)
    }
    
    # optimum / breadth50 (50% of relative max)）
    optimum_env <- env_grid[which.max(grid_df$fit)]
    thr <- 0.5 * max(grid_df$fit, na.rm = TRUE)
    idx50 <- which(grid_df$fit >= thr)
    env50_min <- if (length(idx50)) min(env_grid[idx50]) else NA_real_
    env50_max <- if (length(idx50)) max(env_grid[idx50]) else NA_real_
    breadth50 <- if (length(idx50)) (env50_max - env50_min) else NA_real_
    
    plot_dat <- data.frame(env = dat$env, y = dat$y)
    
  } else {
    # relative: logit(p*) + Gaussian
    prev <- mean(otu_sub > 0, na.rm = TRUE)
    if (prev < min_prev) stop("This OTU prevalence is below min_prev; plot not fitted.")
    mval <- mean(otu_sub, na.rm = TRUE)
    if (mval < min_mean) stop("This OTU mean abundance is below min_mean; plot not fitted.")
    
    p <- otu_sub
    n <- length(p)
    p_star <- (p * (n - 1) + 0.5) / n
    p_star <- pmin(pmax(p_star, 1e-12), 1 - 1e-12)
    y <- log(p_star / (1 - p_star))
    
    dat <- data.frame(env = env_vec, y = y)
    fit <- mgcv::gam(y ~ s(env, k = k_spline), data = dat, method = "REML")
    
    pred <- predict(fit,
                    newdata = data.frame(env = env_grid),
                    se.fit = add_ci,
                    type = "link")
    eta <- as.numeric(pred$fit)
    se  <- if (add_ci) as.numeric(pred$se.fit) else NULL
    
    inv_logit <- function(z) 1 / (1 + exp(-z))
    
    fit_p <- inv_logit(eta)
    up_p  <- lo_p <- NULL
    if (add_ci) {
      up_p <- inv_logit(eta + 1.96*se)
      lo_p <- inv_logit(eta - 1.96*se)
    }
    
    grid_df <- data.frame(env=env_grid, fit=fit_p,
                          fit_upper=if (add_ci) up_p else NA_real_,
                          fit_lower=if (add_ci) lo_p else NA_real_)
    
    y_label <- "Relative abundance (fitted)"
    
    optimum_env <- env_grid[which.max(grid_df$fit)]
    thr <- 0.5 * max(grid_df$fit, na.rm = TRUE)
    idx50 <- which(grid_df$fit >= thr)
    env50_min <- if (length(idx50)) min(env_grid[idx50]) else NA_real_
    env50_max <- if (length(idx50)) max(env_grid[idx50]) else NA_real_
    breadth50 <- if (length(idx50)) (env50_max - env50_min) else NA_real_
    
    # Suggest plotting "raw relative abundance" for points, more intuitive
    plot_dat <- data.frame(env = env_vec, y = otu_sub)
  }
  
  # Plotting (viridis palette: do not hardcode specific colors; otherwise use schemes above)
  if (palette == "viridis") {
    if (!requireNamespace("viridisLite", quietly = TRUE)) {
      stop("palette='viridis' requires package viridisLite.")
    }
    vv <- viridisLite::viridis(3)
    cols$point <- vv[2]
    cols$line  <- vv[1]
    cols$ci    <- vv[3]
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = plot_dat,
      ggplot2::aes(x = env, y = y),
      shape = 21, size = 3, stroke = 0.7,
      fill = cols$point, colour = cols$edge, alpha = point_alpha
    ) +
    { if (add_ci) ggplot2::geom_ribbon(
      data = grid_df,
      ggplot2::aes(x = env, ymin = fit_lower, ymax = fit_upper),
      alpha = 0.25, fill = cols$ci
    ) } +
    ggplot2::geom_line(
      data = grid_df,
      ggplot2::aes(x = env, y = fit),
      linewidth = 1.2, colour = cols$line
    ) +
    ggplot2::geom_vline(xintercept = optimum_env, linetype = "dashed", colour = cols$vline) +
    ggplot2::labs(
      title = paste0("GAM niche curve - ", otu_id, " (", data_type, ")"),
      x = env_var,
      y = y_label
    ) +
    ggplot2::theme_bw()
  
  list(
    data        = plot_dat,
    grid        = grid_df,
    fit         = fit,
    optimum_env = optimum_env,
    env50_min   = env50_min,
    env50_max   = env50_max,
    breadth50   = breadth50,
    plot        = p
  )
}

