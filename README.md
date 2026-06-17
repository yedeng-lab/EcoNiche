# EcoNiche

**EcoNiche** is an R package for estimating taxon-level and community-level niche position and niche breadth under continuous, multidimensional environmental gradients.

EcoNiche is designed for ecological community data in which environmental gradients are continuous, multivariate and potentially confounded by background variables. The package combines constrained ordination, partial constrained ordination, generalized additive models and classical discrete-state niche breadth indices in a reproducible R workflow.

The core idea is to estimate niche position and conditional niche breadth in a continuous environmental space, rather than relying only on arbitrary discrete sample states.

Package version: **1.0.3**

---

## Availability

EcoNiche can be installed from CRAN:

```r
install.packages("EcoNiche")
```

The stable release is distributed through **CRAN**. The GitHub repository provides the development version, source code and issue tracking.

The development version can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("zhoushuotao/EcoNiche")
```

Source code and issue tracking are maintained on GitHub:

- CRAN: <https://CRAN.R-project.org/package=EcoNiche>
- GitHub: <https://github.com/zhoushuotao/EcoNiche>
- Issues: <https://github.com/zhoushuotao/EcoNiche/issues>

---

## What EcoNiche Does

EcoNiche provides three complementary analytical modules.

### 1. CCA / partial CCA: main ordination-based framework

This is the primary EcoNiche workflow. It uses canonical correspondence analysis (CCA) and partial CCA to:

- construct a composite environmental gradient from selected predictors,
- estimate taxon-level niche position as an abundance-weighted mean along the constrained axis,
- estimate conditional niche breadth as abundance-weighted dispersion after covariate control,
- aggregate taxon-level breadth estimates to sample or group levels.

Use this workflow when the ecological question involves continuous multivariate gradients and covariate control.

### 2. GAM: single-gradient response curves

The GAM module fits taxon-specific nonlinear response curves along one selected environmental gradient. It returns response optima and `breadth50`, the interval over which the fitted response remains at or above 50% of the fitted peak.

GAM results are useful for visualizing nonlinear responses and comparing curve-derived response ranges. However, `breadth50` describes a response along one selected gradient and should not be treated as a direct substitute for ordination-based niche breadth in multivariate environmental space.

### 3. Levins and Shannon breadth: discrete-state benchmarks

EcoNiche also calculates classical Levins and Shannon niche breadth metrics. These treat samples, or bins along a gradient, as discrete states.

These metrics are provided mainly as benchmarks for comparison with the continuous-gradient ordination workflow. Gradient-binned Levins' breadth remains sensitive to bin number, binning strategy and the distribution of samples along the gradient.

---

## Conceptual Structure

<img width="4159" height="2870" alt="EcoNiche analytical framework" src="https://github.com/user-attachments/assets/4fa71bb4-9b56-4918-8df5-d352e5082100" />

*Analytical structure of EcoNiche. CCA and partial CCA define a composite environmental gradient for estimating niche position and conditional niche breadth. GAM provides single-gradient response curves, and Levins' index serves as a discrete-state benchmark.*

---

## Built-in Plant Dataset

EcoNiche v1.0.3 includes a plant community dataset for reproducible examples:

- `plant_otu`: species-by-sample abundance matrix, with taxa in rows and samples in columns.
- `plant_env`: sample-by-environment data frame, with samples in rows and environmental variables in columns.
- `plant_group`: named vector defining the group of each sample.

```r
library(EcoNiche)

data("plant_otu")
data("plant_env")
data("plant_group")

identical(colnames(plant_otu), rownames(plant_env))
identical(names(plant_group), colnames(plant_otu))

dim(plant_otu)
dim(plant_env)
length(plant_group)
```

Most EcoNiche functions assume the same data orientation:

- `otu`: taxa in rows, samples in columns.
- `env`: samples in rows, environmental variables in columns.
- `group`: named vector or factor; names must be sample IDs.

---

## Quick Start: Ordination-Based Niche Estimation

The example below uses precipitation, soil moisture, available phosphorus and soil pH to construct the focal multivariate environmental gradient. Annual mean temperature is used as a covariate and as the plotting gradient.

```r
library(EcoNiche)

data("plant_otu")
data("plant_env")
data("plant_group")

res_grad <- cca_workflow_gradient(
  otu = plant_otu,
  env = plant_env,
  sel = c("AMP", "Moisture", "AP", "pH"),
  covariates = "AMT",
  var = "AMT",
  method = "loess",
  galaxy_colnum = FALSE
)

head(res_grad$step3$species)
res_grad$step3$plot
res_grad$step4$plot
```

For group-level or sample-level aggregation:

```r
res_group <- cca_workflow_group(
  otu = plant_otu,
  env = plant_env,
  group = plant_group,
  sel = c("AMP", "Moisture", "AP", "pH"),
  covariates = "AMT",
  var = "AMT",
  choice = "all",
  method = "loess",
  plot_type = "sample",
  galaxy_colnum = FALSE
)

head(res_group$step5$sample_summary)
head(res_group$step5$group_summary)
res_group$step6$plot
```

The wrapper functions are convenient for routine use. For a step-by-step implementation using the six core functions, see the package vignette and the supplementary example code associated with the EcoNiche manuscript.

---

## Step-by-Step CCA / Partial CCA Workflow

The wrapper workflow above is composed of six exported functions:

1. `cca_prep_env()`: prepare selected environmental variables and covariates.
2. `cca_fit_ordination()`: fit CCA and partial CCA ordination models.
3. `cca_calc_species()`: calculate taxon-level niche position and niche breadth.
4. `cca_calc_gradient()`: map site-level niche position to an environmental gradient.
5. `cca_calc_group()`: aggregate breadth estimates to sample and group levels.
6. `cca_plot_group_env()`: visualize aggregated breadth along an environmental gradient.

Example skeleton:

```r
step1 <- cca_prep_env(
  env = plant_env,
  sel = c("AMP", "Moisture", "AP", "pH"),
  constrain = "AMT",
  galaxy_colnum = FALSE
)

step2 <- cca_fit_ordination(
  otu = plant_otu,
  env = plant_env,
  sel = c("AMP", "Moisture", "AP", "pH"),
  covariates = "AMT",
  galaxy_colnum = FALSE
)

step3 <- cca_calc_species(
  otu = plant_otu,
  site_width = step2$site_width,
  site_pos = step2$site_pos,
  method = "loess",
  make_plot = TRUE,
  top_node = 10000
)

step4 <- cca_calc_gradient(
  env = plant_env,
  site_pos = step2$site_pos,
  var = "AMT",
  make_plot = TRUE,
  galaxy_colnum = FALSE
)

step5 <- cca_calc_group(
  otu = plant_otu,
  site_width = step2$site_width,
  group = plant_group,
  choice = "all"
)

step6 <- cca_plot_group_env(
  env = plant_env,
  sample_summary = step5$sample_summary,
  var = "AMT",
  group = plant_group,
  plot_type = "sample",
  method = "loess",
  make_plot = TRUE,
  galaxy_colnum = FALSE
)
```

---

## Choosing `sel` and `covariates`

EcoNiche separates environmental predictors into:

- `sel`: focal environmental variables used to define the primary constrained axis.
- `covariates`: conditioning variables used to control background structure.

If samples were collected along a known gradient, that gradient or its associated variables should usually be placed in `sel`, with background variables in `covariates`.

```r
sel = "Temperature"
covariates = c("pH", "NDVI")
```

If the dominant environmental structure is unclear, inspect correlations or PCA loadings first, then define `sel` and `covariates` according to the ecological hypothesis. Avoid using highly collinear predictors in ways that make the constrained axis difficult to interpret.

---

## GAM-Based Response Analysis

Use GAMs when the question is: for each taxon, where is the fitted response optimum along a selected gradient, and how wide is the response curve?

```r
otu_gam <- plant_otu[rowSums(plant_otu > 0) >= 20, ]

gam_df <- gam_fit_model(
  otu = otu_gam,
  env = plant_env,
  env_var = "AMT",
  data_type = "count",
  count_family = "nb",
  use_offset = TRUE,
  min_prev = 0.10,
  min_total = 20,
  k_spline = 5,
  n_grid = 100,
  verbose = FALSE
)

head(gam_df)
```

Plot a response curve for one taxon:

```r
gam_plot <- gam_plot_species(
  otu = plant_otu,
  env = plant_env,
  env_var = "AMT",
  otu_id = gam_df$OTU[1],
  data_type = "count",
  count_family = "nb",
  add_ci = TRUE
)

gam_plot$plot
gam_plot$optimum_env
gam_plot$breadth50
```

Aggregate GAM-derived breadth to sample level:

```r
site_gam_width <- gam_calc_sitewidth(
  otu = plant_otu,
  niche_df = gam_df,
  width_col = "breadth50",
  weight_mode = "auto"
)

head(site_gam_width)
```

---

## Levins and Shannon Niche Breadth

Classical Levins and Shannon breadth can be calculated by treating samples as discrete states:

```r
width_both <- niche_width_calc(
  otu = plant_otu,
  env = plant_env,
  method = "both",
  standardize = TRUE,
  min_occ = 3,
  min_abund = 5
)

head(width_both)
```

EcoNiche also provides a gradient-binned Levins implementation:

```r
levins_binned <- levins_calc_binned(
  otu = plant_otu,
  env = plant_env,
  axis_var = "AMT",
  nbin = 8,
  bin_method = "equal_freq",
  agg_fun = "mean",
  otu_mode = "auto",
  min_occ = 3,
  min_abund = 5
)

head(levins_binned)
```

The binned implementation makes classical Levins' breadth more comparable with gradient-based analyses, but it remains a discretized approximation. It should be interpreted as a discrete-state reference metric rather than as a fully continuous estimate of niche breadth.

Community-level binned Levins breadth can be visualized along the same environmental gradient:

```r
levins_group <- levins_calc_group(
  otu = plant_otu,
  env = plant_env,
  levins_df = levins_binned,
  grad = "AMT",
  width_col = "levins_Bstd",
  method = "loess",
  make_plot = TRUE
)

head(levins_group$data)
levins_group$plot
```

---

## Which Workflow Should I Use?

- Use CCA / partial CCA when the goal is multivariate niche estimation with covariate control.
- Use GAM when the goal is to visualize or summarize taxon-specific responses along one selected gradient.
- Use Levins or Shannon breadth when a classical discrete-state benchmark is needed.
- Use gradient-binned Levins only as a discretized reference for continuous-gradient analyses.

---

## Common Pitfalls

1. Check sample alignment before analysis: `colnames(otu)` should match `rownames(env)`.
2. Use column names rather than numeric indices unless Galaxy-style numeric columns are required.
3. Filter very rare taxa before GAM fitting, because response curves for rare taxa can be unstable.
4. For count data in GAMs, use `use_offset = TRUE` to account for sample-level total abundance.
5. Interpret binned Levins' breadth as a discrete-state benchmark, not as a fully continuous niche breadth estimate.

---

## Vignette

A vignette under `vignettes/` demonstrates the built-in plant workflow.

To build vignettes locally:

```r
install.packages(c("knitr", "rmarkdown"))
devtools::build_vignettes()
```

---

## Exported Functions

CCA / partial CCA:

- `cca_workflow()`, `cca_workflow_gradient()`, `cca_workflow_group()`
- `cca_prep_env()`, `cca_fit_ordination()`
- `cca_calc_species()`, `cca_calc_gradient()`
- `cca_calc_group()`, `cca_plot_group_env()`

GAM:

- `gam_fit_model()`, `gam_plot_species()`, `gam_calc_sitewidth()`

Levins and Shannon:

- `niche_width_calc()`
- `levins_calc_binned()`, `levins_calc_group()`, `levins_plot_pos_width()`

---

## Citation

If you use EcoNiche in a publication, please cite the software and include the version:

> Zhou S, Feng K, Deng Y. EcoNiche: Community Niche Position and Width Estimation Tools. R package version 1.0.3. CRAN: <https://CRAN.R-project.org/package=EcoNiche>. GitHub: <https://github.com/zhoushuotao/EcoNiche>.

Please replace or supplement this with the formal paper citation after the EcoNiche manuscript is published.

---

## Support

Issues and pull requests are welcome. When reporting a bug, please include:

- a minimal reproducible example,
- `sessionInfo()`,
- dimensions of `otu` and `env`,
- the exact error message or traceback.

---

## License

MIT License. See `LICENSE`.


