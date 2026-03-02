# EcoNiche

**EcoNiche** is an R package for quantifying niche position and niche breadth under continuous, multidimensional environmental gradients.

Traditional niche breadth indices often treat sampling units as discrete states. However, real ecological gradients are continuous, multidimensional, and frequently confounded by covariates (e.g., climate, spatial structure, background environmental effects). EcoNiche provides a standardized workflow to estimate niche parameters within a constrained ordination framework while explicitly controlling for covariates.

The package is built around one central idea:

> Estimate niche position and conditional niche breadth in a continuous environmental space, rather than across arbitrary discrete states.

To achieve this, EcoNiche implements three analytical modules with distinct roles:

---

## 1️⃣ CCA / partial-CCA (core framework)

**This is the primary method of EcoNiche.**

This workflow is based on canonical correspondence analysis (CCA) and partial CCA (pCCA) following Okie et al. (2015) and Feng et al. (2020).

It allows users to:

- Construct a **composite environmental gradient** from multiple predictors
- Quantify **niche position** as the abundance-weighted mean along the dominant constrained axis
- Quantify **conditional niche breadth** as abundance-weighted dispersion after covariate control

Conceptually:

- CCA defines the dominant environmental gradient associated with focal predictors.
- pCCA removes covariate structure (e.g., temperature, spatial effects) to reduce confounding.
- Species’ spread in this covariate-controlled space represents conditional niche breadth.

This framework is specifically designed for:

- Continuous gradients
- Multivariate environmental predictors
- Explicit covariate control
- Community-scale comparisons

---

## 2️⃣ GAM (nonlinear single-gradient modeling and species response curves)

The GAM module serves both as a validation tool and as a biologically interpretable response modeling framework.

It fits species-specific response curves along the composite environmental axis (typically CCA1), allowing:

- Estimation of gradient optima
- Extraction of response breadth (Breadth50)
- Assessment of nonlinear response structure
- Explicit visualization of species-specific response curves

This module serves three complementary purposes:

1. **Cross-validation of ordination-based niche position**  
   GAM-derived optima can be compared directly with CCA-based abundance-weighted means to evaluate positional consistency.

2. **Quantification of one-dimensional tolerance along the dominant composite gradient**  
   Breadth50 provides a practical measure of effective response range along the main constrained axis.

3. **Visualization of individual species response curves**  
   GAM explicitly models the abundance–environment relationship, allowing users to:
   - Identify unimodal vs. skewed response shapes
   - Detect boundary-constrained optima
   - Compare narrow vs. broad response profiles across taxa
   - Interpret ecological strategies at the single-species level

Importantly, the GAM analysis is conducted within the same gradient coordinate system defined by CCA, ensuring conceptual consistency between ordination-based and response-curve-based niche estimation.

---

## 3️⃣ Levins’ niche breadth (discrete-state benchmark)

Levins’ index is implemented for comparison, not as the primary estimator.

EcoNiche provides two Levins implementations:

- Sample-defined states (classical formulation)
- Gradient-binned states (bins along the composite CCA axis)

These serve as sensitivity benchmarks to evaluate:

- Dependence on state definition
- Sensitivity to sampling structure
- Confounding with species commonness (occurrence frequency and total abundance)

In continuous-gradient contexts, Levins’ breadth may reflect occupancy or prevalence more than environmental tolerance. Therefore, it is included primarily for methodological comparison.

---

## Conceptual structure of EcoNiche

<img width="4159" height="2870" alt="github计算框架" src="https://github.com/user-attachments/assets/4fa71bb4-9b56-4918-8df5-d352e5082100" />

* Analytical structure of EcoNiche. CCA/pCCA define a composite environmental gradient for estimating niche position and conditional niche breadth. GAM provides species response curves along the same axis, and Levins’ index serves as a discrete-state comparison.*

---

EcoNiche is designed as a modular analytical framework rather than a fixed analytical pipeline.

At its core, the framework is built on constrained ordination:

**CCA / partial-CCA define the composite environmental gradient and estimate niche parameters in continuous multivariate space.**

From this core framework, users may choose different complementary modules depending on their research question:

---

### Core module: CCA / pCCA  
Used to estimate:

- Niche position (abundance-weighted mean along the constrained axis)
- Conditional niche breadth (abundance-weighted dispersion after covariate control)

This module is recommended when:
- Environmental predictors are multidimensional
- Covariate control is required
- Continuous-gradient niche estimation is the primary objective

---

### Optional module A: GAM  
Used when users wish to:

- Visualize species-specific response curves
- Quantify nonlinear response shapes
- Extract one-dimensional optima and response breadth (Breadth50)
- Cross-validate ordination-based niche position

GAM operates along the composite gradient defined by CCA, ensuring consistency in coordinate space.

---

### Optional module B: Levins’ breadth  
Used when users wish to:

- Compute classical discrete-state niche breadth
- Compare continuous-gradient estimation with sample-defined states
- Evaluate sensitivity to binning or state definition
- Assess potential commonness-driven effects

Levins’ index is included primarily as a comparative benchmark rather than the primary estimator under continuous gradients.

---

In practice, users typically:

1. Define focal environmental predictors and covariates.
2. Estimate niche position and conditional niche breadth using CCA/pCCA.
3. Optionally apply GAM for species-level response modeling.
4. Optionally compute Levins’ breadth for methodological comparison.

This flexible structure allows researchers to tailor analyses to their ecological question while maintaining a coherent gradient-based coordinate system.

> Package version: **0.9.0**

---

## Installation

### Install the development version from GitHub

```r
# install.packages("remotes")
remotes::install_github("zhoushuotao/EcoNiche")
````

### Dependencies

EcoNiche imports:

* **vegan** (CCA / partial-CCA)
* **mgcv** (GAM)
* **ggplot2** (plots)

Suggested:

* **knitr**, **rmarkdown** (vignettes)
* **testthat**, **viridisLite**

---

## Data format (important)

Most functions assume the following:

### 1) `otu`: taxon-by-sample abundance table

* **Rows = taxa** (OTUs/ASVs/species)
* **Columns = samples**
* Values can be **counts** or **relative abundances** (non-negative)

```r
otu[1:3, 1:3]
#        S1 S2 S3
# OTU1   10  0  2
# OTU2    3  1  0
# OTU3    0  0  5
```

### 2) `env`: sample-by-environment table

* **Rows = samples**
* **Columns = environmental variables**
* `rownames(env)` should match `colnames(otu)` (EcoNiche will align by intersection when relevant)

```r
env[1:3, ]
#      Temp   pH   SOC
# S1   15.2  6.7  2.31
# S2   18.1  6.3  1.90
# S3   12.4  6.8  3.05
```

### 3) `group` (only for group workflow)

* Named vector/factor of group labels
* `names(group)` must be sample IDs

---

## Quick start (toy example)

```r
library(EcoNiche)

set.seed(1)

# OTU table: 50 taxa x 30 samples
otu <- matrix(rpois(50 * 30, lambda = 10), nrow = 50, ncol = 30)
rownames(otu) <- paste0("OTU", 1:50)
colnames(otu) <- paste0("S", 1:30)

# Environment table: 30 samples x 3 variables
env <- data.frame(
  Temp = rnorm(30, 15, 3),
  pH   = rnorm(30, 6.5, 0.4),
  SOC  = rlnorm(30, 2, 0.3)
)
rownames(env) <- colnames(otu)
```

---

## Workflow A: CCA / partial-CCA (multi-dimensional niche)

EcoNiche provides a controller function `cca_workflow()` and two explicit workflows:

* `cca_workflow_gradient()`
* `cca_workflow_group()`

### A1) Gradient workflow

Use when you want:

* a **position axis** from CCA (typically CCA1),
* a **width axis** from partial-CCA (conditioning on covariates),
* plus basic diagnostic plots.

```r
res <- cca_workflow(
  mode = "gradient",
  otu = otu,
  env = env,
  sel = "Temp",                 # constrained variable(s) for the width axis (partial-CCA)
  covariates = c("pH", "SOC")    # variables for position axis (CCA) and as covariates in partial-CCA
)

# Species-level niche traits (position/width)
head(res$species$species)

# Plot: niche width vs niche position
res$species$plot

# Site niche position vs a gradient variable
res$gradient$plot
```

### Choosing `sel` and `covariates` in CCA/pCCA

EcoNiche separates environmental predictors into:

- `sel` — focal environmental variable(s) defining the primary constrained axis
- `covariates` — conditioning variables used to control background structure

The appropriate choice depends on study design and ecological hypothesis.

---

#### Scenario 1: Sampling along a known focal gradient (hypothesis-driven design)

If samples were collected explicitly along a known environmental gradient (e.g., temperature, pH, elevation), that gradient should typically be treated as the focal predictor.

Example:

- Samples span a temperature gradient across latitude.
- Other variables (e.g., pH, NDVI) vary but are not the primary research focus.

Recommended setup:

```r
sel = "Temperature"
covariates = c("pH", "NDVI")
```

Interpretation:

- CCA constructs the dominant composite axis associated with the focal gradient.

- pCCA controls for background covariates.

- Niche position reflects species location along the focal gradient.

- Conditional niche breadth reflects dispersion after removing background effects.

This setup is appropriate when the research question is:

> How do species’ niches vary along temperature after accounting for other environmental structure?

#### Scenario 2: No predefined focal gradient (exploratory design) 

If the dominant environmental structure is unclear, or multiple predictors are strongly correlated, a data-driven approach is recommended.

In this case:

1. Perform PCA on environmental variables.

2. Inspect loadings of PC1 and PC2.

3. Identify the dominant composite environmental structure.

4. Use this information to guide CCA/pCCA specification.

Example:

- PC1 loads strongly on precipitation and moisture.

- PC2 reflects soil chemistry.

You may then:

- Use predictors associated with PC1 as sel

- Use remaining structure as covariates

- Or reduce dimensionality first via cca_prep_env()

This approach avoids arbitrarily selecting a single variable when gradients are inherently multivariate.

### A2) Group workflow

Use when samples belong to groups (e.g., habitat types, treatments, regions) and you want group summaries.

```r
group <- rep(c("A", "B", "C"), length.out = ncol(otu))
names(group) <- colnames(otu)

res_g <- cca_workflow(
  mode = "group",
  otu = otu,
  env = env,
  sel = "Temp",
  covariates = c("pH", "SOC"),
  group = group,
  plot_type = "both"  # "sample" / "group" / "both"
)

head(res_g$group$group_summary)
res_g$plot
```

---

## Workflow B: GAM response curves along a single gradient

Use when your main question is:

**For each taxon, what is the optimum along gradient X, and how wide is its response?**

EcoNiche fits one GAM per taxon and reports (among others):

* `optimum_env` (peak along the gradient)
* `env50_min`, `env50_max`, `breadth50` (50% peak breadth)

### B1) Fit all taxa

```r
gam_df <- gam_fit_model(
  otu = otu,
  env = env,
  env_var = "Temp",
  data_type = "count",      # "auto" / "count" / "relative"
  count_family = "nb",      # "nb" / "poisson"
  use_offset = TRUE,        # offset = log(library size) for count mode
  min_prev = 0.10,
  min_total = 100,
  k_spline = 5,
  n_grid = 200,
  verbose = TRUE
)

head(gam_df)
```

### B2) Plot a single taxon

```r
p <- gam_plot_species(
  otu = otu,
  env = env,
  env_var = "Temp",
  otu_id = "OTU1",
  data_type = "count",
  count_family = "nb",
  add_ci = TRUE
)

p$plot
p$optimum_env
p$breadth50
```

### B3) Sample-level summary (abundance-weighted mean breadth)

```r
site_bw <- gam_calc_sitewidth(
  otu = otu,
  niche_df = gam_df,
  width_col = "breadth50",
  weight_mode = "auto"
)

head(site_bw)
```

---

## Workflow C: Levins niche breadth

Levins breadth treats environments as **discrete states** and measures how evenly a taxon occupies these states.

EcoNiche supports:

### C1) Samples as states (classic)

```r
lev_samples <- niche_width_calc(
  otu = otu,
  env = env,           # optional; used for aligning samples
  method = "levins",   # "levins" / "shannon" / "both"
  standardize = TRUE,
  min_occ = 3,
  min_abund = 5
)

head(lev_samples)
```

### C2) Gradient-binned states (recommended for continuous gradients)

Discretize a continuous gradient into bins, aggregate abundances within bins, then compute Levins breadth across bins.

```r
lev_binned <- levins_calc_binned(
  otu = otu,
  env = env,
  axis_var = "Temp",
  nbin = 8,
  bin_method = "equal_freq",  # "equal_freq" (quantiles) or "equal_width"
  agg_fun = "mean",           # "mean" recommended; or "sum"
  otu_mode = "auto",
  min_occ = 3,
  min_abund = 5
)

head(lev_binned)
```

### C3) Community-mean Levins width vs gradient

```r
lev_comm <- levins_calc_group(
  otu = otu,
  env = env,
  levins_df = lev_binned,
  grad = "Temp",
  width_col = "levins_Bstd",
  method = "loess",
  make_plot = TRUE
)

head(lev_comm$data)
lev_comm$plot
```

### C4) Position–width relationship (e.g., compare CCA position vs Levins width)

```r
pos_df <- res$species$species
pos_df$OTU <- rownames(pos_df)

pos_width <- levins_plot_pos_width(
  pos_df = pos_df,
  levins_df = lev_binned,
  id_col = "OTU",
  pos_col = "NichePosition",
  width_col = "levins_Bstd",
  method = "lm",
  make_plot = TRUE
)

pos_width$plot
```

---

## Practical guidance (what to use when)

* **Need covariate control + multi-dimensional niche** → CCA / partial-CCA workflow
* **Need interpretable optimum + breadth along one gradient** → GAM workflow
* **Need a classic breadth index / discrete states** → Levins workflow
* **Gradient is continuous but you still want Levins** → `levins_calc_binned()` (binned states)

---

## Common pitfalls

1. **Sample alignment**

   * Ensure `colnames(otu)` and `rownames(env)` refer to the same sample IDs.
2. **Rare taxa**

   * GAM fits can be unstable for very rare taxa. Use `min_prev`, `min_total`, and/or `min_mean`.
3. **Counts vs relative abundance**

   * For counts, consider the offset (`use_offset=TRUE`) to account for library size.
4. **Numeric column indices**

   * Some functions support numeric indices (Galaxy-style conventions). To avoid confusion, prefer **column names**.

---

## Vignette

A vignette exists under `vignettes/` (currently titled “CoNiche workflow” due to legacy naming).

To build vignettes locally:

```r
install.packages(c("knitr", "rmarkdown"))
devtools::build_vignettes()
```

---

## API overview (exported functions)

**CCA / partial-CCA**

* `cca_workflow()`, `cca_workflow_gradient()`, `cca_workflow_group()`
* `cca_prep_env()`, `cca_fit_ordination()`
* `cca_calc_species()`, `cca_calc_gradient()`
* `cca_calc_group()`, `cca_plot_group_env()`

**GAM**

* `gam_fit_model()`, `gam_plot_species()`, `gam_calc_sitewidth()`

**Levins**

* `niche_width_calc()`
* `levins_calc_binned()`, `levins_calc_group()`, `levins_plot_pos_width()`

---

## Citation

If you use EcoNiche in a publication, please cite the software (and include the version):

> Zhou S, Feng K, Deng Y. EcoNiche: Community Niche Position and Width Estimation Tools (R package), v0.9.0. GitHub: [https://github.com/yedeng-lab/EcoNiche](https://github.com/yedeng-lab/EcoNiche)

(Replace with a formal paper citation if/when a methods manuscript is published.)

---

## Contributing / Support

Issues and pull requests are welcome.

When reporting bugs, please include:

* a minimal reproducible example,
* `sessionInfo()`,
* dimensions of `otu` and `env`,
* the exact error message / traceback.

---

## License

MIT License (see `LICENSE` / `LICENSE.md`).

```
::contentReference[oaicite:0]{index=0}
```
