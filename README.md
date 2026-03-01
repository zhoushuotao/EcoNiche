````markdown
# EcoNiche

**EcoNiche** is an R package for estimating **niche position** and **niche width** for taxa (OTUs/ASVs/species), and for summarizing these traits at the **sample** or **group** level.

It implements three complementary workflows:

1. **CCA / partial-CCA** (multi-dimensional environmental niche; supports covariate control)
2. **GAM along a single environmental gradient** (taxon response curves; optimum + breadth)
3. **Levins niche breadth** (samples-as-states or gradient-binned states)

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

> Zhou S, Feng K, Deng Y. EcoNiche: Community Niche Position and Width Estimation Tools (R package), v0.9.0. GitHub: [https://github.com/zhoushuotao/EcoNiche](https://github.com/zhoushuotao/EcoNiche)

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
