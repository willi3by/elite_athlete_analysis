# Elite athlete structural brain-profile analysis

This is the minimal post-segmentation analysis package for Frontiers manuscript 1904810. It contains only the deidentified analysis dataset, the complete R analysis script, and the package-version lock file.

## Files

- `data/Athlete_Analysis_deidentified.csv` — deidentified post-segmentation analysis dataset
- `revision_analyses.R` — complete analysis workflow
- `renv.lock` — R and package versions used for the analysis

## Run

From the repository root, restore the recorded R environment:

```r
install.packages("renv")
renv::restore()
```

A working C++/Stan toolchain is required. Then run:

```bash
Rscript revision_analyses.R
```

The script starts with a clean `outputs/` directory, refits all Bayesian models, and regenerates the PCA results, bootstrap and parallel analyses, manuscript and supplementary tables, statistical figures, model summaries, diagnostics, fitted model objects, and `sessionInfo.txt`.

Default manuscript settings are 2,000 PCA bootstrap resamples, 1,000 parallel-analysis permutations, four Stan chains, 4,000 iterations per chain, 2,000 warmup iterations, and seed 123.

Optional environment variables are `ANALYSIS_CORES`, `PCA_BOOTSTRAPS`, `PARALLEL_REPS`, `BRMS_CHAINS`, `BRMS_ITER`, `BRMS_WARMUP`, and `BRMS_BACKEND`.

Raw MRI data and the proprietary BrainKey segmentation and normative-reference software are not included. Manuscript Figures 1 and 4 are anatomical rendering panels and cannot be generated from the tabular post-segmentation dataset.
