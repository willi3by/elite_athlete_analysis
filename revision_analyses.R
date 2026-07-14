#!/usr/bin/env Rscript

# Complete post-segmentation analysis for Frontiers manuscript 1904810.
#
# The primary PCA and Bayesian models preserve the variable order, preprocessing,
# response definitions, and unmodified PCA signs used in the original analysis.
# Editorially requested additions (KINARM bootstrap, parallel analysis, alternate
# performance definitions, and hockey sensitivity analyses) are identified below.
#
# Run from the repository root:
#   Rscript revision_analyses.R
#
# Default manuscript settings:
#   2,000 PCA bootstrap resamples
#   1,000 permutation-based parallel-analysis resamples
#   4 Stan chains, 4,000 iterations, 2,000 warmup iterations
#
# Optional environment variables:
#   ANALYSIS_CORES       Number of computational cores (default: up to 4)
#   PCA_BOOTSTRAPS       Number of bootstrap resamples (default: 2000)
#   PARALLEL_REPS        Number of parallel-analysis permutations (default: 1000)
#   BRMS_CHAINS          Number of chains (default: 4)
#   BRMS_ITER            Iterations per chain (default: 4000)
#   BRMS_WARMUP          Warmup iterations per chain (default: 2000)
#   BRMS_BACKEND         Optional brms backend, e.g. cmdstanr

options(stringsAsFactors = FALSE)

required_packages <- c(
  "brms", "posterior", "clue", "readr", "dplyr", "tidyr", "purrr",
  "stringr", "forcats", "ggplot2", "tibble"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running the analysis. See README.md."
  )
}

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
  library(clue)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(tibble)
})

# -----------------------------------------------------------------------------
# Repository and run configuration
# -----------------------------------------------------------------------------

find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    marker <- file.path(current, "data", "Athlete_Analysis_deidentified.csv")
    if (file.exists(marker)) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop(
        "Could not locate repository root. Run from the repository root or a ",
        "subdirectory containing data/Athlete_Analysis_deidentified.csv."
      )
    }
    current <- parent
  }
}

repo_root <- find_repo_root()
setwd(repo_root)

parse_int <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value) || value < minimum) {
    stop(name, " must be an integer >= ", minimum, ".")
  }
  value
}

available_cores <- parallel::detectCores(logical = TRUE)
if (is.na(available_cores)) available_cores <- 1L
default_cores <- max(1L, min(4L, available_cores - 1L))

analysis_cores <- parse_int("ANALYSIS_CORES", default_cores)
n_boot <- parse_int("PCA_BOOTSTRAPS", 2000L)
n_parallel <- parse_int("PARALLEL_REPS", 1000L)
brms_chains <- parse_int("BRMS_CHAINS", 4L)
brms_iter <- parse_int("BRMS_ITER", 4000L)
brms_warmup <- parse_int("BRMS_WARMUP", 2000L)
if (brms_warmup >= brms_iter) stop("BRMS_WARMUP must be smaller than BRMS_ITER.")
brms_backend <- Sys.getenv("BRMS_BACKEND", unset = "")

options(mc.cores = analysis_cores)
set.seed(123)

output_root <- "outputs"

# Every run starts clean and refits all Bayesian models.
if (dir.exists(output_root)) {
  unlink(output_root, recursive = TRUE, force = TRUE)
}

output_dirs <- c(
  output_root,
  file.path(output_root, "data"),
  file.path(output_root, "tables"),
  file.path(output_root, "figures"),
  file.path(output_root, "models"),
  file.path(output_root, "model_summaries"),
  file.path(output_root, "diagnostics"),
  file.path(output_root, "pca")
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

log_file <- file.path(output_root, "analysis_log.txt")
log_message <- function(...) {
  text <- paste0(...)
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", stamp, "] ", text)
  message(line)
  cat(line, "\n", file = log_file, append = TRUE)
}
cat("", file = log_file)

log_message("Repository root: ", repo_root)
log_message("Cores: ", analysis_cores, "; PCA bootstraps: ", n_boot,
            "; parallel-analysis repetitions: ", n_parallel)
log_message("Stan settings: chains=", brms_chains, ", iter=", brms_iter,
            ", warmup=", brms_warmup)

# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

safe_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, path, na = "")
}

clean_region_label <- function(x) {
  labels <- c(
    "Left.Hippocampus_1" = "Left Hippocampus",
    "Right.Hippocampus_1" = "Right Hippocampus",
    "Left.Lateral.Ventricle_1" = "Left Lateral Ventricle",
    "Right.Lateral.Ventricle_1" = "Right Lateral Ventricle",
    "Third.Ventricle_1" = "Third Ventricle",
    "Fourth.Ventricle_1" = "Fourth Ventricle",
    "Frontal.Gray.Matter_1" = "Frontal Gray Matter",
    "Parietal.Gray.Matter_1" = "Parietal Gray Matter",
    "Temporal.Gray.Matter_1" = "Temporal Gray Matter",
    "Occipital.Gray.Matter_1" = "Occipital Gray Matter",
    "Caudate_1" = "Caudate",
    "Putamen_1" = "Putamen",
    "Amygdala_1" = "Amygdala",
    "Thalamus_1" = "Thalamus",
    "White.Matter_1" = "White Matter",
    "Brainstem_1" = "Brainstem",
    "Globus.Pallidus_1" = "Globus Pallidus",
    "Cerebellum_1" = "Cerebellum"
  )
  out <- unname(labels[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

kinarm_labels <- c(
  "ArmPosRScore" = "Arm Position Matching: right-arm task score",
  "ArmPosRM" = "Arm Position Matching: right-arm metric score",
  "ArmPosLScoreArmPosLM" = "Arm Position Matching: left-arm task/metric score",
  "BallOnBarScore" = "Ball on Bar: overall task score",
  "BallOnBarM" = "Ball on Bar: overall metric score",
  "TotalHit" = "Object Hit: total objects hit",
  "HitsLeft" = "Object Hit: left-hand hits",
  "HitsRight" = "Object Hit: right-hand hits",
  "TaskScore" = "Object Hit: task score",
  "Mscore" = "Object Hit: metric score",
  "ObjectHitAvoid_objectsHit" = "Object Hit and Avoid: objects hit",
  "ObjectHitScore" = "Object Hit and Avoid: task score",
  "ObjectHitM" = "Object Hit and Avoid: metric score",
  "ObjectHitZ" = "Object Hit and Avoid: standardized score",
  "ObjectsHitProcessRate" = "Object Hit and Avoid: processing rate",
  "ObjectHitRHspeed" = "Object Hit and Avoid: right-hand speed",
  "ObjectHitLHSpeed" = "Object Hit and Avoid: left-hand speed",
  "ObjectHitDistractorIEErrorsHit" = "Object Hit and Avoid: distractor errors",
  "TaskScore.1" = "Object Hit and Avoid: composite task score",
  "Mscore.1" = "Object Hit and Avoid: composite metric score",
  "Total score (#)" = "Ball on Bar: total score",
  "Time" = "Ball on Bar Level 1 (fixed ball): completion time",
  "time1/time2" = "Ball on Bar Level 1 (fixed ball): time ratio",
  "dwelltime" = "Ball on Bar Level 1 (fixed ball): dwell time",
  "errors" = "Ball on Bar Level 1 (fixed ball): errors",
  "Time.1" = "Ball on Bar Level 2 (angle-dependent ball): completion time",
  "time1/time2.1" = "Ball on Bar Level 2 (angle-dependent ball): time ratio",
  "dwelltime.1" = "Ball on Bar Level 2 (angle-dependent ball): dwell time",
  "errors.1" = "Ball on Bar Level 2 (angle-dependent ball): errors",
  "Time.2" = "Ball on Bar Level 3 (freely moving ball): completion time",
  "TaskScore.2" = "Ball on Bar Level 3 (freely moving ball): task score",
  "Mscore.2" = "Ball on Bar Level 3 (freely moving ball): metric score",
  "Posture speed (m/s)" = "Visually Guided Reaching: posture speed, limb 1",
  "Reaction time (s)" = "Visually Guided Reaching: reaction time, limb 1",
  "Movement time (s)" = "Visually Guided Reaching: movement time, limb 1",
  "Path length ratio" = "Visually Guided Reaching: path-length ratio, limb 1",
  "Posture speed (m/s).1" = "Visually Guided Reaching: posture speed, limb 2",
  "Reaction time (s).1" = "Visually Guided Reaching: reaction time, limb 2",
  "Movement time (s).1" = "Visually Guided Reaching: movement time, limb 2",
  "Path length ratio.1" = "Visually Guided Reaching: path-length ratio, limb 2",
  "LTaskScore" = "Visually Guided Reaching: left-arm task score",
  "LMScore" = "Visually Guided Reaching: left-arm metric score",
  "RTaskScore" = "Visually Guided Reaching: right-arm task score",
  "RMScore" = "Visually Guided Reaching: right-arm metric score",
  "LeftPostureSpeed" = "Visually Guided Reaching: left-arm posture speed",
  "LeftReactionTime" = "Visually Guided Reaching: left-arm reaction time",
  "LeftPathLength" = "Visually Guided Reaching: left-arm path length",
  "LeftMT" = "Visually Guided Reaching: left-arm movement time",
  "RightPostureSpeed" = "Visually Guided Reaching: right-arm posture speed",
  "RightReactionTime" = "Visually Guided Reaching: right-arm reaction time",
  "RightPathLength" = "Visually Guided Reaching: right-arm path length",
  "RightMT" = "Visually Guided Reaching: right-arm movement time",
  "LTaskScore.1" = "Visually Guided Reaching: left-arm task score, repeated block",
  "LMScore.1" = "Visually Guided Reaching: left-arm metric score, repeated block",
  "RTaskScore.1" = "Visually Guided Reaching: right-arm task score, repeated block",
  "RMScore.1" = "Visually Guided Reaching: right-arm metric score, repeated block"
)

format_p <- function(x) {
  if (is.na(x)) return("")
  if (x < 0.001) return("<0.001")
  sprintf("%.3f", x)
}

format_mean_sd <- function(x) sprintf("%.1f ± %.1f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
format_count_pct <- function(n, denom) sprintf("%d (%.0f%%)", n, 100 * n / denom)

categorical_p <- function(tab) {
  if (any(dim(tab) < 2L)) return(NA_real_)
  expected <- suppressWarnings(chisq.test(tab, correct = FALSE)$expected)
  if (any(expected < 5)) {
    # This mirrors the default categorical test used by gtsummary::add_p().
    tryCatch(
      fisher.test(tab, workspace = 2e8)$p.value,
      error = function(e) {
        warning("Exact Fisher test failed; using a deterministic Monte Carlo fallback: ", conditionMessage(e))
        set.seed(123L)
        fisher.test(tab, simulate.p.value = TRUE, B = 200000)$p.value
      }
    )
  } else {
    suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
  }
}

loading_display_wide <- function(loadings, variance_pct, labels, threshold = 0.20, strict = FALSE) {
  pcs <- colnames(loadings)
  long <- as.data.frame(loadings) %>%
    rownames_to_column("Variable") %>%
    pivot_longer(-Variable, names_to = "Component", values_to = "Loading")
  if (strict) {
    long <- long %>% filter(abs(Loading) > threshold)
  } else {
    long <- long %>% filter(abs(Loading) >= threshold)
  }
  long <- long %>%
    mutate(
      Label = labels[Variable],
      Label = if_else(is.na(Label), Variable, Label)
    ) %>%
    group_by(Component) %>%
    arrange(desc(abs(Loading)), .by_group = TRUE) %>%
    mutate(Row = row_number(), Display = sprintf("%s (%.3f)", Label, Loading)) %>%
    ungroup()

  wide <- long %>%
    select(Row, Component, Display) %>%
    pivot_wider(names_from = Component, values_from = Display) %>%
    arrange(Row) %>%
    select(-Row)
  for (pc in pcs) if (!pc %in% names(wide)) wide[[pc]] <- NA_character_
  wide <- wide[, pcs, drop = FALSE]
  names(wide) <- sprintf("%s (%.1f%%)", pcs, variance_pct[seq_along(pcs)])
  list(wide = wide, long = long)
}

bootstrap_pca_stability <- function(x, reference_loadings, reference_variance, n_boot, seed = 123L) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (anyNA(x)) stop("bootstrap_pca_stability requires complete data.")
  k <- ncol(reference_loadings)
  p <- ncol(x)
  n <- nrow(x)

  boot_loadings <- array(
    NA_real_, dim = c(n_boot, p, k),
    dimnames = list(NULL, rownames(reference_loadings), colnames(reference_loadings))
  )
  boot_congruence <- matrix(NA_real_, nrow = n_boot, ncol = k)
  boot_variance <- matrix(NA_real_, nrow = n_boot, ncol = k)

  set.seed(seed)
  success <- 0L
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    xb <- x[idx, , drop = FALSE]
    if (any(apply(xb, 2, sd) <= .Machine$double.eps)) next

    fit <- tryCatch(prcomp(xb, center = TRUE, scale. = TRUE), error = function(e) NULL)
    if (is.null(fit)) next
    candidate <- fit$rotation[, seq_len(k), drop = FALSE]
    similarity <- abs(cor(reference_loadings, candidate))
    assignment <- as.integer(clue::solve_LSAP(similarity, maximum = TRUE))
    matched <- candidate[, assignment, drop = FALSE]
    matched_variance <- (fit$sdev^2 / sum(fit$sdev^2))[assignment]

    for (j in seq_len(k)) {
      agreement <- cor(reference_loadings[, j], matched[, j])
      if (!is.finite(agreement)) next
      if (agreement < 0) {
        matched[, j] <- -matched[, j]
        agreement <- -agreement
      }
      boot_congruence[b, j] <- agreement
    }
    success <- success + 1L
    boot_loadings[b, , ] <- matched
    boot_variance[b, ] <- matched_variance
  }

  keep <- complete.cases(boot_congruence)
  if (!any(keep)) stop("All PCA bootstrap iterations failed.")
  boot_loadings <- boot_loadings[keep, , , drop = FALSE]
  boot_congruence <- boot_congruence[keep, , drop = FALSE]
  boot_variance <- boot_variance[keep, , drop = FALSE]

  component <- map_dfr(seq_len(k), function(j) {
    q <- quantile(boot_congruence[, j], c(.025, .25, .5, .75, .975), na.rm = TRUE)
    vq <- quantile(boot_variance[, j], c(.025, .5, .975), na.rm = TRUE)
    tibble(
      Component = paste0("PC", j),
      `Original variance (%)` = 100 * reference_variance[j],
      `Median congruence` = unname(q[3]),
      `IQR lower` = unname(q[2]),
      `IQR upper` = unname(q[4]),
      `95% congruence lower` = unname(q[1]),
      `95% congruence upper` = unname(q[5]),
      `Bootstrap variance median (%)` = 100 * unname(vq[2]),
      `Bootstrap variance lower 95 (%)` = 100 * unname(vq[1]),
      `Bootstrap variance upper 95 (%)` = 100 * unname(vq[3])
    )
  })

  loading <- map_dfr(seq_len(k), function(j) {
    map_dfr(seq_len(p), function(i) {
      draws <- boot_loadings[, i, j]
      original <- reference_loadings[i, j]
      q <- quantile(draws, c(.025, .5, .975), na.rm = TRUE)
      tibble(
        Component = paste0("PC", j),
        Variable = rownames(reference_loadings)[i],
        `Original loading` = original,
        `Bootstrap median` = unname(q[2]),
        `95% lower` = unname(q[1]),
        `95% upper` = unname(q[3]),
        `Same-sign resamples (%)` = 100 * mean(sign(draws) == sign(original), na.rm = TRUE)
      )
    })
  })

  list(
    component = component,
    loading = loading,
    successful_iterations = nrow(boot_congruence),
    congruence_draws = boot_congruence
  )
}

parallel_analysis <- function(x, n_rep = 1000L, seed = 123L) {
  x <- as.matrix(x)

  if (nrow(x) < 2L || ncol(x) < 1L) {
    stop("parallel_analysis() requires at least two observations and one variable.")
  }
  if (anyNA(x)) {
    stop("parallel_analysis() received missing values; use the same complete-case matrix as the observed PCA.")
  }
  zero_variance <- apply(x, 2, function(z) !is.finite(sd(z)) || sd(z) == 0)
  if (any(zero_variance)) {
    stop(
      "parallel_analysis() received zero-variance variables: ",
      paste(colnames(x)[zero_variance], collapse = ", ")
    )
  }

  # prcomp() returns min(n observations, p variables) singular values.  This is
  # smaller than p when p > n, as in the 44-participant × 56-variable KINARM
  # analysis.  Allocate to the number of eigenvalues prcomp actually returns,
  # rather than to the number of input variables.
  reference_fit <- prcomp(x, center = TRUE, scale. = TRUE)
  n_components <- length(reference_fit$sdev)
  eigenvalues <- matrix(NA_real_, nrow = n_rep, ncol = n_components)

  set.seed(seed)
  for (b in seq_len(n_rep)) {
    permuted <- vapply(
      seq_len(ncol(x)),
      function(j) sample(x[, j], size = nrow(x), replace = FALSE),
      numeric(nrow(x))
    )
    fit <- prcomp(permuted, center = TRUE, scale. = TRUE)
    permuted_eigenvalues <- fit$sdev^2

    if (length(permuted_eigenvalues) != n_components) {
      stop(
        "Unexpected PCA dimension during parallel analysis: expected ",
        n_components, " eigenvalues but obtained ",
        length(permuted_eigenvalues), "."
      )
    }
    eigenvalues[b, ] <- permuted_eigenvalues
  }

  out <- apply(eigenvalues, 2, quantile, probs = 0.95, na.rm = TRUE)
  names(out) <- paste0("PC", seq_along(out))
  out
}

make_scree_plot <- function(observed, parallel95, path, max_components = length(observed)) {
  p <- length(observed)
  broken_stick <- rev(cumsum(rev(1 / seq_len(p))))
  n_show <- min(max_components, p)
  plot_df <- tibble(
    Component = seq_len(n_show),
    Observed = observed[seq_len(n_show)],
    `Broken-stick` = broken_stick[seq_len(n_show)],
    `Parallel analysis (95th percentile)` = parallel95[seq_len(n_show)]
  ) %>%
    pivot_longer(-Component, names_to = "Series", values_to = "Eigenvalue")

  p_plot <- ggplot(plot_df, aes(Component, Eigenvalue, linetype = Series, shape = Series)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.8) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
    scale_x_continuous(breaks = seq_len(n_show)) +
    scale_linetype_manual(values = c(
      "Observed" = "solid",
      "Broken-stick" = "dotted",
      "Parallel analysis (95th percentile)" = "dotdash"
    )) +
    labs(x = "Principal component", y = "Eigenvalue", linetype = NULL, shape = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(path, p_plot, width = 8, height = 5, dpi = 300)
  invisible(plot_df)
}

# -----------------------------------------------------------------------------
# Load and validate deidentified data
# -----------------------------------------------------------------------------

data_path <- file.path("data", "Athlete_Analysis_deidentified.csv")
dat <- readr::read_csv(data_path, show_col_types = FALSE)
if (nrow(dat) != 51L) warning("Expected 51 participants; found ", nrow(dat), ".")

brain_vars <- c(
  "Amygdala_1",
  "Brainstem_1",
  "Caudate_1",
  "Cerebellum_1",
  "Fourth.Ventricle_1",
  "Frontal.Gray.Matter_1",
  "Globus.Pallidus_1",
  "White.Matter_1",
  "Left.Hippocampus_1",
  "Left.Lateral.Ventricle_1",
  "Occipital.Gray.Matter_1",
  "Parietal.Gray.Matter_1",
  "Putamen_1",
  "Right.Hippocampus_1",
  "Right.Lateral.Ventricle_1",
  "Temporal.Gray.Matter_1",
  "Thalamus_1",
  "Third.Ventricle_1"
)

kinarm_vars <- c(
  "ArmPosRScore", "ArmPosRM", "ArmPosLScoreArmPosLM", "BallOnBarScore", "BallOnBarM",
  "TotalHit", "HitsLeft", "HitsRight", "TaskScore", "Mscore",
  "ObjectHitAvoid_objectsHit", "ObjectHitScore", "ObjectHitM", "ObjectHitZ",
  "ObjectsHitProcessRate", "ObjectHitRHspeed", "ObjectHitLHSpeed",
  "ObjectHitDistractorIEErrorsHit", "TaskScore.1", "Mscore.1", "Total score (#)",
  "Time", "time1/time2", "dwelltime", "errors", "Time.1", "time1/time2.1",
  "dwelltime.1", "errors.1", "Time.2", "TaskScore.2", "Mscore.2",
  "Posture speed (m/s)", "Reaction time (s)", "Movement time (s)", "Path length ratio",
  "Posture speed (m/s).1", "Reaction time (s).1", "Movement time (s).1", "Path length ratio.1",
  "LTaskScore", "LMScore", "RTaskScore", "RMScore", "LeftPostureSpeed",
  "LeftReactionTime", "LeftPathLength", "LeftMT", "RightPostureSpeed",
  "RightReactionTime", "RightPathLength", "RightMT", "LTaskScore.1", "LMScore.1",
  "RTaskScore.1", "RMScore.1"
)

required_columns <- unique(c(
  "Participant_ID", "Age", "Age_actual", "Rating", "Sex", "Race", "AthleteType",
  "Concussion_Group", brain_vars, kinarm_vars
))
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0) {
  stop("Missing required dataset columns: ", paste(missing_columns, collapse = ", "))
}

# Coerce analysis variables and define outcomes exactly once.
dat <- dat %>%
  mutate(
    across(all_of(c(brain_vars, kinarm_vars)), as.numeric),
    Sex = factor(Sex, levels = c("F", "M")),
    Rating = as.integer(Rating),
    Performance_Group = case_when(
      Rating %in% c(1L, 2L) ~ "High",
      Rating == 3L ~ "Mid",
      Rating %in% c(4L, 5L) ~ "Low",
      TRUE ~ NA_character_
    ),
    Performance_Group = ordered(Performance_Group, levels = c("Low", "Mid", "High")),
    Performance_Rating5 = ordered(as.character(Rating), levels = c("5", "4", "3", "2", "1")),
    High_Performance = if_else(Rating %in% c(1L, 2L), 1L, 0L, missing = NA_integer_),
    Top_Performance = if_else(Rating == 1L, 1L, 0L, missing = NA_integer_),
    Sport_raw = str_squish(str_to_lower(AthleteType)),
    Sport = case_when(
      is.na(Sport_raw) ~ NA_character_,
      str_detect(Sport_raw, "hockey|goalie") ~ "Hockey",
      str_detect(Sport_raw, "rubgy|rugby") ~ "Rugby",
      str_detect(Sport_raw, "football") ~ "Football",
      str_detect(Sport_raw, "alpine|skiing|skier|figure skater|snowboarding") ~ "Skiing/Skating/Snow",
      str_detect(Sport_raw, "sprinter|runner|1500") ~ "Track/Running",
      str_detect(Sport_raw, "swimmer") ~ "Swimming",
      str_detect(Sport_raw, "basketball") ~ "Basketball",
      str_detect(Sport_raw, "lacrosse") ~ "Lacrosse",
      str_detect(Sport_raw, "bobsled") ~ "Bobsled",
      str_detect(Sport_raw, "baseball") ~ "Baseball",
      str_detect(Sport_raw, "soccer") ~ "Soccer",
      str_detect(Sport_raw, "mma") ~ "MMA",
      str_detect(Sport_raw, "trampoline") ~ "Trampoline",
      TRUE ~ "Other"
    ),
    Sport = factor(Sport),
    Hockey = case_when(
      is.na(Sport) ~ NA_integer_,
      Sport == "Hockey" ~ 1L,
      TRUE ~ 0L
    ),
    Concussion_Group = factor(as.integer(Concussion_Group), levels = c(0L, 1L))
  )

concussion_counts <- table(dat$Concussion_Group, useNA = "ifany")
log_message("Concussion-history counts: ", paste(names(concussion_counts), concussion_counts, collapse = "; "))

# -----------------------------------------------------------------------------
# Table 1: participant characteristics
# -----------------------------------------------------------------------------

performance_levels_display <- c("High", "Mid", "Low")

age_p <- kruskal.test(Age_actual ~ Performance_Group, data = dat)$p.value
sex_p <- categorical_p(table(dat$Sex, dat$Performance_Group))
sport_p <- categorical_p(table(dat$Sport, dat$Performance_Group))
concussion_p <- suppressWarnings(
  chisq.test(table(dat$Concussion_Group, dat$Performance_Group), correct = FALSE)$p.value
)

characteristic_rows <- list()
add_characteristic_row <- function(label, overall = "", high = "", mid = "", low = "", p = "") {
  characteristic_rows[[length(characteristic_rows) + 1L]] <<- tibble(
    Characteristic = label, Overall = overall, High = high, Mid = mid, Low = low, p_value = p
  )
}

add_characteristic_row(
  "Age, years",
  format_mean_sd(dat$Age_actual),
  format_mean_sd(dat$Age_actual[dat$Performance_Group == "High"]),
  format_mean_sd(dat$Age_actual[dat$Performance_Group == "Mid"]),
  format_mean_sd(dat$Age_actual[dat$Performance_Group == "Low"]),
  format_p(age_p)
)
add_characteristic_row("Sex", p = format_p(sex_p))
for (level in levels(dat$Sex)) {
  add_characteristic_row(
    paste0("    ", level),
    format_count_pct(sum(dat$Sex == level, na.rm = TRUE), sum(!is.na(dat$Sex))),
    format_count_pct(sum(dat$Sex == level & dat$Performance_Group == "High", na.rm = TRUE), sum(dat$Performance_Group == "High")),
    format_count_pct(sum(dat$Sex == level & dat$Performance_Group == "Mid", na.rm = TRUE), sum(dat$Performance_Group == "Mid")),
    format_count_pct(sum(dat$Sex == level & dat$Performance_Group == "Low", na.rm = TRUE), sum(dat$Performance_Group == "Low"))
  )
}
add_characteristic_row(paste0("Sport (N=", sum(!is.na(dat$Sport)), ")"), p = format_p(sport_p))
for (level in levels(droplevels(dat$Sport))) {
  add_characteristic_row(
    paste0("    ", level),
    format_count_pct(sum(dat$Sport == level, na.rm = TRUE), sum(!is.na(dat$Sport))),
    format_count_pct(sum(dat$Sport == level & dat$Performance_Group == "High", na.rm = TRUE), sum(!is.na(dat$Sport) & dat$Performance_Group == "High")),
    format_count_pct(sum(dat$Sport == level & dat$Performance_Group == "Mid", na.rm = TRUE), sum(!is.na(dat$Sport) & dat$Performance_Group == "Mid")),
    format_count_pct(sum(dat$Sport == level & dat$Performance_Group == "Low", na.rm = TRUE), sum(!is.na(dat$Sport) & dat$Performance_Group == "Low"))
  )
}
add_characteristic_row(
  paste0("Concussion History (N=", sum(!is.na(dat$Concussion_Group)), ")"),
  p = format_p(concussion_p)
)
for (value in c(0L, 1L)) {
  label <- if (value == 0L) "    No Concussion History" else "    Concussion History"
  add_characteristic_row(
    label,
    format_count_pct(sum(as.character(dat$Concussion_Group) == as.character(value), na.rm = TRUE), sum(!is.na(dat$Concussion_Group))),
    format_count_pct(sum(as.character(dat$Concussion_Group) == as.character(value) & dat$Performance_Group == "High", na.rm = TRUE), sum(!is.na(dat$Concussion_Group) & dat$Performance_Group == "High")),
    format_count_pct(sum(as.character(dat$Concussion_Group) == as.character(value) & dat$Performance_Group == "Mid", na.rm = TRUE), sum(!is.na(dat$Concussion_Group) & dat$Performance_Group == "Mid")),
    format_count_pct(sum(as.character(dat$Concussion_Group) == as.character(value) & dat$Performance_Group == "Low", na.rm = TRUE), sum(!is.na(dat$Concussion_Group) & dat$Performance_Group == "Low"))
  )
}

table1 <- bind_rows(characteristic_rows)
safe_write_csv(table1, file.path(output_root, "tables", "table_1_participant_characteristics.csv"))

# -----------------------------------------------------------------------------
# Brain PCA, bootstrap stability, and scree plot
# -----------------------------------------------------------------------------

log_message("Running brain PCA using the original variable order and prcomp sign convention.")
brain_mat <- as.matrix(dat[, brain_vars])
storage.mode(brain_mat) <- "double"

# Preserve the original analysis behavior: mean-impute each regional variable
# before PCA. The released dataset is complete, but this keeps the workflow exact.
brain_means <- colMeans(brain_mat, na.rm = TRUE)
for (j in seq_len(ncol(brain_mat))) {
  brain_mat[is.na(brain_mat[, j]), j] <- brain_means[j]
}
if (anyNA(brain_mat)) stop("Brain PCA variables still contain missing values after mean imputation.")

brain_pca <- prcomp(brain_mat, center = TRUE, scale. = TRUE)
brain_loadings <- brain_pca$rotation[, 1:5, drop = FALSE]
brain_scores <- scale(brain_pca$x[, 1:5, drop = FALSE])
colnames(brain_loadings) <- paste0("PC", 1:5)
colnames(brain_scores) <- paste0("Brain_PC", 1:5)
dat <- bind_cols(dat, as_tibble(brain_scores))

brain_variance <- brain_pca$sdev^2 / sum(brain_pca$sdev^2)
brain_variance_table <- tibble(
  Component = paste0("PC", seq_along(brain_variance)),
  Eigenvalue = brain_pca$sdev^2,
  Variance = brain_variance,
  `Variance (%)` = 100 * brain_variance,
  `Cumulative variance (%)` = 100 * cumsum(brain_variance)
)
safe_write_csv(brain_variance_table, file.path(output_root, "pca", "brain_pca_variance.csv"))
safe_write_csv(
  as.data.frame(brain_loadings) %>% rownames_to_column("Region") %>%
    mutate(Region_label = clean_region_label(Region), .after = Region),
  file.path(output_root, "pca", "brain_pca_loadings.csv")
)

brain_display <- loading_display_wide(
  brain_loadings, 100 * brain_variance[1:5],
  labels = setNames(clean_region_label(rownames(brain_loadings)), rownames(brain_loadings)),
  threshold = 0.20, strict = FALSE
)
safe_write_csv(brain_display$wide, file.path(output_root, "tables", "table_2_brain_pca_loadings_display.csv"))
safe_write_csv(
  brain_display$long %>% select(Component, Variable, Label, Loading),
  file.path(output_root, "tables", "table_2_brain_pca_loadings_long.csv")
)

log_message("Running brain PCA bootstrap stability assessment.")
brain_boot <- bootstrap_pca_stability(
  brain_mat, brain_loadings, brain_variance[1:5], n_boot = n_boot, seed = 123L
)
log_message("Brain PCA bootstrap successful iterations: ", brain_boot$successful_iterations)
safe_write_csv(brain_boot$component, file.path(output_root, "tables", "supplementary_table_s1a_brain_component_stability.csv"))
brain_loading_all <- brain_boot$loading %>%
  mutate(Region = clean_region_label(Variable), .after = Component)
safe_write_csv(brain_loading_all, file.path(output_root, "pca", "brain_pca_bootstrap_loading_stability_all.csv"))
safe_write_csv(
  brain_loading_all %>% filter(abs(`Original loading`) >= 0.20) %>% select(-Variable),
  file.path(output_root, "tables", "supplementary_table_s1b_brain_loading_stability.csv")
)

log_message("Running brain permutation parallel analysis.")
brain_parallel95 <- parallel_analysis(brain_mat, n_rep = n_parallel, seed = 123L)
brain_scree_data <- make_scree_plot(
  brain_pca$sdev^2, brain_parallel95,
  file.path(output_root, "figures", "supplementary_figure_s1_brain_pca_scree.png"),
  max_components = length(brain_pca$sdev)
)
safe_write_csv(brain_scree_data, file.path(output_root, "pca", "brain_pca_scree_data.csv"))

# -----------------------------------------------------------------------------
# KINARM PCA, bootstrap stability, and scree plot
# -----------------------------------------------------------------------------

log_message("Running KINARM PCA using the original variable order and prcomp sign convention.")
kin_complete <- complete.cases(dat[, kinarm_vars])
kin_dat <- as.matrix(dat[kin_complete, kinarm_vars])
storage.mode(kin_dat) <- "double"
if (nrow(kin_dat) != 44L) warning("Expected 44 complete KINARM cases; found ", nrow(kin_dat), ".")

# Preserve the original analysis exactly: explicitly standardize first, then
# call prcomp() without additional centering or scaling.
kinarm_scaled <- scale(kin_dat)
kin_pca <- prcomp(kinarm_scaled, center = FALSE, scale. = FALSE)
kin_loadings <- kin_pca$rotation[, 1:5, drop = FALSE]
kin_scores <- kin_pca$x[, 1:5, drop = FALSE]  # Unstandardized outcome scores, as originally modeled.
colnames(kin_loadings) <- paste0("PC", 1:5)
colnames(kin_scores) <- paste0("PC", 1:5)

for (j in 1:5) {
  dat[[paste0("KINARM_PC", j)]] <- NA_real_
  dat[kin_complete, paste0("KINARM_PC", j)] <- kin_scores[, j]
}

kin_variance <- kin_pca$sdev^2 / sum(kin_pca$sdev^2)
kin_variance_table <- tibble(
  Component = paste0("PC", seq_along(kin_variance)),
  Eigenvalue = kin_pca$sdev^2,
  Variance = kin_variance,
  `Variance (%)` = 100 * kin_variance,
  `Cumulative variance (%)` = 100 * cumsum(kin_variance)
)
safe_write_csv(kin_variance_table, file.path(output_root, "pca", "kinarm_pca_variance.csv"))
safe_write_csv(
  as.data.frame(kin_loadings) %>% rownames_to_column("Variable") %>%
    mutate(Label = unname(kinarm_labels[Variable]), .after = Variable),
  file.path(output_root, "pca", "kinarm_pca_loadings.csv")
)

log_message("Running KINARM PCA bootstrap stability assessment.")
kin_boot <- bootstrap_pca_stability(
  kin_dat, kin_loadings, kin_variance[1:5], n_boot = n_boot, seed = 123L
)
log_message("KINARM PCA bootstrap successful iterations: ", kin_boot$successful_iterations)
safe_write_csv(kin_boot$component, file.path(output_root, "tables", "supplementary_table_s2a_kinarm_component_stability.csv"))
kin_loading_all <- kin_boot$loading %>%
  mutate(Label = unname(kinarm_labels[Variable]), .after = Variable)
safe_write_csv(kin_loading_all, file.path(output_root, "pca", "kinarm_pca_bootstrap_loading_stability_all.csv"))
safe_write_csv(
  kin_loading_all %>% filter(abs(`Original loading`) >= 0.20) %>%
    select(Component, Variable = Label, `Original loading`, `Bootstrap median`, `95% lower`, `95% upper`, `Same-sign resamples (%)`),
  file.path(output_root, "tables", "supplementary_table_s2b_kinarm_loading_stability.csv")
)

kin_display <- loading_display_wide(
  kin_loadings, 100 * kin_variance[1:5], labels = kinarm_labels,
  threshold = 0.20, strict = TRUE
)
safe_write_csv(kin_display$wide, file.path(output_root, "tables", "supplementary_table_s2c_kinarm_loadings_display.csv"))
safe_write_csv(
  kin_display$long %>% select(Component, Variable, Label, Loading),
  file.path(output_root, "tables", "supplementary_table_s2c_kinarm_loadings_long.csv")
)

log_message("Running KINARM permutation parallel analysis.")
kin_parallel95 <- parallel_analysis(kin_dat, n_rep = n_parallel, seed = 123L)
kin_scree_data <- make_scree_plot(
  kin_pca$sdev^2, kin_parallel95,
  file.path(output_root, "figures", "supplementary_figure_s2_kinarm_pca_scree.png"),
  max_components = 20L
)
safe_write_csv(kin_scree_data, file.path(output_root, "pca", "kinarm_pca_scree_data.csv"))

# Save the complete deidentified analysis dataset with reproduced PC scores.
safe_write_csv(dat, file.path(output_root, "data", "analysis_dataset_with_reproduced_pc_scores.csv"))

# -----------------------------------------------------------------------------
# Participant-level supplementary displays
# -----------------------------------------------------------------------------

log_message("Generating participant-level supplementary figures.")
region_labels <- setNames(clean_region_label(brain_vars), brain_vars)

heatmap_df <- dat %>%
  arrange(Rating, Participant_ID) %>%
  mutate(Participant_Order = row_number(), Participant_Label = as.character(Participant_ID)) %>%
  select(Participant_Order, Participant_Label, Rating, all_of(brain_vars)) %>%
  pivot_longer(all_of(brain_vars), names_to = "Region", values_to = "Percentile") %>%
  mutate(
    Region_Label = factor(region_labels[Region], levels = region_labels[brain_vars]),
    Participant_Label = factor(Participant_Label, levels = rev(unique(Participant_Label))),
    Centered_Percentile = Percentile - 50
  )

rating_boundaries <- dat %>%
  arrange(Rating, Participant_ID) %>%
  count(Rating) %>%
  mutate(Cumulative = cumsum(n), Boundary = nrow(dat) - Cumulative + 0.5) %>%
  filter(Cumulative < nrow(dat))

p_heatmap <- ggplot(heatmap_df, aes(x = Region_Label, y = Participant_Label, fill = Centered_Percentile)) +
  geom_tile() +
  geom_hline(data = rating_boundaries, aes(yintercept = Boundary), linewidth = 0.7, inherit.aes = FALSE) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-50, 50), name = "Percentile points\nrelative to 50"
  ) +
  labs(
    title = "Individual-level population-referenced regional brain-volume percentiles",
    x = NULL, y = "Deidentified participant (ordered by performance rating)"
  ) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1),
    axis.text.y = element_text(size = 5),
    panel.grid = element_blank()
  )
ggsave(
  file.path(output_root, "figures", "supplementary_figure_s3_case_level_heatmap.png"),
  p_heatmap, width = 12, height = 9, dpi = 300
)

regional_df <- dat %>%
  select(Participant_ID, Performance_Group, all_of(brain_vars)) %>%
  pivot_longer(all_of(brain_vars), names_to = "Region", values_to = "Percentile") %>%
  mutate(
    Region_Label = factor(region_labels[Region], levels = region_labels[brain_vars]),
    Performance_Group = factor(Performance_Group, levels = c("High", "Mid", "Low"))
  )

p_regional <- ggplot(regional_df, aes(x = Region_Label, y = Percentile)) +
  geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.65) +
  geom_jitter(width = 0.14, height = 0, size = 0.8, alpha = 0.65) +
  facet_grid(Performance_Group ~ ., scales = "free_y") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Population-referenced percentile") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1), panel.grid.minor = element_blank())
ggsave(
  file.path(output_root, "figures", "supplementary_figure_s4_regional_distributions.png"),
  p_regional, width = 12, height = 8, dpi = 300
)

# -----------------------------------------------------------------------------
# Bayesian model fitting
# -----------------------------------------------------------------------------

brain_terms <- paste0("Brain_PC", 1:5)
base_priors_ordinal <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(normal(0, 2), class = "Intercept")
)
base_priors_binary <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(normal(0, 1.5), class = "Intercept")
)
base_priors_gaussian <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(normal(0, 1.5), class = "Intercept")
)

model_registry <- list()

fit_model <- function(name, formula, data, family, priors, control, pp_type) {
  model_path <- file.path(output_root, "models", paste0(name, ".rds"))
  summary_path <- file.path(output_root, "model_summaries", paste0(name, ".txt"))
  pp_path <- file.path(output_root, "figures", paste0("pp_check_", name, ".png"))

  log_message("Fitting model: ", name, " (N=", nrow(model.frame(formula, data = data, na.action = na.omit)), ")")
  brm_args <- list(
    formula = formula,
    data = data,
    family = family,
    prior = priors,
    chains = brms_chains,
    cores = analysis_cores,
    iter = brms_iter,
    warmup = brms_warmup,
    seed = 123,
    control = control,
    save_pars = save_pars(all = TRUE),
    refresh = max(1L, floor(brms_iter / 10L))
  )
  if (nzchar(brms_backend)) brm_args$backend <- brms_backend
  fit <- do.call(brm, brm_args)
  saveRDS(fit, model_path)

  capture.output(summary(fit), file = summary_path)
  tryCatch({
    # Posterior-predictive subsampling is explicitly seeded so diagnostic plots
    # are reproducible rather than depending on parallel Stan RNG state.
    pp_seed <- 10000L + sum(utf8ToInt(name))
    set.seed(pp_seed)
    pp <- pp_check(fit, type = pp_type, ndraws = 100)
    ggsave(pp_path, pp, width = 7, height = 5, dpi = 200)
  }, error = function(e) {
    log_message("Posterior predictive plot failed for ", name, ": ", conditionMessage(e))
  })

  model_registry[[name]] <<- fit
  fit
}

performance_formula <- as.formula(
  paste("Performance_Group ~", paste(c(brain_terms, "Sex"), collapse = " + "))
)
rating5_formula <- as.formula(
  paste("Performance_Rating5 ~", paste(c(brain_terms, "Sex"), collapse = " + "))
)
high_binary_formula <- as.formula(
  paste("High_Performance ~", paste(c(brain_terms, "Sex"), collapse = " + "))
)
top_binary_formula <- as.formula(
  paste("Top_Performance ~", paste(c(brain_terms, "Sex"), collapse = " + "))
)
performance_hockey_formula <- as.formula(
  paste("Performance_Group ~", paste(c(brain_terms, "Sex", "Hockey"), collapse = " + "))
)
concussion_formula <- as.formula(
  paste("Concussion_Group ~", paste(c(brain_terms, "Sex"), collapse = " + "))
)

fit_performance <- fit_model(
  "performance_primary", performance_formula, dat, cumulative("logit"),
  base_priors_ordinal, list(adapt_delta = 0.97), "bars"
)
fit_rating5 <- fit_model(
  "performance_rating5", rating5_formula, dat, cumulative("logit"),
  base_priors_ordinal, list(adapt_delta = 0.97), "bars"
)
fit_high_binary <- fit_model(
  "performance_high_binary", high_binary_formula, dat, bernoulli("logit"),
  base_priors_binary, list(adapt_delta = 0.97), "bars"
)
fit_top_binary <- fit_model(
  "performance_top_binary", top_binary_formula, dat, bernoulli("logit"),
  base_priors_binary, list(adapt_delta = 0.97), "bars"
)
fit_performance_hockey <- fit_model(
  "performance_hockey_adjusted", performance_hockey_formula,
  filter(dat, !is.na(Hockey)), cumulative("logit"),
  base_priors_ordinal, list(adapt_delta = 0.97), "bars"
)
fit_performance_nonhockey <- fit_model(
  "performance_nonhockey", performance_formula,
  filter(dat, Hockey == 0), cumulative("logit"),
  base_priors_ordinal, list(adapt_delta = 0.97), "bars"
)
fit_concussion <- fit_model(
  "concussion", concussion_formula, dat, bernoulli("logit"),
  base_priors_binary, list(adapt_delta = 0.99, max_treedepth = 15), "bars"
)

kinarm_fits <- vector("list", 5)
for (j in 1:5) {
  formula_j <- as.formula(
    paste0("KINARM_PC", j, " ~ ", paste(c(brain_terms, "Sex"), collapse = " + "))
  )
  kinarm_fits[[j]] <- fit_model(
    paste0("kinarm_pc", j), formula_j, dat, gaussian(),
    base_priors_gaussian, list(adapt_delta = 0.99, max_treedepth = 15), "dens_overlay"
  )
}

# -----------------------------------------------------------------------------
# Model coefficient tables and manuscript analytic figures
# -----------------------------------------------------------------------------

term_label <- function(term) {
  dplyr::recode(
    term,
    "Brain_PC1" = "Brain PC1", "Brain_PC2" = "Brain PC2", "Brain_PC3" = "Brain PC3",
    "Brain_PC4" = "Brain PC4", "Brain_PC5" = "Brain PC5",
    "SexM" = "Sex: M", "Hockey" = "Hockey",
    .default = term
  )
}

summarize_model <- function(model, exponentiate = FALSE) {
  fixed <- as.data.frame(fixef(model, probs = c(0.025, 0.975))) %>%
    rownames_to_column("Term") %>%
    filter(!str_detect(Term, "Intercept")) %>%
    transmute(
      Term,
      Predictor = term_label(Term),
      Estimate = Estimate,
      Posterior_SD = Est.Error,
      Lower = Q2.5,
      Upper = Q97.5
    )
  draws <- as_draws_df(model)
  fixed$Pr_direction <- vapply(fixed$Term, function(term) {
    draw_name <- paste0("b_", term)
    if (!draw_name %in% names(draws)) return(NA_real_)
    p_positive <- mean(draws[[draw_name]] > 0)
    max(p_positive, 1 - p_positive)
  }, numeric(1))
  fixed$N <- nobs(model)
  if (exponentiate) {
    fixed <- fixed %>% mutate(OR = exp(Estimate), OR_Lower = exp(Lower), OR_Upper = exp(Upper))
  }
  fixed
}

format_coefficient_table <- function(x, odds_ratio = FALSE) {
  out <- x %>%
    transmute(
      Predictor,
      `Coefficient (95% CrI)` = sprintf("%.2f (%.2f to %.2f)", Estimate, Lower, Upper),
      `Pr(direction)` = sprintf("%.2f", Pr_direction)
    )
  if (odds_ratio) {
    out <- x %>%
      transmute(
        Predictor,
        `Log-odds coefficient (95% CrI)` = sprintf("%.2f (%.2f to %.2f)", Estimate, Lower, Upper),
        `Odds ratio (95% CrI)` = sprintf("%.2f (%.2f to %.2f)", OR, OR_Lower, OR_Upper),
        `Pr(direction)` = sprintf("%.2f", Pr_direction)
      )
  }
  out
}

performance_results <- summarize_model(fit_performance, exponentiate = TRUE)
rating5_results <- summarize_model(fit_rating5, exponentiate = TRUE)
high_results <- summarize_model(fit_high_binary, exponentiate = TRUE)
top_results <- summarize_model(fit_top_binary, exponentiate = TRUE)
hockey_results <- summarize_model(fit_performance_hockey, exponentiate = TRUE)
nonhockey_results <- summarize_model(fit_performance_nonhockey, exponentiate = TRUE)
concussion_results <- summarize_model(fit_concussion, exponentiate = TRUE)
kinarm_results <- map_dfr(seq_along(kinarm_fits), function(j) {
  summarize_model(kinarm_fits[[j]], exponentiate = FALSE) %>%
    mutate(Outcome = paste0("KINARM PC", j), .before = 1)
})

safe_write_csv(performance_results, file.path(output_root, "tables", "supplementary_table_s3_performance_model_raw.csv"))
safe_write_csv(
  format_coefficient_table(performance_results, odds_ratio = TRUE),
  file.path(output_root, "tables", "supplementary_table_s3_performance_model_formatted.csv")
)

s4a_raw <- bind_rows(
  mutate(rating5_results, Model = "Five-level ordered rating", .before = 1),
  mutate(high_results, Model = "High-performance binary rating", .before = 1),
  mutate(top_results, Model = "Top-rating binary rating", .before = 1)
)
safe_write_csv(s4a_raw, file.path(output_root, "tables", "supplementary_table_s4a_performance_grouping_sensitivity_raw.csv"))
safe_write_csv(
  s4a_raw %>% transmute(
    Model, Predictor,
    `Coefficient (95% CrI)` = sprintf("%.2f (%.2f to %.2f)", Estimate, Lower, Upper),
    `OR (95% CrI)` = sprintf("%.2f (%.2f to %.2f)", OR, OR_Lower, OR_Upper),
    `Pr(direction)` = sprintf("%.2f", Pr_direction)
  ),
  file.path(output_root, "tables", "supplementary_table_s4a_performance_grouping_sensitivity_formatted.csv")
)

s4b_raw <- bind_rows(
  mutate(hockey_results, Model = "Hockey-adjusted primary model (known sport; N=49)", .before = 1),
  mutate(nonhockey_results, Model = "Non-hockey subset (N=27)", .before = 1)
)
safe_write_csv(s4b_raw, file.path(output_root, "tables", "supplementary_table_s4b_sport_sensitivity_raw.csv"))
safe_write_csv(
  s4b_raw %>% transmute(
    Model, Predictor,
    `Coefficient (95% CrI)` = sprintf("%.2f (%.2f to %.2f)", Estimate, Lower, Upper),
    `OR (95% CrI)` = sprintf("%.2f (%.2f to %.2f)", OR, OR_Lower, OR_Upper),
    `Pr(direction)` = sprintf("%.2f", Pr_direction)
  ),
  file.path(output_root, "tables", "supplementary_table_s4b_sport_sensitivity_formatted.csv")
)

safe_write_csv(kinarm_results, file.path(output_root, "tables", "supplementary_table_s5_kinarm_models_raw.csv"))
safe_write_csv(
  kinarm_results %>% transmute(
    Outcome, Predictor,
    `Coefficient (95% CrI)` = sprintf("%.2f (%.2f to %.2f)", Estimate, Lower, Upper),
    `Pr(direction)` = sprintf("%.2f", Pr_direction)
  ),
  file.path(output_root, "tables", "supplementary_table_s5_kinarm_models_formatted.csv")
)

safe_write_csv(concussion_results, file.path(output_root, "tables", "supplementary_table_s6_concussion_model_raw.csv"))
safe_write_csv(
  format_coefficient_table(concussion_results, odds_ratio = TRUE),
  file.path(output_root, "tables", "supplementary_table_s6_concussion_model_formatted.csv")
)

# Figure 2: performance forest plot.
plot_order <- rev(c(paste0("Brain PC", 1:5), "Sex: M"))
p_perf <- performance_results %>%
  mutate(Predictor = factor(Predictor, levels = plot_order)) %>%
  ggplot(aes(x = Predictor, y = OR, ymin = OR_Lower, ymax = OR_Upper)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(width = 0.18) +
  geom_point(size = 2) +
  coord_flip() +
  scale_y_log10() +
  labs(x = NULL, y = "Odds ratio (95% credible interval)") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(output_root, "figures", "figure_2_performance_forest.png"), p_perf, width = 7, height = 4.5, dpi = 300)

# Figure 3: KINARM coefficient forest plots.
p_kin <- kinarm_results %>%
  mutate(
    Predictor = factor(Predictor, levels = plot_order),
    Outcome = factor(Outcome, levels = paste0("KINARM PC", 1:5))
  ) %>%
  ggplot(aes(x = Predictor, y = Estimate, ymin = Lower, ymax = Upper)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(width = 0.18) +
  geom_point(size = 1.8) +
  coord_flip() +
  facet_wrap(~ Outcome, ncol = 2, scales = "free_y") +
  labs(x = NULL, y = "Posterior mean coefficient (95% credible interval)") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(output_root, "figures", "figure_3_kinarm_forest.png"), p_kin, width = 9, height = 10, dpi = 300)

# Figure 5: concussion-history forest plot.
p_conc <- concussion_results %>%
  mutate(Predictor = factor(Predictor, levels = plot_order)) %>%
  ggplot(aes(x = Predictor, y = OR, ymin = OR_Lower, ymax = OR_Upper)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(width = 0.18) +
  geom_point(size = 2) +
  coord_flip() +
  scale_y_log10() +
  labs(x = NULL, y = "Odds ratio (95% credible interval)") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(output_root, "figures", "figure_5_concussion_forest.png"), p_conc, width = 7, height = 4.5, dpi = 300)

# -----------------------------------------------------------------------------
# Model diagnostics and run manifest
# -----------------------------------------------------------------------------

log_message("Exporting model diagnostics.")

# Use posterior's default summary measures so the diagnostic columns have
# stable names across posterior package versions. Passing namespaced functions
# directly through ... can produce columns such as `posterior::rhat`; dplyr
# then interprets bare `rhat` as the function rather than a numeric column.
diagnostics <- imap_dfr(model_registry, function(model, name) {
  model_diagnostics <- posterior::summarise_draws(
    posterior::as_draws_array(model)
  )

  required_diagnostic_columns <- c(
    "variable", "mean", "sd", "rhat", "ess_bulk", "ess_tail"
  )
  missing_diagnostic_columns <- setdiff(
    required_diagnostic_columns, names(model_diagnostics)
  )
  if (length(missing_diagnostic_columns) > 0L) {
    stop(
      "Model diagnostic summary for '", name,
      "' is missing required columns: ",
      paste(missing_diagnostic_columns, collapse = ", "),
      ". Available columns: ",
      paste(names(model_diagnostics), collapse = ", ")
    )
  }

  model_diagnostics %>%
    select(all_of(required_diagnostic_columns)) %>%
    mutate(Model = name, .before = 1)
})
safe_write_csv(diagnostics, file.path(output_root, "diagnostics", "all_model_diagnostics.csv"))

finite_max_or_na <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

finite_min_or_na <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

diagnostic_summary <- diagnostics %>%
  summarise(
    Maximum_Rhat = finite_max_or_na(.data$rhat),
    Minimum_Bulk_ESS = finite_min_or_na(.data$ess_bulk),
    Minimum_Tail_ESS = finite_min_or_na(.data$ess_tail)
  )
safe_write_csv(diagnostic_summary, file.path(output_root, "diagnostics", "model_diagnostic_summary.csv"))

sample_sizes <- tibble(
  Analysis = c(
    "Brain PCA", "KINARM PCA", "Primary performance model", "Five-level performance model",
    "High-performance binary model", "Top-rating binary model", "Hockey-adjusted performance model",
    "Non-hockey performance model", "Concussion-history model",
    paste0("KINARM PC", 1:5, " model")
  ),
  N = c(
    nrow(brain_mat), nrow(kin_dat), nobs(fit_performance), nobs(fit_rating5),
    nobs(fit_high_binary), nobs(fit_top_binary), nobs(fit_performance_hockey),
    nobs(fit_performance_nonhockey), nobs(fit_concussion),
    vapply(kinarm_fits, nobs, numeric(1))
  )
)
safe_write_csv(sample_sizes, file.path(output_root, "diagnostics", "analysis_sample_sizes.csv"))

capture.output(sessionInfo(), file = file.path(output_root, "sessionInfo.txt"))

manifest_lines <- c(
  paste0("Run completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Repository root: ", repo_root),
  paste0("Input dataset: ", data_path),
  paste0("Input rows: ", nrow(dat)),
  paste0("PCA bootstrap resamples requested: ", n_boot),
  paste0("Brain PCA bootstrap successful: ", brain_boot$successful_iterations),
  paste0("KINARM PCA bootstrap successful: ", kin_boot$successful_iterations),
  paste0("Parallel-analysis repetitions: ", n_parallel),
  paste0("Stan chains: ", brms_chains),
  paste0("Stan iterations: ", brms_iter),
  paste0("Stan warmup: ", brms_warmup),
  paste0("Stan backend: ", ifelse(nzchar(brms_backend), brms_backend, "brms default")),
  "Prespecified descriptive Pr(direction) threshold: 0.95",
  "Concussion-history definition: original validated Concussion_Group variable (0/1) supplied in the deidentified dataset.",
  "PCA signs are the unmodified signs returned by the original prcomp workflows; bootstrap components are sign-aligned only to their corresponding original components.",
  "KINARM model outcomes use unstandardized PCA scores; brain predictors are standardized PCA scores.",
  "Figures 1 and 4 are anatomical renderings and are not generated from the tabular post-segmentation dataset."
)
writeLines(manifest_lines, file.path(output_root, "run_manifest.txt"))

log_message("Analysis complete. Generated outputs are in ", output_root, ".")
log_message("Review sessionInfo.txt and model diagnostics before release.")
