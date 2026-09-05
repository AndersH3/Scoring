#!/usr/bin/env Rscript
# ECPE stress tests for Anders Hellstrom's current two-pad scoring system.
# Version 1.0.0, 2026-09-01. Only required add-on package: edmdata.
# Source specification: https://github.com/AndersH3/Scoring/blob/main/item_analysis.r
# Source Git blob: a336d1b3d9964bd5f2e0f90f214e15d540be827d
# This implements the current Somers Dxy version, NOT the older rank-cut version.
# See --help and ECPE_STRESS_TEST_README.md for design and interpretation.

fail <- function(...) stop(sprintf(...), call. = FALSE)
assert <- function(ok, msg) if (!isTRUE(ok)) fail("%s", msg)

defaults <- function() list(
  out = "ecpe_stress_results", seed = 20260901L, bootstrap = 500L,
  repeats = 100L, permutations = 500L, folds = 5L, train_fraction = 0.70,
  sizes = c(50L, 100L, 250L, 500L, 1000L, 2000L),
  score_type = "rest", quick = FALSE, self_test = FALSE, help = FALSE
)

parse_args <- function(args) {
  o <- defaults(); i <- 1L
  flags <- c("quick", "self-test", "help")
  while (i <= length(args)) {
    a <- args[i]
    if (!startsWith(a, "--")) fail("Expected an option, got: %s", a)
    a <- substring(a, 3L)
    if (a %in% flags) {
      o[[gsub("-", "_", a)]] <- TRUE; i <- i + 1L; next
    }
    key <- gsub("-", "_", a)
    if (!(key %in% names(o))) fail("Unknown option --%s", a)
    if (i == length(args) || startsWith(args[i + 1L], "--"))
      fail("Missing value for --%s", a)
    value <- args[i + 1L]
    if (key == "sizes") value <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]]))
    else if (is.numeric(o[[key]])) value <- suppressWarnings(as.numeric(value))
    o[[key]] <- value; i <- i + 2L
  }
  if (o$quick) { o$bootstrap <- 30L; o$repeats <- 10L; o$permutations <- 30L }
  for (k in c("seed", "bootstrap", "repeats", "permutations", "folds", "sizes")) {
    v <- o[[k]]
    assert(length(v) > 0L && all(is.finite(v)) && all(v == floor(v)) &&
             all(v >= if (k == "seed") 0 else 2) && all(v <= .Machine$integer.max),
           paste("Invalid integer option:", k))
    o[[k]] <- as.integer(v)
  }
  assert(o$train_fraction >= 0.2 && o$train_fraction <= 0.9 &&
           is.finite(o$train_fraction), "--train-fraction must be between 0.2 and 0.9")
  assert(o$score_type %in% c("rest", "total"), "--score-type must be rest or total")
  assert(nzchar(o$out), "--out must be nonempty")
  o$sizes <- sort(unique(o$sizes))
  o
}

help <- function() cat(paste0(
  "Usage: Rscript ecpe_scoring_stress_test.R [options]\n",
  "  --out DIR              New or empty output directory [ecpe_stress_results]\n",
  "  --seed N               Reproducible RNG seed [20260901]\n",
  "  --bootstrap N          Calibration-person bootstrap replicates [500]\n",
  "  --repeats N            Replicates per sensitivity condition [100]\n",
  "  --permutations N       Independently shuffled-column null replicates [500]\n",
  "  --sizes 50,100,...      Calibration sample sizes [50,100,250,500,1000,2000]\n",
  "  --folds N              Person-level cross-fitting folds [5]\n",
  "  --train-fraction P     Fraction for fixed calibration/audit split [0.70]\n",
  "  --score-type TYPE      rest (current default) or total\n",
  "  --quick                Override counts to 30/10/30; smoke test only\n",
  "  --self-test            Run mathematical/invariance tests and exit\n",
  "  --help                 Show this help\n\n",
  "Install: install.packages('edmdata', repos='https://cloud.r-project.org')\n",
  "Normal execution also runs self-tests. No packages install automatically.\n",
  "Only complete binary data are accepted. Constant items are retained.\n"
))

validate_x <- function(x) {
  assert(is.matrix(x) && (is.numeric(x) || is.logical(x)), "Responses must be a numeric matrix")
  assert(nrow(x) >= 2L && ncol(x) >= 2L, "Need at least two persons and two items")
  assert(!anyNA(x) && all(is.finite(x)) && all(x %in% c(0, 1)),
         "Responses must contain only 0/1; missing values are not imputed or dropped")
  if (is.null(colnames(x))) colnames(x) <- sprintf("item%02d", seq_len(ncol(x)))
  assert(!anyDuplicated(colnames(x)), "Item names must be unique")
  storage.mode(x) <- "integer"
  x
}

# Exact pair count using the small integer score support. Half credit for ties.
# C = P(score_of_correct > score_of_incorrect) + .5 P(equal), Dxy = 2*C - 1.
# This equals Hmisc::somers2(score, response), and is checked against independent
# pair enumeration and the average-rank identity in the built-in self-tests.
association <- function(y, score) {
  n1 <- sum(y); n0 <- length(y) - n1
  if (n1 == 0L || n0 == 0L)
    return(c(C = NA_real_, Dxy = NA_real_, tie_fraction = NA_real_))
  if (all(score >= 0 & score == floor(score)) && max(score) <= 10000) {
    bins <- as.integer(max(score) + 1L)
    z <- tabulate(as.integer(score[y == 0L]) + 1L, nbins = bins)
    p <- tabulate(as.integer(score[y == 1L]) + 1L, nbins = bins)
    pairs <- as.double(n0) * n1
    tied <- sum(as.double(p) * z)
    C <- sum(p * (cumsum(z) - z + 0.5 * z)) / pairs
    tie_fraction <- tied / pairs
  } else {
    C <- (sum(rank(score, ties.method = "average")[y == 1L]) - n1 * (n1 + 1) / 2) / (as.double(n0) * n1)
    levels <- sort(unique(score)); g <- match(score, levels)
    tie_fraction <- sum(as.double(tabulate(g[y == 0L], length(levels))) *
                          tabulate(g[y == 1L], length(levels))) / (as.double(n0) * n1)
  }
  C <- min(1, max(0, C))
  c(C = C, Dxy = 2 * C - 1, tie_fraction = tie_fraction)
}

fit_weights <- function(x, score_type = "rest") {
  n <- nrow(x); j <- ncol(x); correct <- colSums(x); total <- rowSums(x)
  a <- t(vapply(seq_len(j), function(k)
    association(x[, k], if (score_type == "rest") total - x[, k] else total), numeric(3)))
  dplus <- pmax(0, ifelse(is.na(a[, "Dxy"]), 0, a[, "Dxy"]))
  pad1 <- (correct + 1) / (n + 1)
  pad2 <- (n - n * dplus + 1) / (n + 1)
  rarity <- 1 - log(pad1); discrimination <- 1 - log(pad2)
  raw <- rarity * discrimination
  assert(all(is.finite(raw)) && all(raw > 0), "Invalid item weights")
  data.frame(item = colnames(x), n = n, correct = correct, incorrect = n - correct,
             proportion_correct = correct / n, C = a[, "C"], Dxy = a[, "Dxy"],
             tied_pair_fraction = a[, "tie_fraction"], Cplus = dplus,
             concordance_pair_numerator = n * dplus, pad1 = pad1, pad2 = pad2,
             rarity_factor = rarity, discrimination_factor = discrimination,
             raw_weight = raw, operational_weight = raw / min(raw), share = raw / sum(raw),
             constant = correct == 0 | correct == n,
             small_group = pmin(correct, n - correct) < 5, score_type = score_type,
             row.names = NULL, check.names = FALSE)
}

score_with <- function(x, fit) {
  assert(identical(colnames(x), fit$item), "Scoring item order does not match calibration item order")
  as.vector(x %*% fit$share)
}
safe_cor <- function(x, y, method = "pearson") {
  if (length(x) < 2L || anyNA(x) || anyNA(y) || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  unname(cor(x, y, method = method))
}
percentiles <- function(x) (rank(x, ties.method = "average") - 0.5) / length(x)
# Fractional membership at the top-decile boundary avoids arbitrary ID tie breaks.
top_membership <- function(x, fraction = 0.10) {
  k <- max(1L, ceiling(length(x) * fraction)); cutoff <- sort(x, decreasing = TRUE)[k]
  high <- x > cutoff; tied <- x == cutoff
  as.numeric(high) + as.numeric(tied) * (k - sum(high)) / sum(tied)
}
compare_scores <- function(reference, candidate) {
  shift <- abs(percentiles(candidate) - percentiles(reference))
  a <- top_membership(reference); b <- top_membership(candidate)
  c(pearson = safe_cor(reference, candidate), spearman = safe_cor(reference, candidate, "spearman"),
    rmse = sqrt(mean((candidate - reference)^2)), mean_change = mean(candidate - reference),
    max_abs_change = max(abs(candidate - reference)), mean_rank_shift_pp = 100 * mean(shift),
    p95_rank_shift_pp = 100 * unname(quantile(shift, 0.95)),
    max_rank_shift_pp = 100 * max(shift), top_decile_overlap = sum(pmin(a, b)) / sum(a))
}
weight_metrics <- function(fit, reference = NULL) {
  x <- c(max_share = max(fit$share), effective_items = 1 / sum(fit$share^2),
         operational_max = max(fit$operational_weight), constants = sum(fit$constant),
         small_groups = sum(fit$small_group), negative_items = sum(fit$Dxy < 0, na.rm = TRUE))
  if (!is.null(reference)) x <- c(x, share_L1 = sum(abs(fit$share - reference$share)),
                                 share_max_abs = max(abs(fit$share - reference$share)))
  x
}
quantiles <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(mean = NA, sd = NA, p025 = NA, median = NA, p975 = NA))
  c(mean = mean(x), sd = if (length(x) > 1L) sd(x) else NA_real_,
    setNames(as.numeric(quantile(x, c(.025, .5, .975))), c("p025", "median", "p975")))
}
bind <- function(x) { ans <- do.call(rbind, x); rownames(ans) <- NULL; ans }
row <- function(...) as.data.frame(list(...), stringsAsFactors = FALSE, check.names = FALSE)
summarise_by <- function(d, keys) {
  cols <- setdiff(names(d)[vapply(d, is.numeric, logical(1))], c(keys, "replicate"))
  groups <- split(seq_len(nrow(d)), interaction(d[keys], drop = TRUE, lex.order = TRUE))
  bind(lapply(groups, function(ii) bind(lapply(cols, function(k)
    cbind(d[ii[1L], keys, drop = FALSE], metric = k, as.data.frame(as.list(quantiles(d[ii, k]))))))))
}
write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "NA")

run_self_tests <- function() {
  eq <- function(a, b) assert(isTRUE(all.equal(a, b, tolerance = 1e-11, check.attributes = FALSE)),
                             "Mathematical self-test failed")
  eq(association(c(0, 0, 1, 1), 1:4)[1:2], c(1, 1))
  eq(association(c(0, 0, 1, 1), 4:1)[1:2], c(0, -1))
  eq(association(c(0, 1, 0, 1), rep(1, 4))[1:2], c(.5, 0))
  eq(association(c(0, 1, 0, 1), c(10, 12, 12, 14))[1:2], c(.875, .75))
  assert(is.na(association(rep(0, 4), 1:4)["C"]), "Constant C must be NA")
  # Exhaustive small tied-score cases, compared with actual cross-group pairs.
  yy <- as.matrix(expand.grid(rep(list(0:1), 5)))
  ss <- list(c(0, 1, 1, 2, 3), rep(2, 5), 4:0, c(.1, .2, .2, .7, .9))
  checks <- 0L
  for (i in seq_len(nrow(yy))) for (s in ss) {
    y <- yy[i, ]; if (all(y == y[1L])) next
    dif <- outer(s[y == 1], s[y == 0], "-")
    expected <- mean((dif > 0) + .5 * (dif == 0))
    eq(association(y, s)["C"], expected)
    nr <- sum(y); nz <- length(y) - nr
    eq(expected, (sum(rank(s)[y == 1]) - nr * (nr + 1) / 2) / (nr * nz))
    if (requireNamespace("Hmisc", quietly = TRUE))
      eq(association(y, s)[1:2], Hmisc::somers2(s, y)[c("C", "Dxy")])
    checks <- checks + 1L
  }
  x <- validate_x(cbind(a = c(0, 0, 1, 1), b = c(0, 1, 0, 1),
                       zero = 0L, one = 1L))
  f <- fit_weights(x)
  eq(min(f$operational_weight), 1)
  eq(f$raw_weight[f$item == "zero"], 1 + log(5))
  eq(f$raw_weight[f$item == "one"], 1)
  eq(f$discrimination_factor[f$constant], c(1, 1))
  eq(score_with(x, f), as.vector(x %*% f$operational_weight / sum(f$operational_weight)))
  eq(fit_weights(x[4:1, ])$share, f$share)
  eq(score_with(x[, 4:1], fit_weights(x[, 4:1])), score_with(x, f))
  eq(score_with(x * 0L, f), rep(0, 4)); eq(score_with(x * 0L + 1L, f), rep(1, 4))
  reverse <- validate_x(cbind(a = c(0, 0, 1, 1), b = c(1, 1, 0, 0)))
  eq(fit_weights(reverse)$Dxy, c(-1, -1))
  eq(fit_weights(reverse)$discrimination_factor, c(1, 1))
  perfect <- validate_x(cbind(a = c(0, 0, 1, 1), b = c(0, 0, 1, 1)))
  eq(fit_weights(perfect)$discrimination_factor, rep(1 + log(5), 2))
  assert(all(f$share > 0), "Fixed-weight scores must be increasing in each response")
  bad <- x; bad[1, 1] <- NA
  assert(inherits(try(validate_x(bad), silent = TRUE), "try-error"), "Missing data must be rejected")
  bad[1, 1] <- 2
  assert(inherits(try(validate_x(bad), silent = TRUE), "try-error"), "Nonbinary data must be rejected")
  eq(sum(top_membership(rep(1, 10))), 1)
  eq(compare_scores(1:10, 1:10)["top_decile_overlap"], 1)
  msg <- sprintf("PASS: %d exhaustive pair/rank comparisons plus scoring, ties, constants, reversal, normalization, ordering and input checks. Hmisc comparison: %s.",
                 checks, if (requireNamespace("Hmisc", quietly = TRUE)) "passed" else "not installed (independent pair/rank checks passed)")
  message(msg); msg
}

load_ecpe <- function() {
  if (!requireNamespace("edmdata", quietly = TRUE))
    fail("Install edmdata first: Rscript -e 'install.packages(\"edmdata\", repos=\"https://cloud.r-project.org\")'")
  e <- new.env(parent = baseenv())
  utils::data("items_ecpe", package = "edmdata", envir = e)
  utils::data("qmatrix_ecpe", package = "edmdata", envir = e)
  x <- validate_x(as.matrix(e$items_ecpe)); q <- as.matrix(e$qmatrix_ecpe)
  assert(identical(dim(x), c(2922L, 28L)), "Unexpected ECPE dimensions: expected 2922 x 28")
  assert(identical(dim(q), c(28L, 3L)) && !anyNA(q) && all(q %in% c(0, 1)), "Unexpected ECPE Q matrix")
  # The package supplies matching names. Do not silently assume a row alignment.
  assert(!is.null(rownames(q)) && setequal(rownames(q), colnames(x)),
         "Q-matrix row names must match response item names")
  q <- q[colnames(x), , drop = FALSE]
  list(x = x, q = q)
}

raw_reversals <- function(raw, weighted, tolerance = 1e-12) {
  reverse <- 0; comparable <- 0
  for (v in sort(unique(raw))) {
    low <- weighted[raw == v]; high <- sort(weighted[raw > v])
    comparable <- comparable + length(low) * length(high)
    if (length(high)) reverse <- reverse + sum(findInterval(low - tolerance, high))
  }
  row(comparable_pairs = comparable, strict_reversals = reverse,
      reversal_fraction = if (comparable) reverse / comparable else NA_real_)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  o <- parse_args(args)
  if (o$help) { help(); return(invisible(NULL)) }
  test_result <- run_self_tests()
  if (o$self_test) return(invisible(test_result))
  started <- Sys.time(); dat <- load_ecpe(); x <- dat$x; q <- dat$q
  n <- nrow(x); j <- ncol(x)
  assert(o$folds <= n / 2, "Too many folds")
  assert(!dir.exists(o$out) || length(list.files(o$out, all.files = TRUE, no.. = TRUE)) == 0,
         "Output directory exists and is not empty; choose a new --out directory")
  dir.create(o$out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(o$out, mustWork = TRUE)
  save_table <- function(d, name) { write_csv(d, file.path(out, paste0(name, ".csv"))); invisible(d) }
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(o$seed)
  initial_rng <- .Random.seed
  writeLines(capture.output(dput(o)), file.path(out, "settings.R"))
  writeLines(test_result, file.path(out, "self_tests.txt"))
  save_table(data.frame(person = seq_len(n), x, check.names = FALSE), "ecpe_responses")
  save_table(data.frame(item = rownames(q), q, check.names = FALSE), "ecpe_qmatrix")
  full <- fit_weights(x, o$score_type); rest <- fit_weights(x, "rest"); total <- fit_weights(x, "total")
  save_table(full, "item_weights_full"); save_table(rest, "item_weights_rest"); save_table(total, "item_weights_total")
  raw <- rowMeans(x); full_score <- score_with(x, full)
  full_comparison <- compare_scores(raw, full_score)
  save_table(as.data.frame(as.list(full_comparison)), "full_weighted_vs_raw")
  reversals <- save_table(raw_reversals(rowSums(x), full_score), "raw_score_reversals")
  grouped <- bind(lapply(sort(unique(rowSums(x))), function(k) {
    s <- full_score[rowSums(x) == k]
    row(number_correct = k, persons = length(s), min = min(s), mean = mean(s),
        max = max(s), sd = if (length(s) > 1) sd(s) else NA_real_)
  }))
  save_table(grouped, "score_by_number_correct")
  save_table(row(item = full$item, rest_Dxy = rest$Dxy, total_Dxy = total$Dxy,
                 rest_share = rest$share, total_share = total$share,
                 Dxy_increase = total$Dxy - rest$Dxy,
                 share_change = total$share - rest$share), "rest_vs_total")
  # Q skills overlap; shares in this table need not add to 1 across skills.
  skill <- bind(lapply(seq_len(ncol(q)), function(k) row(skill = colnames(q)[k],
    label = c("Morphosyntactic rules", "Cohesive rules", "Lexical rules")[k],
    items = sum(q[, k]), raw_share_coverage = mean(q[, k]),
    weighted_share_coverage = sum(full$share[q[, k] == 1]))))
  save_table(skill, "q_skill_coverage")

  message("Creating fixed calibration/audit split and person-level cross-fits...")
  train_id <- sort(sample.int(n, floor(n * o$train_fraction)))
  audit_id <- setdiff(seq_len(n), train_id); cal <- x[train_id, ]; audit <- x[audit_id, ]
  nc <- nrow(cal); na <- nrow(audit)
  ref <- fit_weights(cal, o$score_type); ref_score <- score_with(audit, ref)
  save_table(ref, "item_weights_calibration")
  folds <- sample(rep(seq_len(o$folds), length.out = n))
  cross <- numeric(n); cross_fits <- vector("list", o$folds); cross_metrics <- vector("list", o$folds)
  for (k in seq_len(o$folds)) {
    train <- folds != k; f <- fit_weights(x[train, ], o$score_type)
    cross[!train] <- score_with(x[!train, ], f)
    cross_fits[[k]] <- cbind(fold = k, f)
    cross_metrics[[k]] <- cbind(fold = k, as.data.frame(as.list(compare_scores(raw[!train], cross[!train]))))
  }
  save_table(bind(cross_fits), "crossfit_item_weights")
  save_table(bind(cross_metrics), "crossfit_within_fold_metrics")
  save_table(row(person = seq_len(n), calibration = seq_len(n) %in% train_id,
                 audit = seq_len(n) %in% audit_id, fold = folds, number_correct = rowSums(x),
                 raw_proportion = raw, full_weighted = full_score, crossfit_weighted = cross,
                 raw_midrank_percentile = 100 * percentiles(raw),
                 weighted_midrank_percentile = 100 * percentiles(full_score)), "person_scores")

  message(sprintf("Bootstrap: %d calibration-person resamples; scoring %d untouched audit persons...", o$bootstrap, na))
  bshare <- braw <- boper <- matrix(NA_real_, o$bootstrap, j)
  bscore <- matrix(NA_real_, o$bootstrap, na); bdxy <- matrix(NA_real_, o$bootstrap, j)
  bmetrics <- vector("list", o$bootstrap)
  for (b in seq_len(o$bootstrap)) {
    f <- fit_weights(cal[sample.int(nc, nc, replace = TRUE), ], o$score_type)
    bshare[b, ] <- f$share; braw[b, ] <- f$raw_weight; boper[b, ] <- f$operational_weight
    bdxy[b, ] <- f$Dxy; bscore[b, ] <- score_with(audit, f)
    bmetrics[[b]] <- cbind(replicate = b, as.data.frame(as.list(c(compare_scores(ref_score, bscore[b, ]), weight_metrics(f, ref)))))
  }
  boot_summary <- bind(lapply(seq_len(j), function(k) {
    bind(lapply(c("share", "raw_weight", "operational_weight", "Dxy"), function(v) {
      values <- switch(v, share = bshare[, k], raw_weight = braw[, k], operational_weight = boper[, k], Dxy = bdxy[, k])
      cbind(item = ref$item[k], quantity = v, estimate = ref[[v]][k], valid_replicates = sum(is.finite(values)),
            as.data.frame(as.list(quantiles(values))))
    }))
  }))
  save_table(boot_summary, "bootstrap_item_intervals")
  save_table(bind(bmetrics), "bootstrap_audit_metrics")
  person_intervals <- cbind(person = audit_id, reference_score = ref_score,
                            as.data.frame(t(apply(bscore, 2, quantiles))))
  save_table(person_intervals, "bootstrap_audit_score_intervals")

  message("Testing calibration sample sizes and selected-score populations...")
  sizes <- sort(unique(c(o$sizes[o$sizes <= nc], nc))); ss <- list(); ii <- 0L
  for (size in sizes) for (b in seq_len(o$repeats)) {
    f <- fit_weights(cal[sample.int(nc, size), , drop = FALSE], o$score_type)
    ii <- ii + 1L
    ss[[ii]] <- cbind(calibration_n = size, replicate = b,
      as.data.frame(as.list(c(compare_scores(ref_score, score_with(audit, f)), weight_metrics(f, ref)))))
  }
  sample_results <- bind(ss); save_table(sample_results, "sample_size_replicates")
  sample_summary <- save_table(summarise_by(sample_results, "calibration_n"), "sample_size_summary")
  cal_total <- rowSums(cal); cut_low <- unname(quantile(cal_total, .3)); cut_high <- unname(quantile(cal_total, .7))
  pools <- list(unselected = seq_len(nc), lower_30pct_with_ties = which(cal_total <= cut_low),
                upper_30pct_with_ties = which(cal_total >= cut_high))
  shift_n <- min(500L, nc); shifted <- list(); ii <- 0L
  for (label in names(pools)) for (b in seq_len(o$repeats)) {
    pool <- pools[[label]]; ids <- pool[sample.int(length(pool), shift_n, replace = TRUE)]
    f <- fit_weights(cal[ids, ], o$score_type); ii <- ii + 1L
    shifted[[ii]] <- cbind(selection = label, calibration_n = shift_n, pool_n = length(pool), replicate = b,
      as.data.frame(as.list(c(compare_scores(ref_score, score_with(audit, f)), weight_metrics(f, ref)))))
  }
  shifted <- bind(shifted); save_table(shifted, "population_shift_replicates")
  save_table(summarise_by(shifted, "selection"), "population_shift_summary")

  message("Testing response corruption in calibration only...")
  contaminated <- list(); ii <- 0L
  for (kind in c("cell_flip", "random_person_50pct", "all_zero_person", "all_one_person"))
    for (rate in c(.01, .05, .10)) for (b in seq_len(o$repeats)) {
      y <- cal
      if (kind == "cell_flip") {
        at <- sample.int(length(y), max(1L, round(rate * length(y)))); y[at] <- 1L - y[at]
      } else {
        at <- sample.int(nc, max(1L, round(rate * nc)))
        y[at, ] <- switch(kind, random_person_50pct = rbinom(length(at) * j, 1, .5),
                          all_zero_person = 0L, all_one_person = 1L)
      }
      f <- fit_weights(y, o$score_type); ii <- ii + 1L
      contaminated[[ii]] <- cbind(kind = kind, rate = rate, replicate = b,
        changed_cell_fraction = mean(y != cal),
        as.data.frame(as.list(c(compare_scores(ref_score, score_with(audit, f)), weight_metrics(f, ref)))))
    }
  contaminated <- bind(contaminated); save_table(contaminated, "contamination_replicates")
  save_table(summarise_by(contaminated, c("kind", "rate")), "contamination_summary")

  message("Permutation null: preserving each item's correct count, breaking cross-item dependence...")
  null_rest <- null_total <- matrix(NA_real_, o$permutations, j); null_metrics <- list()
  for (b in seq_len(o$permutations)) {
    y <- cal
    for (k in seq_len(j)) y[, k] <- cal[sample.int(nc), k]
    fr <- fit_weights(y, "rest"); ft <- fit_weights(y, "total")
    null_rest[b, ] <- fr$Dxy; null_total[b, ] <- ft$Dxy
    null_metrics[[b]] <- row(replicate = b,
      mean_Dxy_rest = mean(fr$Dxy, na.rm = TRUE), mean_Dxy_total = mean(ft$Dxy, na.rm = TRUE),
      mean_discrimination_factor_rest = mean(fr$discrimination_factor),
      mean_discrimination_factor_total = mean(ft$discrimination_factor),
      max_share_rest = max(fr$share), max_share_total = max(ft$share))
  }
  observed_rest <- fit_weights(cal, "rest"); observed_total <- fit_weights(cal, "total")
  null_table <- bind(lapply(c("rest", "total"), function(mode) {
    a <- if (mode == "rest") null_rest else null_total
    observed <- if (mode == "rest") observed_rest$Dxy else observed_total$Dxy
    p <- vapply(seq_len(j), function(k) if (is.na(observed[k])) NA_real_ else
      (1 + sum(a[, k] >= observed[k], na.rm = TRUE)) / (1 + sum(is.finite(a[, k]))), numeric(1))
    cbind(item = colnames(cal), score_type = mode, observed_Dxy = observed,
          one_sided_p = p, BH_adjusted_p = p.adjust(p, method = "BH"),
          as.data.frame(t(apply(a, 2, quantiles))))
  }))
  save_table(null_table, "permutation_item_tests"); save_table(bind(null_metrics), "permutation_null_metrics")

  message("Testing item deletion, constant/rare/reversed items, duplicates and sample replication...")
  deletion <- bind(lapply(seq_len(j), function(k) {
    ca <- cal[, -k, drop = FALSE]; au <- audit[, -k, drop = FALSE]
    f <- fit_weights(ca, o$score_type); new_score <- score_with(au, f)
    kept <- ref$share[-k]; kept <- kept / sum(kept); fixed_deleted <- as.vector(au %*% kept)
    cbind(deleted_item = colnames(cal)[k], deleted_share = ref$share[k],
          recalibration_only_rmse = sqrt(mean((new_score - fixed_deleted)^2)),
          share_L1_remaining = sum(abs(f$share - kept)),
          as.data.frame(as.list(compare_scores(ref_score, new_score))))
  }))
  save_table(deletion, "item_deletion")
  k <- which.max(ref$share); altered <- list(); altered_weights <- list()
  for (kind in c("append_zero", "append_one", "append_rare_1pct", "append_random_50pct",
                  "append_reversed", "append_one_duplicate", "append_five_duplicates")) {
    copies <- if (kind == "append_five_duplicates") 5L else 1L
    make_extra <- function(y) {
      m <- nrow(y)
      z <- switch(kind, append_zero = rep(0L, m), append_one = rep(1L, m),
        append_rare_1pct = rbinom(m, 1, .01), append_random_50pct = rbinom(m, 1, .5),
        append_reversed = 1L - y[, k], append_one_duplicate = y[, k], append_five_duplicates = y[, k])
      extra <- matrix(rep(z, copies), nrow = m, ncol = copies)
      colnames(extra) <- paste0("synthetic_", seq_len(copies)); cbind(y, extra)
    }
    ca <- make_extra(cal); au <- make_extra(audit); f <- fit_weights(ca, o$score_type)
    altered[[kind]] <- cbind(scenario = kind, copied_or_reversed_item = colnames(cal)[k],
      extra_share = sum(tail(f$share, copies)),
      as.data.frame(as.list(c(compare_scores(ref_score, score_with(au, f)), weight_metrics(f)))))
    altered_weights[[kind]] <- cbind(scenario = kind, f)
  }
  altered <- bind(altered); save_table(altered, "augmented_item_stress")
  save_table(bind(altered_weights), "augmented_item_weights")
  replication <- bind(lapply(c(1L, 2L, 5L, 10L, 20L), function(mult) {
    f <- fit_weights(cal[rep(seq_len(nc), mult), ], o$score_type)
    finite <- is.finite(f$Dxy) & is.finite(ref$Dxy)
    cbind(copies = mult, nominal_n = nc * mult,
      max_abs_Dxy_change = if (any(finite)) max(abs(f$Dxy[finite] - ref$Dxy[finite])) else NA_real_,
      as.data.frame(as.list(c(compare_scores(ref_score, score_with(audit, f)), weight_metrics(f, ref)))))
  }))
  save_table(replication, "sample_replication")

  message("Estimating random split-half consistency on the untouched audit sample...")
  split_half <- bind(lapply(seq_len(o$repeats), function(b) {
    a <- sort(sample.int(j, floor(j / 2))); z <- setdiff(seq_len(j), a)
    fa <- fit_weights(cal[, a], o$score_type); fz <- fit_weights(cal[, z], o$score_type)
    rw <- safe_cor(score_with(audit[, a], fa), score_with(audit[, z], fz))
    rr <- safe_cor(rowMeans(audit[, a]), rowMeans(audit[, z]))
    sb <- function(r) if (is.finite(r) && r > -1) 2 * r / (1 + r) else NA_real_
    row(replicate = b, weighted_half_r = rw, raw_half_r = rr,
        weighted_spearman_brown = sb(rw), raw_spearman_brown = sb(rr),
        corrected_difference = sb(rw) - sb(rr))
  }))
  save_table(split_half, "split_half_consistency")

  # Base-R multi-page vector PDF; no graphics add-on packages are required.
  pdf(file.path(out, "diagnostic_plots.pdf"), width = 10, height = 7, onefile = TRUE)
  tryCatch({
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    plot(rowSums(x), full_score, pch = 16, cex = .35, col = adjustcolor("#246B8E", .2),
         xlab = "Number correct", ylab = "Full-sample weighted score", main = "Raw versus weighted scores")
    lines(grouped$number_correct, grouped$mean, col = "#C65229", lwd = 2)
    abline(a = 0, b = 1 / j, lty = 2, col = "grey50")
    plot(full$proportion_correct, full$share, pch = 19, col = "#246B8E",
         xlab = "Proportion correct", ylab = "Normalized item weight", main = "Difficulty and weight")
    top_items <- order(full$share, decreasing = TRUE)[1:3]
    text(full$proportion_correct[top_items], full$share[top_items],
         labels = full$item[top_items], pos = 4, cex = .75)
    par(mfrow = c(1, 1), mar = c(7, 4, 3, 1))
    bs <- boot_summary[boot_summary$quantity == "share", ]
    plot(seq_len(j), bs$estimate, ylim = range(bs$p025, bs$p975, bs$estimate), xaxt = "n",
         pch = 19, col = "#246B8E", xlab = "", ylab = "Item share",
         main = "Calibration bootstrap: pointwise 95% percentile intervals")
    arrows(seq_len(j), bs$p025, seq_len(j), bs$p975, angle = 90, code = 3, length = .03, col = "#246B8E")
    axis(1, seq_len(j), bs$item, las = 2, cex.axis = .7)
    abline(h = 1 / j, lty = 2, col = "grey50")
    par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
    boxplot(rmse ~ factor(calibration_n), data = sample_results, col = "#D9E8ED",
            xlab = "Calibration sample size", ylab = "Audit score RMSE", main = "Calibration-size sensitivity")
    boxplot(p95_rank_shift_pp ~ factor(calibration_n), data = sample_results, col = "#D9E8ED",
            xlab = "Calibration sample size", ylab = "95th percentile rank shift (pp)", main = "Audit rank sensitivity")
    par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
    nm <- bind(null_metrics)
    boxplot(list(Rest = nm$mean_Dxy_rest, Total = nm$mean_Dxy_total), col = c("#D9E8ED", "#F1D4C9"),
            ylab = "Mean item Dxy", main = "Independent-column null")
    abline(h = 0, lty = 2, col = "grey50")
    boxplot(list(Raw = split_half$raw_spearman_brown, Weighted = split_half$weighted_spearman_brown),
            ylab = "Spearman-Brown corrected correlation", main = "Audit split-half consistency", col = "#D9E8ED")
    par(mfrow = c(1, 1), mar = c(8, 4, 3, 1))
    kind_label <- c(cell_flip = "Cell flips", random_person_50pct = "Random",
                    all_zero_person = "All zero", all_one_person = "All one")
    condition_label <- paste(kind_label[contaminated$kind],
                             paste0(100 * contaminated$rate, "%"), sep = " / ")
    condition <- factor(condition_label, levels = unique(condition_label))
    boxplot(contaminated$rmse ~ condition, las = 2,
            cex.axis = .8, col = "#D9E8ED", xlab = "", ylab = "Audit score RMSE",
            main = "Calibration corruption, original audit responses unchanged")
  }, finally = dev.off())

  median_boot_rmse <- median(bind(bmetrics)$rmse)
  worst <- sample_results[sample_results$calibration_n == min(sizes), ]
  report <- c(
    "# ECPE scoring stress-test run", "",
    sprintf("Completed: %s UTC. Program 1.0.0. Seed %d. Quick mode: %s.", format(Sys.time(), tz = "UTC"), o$seed, o$quick),
    sprintf("Data: %d examinees x %d binary ECPE grammar items; edmdata %s. Calibration: %d, audit: %d.", n, j, as.character(packageVersion("edmdata")), nc, na),
    sprintf("Score reference: %s. Bootstrap %d; sensitivity replicates %d; permutations %d; cross-fitting folds %d.",
            o$score_type, o$bootstrap, o$repeats, o$permutations, o$folds), "",
    "## Exact scoring specification", "",
    "C is cross-group concordance with half credit for tied corresponding scores; Dxy = 2*C - 1. Cplus = max(0,Dxy).",
    "For constant items empirical C and Dxy stay NA, while Cplus=0 is used for weighting.",
    "pad1=(correct+1)/(n+1); pad2=(n-n*Cplus+1)/(n+1); raw_weight=(1-log(pad1))*(1-log(pad2)).",
    "Operational weights = raw_weight/min(raw_weight); normalized score = sum(x*weight)/sum(weight). No rounding of n*Cplus.",
    "GitHub source blob: a336d1b3d9964bd5f2e0f90f214e15d540be827d.", "",
    "## Descriptive results", "",
    sprintf("- Full-sample weighted versus raw Spearman correlation: %.6f.", full_comparison["spearman"]),
    sprintf("- Strict reversals between persons with different raw totals: %.4f%% of %.0f comparable pairs.", 100 * reversals$reversal_fraction, reversals$comparable_pairs),
    sprintf("- Maximum full-sample item share: %.4f%%; effective item count 1/sum(share^2): %.3f.", 100 * max(full$share), 1 / sum(full$share^2)),
    sprintf("- Median audit score RMSE under calibration bootstrap: %.6f (0-1 score units).", median_boot_rmse),
    sprintf("- Median audit score RMSE with n=%d calibration persons: %.6f.", min(sizes), median(worst$rmse)),
    sprintf("- Median split-half corrected correlation: raw %.6f, weighted %.6f; median paired difference %.6f.",
            median(split_half$raw_spearman_brown, na.rm = TRUE), median(split_half$weighted_spearman_brown, na.rm = TRUE),
            median(split_half$corrected_difference, na.rm = TRUE)), "",
    "## Interpretation boundaries", "",
    "The full-sample scores are descriptive. Audit scores always use calibration-only weights. Bootstrap intervals measure calibration uncertainty conditional on these responses and this calibration population; they are not latent-ability or test-retest confidence intervals.",
    "The frozen audit reference is the original calibration fit, not a known true score. Agreement with it measures stability, not scoring accuracy. Neither high raw-score correlation nor high consistency establishes superiority or external validity.",
    "Percentile intervals are pointwise and can be irregular at the Dxy=0 truncation boundary. Permutation tests are one-sided for positive association, with +1 correction and BH adjustment separately by score type. Total-score nulls retain self-inclusion effects by design.",
    "Sample-size draws use no replacement within a draw. Population-shift draws use replacement from raw-score-selected pools; this selection itself changes association. It is not demographic DIF evidence. Contamination changes calibration responses only; 50% random responding is a specified artificial condition, not an ECPE guessing-rate estimate.",
    "Augmented-item and deletion scenarios change test content and the score denominator. Their score shifts cannot be interpreted solely as estimation errors. Duplicate rows in sample_replication.csv add no information: any change isolates the +1 padding's sample-size dependence.",
    "Random split halves are content-overlapping measurements of a multidimensional test. Spearman-Brown values are descriptive under its parallel-form assumptions, not a validated reliability coefficient for ECPE. Q-matrix skills overlap, so skill coverage values do not partition the score.",
    "Top-decile overlap uses fractional membership for boundary ties; rank shifts use average ranks, in percentage points. Exact raw reversals exclude raw-score ties and weighted differences <=1e-12. Cross-fitted scores from different folds use different weights; inspect within-fold metrics before comparing pooled scores.",
    "Constant and negative-Dxy items are retained exactly as in the source. A rare uninformative item still receives the rarity factor. All-zero items can consume denominator weight. These are properties exposed by the tests, not silently changed by this program.", "",
    "## Validation and provenance", "", test_result, "",
    "Input dimensions, complete 0/1 coding and Q-matrix alignment were checked. settings.R, sessionInfo.txt, package_citation.txt, input_md5.csv and run_state.rds record reproducibility information.",
    "The RDS contains audit bootstrap scores and weight draws, split IDs and RNG states. It does not change or refit the original GitHub program.",
    "Implementation generated by OpenAI ChatGPT from Anders Hellstrom's scoring specification and current source. The scoring-method definition remains the user's; the stress-test design and implementation were generated for this request.", "",
    "## Sources", "",
    "- https://tmsalab.github.io/edmdata/reference/items_ecpe.html",
    "- https://tmsalab.github.io/edmdata/reference/qmatrix_ecpe.html",
    "- https://github.com/AndersH3/Scoring/blob/main/item_analysis.r",
    "- https://github.com/harrelfe/Hmisc/blob/master/R/somers2.s", "",
    sprintf("Elapsed time: %.1f seconds.", as.numeric(difftime(Sys.time(), started, units = "secs")))
  )
  writeLines(report, file.path(out, "REPORT.md"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  writeLines(capture.output(print(citation("edmdata"))), file.path(out, "package_citation.txt"))
  input_paths <- file.path(out, c("ecpe_responses.csv", "ecpe_qmatrix.csv"))
  save_table(row(file = basename(input_paths), md5 = unname(tools::md5sum(input_paths))), "input_md5")
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg)) {
    script_path <- sub("^--file=", "", script_arg[1L])
    if (file.exists(script_path)) file.copy(script_path, file.path(out, "program_used.R"), overwrite = FALSE)
  }
  run_state <- list(version = "1.0.0", settings = o, initial_rng = initial_rng, final_rng = .Random.seed,
               train_ids = train_id, audit_ids = audit_id, folds = folds,
               bootstrap = list(share = bshare, raw = braw, operational = boper, Dxy = bdxy, audit_scores = bscore),
               bootstrap_item_order = colnames(x), bootstrap_audit_person_order = audit_id,
               null_rest = null_rest, null_total = null_total)
  # Uncompressed RDS avoids an additional compression-library dependency;
  # verify reopening before declaring the run complete or hashing outputs.
  state_path <- file.path(out, "run_state.rds")
  saveRDS(run_state, state_path, compress = FALSE)
  assert(identical(readRDS(state_path), run_state), "Saved run state failed its read-back check")
  files <- list.files(out, full.names = TRUE)
  save_table(row(file = basename(files), bytes = file.info(files)$size,
                 md5 = unname(tools::md5sum(files))), "output_manifest")
  message(sprintf("Done. Read %s", file.path(out, "REPORT.md")))
  invisible(list(output = out, full = full, calibration = ref))
}

if (sys.nframe() == 0L) {
  tryCatch(main(), error = function(e) {
    message("ERROR: ", conditionMessage(e)); quit(save = "no", status = 1L, runLast = FALSE)
  })
}
