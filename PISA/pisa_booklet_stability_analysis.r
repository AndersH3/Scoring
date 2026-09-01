#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edmdata)
  library(Hmisc)
})

# Current project formula:
#   item-rest Somers Dxy
#   remapped_concordance = max(0, Dxy)
#   raw_weight = [1 + log((n+1)/(n1+1))] *
#                [1 + log((n+1)/(n*(1-remapped)+1))]
#   operational weight = raw_weight / min(raw_weight within form)
#
# PISA edmdata coding:
#   Score 7 (N/A / not administered) -> NA
#   Score 8 (Not Attempted)          -> 0
# Structural NA patterns are used here to reconstruct the principal booklet groups.

x <- as.data.frame(edmdata::items_pisa12_us_math)
items <- names(x)

pattern_key <- apply(is.na(x), 1, function(z) paste0(as.integer(z), collapse=""))
tab <- sort(table(pattern_key), decreasing=TRUE)

# PISA 2012 has 13 principal standard booklet forms. In this U.S. matrix,
# the 13 largest non-empty structural patterns account for almost all students.
pattern_keys <- names(tab)[seq_len(min(13, length(tab)))]

patterns <- lapply(seq_along(pattern_keys), function(k) {
  key <- pattern_keys[k]
  idx <- which(pattern_key == key)
  administered <- names(x)[!is.na(x[idx[1], ])]
  list(
    pattern = sprintf("PATTERN%02d", k),
    key = key,
    idx = idx,
    n_students = length(idx),
    items = administered
  )
})

major_pattern_summary <- do.call(rbind, lapply(patterns, function(p) {
  data.frame(
    pattern = p$pattern,
    n_students = p$n_students,
    n_items_administered = length(p$items),
    stringsAsFactors = FALSE
  )
}))
write.csv(major_pattern_summary, "pisa_major_pattern_summary.csv", row.names=FALSE)

somers_dxy <- function(rest, y) {
  if (length(unique(y)) < 2L || length(unique(rest)) < 1L) return(NA_real_)
  as.numeric(Hmisc::somers2(rest, y)["Dxy"])
}

raw_weight <- function(n, n1, dxy) {
  d_for_weight <- if (is.na(dxy)) 0 else max(0, dxy)
  f1 <- 1 + log((n + 1) / (n1 + 1))
  f2 <- 1 + log((n + 1) / (n * (1 - d_for_weight) + 1))
  f1 * f2
}

# For each item, determine all major patterns that contain it and the
# intersection of administered items across those patterns. Under the PISA
# balanced incomplete-block design, this intersection is a useful empirical
# way to recover a common within-cluster anchor set for stability comparisons.
patterns_containing_item <- lapply(items, function(item) {
  which(vapply(patterns, function(p) item %in% p$items, logical(1)))
})
names(patterns_containing_item) <- items

shared_item_sets <- lapply(items, function(item) {
  pp <- patterns_containing_item[[item]]
  if (length(pp) < 2L) return(character())
  Reduce(intersect, lapply(patterns[pp], `[[`, "items"))
})
names(shared_item_sets) <- items

all_rows <- list()
row_id <- 1L

for (p in patterns) {
  idx <- p$idx
  form <- x[idx, p$items, drop=FALSE]
  stopifnot(!anyNA(form))
  n_form <- nrow(form)
  form_total <- rowSums(form)

  # First pass: form-specific item-rest calibrations and raw weights.
  form_rows <- vector("list", length(p$items))
  for (jj in seq_along(p$items)) {
    item <- p$items[jj]
    y <- as.numeric(form[[item]])
    n1 <- sum(y == 1)
    n0 <- n_form - n1
    rest_form <- form_total - y
    dxy_form <- somers_dxy(rest_form, y)
    remap_form <- if (is.na(dxy_form)) 0 else max(0, dxy_form)
    rw_form <- raw_weight(n_form, n1, dxy_form)

    # Common-anchor rest score: use only items administered in every major
    # pattern containing this target item, excluding the target itself.
    shared <- setdiff(shared_item_sets[[item]], item)
    shared <- intersect(shared, names(form))
    if (length(shared) > 0L) {
      rest_common <- rowSums(form[, shared, drop=FALSE])
      dxy_common <- somers_dxy(rest_common, y)
      rw_common <- raw_weight(n_form, n1, dxy_common)
    } else {
      dxy_common <- NA_real_
      rw_common <- NA_real_
    }

    form_rows[[jj]] <- data.frame(
      pattern = p$pattern,
      n_students = n_form,
      n_items_form = ncol(form),
      item = item,
      n1 = n1,
      n0 = n0,
      p_correct = n1 / n_form,
      dxy_form_rest = dxy_form,
      remapped_form = remap_form,
      raw_weight_form = rw_form,
      shared_rest_item_count = length(shared),
      dxy_common_rest = dxy_common,
      raw_weight_common_rest = rw_common,
      stringsAsFactors = FALSE
    )
  }

  dfp <- do.call(rbind, form_rows)
  divisor <- min(dfp$raw_weight_form, na.rm=TRUE)
  dfp$pattern_divisor <- divisor
  dfp$min1_form_weight <- dfp$raw_weight_form / divisor
  all_rows[[row_id]] <- dfp
  row_id <- row_id + 1L
}

cal <- do.call(rbind, all_rows)
write.csv(cal, "pisa_all_pattern_item_calibrations.csv", row.names=FALSE)

# Repeated-item stability summary
split_item <- split(cal, cal$item)
stability <- do.call(rbind, lapply(split_item, function(d) {
  if (nrow(d) < 2L) return(NULL)

  cv <- function(z) {
    z <- z[is.finite(z)]
    if (length(z) < 2L || mean(z) == 0) return(NA_real_)
    sd(z) / mean(z)
  }

  data.frame(
    item = d$item[1],
    n_patterns = nrow(d),
    patterns = paste(d$pattern, collapse=";"),
    p_correct_mean = mean(d$p_correct),
    p_correct_sd = sd(d$p_correct),
    p_correct_min = min(d$p_correct),
    p_correct_max = max(d$p_correct),
    dxy_form_mean = mean(d$dxy_form_rest, na.rm=TRUE),
    dxy_form_sd = sd(d$dxy_form_rest, na.rm=TRUE),
    dxy_form_min = min(d$dxy_form_rest, na.rm=TRUE),
    dxy_form_max = max(d$dxy_form_rest, na.rm=TRUE),
    raw_weight_form_mean = mean(d$raw_weight_form, na.rm=TRUE),
    raw_weight_form_sd = sd(d$raw_weight_form, na.rm=TRUE),
    raw_weight_form_cv = cv(d$raw_weight_form),
    raw_weight_form_min = min(d$raw_weight_form, na.rm=TRUE),
    raw_weight_form_max = max(d$raw_weight_form, na.rm=TRUE),
    min1_weight_mean = mean(d$min1_form_weight, na.rm=TRUE),
    min1_weight_sd = sd(d$min1_form_weight, na.rm=TRUE),
    common_rest_item_count_min = min(d$shared_rest_item_count),
    dxy_common_mean = mean(d$dxy_common_rest, na.rm=TRUE),
    dxy_common_sd = sd(d$dxy_common_rest, na.rm=TRUE),
    raw_weight_common_mean = mean(d$raw_weight_common_rest, na.rm=TRUE),
    raw_weight_common_sd = sd(d$raw_weight_common_rest, na.rm=TRUE),
    raw_weight_common_cv = cv(d$raw_weight_common_rest),
    stringsAsFactors = FALSE
  )
}))
write.csv(stability, "pisa_repeated_item_stability.csv", row.names=FALSE)

# Pairwise comparisons for each repeated item.
pair_rows <- list()
kk <- 1L
for (item in names(split_item)) {
  d <- split_item[[item]]
  if (nrow(d) < 2L) next
  cmb <- combn(seq_len(nrow(d)), 2)
  for (cc in seq_len(ncol(cmb))) {
    a <- d[cmb[1,cc],]
    b <- d[cmb[2,cc],]
    pair_rows[[kk]] <- data.frame(
      item = item,
      pattern_a = a$pattern,
      pattern_b = b$pattern,
      n_a = a$n_students,
      n_b = b$n_students,
      p_correct_a = a$p_correct,
      p_correct_b = b$p_correct,
      delta_p_correct = b$p_correct - a$p_correct,
      dxy_form_a = a$dxy_form_rest,
      dxy_form_b = b$dxy_form_rest,
      delta_dxy_form = b$dxy_form_rest - a$dxy_form_rest,
      raw_weight_form_a = a$raw_weight_form,
      raw_weight_form_b = b$raw_weight_form,
      ratio_raw_weight_form = b$raw_weight_form / a$raw_weight_form,
      dxy_common_a = a$dxy_common_rest,
      dxy_common_b = b$dxy_common_rest,
      delta_dxy_common = b$dxy_common_rest - a$dxy_common_rest,
      raw_weight_common_a = a$raw_weight_common_rest,
      raw_weight_common_b = b$raw_weight_common_rest,
      ratio_raw_weight_common = b$raw_weight_common_rest / a$raw_weight_common_rest,
      stringsAsFactors = FALSE
    )
    kk <- kk + 1L
  }
}
pairwise <- if (length(pair_rows)) do.call(rbind, pair_rows) else data.frame()
write.csv(pairwise, "pisa_repeated_item_pairwise.csv", row.names=FALSE)

pm <- cal[cal$item == "pm995q02", ]
write.csv(pm, "pisa_pm995q02_across_patterns.csv", row.names=FALSE)

# Simple plots
if (nrow(pm) > 0) {
  png("pisa_pm995q02_across_patterns.png", width=1400, height=900, res=160)
  op <- par(mar=c(8,5,3,1))
  plot(seq_len(nrow(pm)), pm$raw_weight_form,
       xaxt="n", xlab="", ylab="Raw item weight",
       main="pm995q02 across PISA booklet-pattern groups", pch=19)
  axis(1, at=seq_len(nrow(pm)), labels=pm$pattern, las=2)
  par(op)
  dev.off()
}

if (nrow(stability) > 0) {
  s <- stability[order(stability$raw_weight_form_cv, decreasing=TRUE), ]
  png("pisa_repeated_item_weight_stability.png", width=1600, height=1000, res=160)
  op <- par(mar=c(9,5,3,1))
  barplot(s$raw_weight_form_cv,
          names.arg=s$item, las=2,
          ylab="CV of raw item weight across patterns",
          main="Repeated PISA item-weight stability across booklet-pattern groups")
  par(op)
  dev.off()
}

cat("\nMajor booklet-pattern groups:\n")
print(major_pattern_summary)

cat("\nRepeated-item stability (sorted by raw-weight CV):\n")
print(stability[order(stability$raw_weight_form_cv, decreasing=TRUE), ])

cat("\npm995q02 across major patterns:\n")
print(pm)

cat("\nWrote:\n",
    "  pisa_major_pattern_summary.csv\n",
    "  pisa_all_pattern_item_calibrations.csv\n",
    "  pisa_repeated_item_stability.csv\n",
    "  pisa_repeated_item_pairwise.csv\n",
    "  pisa_pm995q02_across_patterns.csv\n",
    "  pisa_pm995q02_across_patterns.png\n",
    "  pisa_repeated_item_weight_stability.png\n", sep="")
