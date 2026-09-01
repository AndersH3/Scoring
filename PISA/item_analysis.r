#!/usr/bin/env Rscript

# Two-pad weighted item scoring using Somers' Dxy and remapped positive
# concordance/discrimination.
#
# Current method
# --------------
# For item j, with n test takers and c_j correct responses:
#
#   rarity_pair_j = (c_j, n)
#
# Item discrimination is calculated from the corresponding score vector.
# The DEFAULT corresponding score is the item-rest score:
#
#   T_{i,-j} = sum_{k != j} x_{ik}
#
# Full total score remains available with --score-type total.
#
# Harrell's Somers Dxy is obtained from Hmisc::somers2():
#
#   Dxy_j = 2*C_j - 1
#
# where C_j is the empirical concordance probability in [0,1].  For weighting,
# C_j is remapped so that every value <= 0.5 maps to 0, and values from 0.5 to
# 1 map linearly to 0..1:
#
#   Cplus_j = max(0, 2*C_j - 1) = max(0, Dxy_j)
#
# For constant all-zero/all-one items, empirical C_j and Dxy_j are undefined.
# The item is retained.  A neutral pre-remap value C_weight_base=0.5 is used
# only for weighting, and therefore maps to Cplus_j=0.
#
# The second pair is:
#
#   concordance_pair_j = (n*Cplus_j, n)
#
# n*Cplus_j is intentionally NOT rounded.
#
# Padding:
#
#   pad1((c_j,n))          = (c_j + 1) / (n + 1)
#   pad2((n*Cplus_j,n))    = (n - n*Cplus_j + 1) / (n + 1)
#
# Raw formula weight:
#
#   raw_w_j = [1 - log(pad1)] * [1 - log(pad2)]
#
# Thus C_j <= 0.5 gives Cplus_j=0, pad2=1, and a concordance factor of 1:
# no positive-discrimination reward and no negative-discrimination penalty.
# Perfect concordance C_j=1 gives Cplus_j=1 and the maximum concordance reward.
#
# Operational item weight:
#
#   w_j = raw_w_j / min_k(raw_w_k)
#
# Person score:
#
#   S_i = sum_j x_ij*w_j / sum_j w_j
#
# After all raw formula weights are calculated, every item weight is divided by
# the minimum raw item weight.  Therefore the lowest operational item weight is
# always 1.  This is a pure multiplicative rescaling, so item-weight ratios and
# normalized person scores are unchanged by this final normalization.
#
# Required packages:
#   data.table, optparse, Hmisc
#
# Research-only packages:
#   boot  (only when --bootstrap-reps > 0)
#   coin  (only when --permutation-reps > 0)

fail <- function(format, ...) {
  stop(sprintf(format, ...), call. = FALSE)
}

near <- function(a, b, tolerance = 1e-10) {
  if (length(a) != length(b)) return(FALSE)
  ok_na <- is.na(a) & is.na(b)
  ok_num <- !is.na(a) & !is.na(b) &
    abs(a - b) <= tolerance * pmax(1, abs(a), abs(b))
  all(ok_na | ok_num)
}

check_packages <- function(packages) {
  installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  if (!all(installed)) {
    missing <- packages[!installed]
    fail(
      paste0(
        "Missing R package(s): %s.\n",
        "Install the required packages with:\n",
        "Rscript -e 'install.packages(c(\"data.table\",\"optparse\",\"Hmisc\",",
        "\"boot\",\"coin\"), repos=\"https://cloud.r-project.org\")'"
      ),
      paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}

near_integer <- function(x, tolerance = 1e-12) {
  length(x) == 1L &&
    !is.na(x) &&
    is.finite(x) &&
    abs(x - round(x)) <= tolerance
}

format_pair_number <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) {
    return(NA_character_)
  }

  if (near_integer(x)) {
    sprintf("%.0f", round(x))
  } else {
    formatC(x, digits = 12L, format = "fg", flag = "#")
  }
}

pair_label <- function(pair) {
  if (!is.numeric(pair) || length(pair) != 2L ||
      anyNA(pair) || any(!is.finite(pair))) {
    return(NA_character_)
  }

  sprintf(
    "(%s,%s)",
    format_pair_number(pair[[1L]]),
    format_pair_number(pair[[2L]])
  )
}

validate_rarity_pair <- function(pair, function_name = "pad1") {
  if (!is.numeric(pair) || length(pair) != 2L || anyNA(pair)) {
    fail("%s() requires a numeric pair (a,b).", function_name)
  }

  if (any(!is.finite(pair)) || any(pair != floor(pair))) {
    fail("%s() requires finite integer-valued components.", function_name)
  }

  if (pair[[1L]] < 0 || pair[[2L]] < 1 || pair[[1L]] > pair[[2L]]) {
    fail("%s() requires 0 <= a <= b and b >= 1.", function_name)
  }

  as.integer(pair)
}

validate_concordance_pair <- function(pair, function_name = "pad2") {
  if (!is.numeric(pair) || length(pair) != 2L || anyNA(pair)) {
    fail("%s() requires a numeric pair (a,b).", function_name)
  }

  if (any(!is.finite(pair))) {
    fail("%s() requires finite pair components.", function_name)
  }

  numerator <- as.double(pair[[1L]])
  denominator <- as.double(pair[[2L]])

  if (!near_integer(denominator) || denominator < 1) {
    fail("%s() requires an integer-valued denominator b >= 1.", function_name)
  }

  tolerance <- 1e-10 * max(1, denominator)

  if (numerator < -tolerance || numerator > denominator + tolerance) {
    fail("%s() requires 0 <= a <= b.", function_name)
  }

  numerator <- min(max(numerator, 0), denominator)
  c(numerator, round(denominator))
}

difficulty_padded_pair <- function(pair) {
  validated <- validate_rarity_pair(pair, "pad1")
  c(validated[[1L]] + 1L, validated[[2L]] + 1L)
}

concordance_padded_pair <- function(pair) {
  validated <- validate_concordance_pair(pair, "pad2")
  c(
    validated[[2L]] + 1 - validated[[1L]],
    validated[[2L]] + 1
  )
}

pad1 <- function(pair) {
  transformed <- difficulty_padded_pair(pair)
  transformed[[1L]] / transformed[[2L]]
}

pad2 <- function(pair) {
  transformed <- concordance_padded_pair(pair)
  transformed[[1L]] / transformed[[2L]]
}

difficulty_weight_factor <- function(pair) {
  1 - log(pad1(pair))
}

concordance_weight_factor <- function(pair) {
  1 - log(pad2(pair))
}

item_weight <- function(rarity_pair, concordance_pair) {
  difficulty_weight_factor(rarity_pair) *
    concordance_weight_factor(concordance_pair)
}

clamp_probability <- function(value, tolerance = 1e-12) {
  if (length(value) != 1L || is.na(value) || !is.finite(value)) {
    return(NA_real_)
  }

  if (value < -tolerance || value > 1 + tolerance) {
    fail("Probability %.17g is outside [0,1].", value)
  }

  min(max(as.double(value), 0), 1)
}

clamp_dxy <- function(value, tolerance = 1e-12) {
  if (length(value) != 1L || is.na(value) || !is.finite(value)) {
    return(NA_real_)
  }

  if (value < -1 - tolerance || value > 1 + tolerance) {
    fail("Somers Dxy %.17g is outside [-1,1].", value)
  }

  min(max(as.double(value), -1), 1)
}


remap_concordance <- function(concordance) {
  unit_value <- clamp_probability(concordance)

  if (is.na(unit_value)) {
    return(NA_real_)
  }

  # Piecewise-linear remapping requested for the weighting stage:
  #   C <= 0.5 -> 0
  #   0.5 < C <= 1 -> 2*C - 1
  # This is exactly max(0, Dxy), because Dxy = 2*C - 1.
  as.double(pmax(0, 2 * unit_value - 1))
}

unit_interval_pair <- function(value, n) {
  if (length(n) != 1L || is.na(n) || !is.finite(n) ||
      n != floor(n) || n < 1) {
    fail("The number of test takers n must be a positive integer.")
  }

  unit_value <- clamp_probability(value)

  if (is.na(unit_value)) {
    return(c(NA_real_, as.double(n)))
  }

  # Deliberately preserve the fractional numerator n*value.
  c(as.double(n) * unit_value, as.double(n))
}

# Backward-compatible helper name for code that supplied an ordinary C index.
concordance_pair_from_probability <- unit_interval_pair

pair_summary <- function(
  numerator,
  denominator,
  prefix,
  padded_pair_function,
  pad_function
) {
  raw_pair <- c(as.double(numerator), as.double(denominator))

  if (is.na(raw_pair[[1L]])) {
    values <- list(
      NA_character_,
      NA_real_,
      NA_real_,
      NA_character_,
      NA_real_,
      NA_real_
    )
  } else {
    transformed_pair <- padded_pair_function(raw_pair)

    values <- list(
      pair_label(raw_pair),
      as.double(raw_pair[[1L]]),
      as.double(raw_pair[[1L]] / raw_pair[[2L]]),
      pair_label(transformed_pair),
      as.double(transformed_pair[[1L]]),
      as.double(pad_function(raw_pair))
    )
  }

  stats::setNames(
    values,
    paste0(
      prefix,
      c(
        "_pair",
        "_m",
        "_proportion",
        "_padded_pair",
        "_padded_m",
        "_padded_proportion"
      )
    )
  )
}

cross_group_pair_counts <- function(item_responses, corresponding_scores) {
  y <- as.integer(item_responses)
  x <- as.double(corresponding_scores)

  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)

  if (n1 == 0L || n0 == 0L) {
    return(list(
      P = NA_real_,
      Q = NA_real_,
      T = NA_real_,
      U = NA_real_,
      cross_group_pairs = as.double(n1 * n0)
    ))
  }

  score_values <- sort(unique(x))
  neg_below <- 0
  P <- 0
  Q <- 0
  T <- 0

  for (score in score_values) {
    positive_here <- sum(y == 1L & x == score)
    negative_here <- sum(y == 0L & x == score)
    negative_above <- n0 - neg_below - negative_here

    P <- P + positive_here * neg_below
    Q <- Q + positive_here * negative_above
    T <- T + positive_here * negative_here

    neg_below <- neg_below + negative_here
  }

  list(
    P = as.double(P),
    Q = as.double(Q),
    T = as.double(T),
    U = as.double(P + 0.5 * T),
    cross_group_pairs = as.double(n0 * n1)
  )
}

somers_hmisc <- function(item_responses, corresponding_scores) {
  y <- as.integer(item_responses)
  x <- as.double(corresponding_scores)

  if (length(y) != length(x) || length(y) < 2L) {
    fail("Somers Dxy requires equally sized vectors of length at least 2.")
  }

  if (anyNA(y) || any(!(y %in% c(0L, 1L)))) {
    fail("Item responses must contain only 0 and 1.")
  }

  if (anyNA(x) || any(!is.finite(x))) {
    fail("Corresponding scores cannot contain missing or non-finite values.")
  }

  n <- length(y)
  n1 <- sum(y == 1L)
  n0 <- n - n1

  if (n1 == 0L || n0 == 0L) {
    return(list(
      C = NA_real_,
      Dxy = NA_real_,
      n = as.integer(n),
      n1 = as.integer(n1),
      n0 = as.integer(n0),
      mean_positive_rank = NA_real_,
      score_levels = as.integer(length(unique(x)))
    ))
  }

  result <- Hmisc::somers2(
    x = x,
    y = y,
    na.rm = FALSE
  )

  C <- clamp_probability(unname(result[["C"]]))
  Dxy <- clamp_dxy(unname(result[["Dxy"]]))
  ranks <- rank(x, ties.method = "average")

  list(
    C = as.double(C),
    Dxy = as.double(Dxy),
    n = as.integer(n),
    n1 = as.integer(n1),
    n0 = as.integer(n0),
    mean_positive_rank = as.double(mean(ranks[y == 1L])),
    score_levels = as.integer(length(unique(x)))
  )
}

validate_item_statistics <- function(
  item_name,
  item_responses,
  corresponding_scores,
  somers,
  pair_counts,
  extended_validation = FALSE,
  tolerance = 1e-10
) {
  n1 <- somers$n1
  n0 <- somers$n0

  if (n1 == 0L || n0 == 0L) {
    return("not_applicable_constant_item")
  }

  C <- somers$C
  Dxy <- somers$Dxy
  cross_pairs <- n0 * n1

  checks <- c(
    C_range = C >= -tolerance && C <= 1 + tolerance,
    Dxy_range = Dxy >= -1 - tolerance && Dxy <= 1 + tolerance,
    Dxy_identity = near(Dxy, 2 * C - 1, tolerance),
    pair_total = near(
      pair_counts$P + pair_counts$Q + pair_counts$T,
      cross_pairs,
      tolerance
    ),
    U_identity = near(
      pair_counts$U,
      pair_counts$P + 0.5 * pair_counts$T,
      tolerance
    ),
    C_from_U = near(
      C,
      pair_counts$U / cross_pairs,
      tolerance
    ),
    Dxy_from_pairs = near(
      Dxy,
      (pair_counts$P - pair_counts$Q) / cross_pairs,
      tolerance
    )
  )

  if (isTRUE(extended_validation)) {
    y <- as.integer(item_responses)
    x <- as.double(corresponding_scores)
    ranks <- rank(x, ties.method = "average")

    C_rank <- (
      mean(ranks[y == 1L]) - (n1 + 1) / 2
    ) / n0

    checks <- c(
      checks,
      Hmisc_vs_rank_formula = near(C, C_rank, tolerance)
    )
  }

  failed <- names(checks)[!checks]

  if (length(failed) > 0L) {
    fail(
      "Internal consistency check failed for item '%s': %s",
      item_name,
      paste(failed, collapse = ", ")
    )
  }

  if (isTRUE(extended_validation)) {
    "passed_extended_validation"
  } else {
    "passed"
  }
}

bootstrap_uncertainty <- function(
  item_responses,
  corresponding_scores,
  reps,
  conf_level
) {
  if (reps <= 0L) {
    return(list(
      bootstrap_reps = 0L,
      bootstrap_valid_reps = NA_integer_,
      bootstrap_dxy_se = NA_real_,
      bootstrap_dxy_ci_low = NA_real_,
      bootstrap_dxy_ci_high = NA_real_,
      bootstrap_C_se = NA_real_,
      bootstrap_C_ci_low = NA_real_,
      bootstrap_C_ci_high = NA_real_,
      bootstrap_remapped_concordance_se = NA_real_,
      bootstrap_remapped_concordance_ci_low = NA_real_,
      bootstrap_remapped_concordance_ci_high = NA_real_,
      bootstrap_raw_weight_se = NA_real_,
      bootstrap_raw_weight_ci_low = NA_real_,
      bootstrap_raw_weight_ci_high = NA_real_
    ))
  }

  data <- data.frame(
    y = as.integer(item_responses),
    x = as.double(corresponding_scores)
  )

  statistic <- function(data, indices) {
    yy <- data$y[indices]
    xx <- data$x[indices]
    n <- length(yy)
    n1 <- sum(yy == 1L)
    n0 <- n - n1

    if (n1 == 0L || n0 == 0L) {
      return(c(
        Dxy = NA_real_,
        C = NA_real_,
        remapped = NA_real_,
        weight = NA_real_
      ))
    }

    result <- Hmisc::somers2(xx, yy, na.rm = FALSE)
    C <- clamp_probability(as.double(unname(result[["C"]])))
    Dxy <- clamp_dxy(as.double(unname(result[["Dxy"]])))
    remapped <- remap_concordance(C)

    rarity <- c(as.double(n1), as.double(n))
    concordance_pair <- unit_interval_pair(remapped, n)
    weight <- item_weight(rarity, concordance_pair)

    c(Dxy = Dxy, C = C, remapped = remapped, weight = weight)
  }

  boot_object <- boot::boot(
    data = data,
    statistic = statistic,
    R = as.integer(reps)
  )

  estimates <- as.matrix(boot_object$t)
  colnames(estimates) <- c("Dxy", "C", "remapped", "weight")

  alpha <- 1 - conf_level

  summarize_column <- function(values) {
    values <- values[is.finite(values)]

    if (length(values) < 2L) {
      return(c(se = NA_real_, low = NA_real_, high = NA_real_))
    }

    ci <- stats::quantile(
      values,
      probs = c(alpha / 2, 1 - alpha / 2),
      names = FALSE,
      type = 7
    )

    c(
      se = stats::sd(values),
      low = ci[[1L]],
      high = ci[[2L]]
    )
  }

  dxy_summary <- summarize_column(estimates[, "Dxy"])
  C_summary <- summarize_column(estimates[, "C"])
  remapped_summary <- summarize_column(estimates[, "remapped"])
  weight_summary <- summarize_column(estimates[, "weight"])

  valid_reps <- sum(
    is.finite(estimates[, "Dxy"]) &
      is.finite(estimates[, "C"]) &
      is.finite(estimates[, "remapped"]) &
      is.finite(estimates[, "weight"])
  )

  list(
    bootstrap_reps = as.integer(reps),
    bootstrap_valid_reps = as.integer(valid_reps),
    bootstrap_dxy_se = as.double(dxy_summary[["se"]]),
    bootstrap_dxy_ci_low = as.double(dxy_summary[["low"]]),
    bootstrap_dxy_ci_high = as.double(dxy_summary[["high"]]),
    bootstrap_C_se = as.double(C_summary[["se"]]),
    bootstrap_C_ci_low = as.double(C_summary[["low"]]),
    bootstrap_C_ci_high = as.double(C_summary[["high"]]),
    bootstrap_remapped_concordance_se =
      as.double(remapped_summary[["se"]]),
    bootstrap_remapped_concordance_ci_low =
      as.double(remapped_summary[["low"]]),
    bootstrap_remapped_concordance_ci_high =
      as.double(remapped_summary[["high"]]),
    bootstrap_raw_weight_se = as.double(weight_summary[["se"]]),
    bootstrap_raw_weight_ci_low = as.double(weight_summary[["low"]]),
    bootstrap_raw_weight_ci_high = as.double(weight_summary[["high"]])
  )
}

permutation_inference <- function(
  item_responses,
  corresponding_scores,
  reps
) {
  if (reps <= 0L) {
    return(list(
      permutation_reps = 0L,
      permutation_p_value = NA_real_
    ))
  }

  y <- as.integer(item_responses)
  x <- as.double(corresponding_scores)
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)

  if (n1 == 0L || n0 == 0L) {
    return(list(
      permutation_reps = as.integer(reps),
      permutation_p_value = NA_real_
    ))
  }

  if (length(unique(x)) == 1L) {
    return(list(
      permutation_reps = as.integer(reps),
      permutation_p_value = 1
    ))
  }

  data <- data.frame(
    score = x,
    response = factor(y, levels = c(0, 1))
  )

  test <- coin::wilcox_test(
    score ~ response,
    data = data,
    distribution = coin::approximate(
      nresample = as.integer(reps)
    ),
    alternative = "two.sided"
  )

  list(
    permutation_reps = as.integer(reps),
    permutation_p_value = as.double(coin::pvalue(test))
  )
}

make_option_parser <- function() {
  options <- list(
    optparse::make_option(
      c("-o", "--output"),
      type = "character",
      default = NULL,
      help = "Alias for --items-output FILE.",
      metavar = "FILE"
    ),
    optparse::make_option(
      "--items-output",
      type = "character",
      default = NULL,
      dest = "items_output",
      help = "Write the full item audit table [default: item_weights.csv].",
      metavar = "FILE"
    ),
    optparse::make_option(
      "--persons-output",
      type = "character",
      default = NULL,
      dest = "persons_output",
      help = "Write test-taker scores [default: person_scores.csv].",
      metavar = "FILE"
    ),
    optparse::make_option(
      "--id-column",
      type = "character",
      default = "auto",
      dest = "id_column",
      help = "Person-ID column name or one-based number [default: auto].",
      metavar = "COLUMN"
    ),
    optparse::make_option(
      "--no-id-column",
      action = "store_true",
      default = FALSE,
      dest = "no_id_column",
      help = "Treat every input column as an item."
    ),
    optparse::make_option(
      "--sep",
      type = "character",
      default = "auto",
      help = "auto, comma, semicolon, tab, or one literal character.",
      metavar = "SEPARATOR"
    ),
    optparse::make_option(
      "--no-header",
      action = "store_true",
      default = FALSE,
      dest = "no_header",
      help = "The input has no header row."
    ),
    optparse::make_option(
      "--score-type",
      type = "character",
      default = "rest",
      dest = "score_type",
      help = "Corresponding score: rest or total [default: rest].",
      metavar = "TYPE"
    ),
    optparse::make_option(
      "--rest-score",
      action = "store_true",
      default = FALSE,
      dest = "rest_score",
      help = "Deprecated alias for --score-type rest."
    ),
    optparse::make_option(
      "--undefined-somers",
      type = "character",
      default = "neutral",
      dest = "undefined_somers",
      help = paste(
        "Constant-item policy: neutral or error [default: neutral].",
        paste(
          "neutral keeps the item, assigns pre-remap C=0.5 for weighting,",
          "then remaps it to 0."
        )
      ),
      metavar = "POLICY"
    ),
    optparse::make_option(
      "--min-group-warning",
      type = "integer",
      default = 5L,
      dest = "min_group_warning",
      help = "Warn when min(n0,n1) is below this threshold [default: 5].",
      metavar = "N"
    ),
    optparse::make_option(
      "--diagnostics",
      action = "store_true",
      default = FALSE,
      help = "Print diagnostic columns, including Mann-Whitney U."
    ),
    optparse::make_option(
      "--research-diagnostics",
      action = "store_true",
      default = FALSE,
      dest = "research_diagnostics",
      help = "Print P/Q/T pair counts and research-inference columns."
    ),
    optparse::make_option(
      "--bootstrap-reps",
      type = "integer",
      default = 0L,
      dest = "bootstrap_reps",
      help = "Research bootstrap replicates; 0 disables [default: 0].",
      metavar = "R"
    ),
    optparse::make_option(
      "--bootstrap-conf-level",
      type = "double",
      default = 0.95,
      dest = "bootstrap_conf_level",
      help = "Bootstrap percentile CI confidence level [default: 0.95].",
      metavar = "P"
    ),
    optparse::make_option(
      "--permutation-reps",
      type = "integer",
      default = 0L,
      dest = "permutation_reps",
      help = "Monte Carlo permutation replicates; 0 disables [default: 0].",
      metavar = "R"
    ),
    optparse::make_option(
      "--seed",
      type = "integer",
      default = 20260827L,
      help = "RNG seed for research resampling [default: 20260827].",
      metavar = "INTEGER"
    ),
    optparse::make_option(
      "--validation",
      action = "store_true",
      default = FALSE,
      help = paste(
        "Run extended validation: Hmisc results are cross-checked against",
        "independent rank/pair-count identities."
      )
    ),
    optparse::make_option(
      "--regression-expected",
      type = "character",
      default = NULL,
      dest = "regression_expected",
      help = "Compare core item results with a frozen expected CSV.",
      metavar = "FILE"
    ),
    optparse::make_option(
      "--self-test",
      action = "store_true",
      default = FALSE,
      dest = "self_test",
      help = "Run built-in mathematical/integration tests."
    )
  )

  optparse::OptionParser(
    usage = "%prog INPUT.csv [options]\n       %prog [options] INPUT.csv",
    option_list = options,
    description = paste(
      paste(
        "Two-pad item weighting using Hmisc::somers2() with",
        "C <= 0.5 mapped to 0 and C in (0.5,1] mapped linearly to (0,1]."
      ),
      "Item-rest scores are the default corresponding scores."
    )
  )
}

parse_arguments <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  parser <- make_option_parser()

  input_first <- length(arguments) > 0L &&
    !startsWith(arguments[[1L]], "-")

  leading_input <- if (input_first) arguments[[1L]] else NULL
  option_arguments <- if (input_first) arguments[-1L] else arguments

  parsed <- optparse::parse_args(
    parser,
    args = option_arguments,
    positional_arguments = c(0L, 1L),
    convert_hyphens_to_underscores = TRUE
  )

  if (!is.null(leading_input) && length(parsed$args) > 0L) {
    fail("Only one input CSV file can be supplied.")
  }

  parsed$options$input <- if (!is.null(leading_input)) {
    leading_input
  } else if (length(parsed$args) == 1L) {
    parsed$args[[1L]]
  } else {
    NULL
  }

  if (isTRUE(parsed$options$no_id_column)) {
    if (!identical(parsed$options$id_column, "auto")) {
      fail("Use either --id-column or --no-id-column, not both.")
    }
    parsed$options$id_column <- "none"
  }

  if (!is.null(parsed$options$output) &&
      !is.null(parsed$options$items_output) &&
      !identical(parsed$options$output, parsed$options$items_output)) {
    fail("Use either --output or --items-output, or give both the same file.")
  }

  parsed$options$items_output <- if (!is.null(parsed$options$items_output)) {
    parsed$options$items_output
  } else if (!is.null(parsed$options$output)) {
    parsed$options$output
  } else {
    "item_weights.csv"
  }

  if (is.null(parsed$options$persons_output)) {
    parsed$options$persons_output <- "person_scores.csv"
  }

  if (identical(parsed$options$items_output, parsed$options$persons_output)) {
    fail("Item and person outputs must be different files.")
  }

  score_type <- tolower(trimws(parsed$options$score_type))
  if (!(score_type %in% c("rest", "total"))) {
    fail("--score-type must be 'rest' or 'total'.")
  }

  if (isTRUE(parsed$options$rest_score) && identical(score_type, "total")) {
    fail("--rest-score conflicts with --score-type total.")
  }

  if (isTRUE(parsed$options$rest_score)) {
    score_type <- "rest"
  }

  parsed$options$score_type <- score_type

  policy <- tolower(trimws(parsed$options$undefined_somers))
  if (!(policy %in% c("neutral", "error"))) {
    fail("--undefined-somers must be 'neutral' or 'error'.")
  }
  parsed$options$undefined_somers <- policy

  if (is.na(parsed$options$min_group_warning) ||
      parsed$options$min_group_warning < 0L) {
    fail("--min-group-warning must be a nonnegative integer.")
  }

  for (field in c("bootstrap_reps", "permutation_reps")) {
    value <- parsed$options[[field]]
    if (is.na(value) || value < 0L) {
      fail("--%s must be a nonnegative integer.", gsub("_", "-", field))
    }
    if (value == 1L) {
      fail("--%s must be 0 or at least 2.", gsub("_", "-", field))
    }
  }

  if (!is.finite(parsed$options$bootstrap_conf_level) ||
      parsed$options$bootstrap_conf_level <= 0 ||
      parsed$options$bootstrap_conf_level >= 1) {
    fail("--bootstrap-conf-level must be strictly between 0 and 1.")
  }

  if (is.na(parsed$options$seed)) {
    fail("--seed must be an integer.")
  }

  parsed$options
}

resolve_separator <- function(value) {
  separator <- switch(
    value,
    auto = "auto",
    comma = ",",
    semicolon = ";",
    tab = "\t",
    value
  )

  if (!identical(separator, "auto") &&
      nchar(separator, type = "chars") != 1L) {
    fail("Separator must be auto, comma, semicolon, tab, or one character.")
  }

  separator
}

resolve_id_column <- function(data, specification = "auto") {
  if (identical(specification, "none")) {
    return(NULL)
  }

  if (!identical(specification, "auto")) {
    if (grepl("^[1-9][0-9]*$", specification)) {
      column <- as.integer(specification)
      if (is.na(column) || column > ncol(data)) {
        fail("Person-ID column %s is outside the CSV range.", specification)
      }
      return(column)
    }

    column <- match(specification, names(data))
    if (is.na(column)) {
      fail("Person-ID column '%s' was not found.", specification)
    }
    return(column)
  }

  first_name <- tolower(trimws(names(data)[[1L]]))
  id_named <- grepl(
    paste0(
      "^(id|person|person_?id|person_?number|participant|",
      "participant_?id|respondent|respondent_?id)$"
    ),
    first_name
  )

  first_is_binary <- all(trimws(data[[1L]]) %in% c("0", "1"))

  if (id_named || !first_is_binary) 1L else NULL
}

read_response_matrix <- function(
  path,
  id_column = "auto",
  separator = "auto",
  header = TRUE
) {
  if (!file.exists(path)) {
    fail("Input CSV does not exist: %s", path)
  }

  input <- data.table::fread(
    file = path,
    header = header,
    sep = resolve_separator(separator),
    colClasses = "character",
    check.names = FALSE,
    strip.white = TRUE,
    fill = FALSE,
    blank.lines.skip = TRUE,
    na.strings = c("", "NA"),
    showProgress = FALSE
  )

  if (nrow(input) < 1L || ncol(input) < 1L) {
    fail("The CSV must contain at least one person and one item.")
  }

  if (anyDuplicated(names(input))) {
    fail("CSV column names must be unique.")
  }

  id_index <- resolve_id_column(input, id_column)

  if (is.null(id_index)) {
    person_ids <- as.character(seq_len(nrow(input)))
    item_data <- input
  } else {
    person_ids <- trimws(input[[id_index]])

    if (anyNA(person_ids) ||
        any(!nzchar(person_ids)) ||
        anyDuplicated(person_ids)) {
      fail("Person identifiers must be nonempty, nonmissing, and unique.")
    }

    item_data <- data.table::copy(input)
    data.table::set(item_data, j = id_index, value = NULL)
  }

  if (ncol(item_data) < 1L || any(!nzchar(names(item_data)))) {
    fail("The CSV must contain at least one uniquely named item column.")
  }

  values <- as.matrix(item_data)
  invalid <- which(
    is.na(values) | !(values %in% c("0", "1")),
    arr.ind = TRUE
  )

  if (nrow(invalid) > 0L) {
    row <- invalid[1L, 1L]
    column <- invalid[1L, 2L]
    value <- values[row, column]

    if (is.na(value) || !nzchar(value)) {
      value <- "<missing>"
    }

    fail(
      "Item '%s' has invalid response '%s' for person '%s'; only 0/1 allowed.",
      colnames(values)[[column]],
      value,
      person_ids[[row]]
    )
  }

  matrix(
    as.integer(values),
    nrow = nrow(values),
    ncol = ncol(values),
    dimnames = list(person_ids, colnames(values))
  )
}

analyse_item <- function(
  item_name,
  item_responses,
  corresponding_scores,
  undefined_somers = "neutral",
  min_group_warning = 5L,
  research_diagnostics = FALSE,
  bootstrap_reps = 0L,
  bootstrap_conf_level = 0.95,
  permutation_reps = 0L,
  validation = FALSE
) {
  people <- length(item_responses)

  if (people == 0L || length(corresponding_scores) != people) {
    fail("Item '%s' requires equally sized nonempty vectors.", item_name)
  }

  responses <- as.integer(item_responses)

  if (anyNA(responses) || any(!(responses %in% c(0L, 1L)))) {
    fail("Item '%s' contains values other than 0 and 1.", item_name)
  }

  if (!is.numeric(corresponding_scores) ||
      anyNA(corresponding_scores) ||
      any(!is.finite(corresponding_scores))) {
    fail("Item '%s' has invalid corresponding scores.", item_name)
  }

  somers <- somers_hmisc(responses, corresponding_scores)
  n1 <- somers$n1
  n0 <- somers$n0
  pair_counts <- cross_group_pair_counts(responses, corresponding_scores)

  status <- if (n1 == 0L || n0 == 0L) {
    "constant_item"
  } else if (somers$score_levels == 1L) {
    "constant_corresponding_score"
  } else {
    "ok"
  }

  if (is.na(somers$Dxy) && identical(undefined_somers, "error")) {
    fail("Item '%s' has undefined Somers Dxy (%s).", item_name, status)
  }

  validation_status <- validate_item_statistics(
    item_name = item_name,
    item_responses = responses,
    corresponding_scores = corresponding_scores,
    somers = somers,
    pair_counts = pair_counts,
    extended_validation = validation
  )

  rarity <- c(as.double(n1), as.double(people))

  # Somers Dxy/C are empirically undefined for constant binary items.  The item
  # is nevertheless retained.  For weighting only, assign a neutral PRE-REMAP
  # concordance C=0.5 (Dxy=0).  Then apply the requested piecewise-linear map:
  #
  #   remapped_concordance = max(0, 2*C - 1) = max(0, Dxy)
  #
  # Thus all C <= 0.5 receive remapped_concordance=0 and a concordance factor
  # of exactly 1.  Empirical somers_dxy/concordance_probability remain NA for
  # constant items so the output does not misrepresent an estimated association.
  concordance_imputed <- is.na(somers$C)
  weight_base_concordance_probability <- if (concordance_imputed) {
    0.5
  } else {
    somers$C
  }
  weight_base_somers_dxy <- 2 * weight_base_concordance_probability - 1
  remapped_concordance <- remap_concordance(
    weight_base_concordance_probability
  )

  # Algebraic validation of the remapping itself.
  expected_remapped <- max(0, weight_base_somers_dxy)
  if (!near(remapped_concordance, expected_remapped, tolerance = 1e-12)) {
    fail(
      "Remapped-concordance identity failed for item '%s'.",
      item_name
    )
  }

  concordance_pair <- unit_interval_pair(
    remapped_concordance,
    people
  )

  rarity_factor <- difficulty_weight_factor(rarity)
  concordance_factor <- concordance_weight_factor(concordance_pair)
  calculated_weight <- rarity_factor * concordance_factor

  weight_status <- if (!is.finite(calculated_weight)) {
    "undefined_weight"
  } else if (calculated_weight <= 0) {
    "nonpositive_weight"
  } else if (concordance_imputed) {
    "ok_constant_neutral_pre_remap_to_zero"
  } else {
    "ok"
  }

  small_group <- n1 > 0L &&
    n0 > 0L &&
    min(n0, n1) < min_group_warning

  small_group_message <- if (small_group) {
    sprintf(
      "min(n0,n1)=%d is below warning threshold %d",
      min(n0, n1),
      min_group_warning
    )
  } else {
    ""
  }

  bootstrap <- if (n1 > 0L && n0 > 0L) {
    bootstrap_uncertainty(
      item_responses = responses,
      corresponding_scores = corresponding_scores,
      reps = bootstrap_reps,
      conf_level = bootstrap_conf_level
    )
  } else {
    bootstrap_uncertainty(
      item_responses = c(0L, 1L),
      corresponding_scores = c(0, 1),
      reps = 0L,
      conf_level = bootstrap_conf_level
    )
  }

  permutation <- if (n1 > 0L && n0 > 0L) {
    permutation_inference(
      item_responses = responses,
      corresponding_scores = corresponding_scores,
      reps = permutation_reps
    )
  } else {
    permutation_inference(
      item_responses = c(0L, 1L),
      corresponding_scores = c(0, 1),
      reps = 0L
    )
  }

  research_values <- if (isTRUE(research_diagnostics)) {
    list(
      concordant_pairs = as.double(pair_counts$P),
      discordant_pairs = as.double(pair_counts$Q),
      tied_cross_pairs = as.double(pair_counts$T)
    )
  } else {
    list(
      concordant_pairs = NA_real_,
      discordant_pairs = NA_real_,
      tied_cross_pairs = NA_real_
    )
  }

  values <- c(
    list(
      item = as.character(item_name),
      number_of_test_takers = as.integer(people),
      n1 = as.integer(n1),
      n0 = as.integer(n0),
      cross_group_pairs = as.double(n0 * n1),
      corresponding_score_levels = as.integer(somers$score_levels),
      mean_positive_rank = as.double(somers$mean_positive_rank),
      somers_dxy = as.double(somers$Dxy),
      concordance_probability = as.double(somers$C),
      weight_base_somers_dxy = as.double(weight_base_somers_dxy),
      weight_base_concordance_probability =
        as.double(weight_base_concordance_probability),
      remapped_concordance = as.double(remapped_concordance),
      positive_somers_dxy = as.double(remapped_concordance),
      # Backward-compatible aliases for the pre-remap weighting base.
      weight_somers_dxy = as.double(weight_base_somers_dxy),
      weight_concordance_probability =
        as.double(weight_base_concordance_probability),
      concordance_imputed_for_weight = isTRUE(concordance_imputed),
      concordance_weight_source = if (concordance_imputed) {
        "neutral_C_0.5_then_remap_to_0_constant_item"
      } else if (somers$C <= 0.5) {
        "empirical_Hmisc_somers2_then_remap_to_0"
      } else {
        "empirical_Hmisc_somers2_linear_remap_above_0.5"
      },
      remapping_validation_status = "passed",
      mann_whitney_u = as.double(pair_counts$U)
    ),
    research_values,
    pair_summary(
      n1,
      people,
      "rarity",
      difficulty_padded_pair,
      pad1
    ),
    pair_summary(
      concordance_pair[[1L]],
      people,
      "concordance",
      concordance_padded_pair,
      pad2
    ),
    list(
      rarity_weight_factor = as.double(rarity_factor),
      concordance_weight_factor = as.double(concordance_factor),
      raw_item_weight = as.double(calculated_weight),
      status = status,
      weight_status = weight_status,
      small_group_warning = isTRUE(small_group),
      small_group_message = small_group_message,
      validation_status = validation_status
    ),
    bootstrap,
    permutation
  )

  data.table::as.data.table(values)
}

normalize_item_weights <- function(results) {
  if (!data.table::is.data.table(results)) {
    results <- data.table::as.data.table(results)
  }

  raw_weights <- as.double(results$raw_item_weight)

  if (length(raw_weights) < 1L ||
      anyNA(raw_weights) ||
      any(!is.finite(raw_weights))) {
    fail("All raw item weights must be finite before normalization.")
  }

  if (any(raw_weights <= 0)) {
    fail(
      paste(
        "Minimum-to-one normalization requires every raw item weight",
        "to be strictly positive."
      )
    )
  }

  divisor <- min(raw_weights)
  normalized <- raw_weights / divisor

  # A multiplicative common-factor transformation must preserve all ratios.
  if (!near(min(normalized), 1, tolerance = 1e-12)) {
    fail("Internal weight normalization check failed: minimum is not 1.")
  }

  data.table::set(
    results,
    j = "item_weight",
    value = as.double(normalized)
  )
  data.table::set(
    results,
    j = "weight_scale_divisor",
    value = as.double(divisor)
  )
  data.table::set(
    results,
    j = "minimum_item_weight",
    value = as.double(min(normalized))
  )
  data.table::set(
    results,
    j = "maximum_item_weight",
    value = as.double(max(normalized))
  )

  results
}

analyse_items <- function(
  response_matrix,
  score_type = "rest",
  undefined_somers = "neutral",
  min_group_warning = 5L,
  research_diagnostics = FALSE,
  bootstrap_reps = 0L,
  bootstrap_conf_level = 0.95,
  permutation_reps = 0L,
  validation = FALSE
) {
  if (!is.matrix(response_matrix)) {
    response_matrix <- as.matrix(response_matrix)
  }

  if (!is.numeric(response_matrix) ||
      anyNA(response_matrix) ||
      any(!(response_matrix %in% c(0L, 1L)))) {
    fail("Response matrix must be numeric, complete, and binary 0/1.")
  }

  if (nrow(response_matrix) < 1L || ncol(response_matrix) < 1L) {
    fail("Response matrix must contain at least one person and one item.")
  }

  if (is.null(colnames(response_matrix))) {
    colnames(response_matrix) <- paste0("item_", seq_len(ncol(response_matrix)))
  }

  total_scores <- rowSums(response_matrix)

  per_item <- lapply(
    seq_len(ncol(response_matrix)),
    function(index) {
      responses <- response_matrix[, index]

      corresponding_scores <- if (identical(score_type, "rest")) {
        total_scores - responses
      } else {
        total_scores
      }

      analyse_item(
        item_name = colnames(response_matrix)[[index]],
        item_responses = responses,
        corresponding_scores = corresponding_scores,
        undefined_somers = undefined_somers,
        min_group_warning = min_group_warning,
        research_diagnostics = research_diagnostics,
        bootstrap_reps = bootstrap_reps,
        bootstrap_conf_level = bootstrap_conf_level,
        permutation_reps = permutation_reps,
        validation = validation
      )
    }
  )

  results <- data.table::rbindlist(per_item, use.names = TRUE, fill = TRUE)

  data.table::set(
    results,
    j = "score_type",
    value = if (identical(score_type, "rest")) "item_rest" else "total"
  )

  results <- normalize_item_weights(results)

  summary <- summarise_weights(results$item_weight)

  for (field in names(summary)) {
    data.table::set(results, j = field, value = summary[[field]])
  }

  results
}

summarise_weights <- function(weights) {
  if (!is.numeric(weights) || length(weights) < 1L) {
    fail("At least one numeric item weight is required.")
  }

  maximum <- sum(weights)
  minimum_achievable <- sum(pmin(weights, 0))
  maximum_achievable <- sum(pmax(weights, 0))
  flags <- character()

  if (any(!is.finite(weights))) {
    flags <- c(flags, "non_finite_item_weights")
  }

  if (any(weights < 0, na.rm = TRUE)) {
    flags <- c(flags, "negative_item_weights")
  }

  if (!is.finite(maximum)) {
    flags <- c(flags, "non_finite_max_total_score")
  } else if (maximum == 0) {
    flags <- c(flags, "zero_max_total_score")
  } else if (maximum < 0) {
    flags <- c(flags, "negative_max_total_score")
  }

  list(
    max_total_score = as.double(maximum),
    minimum_achievable_total_score = as.double(minimum_achievable),
    maximum_achievable_total_score = as.double(maximum_achievable),
    maximum_is_true_maximum = all(is.finite(weights)) && all(weights >= 0),
    zero_weight_count = sum(weights == 0, na.rm = TRUE),
    negative_weight_count = sum(weights < 0, na.rm = TRUE),
    non_finite_weight_count = sum(!is.finite(weights)),
    scoring_model_status = if (length(flags) == 0L) {
      "ok"
    } else {
      paste(flags, collapse = ";")
    }
  )
}

score_test_takers <- function(response_matrix, item_results) {
  if (!is.matrix(response_matrix)) {
    response_matrix <- as.matrix(response_matrix)
  }

  if (ncol(response_matrix) != nrow(item_results)) {
    fail("Number of item weights must match response columns.")
  }

  if (!is.null(colnames(response_matrix)) &&
      !identical(colnames(response_matrix), as.character(item_results$item))) {
    fail("Response columns and item weights must have identical names/order.")
  }

  weights <- as.double(item_results$item_weight)
  formula_weights <- as.double(item_results$raw_item_weight)
  summary <- summarise_weights(weights)
  formula_summary <- summarise_weights(formula_weights)

  contributions <- sweep(response_matrix, 2L, weights, FUN = "*")
  contributions[response_matrix == 0L] <- 0

  raw_scores <- as.vector(response_matrix %*% weights)
  formula_raw_scores <- as.vector(response_matrix %*% formula_weights)
  maximum <- summary$max_total_score
  formula_maximum <- formula_summary$max_total_score

  normalized <- if (is.finite(maximum) && maximum != 0) {
    raw_scores / maximum
  } else {
    rep(NA_real_, nrow(response_matrix))
  }

  formula_normalized <- if (is.finite(formula_maximum) &&
                            formula_maximum != 0) {
    formula_raw_scores / formula_maximum
  } else {
    rep(NA_real_, nrow(response_matrix))
  }

  if (!near(normalized, formula_normalized, tolerance = 1e-10)) {
    fail(
      paste(
        "Multiplicative item-weight normalization unexpectedly changed",
        "normalized person scores."
      )
    )
  }

  normalization_status <- if (!is.finite(maximum)) {
    rep("undefined_non_finite_maximum", nrow(response_matrix))
  } else if (maximum == 0) {
    rep("undefined_zero_maximum", nrow(response_matrix))
  } else {
    ifelse(
      !is.finite(normalized),
      "undefined_non_finite_person_score",
      ifelse(
        normalized < 0 | normalized > 1,
        "outside_unit_interval",
        "ok"
      )
    )
  }

  people <- rownames(response_matrix)
  if (is.null(people)) {
    people <- as.character(seq_len(nrow(response_matrix)))
  }

  correct_counts <- rowSums(response_matrix)

  correct_items <- vapply(
    seq_len(nrow(response_matrix)),
    function(index) {
      paste(
        item_results$item[response_matrix[index, ] == 1L],
        collapse = "|"
      )
    },
    character(1)
  )

  scores <- data.table::data.table(
    person = people,
    number_correct = as.integer(correct_counts),
    number_of_items = as.integer(ncol(response_matrix)),
    correct_item_pair = sprintf(
      "(%d,%d)",
      correct_counts,
      ncol(response_matrix)
    ),
    unweighted_proportion_correct = correct_counts / ncol(response_matrix),
    correct_items = correct_items,
    formula_weighted_raw_score = as.double(formula_raw_scores),
    formula_max_total_score = as.double(formula_maximum),
    weighted_raw_score = as.double(raw_scores),
    max_total_score = as.double(maximum),
    normalized_weighted_score = as.double(normalized),
    normalized_weighted_percent = as.double(100 * normalized),
    minimum_achievable_total_score =
      summary$minimum_achievable_total_score,
    maximum_achievable_total_score =
      summary$maximum_achievable_total_score,
    maximum_is_true_maximum =
      summary$maximum_is_true_maximum,
    normalization_status = normalization_status,
    scoring_model_status = summary$scoring_model_status,
    score_type = item_results$score_type[[1L]]
  )

  response_columns <- data.table::as.data.table(response_matrix)
  contribution_columns <- data.table::as.data.table(contributions)

  data.table::setnames(
    response_columns,
    make.names(
      paste0("response_item_", item_results$item),
      unique = TRUE
    )
  )

  data.table::setnames(
    contribution_columns,
    make.names(
      paste0("contribution_item_", item_results$item),
      unique = TRUE
    )
  )

  cbind(scores, response_columns, contribution_columns)
}

analyse_test <- function(response_matrix, options) {
  research_enabled <- isTRUE(options$research_diagnostics) ||
    options$bootstrap_reps > 0L ||
    options$permutation_reps > 0L

  items <- analyse_items(
    response_matrix = response_matrix,
    score_type = options$score_type,
    undefined_somers = options$undefined_somers,
    min_group_warning = options$min_group_warning,
    research_diagnostics = research_enabled,
    bootstrap_reps = options$bootstrap_reps,
    bootstrap_conf_level = options$bootstrap_conf_level,
    permutation_reps = options$permutation_reps,
    validation = options$validation
  )

  list(
    items = items,
    persons = score_test_takers(response_matrix, items),
    summary = summarise_weights(items$item_weight)
  )
}

compare_regression_results <- function(actual, expected_path) {
  if (is.null(expected_path)) {
    return(invisible(TRUE))
  }

  if (!file.exists(expected_path)) {
    fail("Regression expected-results CSV does not exist: %s", expected_path)
  }

  expected <- data.table::fread(
    expected_path,
    na.strings = c("", "NA"),
    showProgress = FALSE
  )

  required <- c(
    "item",
    "n1",
    "n0",
    "cross_group_pairs",
    "somers_dxy",
    "concordance_probability",
    "weight_somers_dxy",
    "weight_concordance_probability",
    "remapped_concordance",
    "mann_whitney_u",
    "raw_item_weight",
    "item_weight",
    "status",
    "weight_status",
    "score_type"
  )

  missing <- setdiff(required, names(expected))
  if (length(missing) > 0L) {
    fail(
      "Regression expected CSV is missing: %s",
      paste(missing, collapse = ", ")
    )
  }

  if (!identical(as.character(actual$item), as.character(expected$item))) {
    fail("Regression check failed: item names/order differ.")
  }

  numeric_columns <- c(
    "n1",
    "n0",
    "cross_group_pairs",
    "somers_dxy",
    "concordance_probability",
    "weight_somers_dxy",
    "weight_concordance_probability",
    "remapped_concordance",
    "mann_whitney_u",
    "raw_item_weight",
    "item_weight"
  )

  for (column in numeric_columns) {
    if (!near(
      as.double(actual[[column]]),
      as.double(expected[[column]]),
      tolerance = 1e-9
    )) {
      fail("Regression check failed for numeric column '%s'.", column)
    }
  }

  character_columns <- c("status", "weight_status", "score_type")

  for (column in character_columns) {
    if (!identical(
      as.character(actual[[column]]),
      as.character(expected[[column]])
    )) {
      fail("Regression check failed for text column '%s'.", column)
    }
  }

  message("Frozen regression-output check passed.")
  invisible(TRUE)
}

run_self_tests <- function(validation = TRUE) {
  # Padding.
  stopifnot(identical(pad1(c(2L, 4L)), 3 / 5))
  stopifnot(isTRUE(all.equal(pad2(c(0, 4)), 1)))
  stopifnot(isTRUE(all.equal(pad2(c(2, 4)), 3 / 5)))
  stopifnot(isTRUE(all.equal(pad2(c(4, 4)), 1 / 5)))

  # Requested piecewise-linear concordance remapping.
  stopifnot(isTRUE(all.equal(remap_concordance(0.00), 0.00)))
  stopifnot(isTRUE(all.equal(remap_concordance(0.25), 0.00)))
  stopifnot(isTRUE(all.equal(remap_concordance(0.50), 0.00)))
  stopifnot(isTRUE(all.equal(remap_concordance(0.60), 0.20)))
  stopifnot(isTRUE(all.equal(remap_concordance(0.75), 0.50)))
  stopifnot(isTRUE(all.equal(remap_concordance(0.90), 0.80)))
  stopifnot(isTRUE(all.equal(remap_concordance(1.00), 1.00)))

  # Perfect concordance.
  perfect_y <- c(0, 0, 1, 1)
  perfect_x <- 1:4
  perfect <- somers_hmisc(perfect_y, perfect_x)
  perfect_pairs <- cross_group_pair_counts(perfect_y, perfect_x)

  stopifnot(isTRUE(all.equal(perfect$C, 1)))
  stopifnot(isTRUE(all.equal(perfect$Dxy, 1)))
  stopifnot(isTRUE(all.equal(perfect_pairs$U, 4)))

  validate_item_statistics(
    "perfect",
    perfect_y,
    perfect_x,
    perfect,
    perfect_pairs,
    extended_validation = validation
  )

  # Complete reversal.
  reverse <- somers_hmisc(perfect_y, 4:1)
  stopifnot(isTRUE(all.equal(reverse$C, 0)))
  stopifnot(isTRUE(all.equal(reverse$Dxy, -1)))

  # All scores tied -> C=.5, Dxy=0.
  tied <- somers_hmisc(
    c(0, 1, 0, 1),
    c(10, 10, 10, 10)
  )
  stopifnot(isTRUE(all.equal(tied$C, 0.5)))
  stopifnot(isTRUE(all.equal(tied$Dxy, 0)))

  # Mixed/tied example.
  mixed_y <- c(0, 1, 0, 1)
  mixed_x <- c(10, 12, 12, 14)
  mixed <- somers_hmisc(mixed_y, mixed_x)
  mixed_pairs <- cross_group_pair_counts(mixed_y, mixed_x)

  stopifnot(isTRUE(all.equal(mixed$C, 0.875)))
  stopifnot(isTRUE(all.equal(mixed$Dxy, 0.75)))
  stopifnot(isTRUE(all.equal(mixed_pairs$U, 3.5)))

  validate_item_statistics(
    "mixed",
    mixed_y,
    mixed_x,
    mixed,
    mixed_pairs,
    extended_validation = validation
  )

  # Constant empirical Somers remains undefined, but constant items are retained.
  # Neutral pre-remap C=0.5 is used for weighting, and maps to 0.
  constant <- somers_hmisc(c(0, 0, 0, 0), 1:4)
  stopifnot(is.na(constant$C))
  stopifnot(is.na(constant$Dxy))

  constant_item <- analyse_item(
    "constant_zero",
    c(0, 0, 0, 0),
    c(0, 1, 2, 3),
    undefined_somers = "neutral"
  )
  stopifnot(is.na(constant_item$somers_dxy))
  stopifnot(is.na(constant_item$concordance_probability))
  stopifnot(isTRUE(all.equal(
    constant_item$weight_base_concordance_probability, 0.5
  )))
  stopifnot(isTRUE(all.equal(constant_item$remapped_concordance, 0)))
  stopifnot(isTRUE(all.equal(constant_item$concordance_weight_factor, 1)))
  stopifnot(isTRUE(constant_item$concordance_imputed_for_weight))
  stopifnot(constant_item$raw_item_weight > 0)

  # Complete reversal and all-tied scores both receive no concordance reward.
  reverse_item <- analyse_item(
    "reverse",
    perfect_y,
    4:1,
    undefined_somers = "neutral"
  )
  stopifnot(isTRUE(all.equal(reverse_item$remapped_concordance, 0)))
  stopifnot(isTRUE(all.equal(reverse_item$concordance_weight_factor, 1)))

  tied_item <- analyse_item(
    "tied",
    c(0, 1, 0, 1),
    c(10, 10, 10, 10),
    undefined_somers = "neutral"
  )
  stopifnot(isTRUE(all.equal(tied_item$remapped_concordance, 0)))
  stopifnot(isTRUE(all.equal(tied_item$concordance_weight_factor, 1)))

  # Minimum-to-one normalization preserves ratios.
  test_weights <- data.table::data.table(
    raw_item_weight = c(2, 3, 5)
  )
  test_weights <- normalize_item_weights(test_weights)
  stopifnot(isTRUE(all.equal(test_weights$item_weight, c(1, 1.5, 2.5))))
  stopifnot(isTRUE(all.equal(test_weights$weight_scale_divisor, c(2, 2, 2))))

  message(paste(
    "All built-in Hmisc/Somers/remapped-concordance/normalization",
    "self-tests passed."
  ))
  invisible(TRUE)
}

display_results <- function(items, persons, options) {
  base_columns <- c(
    "item",
    "n1",
    "n0",
    "rarity_pair",
    "somers_dxy",
    "concordance_probability",
    "remapped_concordance",
    "concordance_pair",
    "raw_item_weight",
    "item_weight",
    "status",
    "score_type"
  )

  diagnostic_columns <- c(
    base_columns,
    "cross_group_pairs",
    "corresponding_score_levels",
    "mean_positive_rank",
    "mann_whitney_u",
    "rarity_padded_pair",
    "rarity_padded_proportion",
    "concordance_padded_pair",
    "concordance_padded_proportion",
    "rarity_weight_factor",
    "concordance_weight_factor",
    "weight_base_concordance_probability",
    "concordance_imputed_for_weight",
    "concordance_weight_source",
    "remapping_validation_status",
    "small_group_warning",
    "validation_status"
  )

  research_columns <- c(
    diagnostic_columns,
    "concordant_pairs",
    "discordant_pairs",
    "tied_cross_pairs",
    "bootstrap_reps",
    "bootstrap_valid_reps",
    "bootstrap_dxy_se",
    "bootstrap_dxy_ci_low",
    "bootstrap_dxy_ci_high",
    "bootstrap_C_se",
    "bootstrap_C_ci_low",
    "bootstrap_C_ci_high",
    "bootstrap_remapped_concordance_se",
    "bootstrap_remapped_concordance_ci_low",
    "bootstrap_remapped_concordance_ci_high",
    "bootstrap_raw_weight_se",
    "bootstrap_raw_weight_ci_low",
    "bootstrap_raw_weight_ci_high",
    "permutation_reps",
    "permutation_p_value"
  )

  item_columns <- if (isTRUE(options$research_diagnostics) ||
                      options$bootstrap_reps > 0L ||
                      options$permutation_reps > 0L) {
    research_columns
  } else if (isTRUE(options$diagnostics)) {
    diagnostic_columns
  } else {
    base_columns
  }

  person_columns <- if (isTRUE(options$diagnostics)) {
    names(persons)
  } else {
    c(
      "person",
      "number_correct",
      "weighted_raw_score",
      "max_total_score",
      "normalized_weighted_score",
      "normalization_status"
    )
  }

  message("Item results:")
  print(items[, ..item_columns], row.names = FALSE)

  message("Test-taker weighted scores:")
  print(persons[, ..person_columns], row.names = FALSE)
}

main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  check_packages(c("data.table", "optparse", "Hmisc"))
  options <- parse_arguments(arguments)

  if (options$bootstrap_reps > 0L) {
    check_packages("boot")
  }

  if (options$permutation_reps > 0L) {
    check_packages("coin")
  }

  set.seed(options$seed)

  if (isTRUE(options$self_test)) {
    return(invisible(run_self_tests(validation = TRUE)))
  }

  if (is.null(options$input)) {
    fail("Please supply a person-by-item CSV; use --help for usage.")
  }

  responses <- read_response_matrix(
    path = options$input,
    id_column = options$id_column,
    separator = options$sep,
    header = !isTRUE(options$no_header)
  )

  results <- analyse_test(responses, options)
  items <- results$items
  persons <- results$persons
  summary <- results$summary

  compare_regression_results(
    actual = items,
    expected_path = options$regression_expected
  )

  display_results(items, persons, options)

  data.table::fwrite(
    items,
    file = options$items_output,
    na = "NA"
  )

  data.table::fwrite(
    persons,
    file = options$persons_output,
    na = "NA"
  )

  message(sprintf("Item audit table written to %s", options$items_output))
  message(sprintf("Person scores written to %s", options$persons_output))
  message(sprintf(
    "Maximum total score, sum(item weights): %.15g",
    summary$max_total_score
  ))

  small_items <- items$item[items$small_group_warning]
  if (length(small_items) > 0L) {
    message(sprintf(
      "Warning: small response group for item(s): %s",
      paste(small_items, collapse = ", ")
    ))
  }

  constant_items <- items$item[items$status == "constant_item"]

  if (length(constant_items) > 0L) {
    message(sprintf(
      paste(
        "Constant item(s) retained; empirical Somers Dxy is undefined.",
        paste(
          "Neutral pre-remap C=0.5 is used for weighting and maps to",
          "remapped concordance 0: %s"
        )
      ),
      paste(constant_items, collapse = ", ")
    ))
  }

  message(sprintf(
    "Raw-weight scale divisor (minimum raw item weight): %.15g",
    items$weight_scale_divisor[[1L]]
  ))
  message(sprintf(
    "Operational item-weight range: %.15g to %.15g",
    min(items$item_weight),
    max(items$item_weight)
  ))

  invisible(results)
}

if (sys.nframe() == 0L) {
  tryCatch(
    main(),
    error = function(error) {
      message(sprintf("Error: %s", conditionMessage(error)))
      quit(save = "no", status = 1L, runLast = FALSE)
    }
  )
}
