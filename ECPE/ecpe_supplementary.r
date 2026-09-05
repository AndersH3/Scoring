#!/usr/bin/env Rscript
# Supplementary analyses for the ECPE report, 2026-09-01.
# Run from the source-bundle root: Rscript code/ecpe_supplementary.R
# These are NEW report analyses, separate from the user's original stress run.
source("code/ecpe_scoring_stress_test.R")
dir.create("supplementary", showWarnings = FALSE)
dat <- read.csv("data/ecpe_responses.csv", check.names = FALSE)
x <- validate_x(as.matrix(dat[, -1L])); state <- readRDS("data/run_state.rds")
cal <- x[state$train_ids, ]; audit <- x[state$audit_ids, ]
f <- fit_weights(cal); s <- score_with(cal, f); a <- score_with(audit, f)
q <- rowMeans(cal); qa <- rowMeans(audit)
z <- q * (1 - q)
kappa_unconstrained <- sum(z * (s - q)) / sum(z^2)
kappa <- max(-1, min(1, kappa_unconstrained))
linear <- lm(s ~ q)
quadratic <- lm(s ~ q + I(q^2))
predictions <- list(
  raw_proportion = qa,
  linear = as.vector(cbind(1, qa) %*% coef(linear)),
  ordinary_quadratic = as.vector(cbind(1, qa, qa^2) %*% coef(quadratic)),
  endpoint_monotone_quadratic = qa + kappa * qa * (1 - qa)
)
regression <- do.call(rbind, lapply(names(predictions), function(name) {
  y <- predictions[[name]]
  data.frame(model = name, audit_rmse = sqrt(mean((y - a)^2)),
             audit_mae = mean(abs(y - a)), audit_R2 = 1 - sum((y - a)^2) / sum((a - mean(a))^2))
}))
write.csv(regression, "supplementary/regression_comparison.csv", row.names = FALSE)
write.csv(data.frame(kappa = kappa, kappa_unconstrained = kappa_unconstrained,
                     linear_intercept = coef(linear)[1], linear_slope = coef(linear)[2],
                     quadratic_intercept = coef(quadratic)[1],
                     quadratic_linear = coef(quadratic)[2], quadratic_squared = coef(quadratic)[3]),
          "supplementary/regression_coefficients.csv", row.names = FALSE)
write.csv(data.frame(person = state$audit_ids, q = qa, frozen_weighted_score = a,
                     constrained_prediction = predictions$endpoint_monotone_quadratic),
          "supplementary/regression_audit_predictions.csv", row.names = FALSE)
# Known-null appended items, generated independently of the ECPE response rows.
# All existing weights are re-estimated on the augmented calibration matrix.
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20260902L)
reps <- 1000L; rates <- c(.001, .01, .05, .5)
out <- vector("list", length(rates) * reps); idx <- 0L
for (p in rates) {
  message(sprintf("Null appended items: p=%g, %d replicates", p, reps))
  for (b in seq_len(reps)) {
    y <- cbind(cal, synthetic_null = rbinom(nrow(cal), 1, p))
    w <- fit_weights(y); v <- w[nrow(w), ]; idx <- idx + 1L
    out[[idx]] <- data.frame(probability = p, replicate = b, correct = v$correct,
      Dxy = v$Dxy, raw_weight = v$raw_weight, share = v$share,
      constant = v$constant, share_exceeds_all_original = v$share > max(head(w$share, -1L)))
  }
}
draws <- do.call(rbind, out)
write.csv(draws, "supplementary/null_item_replicates.csv", row.names = FALSE)
summaries <- do.call(rbind, lapply(rates, function(p) {
  d <- draws[draws$probability == p, ]
  data.frame(probability = p, replicates = reps, expected_correct = p * nrow(cal),
    share_p025 = quantile(d$share, .025), share_median = median(d$share), share_p975 = quantile(d$share, .975),
    max_weight_rate = mean(d$share_exceeds_all_original), constant_rate = mean(d$constant),
    positive_Dxy_rate_among_nonconstant = mean(d$Dxy[!d$constant] > 0),
    Dxy_mean_nonconstant = mean(d$Dxy, na.rm = TRUE))
}))
write.csv(summaries, "supplementary/null_item_summary.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "supplementary/sessionInfo.txt")
message("Supplementary analyses completed.")
