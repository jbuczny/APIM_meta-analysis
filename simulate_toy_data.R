###############################################################################
##  simulate_toy_data.R
##  Generates two toy dyadic datasets (N = 130 and N = 174 dyads) for testing
##  apim_meta_analysis.R.
##
##  Design
##  - Distinguishable dyads (p1 = Partner 1, p2 = Partner 2), wide format:
##    one row per dyad.
##  - Three constructs per partner: X (predictor), M (mediator), Y (outcome),
##    each measured with 4 Likert items (1-5) -> 24 items per dyad.
##  - True-score (latent) APIMeM with KNOWN standardized paths, so the
##    reliability-corrected meta-analysis can be checked against the truth.
##  - Study 2 uses less reliable scales (lower loadings) than Study 1, so the
##    uncorrected analysis is attenuated more in Study 2.
##  - ~1.5 % MCAR item-level missingness to exercise the missing-data code.
##
##  Population (latent, standardized) paths
##    actor a  = .40 (P1), .35 (P2)      partner a  = .15 (X1->M2), .10 (X2->M1)
##    actor b  = .30 (P1), .35 (P2)      partner b  = .12 (M1->Y2), .10 (M2->Y1)
##    actor c' = .15 (P1), .10 (P2)      partner c' = .08 (X1->Y2), .05 (X2->Y1)
###############################################################################

set.seed(20260911)

vars <- c("X1", "X2", "M1", "M2", "Y1", "Y2")

## ---- 1. Structural (true-score) model ----------------------------------------
## B[i, j] = effect of variable j on variable i
B <- matrix(0, 6, 6, dimnames = list(vars, vars))
B["M1", "X1"] <- .40; B["M1", "X2"] <- .10    # actor a (P1), partner a (X2->M1)
B["M2", "X2"] <- .35; B["M2", "X1"] <- .15    # actor a (P2), partner a (X1->M2)
B["Y1", "M1"] <- .30; B["Y1", "M2"] <- .10    # actor b (P1), partner b (M2->Y1)
B["Y2", "M2"] <- .35; B["Y2", "M1"] <- .12    # actor b (P2), partner b (M1->Y2)
B["Y1", "X1"] <- .15; B["Y1", "X2"] <- .05    # actor c' (P1), partner c' (X2->Y1)
B["Y2", "X2"] <- .10; B["Y2", "X1"] <- .08    # actor c' (P2), partner c' (X1->Y2)

r_x   <- .30   # correlation between partners' X (non-independence)
r_eM  <- .20   # residual correlation between partners' M
r_eY  <- .25   # residual correlation between partners' Y

## Residual variances chosen so that every latent variable has variance 1
implied_sigma <- function(B, r_x, r_eM, r_eY) {
  Psi <- diag(6); dimnames(Psi) <- dimnames(B)
  Psi["X1", "X2"] <- Psi["X2", "X1"] <- r_x
  IB  <- solve(diag(6) - B)
  for (it in 1:200) {                       # fixed-point iteration (recursive model)
    S <- IB %*% Psi %*% t(IB)
    for (v in c("M1", "M2", "Y1", "Y2")) Psi[v, v] <- Psi[v, v] + (1 - S[v, v])
    Psi["M1", "M2"] <- Psi["M2", "M1"] <- r_eM * sqrt(Psi["M1", "M1"] * Psi["M2", "M2"])
    Psi["Y1", "Y2"] <- Psi["Y2", "Y1"] <- r_eY * sqrt(Psi["Y1", "Y1"] * Psi["Y2", "Y2"])
  }
  S <- IB %*% Psi %*% t(IB)
  stopifnot(max(abs(diag(S) - 1)) < 1e-10)
  S
}
Sigma_T <- implied_sigma(B, r_x, r_eM, r_eY)

## ---- 2. Measurement model + Likert discretization ---------------------------
make_study <- function(n_dyads, loadings, thresholds, p_miss = .015) {
  Tscores <- MASS::mvrnorm(n_dyads, mu = rep(0, 6), Sigma = Sigma_T)
  colnames(Tscores) <- vars
  out <- data.frame(dyad_id = sprintf("D%03d", seq_len(n_dyads)))
  for (v in vars) {
    construct <- tolower(substr(v, 1, 1))
    member    <- paste0("p", substr(v, 2, 2))
    for (i in seq_along(loadings)) {
      lam  <- loadings[i]
      ystar <- lam * Tscores[, v] + rnorm(n_dyads, 0, sqrt(1 - lam^2))
      item  <- as.integer(cut(ystar, breaks = c(-Inf, thresholds, Inf)))  # 1..5
      item[runif(n_dyads) < p_miss] <- NA
      out[[sprintf("%s_%s_%d", construct, member, i)]] <- item
    }
  }
  out
}

## Slightly asymmetric thresholds -> mildly skewed items (non-normality)
thr <- c(-1.6, -0.7, 0.3, 1.2)

study1 <- make_study(130, loadings = c(.80, .75, .78, .72), thresholds = thr)
study2 <- make_study(174, loadings = c(.62, .58, .66, .55), thresholds = thr)

write.csv(study1, "toy_study1.csv", row.names = FALSE, na = "")
write.csv(study2, "toy_study2.csv", row.names = FALSE, na = "")

## Population values saved for checking recovery
pop <- data.frame(
  path  = c("a1", "a2", "a12", "a21", "b1", "b2", "b12", "b21",
            "c1", "c2", "c12", "c21"),
  true_std = c(B["M1","X1"], B["M2","X2"], B["M2","X1"], B["M1","X2"],
               B["Y1","M1"], B["Y2","M2"], B["Y2","M1"], B["Y1","M2"],
               B["Y1","X1"], B["Y2","X2"], B["Y2","X1"], B["Y1","X2"]))
write.csv(pop, "toy_population_values.csv", row.names = FALSE)

cat("Wrote toy_study1.csv (", nrow(study1), "dyads ) and toy_study2.csv (",
    nrow(study2), "dyads )\n")
