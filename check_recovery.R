## Structural check: large-N recovery of true latent paths by the corrected model
src <- readLines("apim_meta_analysis.R")
stop_at <- grep("^## 5. RUN", src) - 2
eval(parse(text = src[1:stop_at]))
set.seed(1)
source_lines <- readLines("simulate_toy_data.R")
eval(parse(text = source_lines[grep("^vars <-", source_lines):(grep("^## ---- 2", source_lines)-1)]))
big <- function(n, lam, discrete) {
  T <- MASS::mvrnorm(n, rep(0,6), Sigma_T); colnames(T) <- vars
  out <- data.frame(id = 1:n)
  for (v in vars) for (i in seq_along(lam)) {
    y <- lam[i]*T[,v] + rnorm(n, 0, sqrt(1-lam[i]^2))
    if (discrete) y <- as.integer(cut(y, c(-Inf, -1.6,-.7,.3,1.2, Inf)))
    out[[sprintf("%s_p%s_%d", tolower(substr(v,1,1)), substr(v,2,2), i)]] <- y
  }
  out
}
truth <- c(a1=.40,a2=.35,a12=.15,a21=.10,b1=.30,b2=.35,b12=.12,b21=.10,c1=.15,c2=.10,c12=.08,c21=.05)
for (disc in c(FALSE, TRUE)) {
  d  <- big(20000, c(.62,.58,.66,.55), disc)
  ps <- prep_study(d, "big", CFG)
  r  <- analyse_study(ps)
  cat("\nDiscretized:", disc, " omega:", round(ps$rel,3), "\n")
  print(round(rbind(truth = truth[PATHS$label], uncorrected = r$uncorrected$est, corrected = r$corrected$est), 3))
}
