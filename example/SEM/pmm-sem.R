library(mice)
library(lavaan)

SEM_ec <- '
  F1 =~ e11 + e12 + e13 + e14 + e15 + e16 + e17 + e18
  F2 =~ c81 + c82 + c83 + c84 + c85 + c86 + c87
  F2 =~ F1
'

df <- dfec

M    <- 5L
SEED <- 20250428L
set.seed(SEED)

rubin_pool <- function(q_vec, se_vec, alpha = 0.05) {
  m <- length(q_vec)
  Ubar <- mean(se_vec^2, na.rm = TRUE)
  Qbar <- mean(q_vec, na.rm = TRUE)
  B    <- stats::var(q_vec, na.rm = TRUE)
  if (is.na(B)) B <- 0
  
  Tvar   <- Ubar + (1 + 1/m) * B
  se_tot <- sqrt(Tvar)
  
  eps <- .Machine$double.eps
  df  <- ifelse(B < eps, Inf,
                (m - 1) * (1 + Ubar / ((1 + 1/m) * B))^2)
  
  tval <- Qbar / se_tot
  pval <- if (is.infinite(df)) 2 * (1 - pnorm(abs(tval))) else 2 * (1 - pt(abs(tval), df))
  qcrit <- if (is.infinite(df)) qnorm(1 - alpha/2) else qt(1 - alpha/2, df)
  ciL <- Qbar - qcrit * se_tot
  ciU <- Qbar + qcrit * se_tot
  
  c(est = Qbar, se = se_tot, t = tval, df = df, p = pval, ci.low = ciL, ci.up = ciU)
}

extract_std_path <- function(fit, lhs = "F2", rhs = "F1") {
  ss <- standardizedSolution(fit, se = TRUE)
  ss <- ss[ss$op == "=~" & ss$lhs == lhs & ss$rhs == rhs, , drop = FALSE]
  if (nrow(ss) < 1) stop("can not find: ", lhs, " =~ ", rhs)
  
  std_cols <- c("std.all", "est.std", "est.std.all")
  col_hit  <- std_cols[std_cols %in% names(ss)]
  if (length(col_hit) == 0) stop("can not find（std.all / est.std）")
  est_std <- as.numeric(ss[[col_hit[1]]])
  
  if (!("se" %in% names(ss))) stop("In standardizedSolution, can not find se")
  se_std <- as.numeric(ss$se)
  
  c(est = est_std, se = se_std)
}

imp <- mice(df, m = M, method = "pmm", seed = SEED, print = FALSE)

q <- numeric(0)
s <- numeric(0)
ok <- 0L

for (j in seq_len(M)) {
  dj <- try(mice::complete(imp, j), silent = TRUE)
  if (inherits(dj, "try-error")) next
  
  fit <- try(sem(SEM_ec, data = dj, std.lv = TRUE), silent = TRUE)
  if (inherits(fit, "try-error")) next
  
  es <- try(extract_std_path(fit, lhs = "F2", rhs = "F1"), silent = TRUE)
  if (inherits(es, "try-error")) next
  
  ok <- ok + 1L
  q[ok] <- es["est"]
  s[ok] <- es["se"]
}

if (ok == 0L) stop("All imputations failed during model fitting or parameter extraction. 
Please check the data and the model specification.")
if (ok < M) warning(sprintf("Only %d out of %d imputations were successful; pooling 
was performed using the successful imputations only.", ok, M))

pool <- rubin_pool(q, s)

out <- data.frame(
  path   = "F2 =~ F1",
  m_used = ok,
  est    = pool["est"],
  se     = pool["se"],
  t      = pool["t"],
  df     = pool["df"],
  p      = pool["p"],
  ci.low = pool["ci.low"],
  ci.up  = pool["ci.up"]
)

out
