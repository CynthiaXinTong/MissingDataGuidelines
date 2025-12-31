library(mice)
library(lavaan)

gcmodel <- '
  i =~ 1*V1 + 1*V2 + 1*V3 + 1*V4
  s =~ 0*V1 + 1*V2 + 2*V3 + 3*V4
'

.param_est_se <- function(fit) {
  pe <- parameterEstimates(fit)
  
  idx <- c(13:15, 20:21)
  if (max(idx) <= nrow(pe)) {
    sub <- pe[idx, c("est", "se")]
  } else {
    want <- rbind(
      c("i","=~","V1"),
      c("i","=~","V2"),
      c("i","=~","V3"),
      c("s","=~","V2"),
      c("s","=~","V4")
    )
    key  <- paste(pe$lhs, pe$op, pe$rhs)
    pick <- match(apply(want, 1, paste, collapse = " "), key)
    if (any(is.na(pick))) stop("can not find five parameters")
    sub <- pe[pick, c("est", "se")]
  }
  
  as.matrix(sub)
}

.rubin_pool <- function(q_mat, se_mat, alpha = 0.05) {
  m <- nrow(q_mat)
  
  Ubar <- colMeans(se_mat^2, na.rm = TRUE)                     # within-imputation var
  Qbar <- colMeans(q_mat,    na.rm = TRUE)                     # pooled estimate
  B    <- apply(q_mat, 2, function(x) stats::var(x, na.rm = TRUE))  # between var
  
  if (m <= 1 || any(is.na(B))) B[is.na(B)] <- 0
  
  Tvar   <- Ubar + (1 + 1/m) * B
  eps    <- .Machine$double.eps
  df     <- ifelse(B < eps, Inf, (m - 1) * (1 + Ubar / ((1 + 1/m) * B))^2)
  se_tot <- sqrt(Tvar)
  tval   <- Qbar / se_tot
  
  pval <- ifelse(is.infinite(df),
                 2 * (1 - pnorm(abs(tval))),
                 2 * (1 - pt(abs(tval), df)))
  
  q_crit <- ifelse(is.infinite(df),
                   qnorm(1 - alpha/2),
                   qt(1 - alpha/2, df))
  
  ci_low <- Qbar - q_crit * se_tot
  ci_up  <- Qbar + q_crit * se_tot
  
  cbind(est = Qbar, se = se_tot, t = tval, df = df, p = pval,
        ci.low = ci_low, ci.up = ci_up)
}

# data：dataagree
df <- dataagree

M    <- 5L
SEED <- 20250428L

set.seed(SEED)
imp <- mice(df, m = M, method = "cart", seed = SEED, print = FALSE)

q_list <- list()
se_list <- list()
ok <- 0L

for (j in seq_len(M)) {
  dj <- try(mice::complete(imp, j), silent = TRUE)
  if (inherits(dj, "try-error")) next
  
  fitj <- try(growth(gcmodel, data = dj), silent = TRUE)
  if (inherits(fitj, "try-error")) next
  
  ok <- ok + 1L
  estse <- .param_est_se(fitj)
  q_list[[ok]]  <- estse[, 1]
  se_list[[ok]] <- estse[, 2]
}

if (ok == 0L) stop("All imputations failed to fit; please check the data and model settings.")

q_mat  <- do.call(rbind, q_list)   # ok × 5
se_mat <- do.call(rbind, se_list)  # ok × 5

tab <- .rubin_pool(q_mat, se_mat)  # 5 × 7

param_names <- paste0("par", 1:5)
res_out <- cbind(param = param_names, round(tab, 6))
res_out



