library(lavaan)
library(dplyr)

## ---------- Utility functions ----------
calibrate_alpha <- function(eta, target){
  if (target <= 0) return(-Inf)
  if (target >= 1) return( Inf)
  f <- function(a) mean(plogis(a + eta)) - target
  uniroot(f, c(-50, 50))$root
}

# MAR missingness with split/grouped anchor variables
make_mar_split_anchors <- function(df6,
                                   target_rates,
                                   anchor_map,
                                   gamma_vec = c(0,0,1.5,0,0,1.5)) {
  stopifnot(ncol(df6) == 6, all(names(df6) == paste0("x",1:6)))
  Z1 <- scale(df6[["x1"]])[,1]
  Z4 <- scale(df6[["x4"]])[,1]
  out <- df6
  n <- nrow(df6)
  
  for (j in 1:6) {
    r <- target_rates[j]
    if (r <= 0) next
    anc <- anchor_map[j]
    eta <- switch(anc,
                  "x1" = gamma_vec[j] * Z1,
                  "x4" = gamma_vec[j] * Z4,
                  stop(sprintf("Column x%d requires an anchor variable, but anchor_map is NA", j)))
    a <- calibrate_alpha(eta, r)
    p <- plogis(a + eta)
    miss <- runif(n) < p
    out[miss, j] <- NA
  }
  out
}

# ---------- MNAR: Auxiliary-threshold approach ----------
make_mnar_aux_all <- function(df6,
                              target_rates = rep(0.3, 6),
                              corAb = 0.8,
                              noise_sd = 1,
                              direction = "high") {
  stopifnot(ncol(df6) == 6, all(names(df6) == paste0("x",1:6)))
  n <- nrow(df6)
  out <- as.data.frame(df6)
  
  to_len6 <- function(x, name) {
    if (length(x) == 1) rep(x, 6)
    else if (length(x) == 6) x
    else stop(sprintf("%s must be a scalar or a length-6 vector", name))
  }
  target_rates <- to_len6(target_rates, "target_rates")
  corAb        <- to_len6(corAb,        "corAb")
  noise_sd     <- to_len6(noise_sd,     "noise_sd")
  direction    <- to_len6(direction,    "direction")
  direction[!direction %in% c("high","low")] <- "high"
  
  for (j in 1:6) {
    r <- target_rates[j]
    if (r <= 0) next
    z <- as.numeric(scale(df6[[j]]))
    a <- (corAb[j] / sqrt(1 - corAb[j]^2)) * noise_sd[j]
    aux <- a * z + rnorm(n, 0, noise_sd[j])
    n_miss <- round(r * n)
    if (n_miss > 0) {
      ord <- if (direction[j] == "high") order(aux, decreasing = TRUE)
      else order(aux, decreasing = FALSE)
      idx <- ord[seq_len(n_miss)]
      out[idx, j] <- NA
    }
  }
  out
}

## ---------- Simulation model ----------
pp <- 1
model_population <- '
  F1 =~ 0.6*x1 + 0.6*x2 + 0.6*x3
  F2 =~ 0.6*x4 + 0.6*x5 + 0.6*x6
  F2 ~ 0.2*F1

  x1~0*1; x2~0*1; x3~0*1; x4~0*1; x5~0*1; x6~0*1
  F1~0*1; F2~0*1

  x1~~0.64*x1; x2~~0.64*x2; x3~~0.64*x3
  x4~~0.64*x4; x5~~0.64*x5; x6~~0.64*x6
  F1~~1*F1; F2~~1*F2
'

# Five missing-rate levels
miss.rate0 <- 0
miss.rate1 <- 0.05
miss.rate2 <- 0.10
miss.rate3 <- 0.30
miss.rate4 <- 0.45

set.seed(202511104)

# Anchor mapping for MAR
anchor_map <- c(NA, NA, "x1", NA, NA, "x4")

# MAR missing-rate template
mk_tr_mar <- function(r) c(0, 0, r, 0, 0, r)
# MNAR missing-rate template (missingness applied to the same target columns here)
mk_tr_mnar <- function(r) c(0, 0, r, 0, 0, r)

for (N in c(100, 200, 500, 1000)) {
  for (rep1 in 1:300) {
    # ---- Complete data ----
    dat <- simulateData(model_population, sample.nobs = N, return.type = 'data.frame')
    df6 <- dplyr::select(dat, x1:x6)
    id <- rep(rep1, N)
    
    # Write complete data
    write.table(cbind(id, df6),
                file = paste('y.','pp',pp,'.N',N,'.txt',sep=''),
                append = TRUE, row.names = FALSE, col.names = FALSE)
    
    # ---- Generate MAR ----
    tr1 <- mk_tr_mar(miss.rate1)
    tr2 <- mk_tr_mar(miss.rate2)
    tr3 <- mk_tr_mar(miss.rate3)
    tr4 <- mk_tr_mar(miss.rate4)
    
    yMAR1 <- make_mar_split_anchors(df6, target_rates = tr1, anchor_map = anchor_map, gamma_vec = gamma_vec)
    yMAR2 <- make_mar_split_anchors(df6, target_rates = tr2, anchor_map = anchor_map, gamma_vec = gamma_vec)
    yMAR3 <- make_mar_split_anchors(df6, target_rates = tr3, anchor_map = anchor_map, gamma_vec = gamma_vec)
    yMAR4 <- make_mar_split_anchors(df6, target_rates = tr4, anchor_map = anchor_map, gamma_vec = gamma_vec)
    
    # ---- Generate MNAR (Aux threshold approach) ----
    tr1_mnar <- mk_tr_mnar(miss.rate1)
    tr2_mnar <- mk_tr_mnar(miss.rate2)
    tr3_mnar <- mk_tr_mnar(miss.rate3)
    tr4_mnar <- mk_tr_mnar(miss.rate4)
    
    yMNAR1 <- make_mnar_aux_all(df6, target_rates = tr1_mnar, corAb = 0.8, noise_sd = 1, direction = "high")
    yMNAR2 <- make_mnar_aux_all(df6, target_rates = tr2_mnar, corAb = 0.8, noise_sd = 1, direction = "high")
    yMNAR3 <- make_mnar_aux_all(df6, target_rates = tr3_mnar, corAb = 0.8, noise_sd = 1, direction = "high")
    yMNAR4 <- make_mnar_aux_all(df6, target_rates = tr4_mnar, corAb = 0.8, noise_sd = 1, direction = "high")
    
    # ---- Write MAR outputs ----
    write.table(cbind(id, yMAR1), file = paste('MAR1.','pp',pp,'.N',N,'.txt',sep=''),
                append = TRUE, row.names = FALSE, col.names = FALSE)
    write.table(cbind(id, yMAR2), file = paste('MAR2.','pp',pp,'.N',N,'.txt',sep=''),
                append = TRUE, row.names = FALSE, col.names = FALSE)
    write.table(cbind(id, yMAR3), file = paste('MAR3.','pp',pp,'.N',N,'.txt',sep=''),
                append = TRUE, row.names = FALSE, col.names = FALSE)
    write.table(cbind(id, yMAR4), file = paste('MAR4.','pp',pp,'.N',N,'.txt',sep=''),
                append = TRUE, row.names = FALSE, col.names = FALSE)
    
    # ---- Write MNAR outputs ----
    write.table(cbind(id, yMNAR1), file = paste('MNAR1.','pp',pp,'.N',N,'.txt',sep=''),
                append = TRUE, row.names = FALSE, col.names = FALSE)
    write.table(cbind(id, yMNAR2), file = paste('MNAR2.','pp',pp,'.N',N,'.txt',sep=''),
                append = TRUE, row.names = FALSE, col.names = FALSE)
    write.table(cbind(id, yMNAR3), file = paste('MNAR3.','pp',pp,'.N',N,'.txt',sep=''),
                append = TRUE, row.names = FALSE, col.names = FALSE)
    write.table(cbind(id, yMNAR4), file = paste('MNAR4.','pp',pp,'.N',N,'.txt',sep=''),
                append = TRUE, row.names = FALSE, col.names = FALSE)
  }
}
