library(dplyr)

# ------------------------------------------------------------
# (1) data generator
# ------------------------------------------------------------
generate_data <- function(N, rep) {
  set.seed(rep)
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  y  <- -0.3 * x1 + 0.14 * x2 + rnorm(N)  # the coefficient settings were obtained from the empirical example by running MissingDataGuidelines/example/regression/missing_data_empirical_example.R
  data.frame(y, x1, x2)
}

# ------------------------------------------------------------
# (2) missing data mechanisms
# ------------------------------------------------------------
MAR_generator_Y_X1 <- function(data, mr) {
  idx <- which(data$x1 > qnorm(1 - mr, mean(data$x1), sd(data$x1)))
  data[idx, "y"] <- NA
  data["y"]
}

MNAR_generator_Y_Self <- function(data, mr) {
  idx <- which(data$y > qnorm(1 - mr, mean(data$y), sd(data$y)))
  data[idx, "y"] <- NA
  data["y"]
}

MNAR_generator_X <- function(data, mr) {
  for (j in 2:ncol(data)) {
    idx <- which(data[[j]] > qnorm(1 - mr, mean(data[[j]]), sd(data[[j]])))
    data[idx, j] <- NA
  }
  data[, 2:3]
}

# ------------------------------------------------------------
# (3) simulation parameters
# ------------------------------------------------------------
N_grid  <- c(100, 200, 500, 1000)
mt_grid <- c("mar", "mnar")
mr_grid <- c(0.05, 0.15, 0.30, 0.45)
R        <- 300
rep_index <- 20251108

# ------------------------------------------------------------
# (4) output directory
# ------------------------------------------------------------
data_dir <- Sys.getenv(
  "MISSDATA_DATA_DIR",
  unset = file.path("data", "missing")
)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# (5) main loop
# ------------------------------------------------------------
pb <- txtProgressBar(min = 0, max = length(N_grid) * R, style = 3)
step <- 0

for (N in N_grid) {
  
  # ---- build the big complete dataset for this N ----
  complete_all <- vector("list", R)
  
  for (rep in seq_len(R)) {
    dat <- generate_data(N, rep_index + rep)
    dat$id <- rep  # add replication ID
    dat <- dat %>% select(id, everything())
    complete_all[[rep]] <- dat
  }
  
  # bind and save
  complete_df <- bind_rows(complete_all)
  write.table(complete_df,
              file = file.path(data_dir, sprintf("complete_data_N_%d.txt", N)),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  # ---- now create the 8 missing variants for this N ----
  for (mt in mt_grid) {
    for (mr in mr_grid) {
      
      missing_all <- vector("list", R)
      
      for (rep in seq_len(R)) {
        data <- complete_all[[rep]] %>% select(-id)
        if (mt == "mar") {
          data_y <- MAR_generator_Y_X1(data, mr)
        } else {
          data_y <- MNAR_generator_Y_Self(data, mr)
        }
        data_x <- MNAR_generator_X(data, mr)
        df_miss <- cbind(id = rep, data_y, data_x)
        missing_all[[rep]] <- df_miss
      }
      
      missing_df <- bind_rows(missing_all)
      fname <- file.path(
        data_dir,
        sprintf("missing_N_%d_mt_%s_mr_%.2f.txt", N, mt, mr)
      )
      write.table(missing_df, file = fname,
                  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    }
  }
  
  step <- step + 1L
  setTxtProgressBar(pb, step)
}
close(pb)

cat("All datasets written in:", data_dir, "\n")
