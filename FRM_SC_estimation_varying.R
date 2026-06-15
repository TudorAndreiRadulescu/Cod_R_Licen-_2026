# =====================================================================
# FRM_SC_estimation_varying.R — Daily FRM estimation (varying universe)
# Reads stage_20, runs quantile regression for each bank each day,
# saves adjacency CSVs and stage_30_varying.RData.
# =====================================================================
source("FRM_SC_config.R")
source("FRM_SC_utils.R")

load(file.path(output_path, "stage_20_loaded.RData"))
# Loaded: ticker, stock_return, macro_return, mktcap_index, M_stock, M_macro, M

FRM_individ <- list()
J_dynamic   <- numeric(0)
idx_out     <- 0L

N0 <- idx_start_at_or_after(ticker, date_start)
N1 <- idx_end_at_or_before( ticker, date_end)
stopifnot(is.finite(N0), is.finite(N1), N1 > N0)

dir.create(file.path(output_path, "Adj_Matrices"), recursive=TRUE, showWarnings=FALSE)

J_eff <- min(J, M_stock)

for (t in N0:N1) {
  # Need at least s days of history before we can form a rolling window
  if (t < (N0 + s)) next

  idx_out <- idx_out + 1L

  # Top-J banks on day t (always 1:7 since mktcap_index is static)
  biggest_index <- as.integer(mktcap_index[t, 2:(J_eff+1), drop=FALSE])

  # Rolling window: the s trading days preceding day t
  rows_ret <- (t - s):(t - 1)
  if (min(rows_ret) < 1) next

  # Build X = [bank returns | macro returns]
  X <- cbind(
    stock_return[rows_ret, biggest_index, drop=FALSE],
    if (M_macro > 0) macro_return[rows_ret, , drop=FALSE] else NULL
  )
  X[!is.finite(X)] <- 0

  # Drop columns that are entirely zero in this window (causes singular matrices)
  X <- X[, colSums(X != 0) > 0, drop=FALSE]

  M_t <- ncol(X)
  J_t <- M_t - M_macro   # number of banks that survived the zero-filter
  J_dynamic[idx_out] <- J_t

  if (J_t <= 0) {
    FRM_individ[[idx_out]] <- matrix(numeric(0), nrow=1)
    next
  }

  adj_matrix  <- matrix(0, M_t, M_t)
  est_lambda_t <- rep(NA_real_, M_t)

  for (k in 1:M_t) {
    est <- tryCatch(
      FRM_Quantile_Regression(as.matrix(X), k, tau, I),
      error = function(e) NULL
    )
    if (is.null(est)) next

    gacv   <- est$Cgacv
    k_best <- if (!any(is.finite(gacv))) length(gacv) else which.min(gacv[is.finite(gacv)])

    est_lambda_t[k]   <- abs(data.matrix(est$lambda[k_best]))
    adj_matrix[k, -k] <- t(as.matrix(est$beta[k_best, ]))
  }

  # Keep only the bank (non-macro) lambdas for the FRM index
  lambda_banks <- t(data.frame(est_lambda_t[1:J_t]))
  colnames(lambda_banks) <- colnames(X)[1:J_t]
  FRM_individ[[idx_out]] <- lambda_banks

  colnames(adj_matrix) <- colnames(X)
  rownames(adj_matrix) <- colnames(X)

  write.csv(
    adj_matrix,
    file.path(output_path, "Adj_Matrices", paste0("adj_matrix_", ticker[t], ".csv")),
    quote=FALSE
  )
}

# Name each snapshot by its ISO date
if (length(FRM_individ)) {
  vary_start <- N0 + s
  snap_dates <- as.Date(as.character(ticker[vary_start:N1]), "%Y%m%d")
  names(FRM_individ) <- format(snap_dates[seq_along(FRM_individ)], "%Y-%m-%d")
}

save(FRM_individ, J_dynamic,
     file=file.path(output_path, "stage_30_varying.RData"))

message("Done. ", length(FRM_individ), " daily snapshots saved -> stage_30_varying.RData")
