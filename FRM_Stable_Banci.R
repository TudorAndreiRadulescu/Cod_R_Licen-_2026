
# ============== FRM_SC_config.R (Block 1) ==============
rm(list = ls(all = TRUE)); graphics.off(); options(digits = 6, warn = -1)
rm(list = setdiff(ls(), c()))
suppressPackageStartupMessages({ library(zoo) })

# --- Edit this single line to your local project folder ---
wdir <- "C:/Users/troiu/Desktop/fotos/licenta/R/Claude Version"
setwd(wdir)

# --- Channel & dates (YYYYMMDD integers) ---
channel <- "Banci"

# Source data span (what exists on disk)
date_start_source <- 20010103
date_end_source   <- 20251231

# Analysis span (what to compute)
date_start <- 20010103
date_end   <- 20251231

# Fixed-window span (for “Fixed” universe products)
date_start_fixed <- 20010103
date_end_fixed   <- 20251231

# --- Model knobs ---
s   <- 90
tau <- 0.05
I   <- 25
L   <- 5
J   <- 7
stock_main <- "BT"
quantiles  <- c(0.99,0.95,0.90,0.85,0.80,0.75,0.70,0.65,0.60,0.55,0.50,0.25)

# --- Paths ---
input_path  <- file.path("Input", channel, paste0(date_start_source, "-", date_end_source))
output_path <- if (tau == 0.05 & s == 90) file.path("Output", channel) else
  file.path("Output", channel, paste0("Sensitivity/tau=", 100*tau, "/s=", s))
website_path <- file.path("Website", channel)
date_banci <- read.csv("Baza_FRM_Banci.csv")
dates_vector <- date_banci$Date
data_matrix <- as.matrix(date_banci[, -which(names(date_banci) == "Date")])
banci_names <- colnames(data_matrix)

# --- Ensure folders exist ---
dirs <- c(
  output_path,
  file.path(output_path, "Adj_Matrices"),
  file.path(output_path, "Adj_Matrices/Fixed"),
  file.path(output_path, "Lambda"),
  file.path(output_path, "Lambda/Fixed"),
  file.path(output_path, "Lambda/Quantiles"),
  file.path(output_path, "Top"),
  file.path(output_path, "Network"),
  file.path(output_path, "Network/Fixed"),
  file.path(output_path, "Macro"),
  file.path(output_path, "Boxplot"),
  website_path,
  file.path(website_path, date_end)
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# --- Packages (quiet) ---
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(lubridate); library(zoo)
  library(ggplot2); library(data.table); library(igraph); library(magick); library(scales)
  library(stringr); library(plotly); library(reshape2); library(quadprog); library(MASS)
})

# --- Algorithm (required later) ---
if (!file.exists("FRM_Statistics_Algorithm.R"))
  stop("Missing FRM_Statistics_Algorithm.R in project root.")
source("FRM_Statistics_Algorithm.R")
# ========================================================

# ============== FRM_SC_utils.R (Block 2) ==============

# --- Transparent plot theme (axes only, no title, legend bottom by default) ---
theme_transparent_bottom <- theme(
  panel.background      = element_rect(fill = "white", colour = NA),
  plot.background       = element_rect(fill = "white", colour = NA),
  legend.background     = element_rect(fill = "white", colour = NA),
  legend.box.background = element_rect(fill = "white", colour = NA),
  legend.position       = "bottom",
  legend.direction      = "horizontal",
  legend.box.just       = "center",
  axis.line             = element_line(colour = "black"),
  axis.text             = element_text(size = 12, colour = "black"),
  axis.title            = element_text(size = 14, colour = "black"),
  panel.grid            = element_blank()
)

# Slightly stricter minimal variant (no legend)
theme_transparent_min <- theme_classic(base_size = 12) +
  theme(
    panel.background      = element_rect(fill = "white", colour = NA),
    plot.background       = element_rect(fill = "white", colour = NA),
    legend.position       = "none",
    axis.text             = element_text(size = 12, colour = "black"),
    axis.title            = element_text(size = 14, colour = "black")
  )

# --- Risk colors (for FRM percentile scatter) ---
risk_colors <- c(
  "1. Low risk"     = "green",
  "2. General risk" = "blue",
  "3. Elevated risk"= "yellow",
  "4. High risk"    = "orange",
  "5. Severe risk"  = "red"
)

# --- Robust, flexible file picker (tries strict regex, then fallback) ---
pick_file_flexible <- function(path, strong_regex, weak_regex, role = "file") {
  f <- list.files(path, pattern = strong_regex, full.names = TRUE, ignore.case = TRUE)
  if (!length(f)) f <- list.files(path, pattern = weak_regex, full.names = TRUE, ignore.case = TRUE)
  if (!length(f)) stop(sprintf("No %s found in %s", role, normalizePath(path)))
  # pick the one with the latest 8-digit date in name, else first
  extract_max_date_from_name <- function(x) {
    b <- basename(x)
    ds <- gregexpr("[0-9]{8}", b, perl = TRUE); vals <- regmatches(b, ds)[[1]]
    if (length(vals)) as.integer(max(vals)) else NA_integer_
  }
  name_dates <- vapply(f, extract_max_date_from_name, integer(1))
  if (any(is.na(name_dates))) return(f[which.max(file.info(f)$mtime)])
  f[which.max(name_dates)]
}

# --- Pretty bank names ---
BANK_NAME_MAP <- c(
  AlphaBank  = "Alpha Bank",
  BNPPariBas = "BNP Paribas",
  BRD        = "BRD",
  BT         = "Banca Transilvania",
  Erste      = "Erste Bank",
  ING        = "ING Bank",
  Unicredit  = "UniCredit"
)

pretty_bank <- function(x) {
  y <- BANK_NAME_MAP[x]
  ifelse(is.na(y), x, y)
}

# --- Date helpers (robust parsing + index finders) ---
parse_date_robust <- function(x) {
  d <- suppressWarnings(as.Date(x))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, "%Y-%m-%d"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, "%m/%d/%Y"))
  if (all(is.na(d))) {
    dt <- suppressWarnings(as.POSIXct(x, tz = "UTC"))
    if (!all(is.na(dt))) d <- as.Date(dt, tz = "UTC")
  }
  d
}

idx_start_at_or_after <- function(key, target) {
  w <- which(key >= target); if (length(w)) w[1] else NA_integer_
}
idx_end_at_or_before <- function(key, target) {
  w <- which(key <= target); if (length(w)) tail(w, 1) else NA_integer_
}

# --- Small numeric utils ---
clamp01 <- function(v) pmin(pmax(v, 0), 1)


# ===============================================================
# Build FIXED adjacency matrices from CSV inputs (no RData needed)
# - Uses FRM_SC_config.R (paths, dates, s/tau/I/J, channel)
# - Reads your input CSVs via FRM_SC_load_data.R (if needed)
# - Writes: Output/<channel>/Adj_Matrices/Fixed/adj_matrix_YYYYMMDD.csv
# - Writes: Output/<channel>/Lambda/Fixed/lambdas_fixed_<start>_<end>.csv
# ===============================================================

rm(list = setdiff(ls(), c()))
suppressPackageStartupMessages({ library(zoo) })

# ---- config + utils + FRM algorithm ----
source("FRM_SC_config.R")
source("FRM_SC_utils.R")
if (!file.exists("FRM_Statistics_Algorithm.R"))
  stop("Missing FRM_Statistics_Algorithm.R")
source("FRM_Statistics_Algorithm.R")

# ---- load stage_20 or build it from CSVs once ----
st20 <- file.path(output_path, "stage_20_loaded.RData")
if (!file.exists(st20)) {
  message("stage_20_loaded.RData not found — running FRM_SC_load_data.R once...")
  source("FRM_SC_load_data.R")
}
load(file.path(output_path, "stage_20_loaded.RData"))
# objects we expect now: dates, ticker, stock_return, macro_return, M_stock, M_macro, mktcap_index, s, tau, I, J

# ---- ensure output folders ----
fixed_dir  <- file.path(output_path, "Adj_Matrices", "Fixed")
lambda_dir <- file.path(output_path, "Lambda", "Fixed")
dir.create(fixed_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(lambda_dir, recursive = TRUE, showWarnings = FALSE)

# ---- helpers for date indexing ----
idx_start_at_or_after <- function(key, target){ w <- which(key >= target); if (length(w)) w[1] else NA_integer_ }
idx_end_at_or_before  <- function(key, target){ w <- which(key <= target); if (length(w)) tail(w,1) else NA_integer_ }

# ---- fixed window from config ----
N0_fixed <- idx_start_at_or_after(ticker, date_start_fixed)
N1_fixed <- idx_end_at_or_before( ticker, date_end_fixed)
if (is.na(N0_fixed) || is.na(N1_fixed) || (N1_fixed - N0_fixed + 1) < s) {
  stop(sprintf("Fixed window [%s, %s] too short for s=%d.", date_start_fixed, date_end_fixed, s))
}

# ---- lock Top-J once at the START of the fixed window ----
J_eff <- min(J, M_stock)
biggest_index_fixed <- as.integer(mktcap_index[N0_fixed, 2:(J_eff + 1), drop = FALSE])
M_J <- J_eff + M_macro

bank_names  <- colnames(stock_return)
macro_names <- if (M_macro > 0) colnames(macro_return) else character(0)
col_full    <- c(bank_names[biggest_index_fixed], macro_names)

# ---- iterate days in fixed window and create CSV adjacencies ----
N_fixed <- N1_fixed - N0_fixed + 1
lambdas_fixed <- matrix(NA_real_, N_fixed, J_eff + 1)
lambdas_fixed[, 1] <- ticker[N0_fixed:N1_fixed]

n_written <- 0L
for (t in N0_fixed:N1_fixed) {
  rows_ret <- (t - s):(t - 1)
  if (min(rows_ret) < 1) next
  
  X <- cbind(
    stock_return[rows_ret, biggest_index_fixed, drop = FALSE],
    if (M_macro > 0) macro_return[rows_ret, , drop = FALSE] else NULL
  )
  mode(X) <- "numeric"; X[!is.finite(X)] <- 0
  
  # 1. Assign names to X BEFORE dropping any columns
  colnames(X) <- col_full 
  
  # Drop only columns that are entirely zero 
  keep <- colSums(X != 0) > 0
  X <- X[, keep, drop = FALSE]
  if (ncol(X) == 0) next
  
  A <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  lam_row <- rep(NA_real_, J_eff)
  
  for (k in 1:ncol(X)) {
    est <- tryCatch(FRM_Quantile_Regression(as.matrix(X), k, tau, I), error = function(e) NULL)
    if (is.null(est)) next
    gacv   <- est$Cgacv
    best   <- if (!any(is.finite(gacv))) length(gacv) else which.min(gacv[is.finite(gacv)])
    est_lambda <- abs(data.matrix(est$lambda[best]))
    est_beta   <- t(as.matrix(est$beta[best, ]))
    A[k, -k] <- est_beta
    
    # Keep track of original column index for lam_row if needed. 
    # (Note: if you dropped columns, 'k' no longer maps 1:1 to J_eff, 
    # so ensure this step behaves as you intend mathematically).
    if (k <= J_eff) lam_row[k] <- est_lambda 
  }
  
  colnames(A) <- colnames(X)
  rownames(A) <- colnames(X)
  
  out_csv <- file.path(fixed_dir, paste0("adj_matrix_", ticker[t], ".csv"))
  write.csv(A, out_csv, quote = FALSE)
  n_written <- n_written + 1L
  
  lambdas_fixed[t - N0_fixed + 1, 2:(J_eff + 1)] <- lam_row
}

# ---- save helper lambdas CSV for the fixed universe ----
colnames(lambdas_fixed) <- c("date", bank_names[biggest_index_fixed])
lam_csv <- file.path(lambda_dir, paste0("lambdas_fixed_", date_start_fixed, "_", date_end_fixed, ".csv"))
write.csv(lambdas_fixed, lam_csv, row.names = FALSE, quote = FALSE)

message("✅ Fixed adjacency CSVs written: ", n_written)
message("📁 Fixed folder: ", normalizePath(fixed_dir))
message("🧾 Lambdas CSV:  ", normalizePath(lam_csv))


# ======================================================

# ============== FRM_SC_load_data.R (Block 3) ==============
source("FRM_SC_config.R")
source("FRM_SC_utils.R")

cat("\n[1] Incarcam datele preparate din Python (Baza_FRM_Banci.csv)...\n")

df_baza <- read.csv("Baza_FRM_Banci.csv")

# Setam data si tickerul
dates <- as.Date(as.character(df_baza$Date), format="%Y%m%d")
ticker <- as.integer(df_baza$Date)
N <- length(ticker)

# Separam bancile de macro (Numele corespund exact cu cele din CSV-ul tau)
banci_cols <- c("AlphaBank", "BNPPariBas", "BRD", "BT", "Erste", "ING", "Unicredit")
macro_cols <- setdiff(colnames(df_baza), c("Date", banci_cols))

cat("[2] Construim matricile de randamente...\n")
stock_return <- as.matrix(df_baza[, banci_cols])
macro_return <- if(length(macro_cols) > 0) as.matrix(df_baza[, macro_cols]) else matrix(, N, 0)

M_stock <- ncol(stock_return)
M_macro <- ncol(macro_return)
M <- M_stock + M_macro

# Simulam indexul de capitalizare (le spunem algoritmului sa foloseasca mereu toate cele 7 banci)
mktcap_index_mat <- matrix(rep(1:M_stock, each=N), nrow=N, ncol=M_stock)
mktcap_index <- cbind(ticker, mktcap_index_mat)

# Cream mediul RData cerut de blocurile urmatoare
cat("[3] Salvam datele in stage_20_loaded.RData...\n")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
save(ticker, stock_return, macro_return, mktcap_index, M_stock, M_macro, M,
     file = file.path(output_path, "stage_20_loaded.RData"))
cat("✅ GATA! Datele au fost incarcate cu succes.\n\n")
# ============== FRM_SC_estimation_varying.R (Block 4) ==============
# Builds daily adjacency matrices (varying Top-J by market cap)
# and a per-day lambda list (FRM_individ), then saves stage_30.

# --- env & deps ---
source("FRM_SC_config.R")
if (file.exists("FRM_SC_utils.R")) source("FRM_SC_utils.R")
stopifnot(file.exists("FRM_Statistics_Algorithm.R"))
source("FRM_Statistics_Algorithm.R")

st20 <- file.path(output_path, "stage_20_loaded.RData")
stopifnot(file.exists(st20))
load(st20)
# Expect: dates, ticker, stock_return, macro_return, M_stock, M_macro, mktcap_index, s, tau, I, J

# --- helpers (local) ---
idx_start_at_or_after <- function(key, target){ w <- which(key >= target); if (length(w)) w[1] else NA_integer_ }
idx_end_at_or_before  <- function(key, target){ w <- which(key <= target); if (length(w)) tail(w,1) else NA_integer_ }

# --- time window ---
N0 <- idx_start_at_or_after(ticker, date_start)
N1 <- idx_end_at_or_before( ticker, date_end)
stopifnot(is.finite(N0), is.finite(N1), N1 > N0)

# --- ensure output dir for matrices ---
dir.create(file.path(output_path, "Adj_Matrices"), recursive = TRUE, showWarnings = FALSE)

FRM_individ <- list()
J_dynamic   <- integer(0)
idx_out     <- 0L

J_eff <- min(J, M_stock)

for (t in N0:N1) {
  # need at least s days of history
  if (t < (N0 + s)) next
  if ((t - N0) %% 100 == 0)
    message(round(100 * (t - N0) / (N1 - N0), 1), "% — ziua ", ticker[t])
  # Top-J (by market cap) *on day t*
  biggest_index <- as.integer(mktcap_index[t, 2:(J_eff + 1), drop = FALSE])
  
  # returns window
  rows_ret <- (t - s):(t - 1)
  if (min(rows_ret) < 1) next
  
  # assemble X = [banks_in_topJ  |  macro (if any)]
  bank_cols_all  <- colnames(stock_return)[biggest_index]
  macro_cols_all <- if (M_macro > 0) colnames(macro_return) else character(0)
  
  X <- cbind(
    stock_return[rows_ret, biggest_index, drop = FALSE],
    if (M_macro > 0) macro_return[rows_ret, , drop = FALSE] else NULL
  )
  X[!is.finite(X)] <- 0
  # drop columns that are all zeros in the window
  keep <- colSums(X != 0) > 0
  X <- X[, keep, drop = FALSE]
  
  if (!ncol(X)) next
  # identify which remaining cols are COINS (exclude macro for J_t)
  bank_pos <- which(colnames(X) %in% bank_cols_all)
  J_t      <- length(bank_pos)
  J_dynamic <- c(J_dynamic, J_t)
  
  # skip if no banks survived
  if (J_t <= 0) {
    idx_out <- idx_out + 1L
    FRM_individ[[idx_out]] <- matrix(numeric(0), nrow = 1)
    next
  }
  
  # estimate per target k
  M_t <- ncol(X)
  A   <- matrix(0, M_t, M_t)
  lam_vec <- rep(NA_real_, M_t)
  
  for (k in 1:M_t) {
    est <- tryCatch(FRM_Quantile_Regression(as.matrix(X), k, tau, I), error = function(e) NULL)
    if (is.null(est)) next
    
    gacv <- est$Cgacv
    k_best <- if (!any(is.finite(gacv))) length(gacv) else which.min(gacv[is.finite(gacv)])
    
    # lambda (scalar per target), betas (row without diagonal)
    est_lambda <- abs(data.matrix(est$lambda[k_best]))
    est_beta   <- t(as.matrix(est$beta[k_best, ]))
    
    A[k, -k] <- est_beta
    lam_vec[k] <- est_lambda
  }
  
  colnames(A) <- colnames(X)
  rownames(A) <- colnames(X)
  
  # keep lambda entries for COINS only, in bank_pos order
  lam_banks <- t(as.data.frame(lam_vec[bank_pos]))
  colnames(lam_banks) <- colnames(X)[bank_pos]
  
  # store snapshot
  idx_out <- idx_out + 1L
  FRM_individ[[idx_out]] <- lam_banks
  
  # write adjacency CSV for this date
  out_csv <- file.path(output_path, "Adj_Matrices", paste0("adj_matrix_", ticker[t], ".csv"))
  write.csv(A, out_csv, quote = FALSE)
}

# --- name snapshots by date (ISO) and save stage_30 ---
if (length(FRM_individ)) {
  vary_start <- (N0 + s)
  snap_dates <- as.Date(as.character(ticker[vary_start:N1]), "%Y%m%d")
  names(FRM_individ) <- format(snap_dates[seq_along(FRM_individ)], "%Y-%m-%d")
}

save(FRM_individ, J_dynamic, file = file.path(output_path, "stage_30_varying.RData"))
message("✓ varying estimation complete: ", length(FRM_individ), " snapshots -> stage_30_varying.RData")
# Dependencies & data
source("FRM_SC_config.R")
if (file.exists("FRM_SC_utils.R")) source("FRM_SC_utils.R")
source("FRM_Statistics_Algorithm.R")
if (!file.exists(file.path(output_path, "stage_20_loaded.RData"))) source("FRM_SC_load_data.R")
load(file.path(output_path, "stage_20_loaded.RData"))

build_fixed_all <- function(){
  fixed_dir <- file.path(output_path, "Adj_Matrices", "Fixed")
  dir.create(fixed_dir, recursive = TRUE, showWarnings = FALSE)
  
  # helpers
  idx_start_at_or_after <- function(key, target){ w <- which(key >= target); if (length(w)) w[1] else NA_integer_ }
  idx_end_at_or_before  <- function(key, target){ w <- which(key <= target); if (length(w)) tail(w,1) else NA_integer_ }
  
  N0 <- idx_start_at_or_after(ticker, date_start_fixed)
  N1 <- idx_end_at_or_before( ticker, date_end_fixed)
  stopifnot(is.finite(N0), is.finite(N1), N1 > N0, s > 0)
  
  J_eff <- min(J, ncol(stock_return))
  biggest_index_fixed <- as.integer(mktcap_index[N0, 2:(J_eff + 1), drop = FALSE])
  
  n_written <- 0L
  for (t in N0:N1) {
    rows <- (t - s):(t - 1); if (min(rows) < 1) next
    X <- cbind(
      stock_return[rows, biggest_index_fixed, drop = FALSE],
      if (exists("macro_return") && ncol(macro_return)) macro_return[rows, , drop = FALSE]
    )
    X[!is.finite(X)] <- 0
    
    keep <- colSums(X != 0) > 0
    X <- X[, keep, drop = FALSE]
    if (ncol(X) == 0) next
    
    # col_full derivat din X dupa filtrare, nu din lista originala
    col_full <- colnames(X)   # ← MUTAT AICI
    
    M <- ncol(X); A <- matrix(0, M, M)
    for (k in 1:M) {
      est <- tryCatch(FRM_Quantile_Regression(as.matrix(X), k, tau, I), error = function(e) NULL)
      if (is.null(est)) next
      gacv <- est$Cgacv
      best <- if (!any(is.finite(gacv))) length(gacv) else which.min(gacv[is.finite(gacv)])
      A[k, -k] <- t(as.matrix(est$beta[best, ]))
    }
    colnames(A) <- rownames(A) <- colnames(X)
    
    out_csv <- file.path(fixed_dir, paste0("adj_matrix_", ticker[t], ".csv"))
    write.csv(A, out_csv, quote = FALSE)
    n_written <- n_written + 1L
  }
  message("✅ Fixed adjacency CSVs written: ", n_written,
          " -> ", normalizePath(fixed_dir, mustWork = FALSE))
}

build_fixed_all()


# ============== FRM_SC_history_outputs.R (Block 5) ==============
# - Appends new daily lambdas to history (RDS)
# - Writes lambdas_wide.csv and FRM_<channel>_index.csv
# - Saves transparent risk-colored scatter PNG (no title, no "Date" label)
# - Saves stage_50_index.RData

# --- env & deps ---
source("FRM_SC_config.R")
source("FRM_SC_utils.R")

st20 <- file.path(output_path, "stage_20_loaded.RData")
st30 <- file.path(output_path, "stage_30_varying.RData")
stopifnot(file.exists(st20), file.exists(st30))

load(st20)  # brings: dates, etc. (not strictly required here)
load(st30)  # brings: FRM_individ (named list yyyy-mm-dd -> 1xJ_t lambdas)

dir.create(file.path(output_path, "Lambda"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(website_path, date_end), recursive = TRUE, showWarnings = FALSE)

# --- append to history RDS (dedupe by date; latest wins) ---
rds_path <- file.path(output_path, "Lambda", paste0("FRM_", channel, ".rds"))
FRM_history <- if (file.exists(rds_path)) c(readRDS(rds_path), FRM_individ) else FRM_individ

# ensure names are Date strings and de-duplicate (keep latest snapshot per date)
nm <- names(FRM_history)
if (is.null(nm) || any(!nzchar(nm))) {
  # fallback: try to derive from list order if names missing
  names(FRM_history) <- nm <- format(seq_along(FRM_history), "%Y-%m-%d")
}
# Keep last occurrence for each date
keep <- !duplicated(nm, fromLast = TRUE)
FRM_history <- FRM_history[keep]

# order by date
iso_dates <- suppressWarnings(as.Date(names(FRM_history)))
ord <- order(iso_dates)
FRM_history <- FRM_history[ord]
saveRDS(FRM_history, rds_path)

# --- lambdas_wide.csv (rows = dates, cols = all unique banks) ---
stock_names <- unique(unlist(lapply(FRM_history, colnames)))
if (length(stock_names)) {
  N_h <- length(FRM_history)
  out <- matrix(0, nrow = N_h, ncol = length(stock_names) + 1)
  for (k in seq_along(stock_names)) {
    ck <- stock_names[k]
    for (t in seq_len(N_h)) {
      mt <- FRM_history[[t]]
      if (!is.null(mt) && ncol(mt) && ck %in% colnames(mt)) out[t, k + 1] <- mt[1, ck]
    }
  }
  out[, 1] <- format(as.Date(names(FRM_history)), "%Y-%m-%d")
  colnames(out) <- c("date", stock_names)
  write.csv(out, file.path(output_path, "Lambda", "lambdas_wide.csv"), row.names = FALSE)
}

# --- FRM index (raw mean of available lambdas each day) ---
FRM_index <- data.frame(
  date = as.Date(names(FRM_history)),
  frm  = vapply(
    FRM_history,
    function(m){
      v <- suppressWarnings(as.numeric(m[1, ]))
      v <- v[is.finite(v)]
      if (!length(v)) NA_real_ else mean(v)
    },
    numeric(1)
  )
)

write.csv(FRM_index, file.path(output_path, "Lambda", paste0("FRM_", channel, "_index.csv")),
          row.names = FALSE)

# --- Risk coloring using ECDF of available FRM values ---
good_frm <- FRM_index$frm[is.finite(FRM_index$frm)]
risk_ecdf <- if (length(good_frm)) ecdf(good_frm) else function(x) NA_real_

FRM_plot <- FRM_index
FRM_plot$risk_pct <- round(100 * risk_ecdf(FRM_plot$frm), 2)
FRM_plot$`Risk level` <- factor(
  ifelse(FRM_plot$risk_pct < 20, "1. Low risk",
         ifelse(FRM_plot$risk_pct < 40, "2. General risk",
                ifelse(FRM_plot$risk_pct < 60, "3. Elevated risk",
                       ifelse(FRM_plot$risk_pct < 80, "4. High risk", "5. Severe risk")))),
  levels = names(risk_colors)
)

# --- Transparent risk scatter (no title, no x label "Date") ---
png(file.path(website_path, date_end, paste0("FRMColor_", channel, ".png")),
    width = 900, height = 600, bg = "white")
print(
  ggplot(FRM_plot, aes(x = date, y = frm)) +
    geom_point(aes(color = `Risk level`), size = 1, na.rm = TRUE) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    scale_color_manual(values = risk_colors) +
    labs(title = NULL, x = NULL, y = paste0("FRM@", channel)) +
    theme_transparent_bottom
)
dev.off()

# --- Standardizare FRM in [0,1] pentru scatter ---
frm_vals <- FRM_plot$frm
frm_min  <- min(frm_vals, na.rm = TRUE)
frm_max  <- max(frm_vals, na.rm = TRUE)
FRM_plot$frm_norm <- (frm_vals - frm_min) / (frm_max - frm_min)

# --- Reface colorarea pe valorile normalizate ---
good_norm  <- FRM_plot$frm_norm[is.finite(FRM_plot$frm_norm)]
risk_ecdf  <- if (length(good_norm)) ecdf(good_norm) else function(x) NA_real_
FRM_plot$risk_pct <- round(100 * risk_ecdf(FRM_plot$frm_norm), 2)
FRM_plot$`Risk level` <- factor(
  ifelse(FRM_plot$risk_pct < 20, "1. Low risk",
         ifelse(FRM_plot$risk_pct < 40, "2. General risk",
                ifelse(FRM_plot$risk_pct < 60, "3. Elevated risk",
                       ifelse(FRM_plot$risk_pct < 80, "4. High risk", "5. Severe risk")))),
  levels = names(risk_colors)
)

png(file.path(website_path, date_end, paste0("FRMColor_", channel, ".png")),
    width=900, height=600, bg="white")
print(
  ggplot(FRM_plot, aes(x=date, y=frm_norm)) +   # ← frm_norm nu frm
    geom_point(aes(color=`Risk level`), size=1, na.rm=TRUE) +
    scale_x_date(
      breaks = seq(as.Date("2001-01-01"), as.Date("2026-01-01"), by = "1 year"),
      date_labels = "%Y"
    ) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) +
    scale_color_manual(values=risk_colors) +
    labs(title=NULL, x=NULL, y=paste0("FRM@", channel, " (normalized)")) +
    theme_transparent_bottom +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
)
dev.off()

# --- Stage 50 for downstream plots ---
save(FRM_index, file = file.path(output_path, "stage_50_index.RData"))
message("✓ history updated: ",
        nrow(FRM_index), " points | wrote lambdas_wide.csv, FRM_", channel, "_index.csv, ",
        "FRMColor_", channel, ".png, and stage_50_index.RData")
# ============== Granger Causality (Block 5.1) ==============

#test for stationarity with adf
library(tseries)
adf.test(na.omit(dat_g$frm))          # p < 0.05 = stationar = ok
adf.test(na.omit(dat_g$BT))           # idem per banca

#daca nu e stationar, facem prima diferenta
#dat_g$frm_diff <- c(NA, diff(dat_g$frm))

# =====================================================================
# CAUZALITATE GRANGER — FRM vs Lambda-uri individuale per banca
# Testeaza: (1) banca -> FRM  si  (2) FRM -> banca
# La lag-uri: 1 zi, 5 zile (1 sapt), 10 zile (2 sapt), 22 zile (1 luna)
# Output: Granger_FRM_vs_Banks.csv + Granger_Heatmap_lag5.png
# =====================================================================

# --- 1. Dependente & context ---
source("FRM_SC_config.R")
source("FRM_SC_utils.R")

if (!requireNamespace("lmtest",   quietly = TRUE)) install.packages("lmtest")
if (!requireNamespace("tseries",  quietly = TRUE)) install.packages("tseries")
if (!requireNamespace("ggplot2",  quietly = TRUE)) install.packages("ggplot2")

suppressPackageStartupMessages({
  library(lmtest)
  library(tseries)
  library(ggplot2)
})

# --- 2. Incarca FRM_index ---
if (!exists("FRM_index") || !is.data.frame(FRM_index)) {
  st50 <- file.path(output_path, "stage_50_index.RData")
  if (file.exists(st50)) {
    load(st50)
    message("FRM_index incarcat din stage_50_index.RData")
  } else {
    frm_csv <- file.path(output_path, "Lambda", paste0("FRM_", channel, "_index.csv"))
    if (!file.exists(frm_csv))
      stop("Nu gasesc FRM_index — ruleaza mai intai Block 5 (FRM_SC_history_outputs.R)")
    FRM_index <- read.csv(frm_csv, stringsAsFactors = FALSE)
    message("FRM_index incarcat din CSV")
  }
}
FRM_index$date <- as.Date(FRM_index$date)
FRM_index <- FRM_index[is.finite(FRM_index$frm), ]

# --- 3. Incarca lambdas_wide ---
lam_csv <- file.path(output_path, "Lambda", "lambdas_wide.csv")
if (!file.exists(lam_csv))
  stop("lambdas_wide.csv nu exista — ruleaza mai intai Block 5")
LAM_df <- read.csv(lam_csv, stringsAsFactors = FALSE)
LAM_df$date <- as.Date(LAM_df$date)
message("lambdas_wide.csv incarcat: ", nrow(LAM_df), " randuri, ", ncol(LAM_df)-1, " banci")

# --- 4. Merge si curatare ---
dat_g <- merge(FRM_index[, c("date","frm")], LAM_df, by = "date", all = FALSE)
dat_g <- dat_g[order(dat_g$date), ]
dat_g <- dat_g[complete.cases(dat_g), ]
message("Date pentru Granger: ", nrow(dat_g), " obs | ",
        as.character(min(dat_g$date)), " -> ", as.character(max(dat_g$date)))

# --- 5. Verifica stationaritatea (ADF test) ---
cat("\n=== TEST STATIONARITATE (ADF) ===\n")
cat("H0: seria are radacina unitara (non-stationara)\n")
cat("p < 0.05 -> stationara -> Granger valid direct\n\n")

adf_frm <- adf.test(dat_g$frm, alternative = "stationary")
cat(sprintf("FRM index:  ADF = %.4f,  p = %.4f  %s\n",
            adf_frm$statistic, adf_frm$p.value,
            ifelse(adf_frm$p.value < 0.05, "✓ stationar", "⚠ NON-stationar")))

banci_cols <- c("AlphaBank","BNPPariBas","BRD","BT","Erste","ING","Unicredit")
banci_in_dat <- intersect(banci_cols, names(dat_g))

non_stationary <- character(0)
for (bank in banci_in_dat) {
  adf_b <- suppressWarnings(adf.test(dat_g[[bank]], alternative = "stationary"))
  status <- ifelse(adf_b$p.value < 0.05, "✓ stationar", "⚠ NON-stationar")
  cat(sprintf("%-12s:  ADF = %.4f,  p = %.4f  %s\n",
              bank, adf_b$statistic, adf_b$p.value, status))
  if (adf_b$p.value >= 0.05) non_stationary <- c(non_stationary, bank)
}

# Daca FRM e non-stationar, folosim prima diferenta
use_diff_frm <- adf_frm$p.value >= 0.05
if (use_diff_frm) {
  dat_g$frm_use <- c(NA, diff(dat_g$frm))
  cat("\n⚠ FRM non-stationar — se foloseste prima diferenta (d.FRM)\n")
} else {
  dat_g$frm_use <- dat_g$frm
  cat("\n✓ FRM stationar — se folosesc valorile originale\n")
}

# Diferenta pentru banci non-stationare
for (bank in banci_in_dat) {
  col_use <- paste0(bank, "_use")
  if (bank %in% non_stationary) {
    dat_g[[col_use]] <- c(NA, diff(dat_g[[bank]]))
  } else {
    dat_g[[col_use]] <- dat_g[[bank]]
  }
}

# Elimina NA-urile introduse de diff
dat_g <- dat_g[complete.cases(dat_g[, c("frm_use",
                                        paste0(banci_in_dat, "_use"))]), ]
message("\nObs dupa diferentiere: ", nrow(dat_g))

# --- 6. Teste Granger ---
lags_to_test <- c(1, 5, 10, 22)
results_list <- list()

cat("\n=== CAUZALITATE GRANGER ===\n")
cat("*** p<0.01  ** p<0.05  * p<0.10\n")

for (lag in lags_to_test) {
  cat(sprintf("\n--- LAG = %d zile ---\n", lag))
  
  for (bank in banci_in_dat) {
    col_use <- paste0(bank, "_use")
    
    # Directia 1: banca -> FRM
    tryCatch({
      gt1 <- grangertest(dat_g$frm_use ~ dat_g[[col_use]], order = lag)
      p1  <- gt1$`Pr(>F)`[2]
      f1  <- gt1$F[2]
      sig1 <- ifelse(p1 < 0.01, "***", ifelse(p1 < 0.05, "**",
                                              ifelse(p1 < 0.10, "*", "")))
      cat(sprintf("  %-12s -> FRM  : F=%.3f  p=%.4f %s\n", bank, f1, p1, sig1))
      results_list[[length(results_list)+1]] <- data.frame(
        lag=lag, direction=paste0(bank," -> FRM"),
        from=bank, to="FRM",
        F_stat=round(f1,4), p_value=round(p1,4), signif=sig1,
        stringsAsFactors=FALSE
      )
    }, error = function(e) cat(sprintf("  %-12s -> FRM  : ERROR - %s\n", bank, e$message)))
    
    # Directia 2: FRM -> banca
    tryCatch({
      gt2 <- grangertest(dat_g[[col_use]] ~ dat_g$frm_use, order = lag)
      p2  <- gt2$`Pr(>F)`[2]
      f2  <- gt2$F[2]
      sig2 <- ifelse(p2 < 0.01, "***", ifelse(p2 < 0.05, "**",
                                              ifelse(p2 < 0.10, "*", "")))
      cat(sprintf("  FRM -> %-12s : F=%.3f  p=%.4f %s\n", bank, f2, p2, sig2))
      results_list[[length(results_list)+1]] <- data.frame(
        lag=lag, direction=paste0("FRM -> ",bank),
        from="FRM", to=bank,
        F_stat=round(f2,4), p_value=round(p2,4), signif=sig2,
        stringsAsFactors=FALSE
      )
    }, error = function(e) cat(sprintf("  FRM -> %-12s : ERROR - %s\n", bank, e$message)))
  }
}

granger_results <- do.call(rbind, results_list)

# --- 7. Salveaza CSV ---
dir.create(file.path(output_path, "Lambda"), recursive=TRUE, showWarnings=FALSE)
out_csv <- file.path(output_path, "Lambda", "Granger_FRM_vs_Banks.csv")
write.csv(granger_results, out_csv, row.names=FALSE)
message("\nGranger CSV salvat: ", normalizePath(out_csv, mustWork=FALSE))

# --- 8. Heatmap p-values pentru lag=5 ---
dir.create(file.path(website_path, date_end), recursive=TRUE, showWarnings=FALSE)

for (lag_plot in c(5, 22)) {
  g_sub <- granger_results[granger_results$lag == lag_plot, ]
  if (!nrow(g_sub)) next
  
  # Ordinea nodurilor: FRM primul, apoi banci alfabetic
  node_order <- c("FRM", sort(banci_in_dat))
  g_sub$from <- factor(g_sub$from, levels = rev(node_order))
  g_sub$to   <- factor(g_sub$to,   levels = node_order)
  
  out_png <- file.path(website_path, date_end,
                       paste0("Granger_Heatmap_lag", lag_plot, ".png"))
  png(out_png, width=1000, height=800, bg="white", res=110)
  print(
    ggplot(g_sub, aes(x=to, y=from, fill=-log10(p_value+1e-10))) +
      geom_tile(color="white", linewidth=0.8) +
      geom_text(aes(label=paste0(signif, "\np=", sprintf("%.3f", p_value))),
                size=3.2, lineheight=1.2) +
      scale_fill_gradient2(
        low="white", mid="#FCD34D", high="#EF4444",
        midpoint=1.3,                # -log10(0.05) = 1.3
        name="-log10(p)",
        limits=c(0, 4)
      ) +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0)) +
      labs(title=NULL,
           x=paste0("Cauzeaza  (lag=", lag_plot, " zile)"),
           y="Este cauzat de") +
      theme_transparent_bottom +
      theme(
        axis.text.x  = element_text(angle=45, hjust=1, size=11),
        axis.text.y  = element_text(size=11),
        panel.border = element_rect(fill=NA, color="black", linewidth=0.8)
      )
  )
  dev.off()
  message("Granger heatmap salvat: ", normalizePath(out_png, mustWork=FALSE))
}

message("\n✓ Bloc Granger complet.")


# ============== FRM_SC_frm_plot_stable.R (Block 6) ==============
# Transparent blue FRM line, no title, no "Date" x-label.

source("FRM_SC_utils.R")
load(file.path(output_path, "stage_50_index.RData"))  # brings FRM_index (date, frm)

dir.create(file.path(website_path, date_end), recursive = TRUE, showWarnings = FALSE)
out_png <- file.path(website_path, date_end, "FRM_Banci_Index.png")

png(out_png, width = 1200, height = 700, bg = "white")
print(
  ggplot(FRM_index, aes(x = date, y = frm)) +
    geom_line(linewidth = 0.9, color = "blue") +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(title = NULL, x = NULL, y = paste0("FRM@", channel)) +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
)
dev.off()
message("✓ wrote: ", normalizePath(out_png))

# ============== FRM_SC_adj_heatmap.R — crash-day helper ==============
# Make a heatmap for a specific date (or auto-pick the crash day).
# Defaults to Fixed adjacencies. Transparent background; no title.

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

source("FRM_SC_config.R")
source("FRM_SC_utils.R")

# ---------- core heatmap (labels in CAPS) ----------
# ---------- core heatmap (labels in CAPS) ----------
suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

source("FRM_SC_config.R")
source("FRM_SC_utils.R")

# ---------- core heatmap (labels in CAPS; zeros blank) ----------
suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

source("FRM_SC_utils.R")
source("FRM_SC_config.R")

plot_adj_matrix <- function(
    adj_file,
    out_dir   = file.path(output_path, "Network"),
    digits    = 2,
    zero_eps  = 1e-3,   # |value| < zero_eps -> blank cell (no color, no label)
    text_size = 3
) {
  stopifnot(file.exists(adj_file))
  
  date_str <- gsub("adj_matrix_|\\.csv", "", basename(adj_file))
  mat <- as.matrix(read.csv(adj_file, row.names = 1, check.names = FALSE))
  mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  
  L <- max(abs(mat), na.rm = TRUE); if (!is.finite(L) || L == 0) L <- 1e-8
  
  rn <- rownames(mat); cn <- colnames(mat)
  label_caps <- function(v) toupper(v)
  rn_cap <- label_caps(rn); cn_cap <- label_caps(cn)
  
  df_long <- reshape2::melt(mat, varnames = c("From","To"), value.name = "Weight")
  df_long$From <- factor(df_long$From, levels = rn, labels = rn_cap)
  df_long$To   <- factor(df_long$To,   levels = cn, labels = cn_cap)
  
  # blank near-zeros
  df_long$Weight_plot <- df_long$Weight
  df_long$Weight_plot[abs(df_long$Weight_plot) < zero_eps] <- NA
  df_long$lab <- ifelse(abs(df_long$Weight) < zero_eps, "", sprintf(paste0("%.", digits, "f"), df_long$Weight))
  df_long$lab_col <- ifelse(abs(df_long$Weight) > 0.6 * L, "white", "black")
  
  p <- ggplot(df_long, aes(x = To, y = From, fill = Weight_plot)) +
    geom_tile(color = NA) +                                   # <-- no borders
    geom_text(aes(label = lab, colour = I(lab_col)), size = text_size) +
    scale_fill_gradient2(
      low  = "#3B82F6", mid = "white", high = "#EF4444",
      midpoint = 0, limits = c(-L, L), oob = scales::squish,
      na.value = "transparent", name = NULL
    ) +
    scale_x_discrete(expand = c(0, 0)) +                      # no outer padding
    scale_y_discrete(expand = c(0, 0)) +
    labs(title = NULL, x = "To", y = "From") +
    coord_fixed() +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),                     # <-- no grid
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background      = element_rect(fill = "white", colour = NA),
      plot.background       = element_rect(fill = "white", colour = NA),
      legend.box.background = element_rect(fill = "white"),
      legend.background     = element_rect(fill = "white")
    )
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_png <- file.path(out_dir, paste0("AdjMatrix_", date_str, ".png"))
  ggsave(out_png, p, width = 8, height = 6, bg = "white", dpi = 300)
  message("✓ saved: ", normalizePath(out_png))
  invisible(p)
}


# ---------- helpers ----------
to_int_yyyymmdd <- function(x){
  if (is.numeric(x)) return(as.integer(x))
  d <- suppressWarnings(as.Date(x))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, "%Y-%m-%d"))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, "%m/%d/%Y"))
  if (is.na(d)) stop("Cannot parse date: ", x)
  as.integer(format(d, "%Y%m%d"))
}

list_adj_files <- function(fixed = TRUE){
  dir <- if (fixed) file.path(output_path, "Adj_Matrices", "Fixed") else file.path(output_path, "Adj_Matrices")
  if (!dir.exists(dir)) return(character(0))
  list.files(dir, pattern = "^[Aa]dj_?matrix_[0-9]{8}\\.csv$", full.names = TRUE)
}

extract_date_from_name <- function(fp){
  b <- basename(fp)
  m <- regexpr("[0-9]{8}", b, perl = TRUE)
  if (m[1] == -1) return(NA_integer_)
  as.integer(substr(b, m[1], m[1] + attr(m, "match.length")[1] - 1))
}

pick_file_for_date <- function(date_int, fixed = TRUE, nearest_ok = TRUE){
  files <- list_adj_files(fixed)
  if (!length(files)) return(NA_character_)
  d <- vapply(files, extract_date_from_name, integer(1))
  exact <- which(d == date_int)
  if (length(exact)) return(files[exact[1]])
  if (!nearest_ok) return(NA_character_)
  ok <- is.finite(d)
  if (!any(ok)) return(NA_character_)
  ii <- which.min(abs(d[ok] - date_int))
  files[ which(ok)[ii] ]
}

# ---------- crash-day finders ----------
# a) Largest daily drop in FRM@Banci (uses CSV from your pipeline)
find_crash_day_frm <- function(){
  frm_csv <- file.path(output_path, "Lambda", paste0("FRM_", channel, "_index.csv"))
  if (!file.exists(frm_csv)) stop("FRM index CSV not found: ", frm_csv)
  X <- read.csv(frm_csv, stringsAsFactors = FALSE, check.names = FALSE)
  names(X) <- tolower(names(X))
  stopifnot(all(c("date","frm") %in% names(X)))
  d <- suppressWarnings(as.Date(X$date))
  if (any(is.na(d))) d <- suppressWarnings(as.Date(X$date, "%Y-%m-%d"))
  X$date <- d
  X <- X[is.finite(X$frm) & !is.na(X$date), ]
  if (nrow(X) < 2) stop("Not enough FRM observations.")
  ch <- diff(X$frm)                  # today - yesterday
  i  <- which.min(ch)                # index in diff
  list(date = X$date[i + 1], yyyymmdd = as.integer(format(X$date[i + 1], "%Y%m%d")), drop = ch[i])
}

# b) Largest daily drop for a given bank (requires stage_20)
find_crash_day_bank <- function(bank = NULL){
  st20 <- file.path(output_path, "stage_20_loaded.RData")
  if (!file.exists(st20)) stop("stage_20_loaded.RData not found at ", st20)
  load(st20)  # brings: dates, stock_return, etc.
  if (is.null(bank)) bank <- if (exists("stock_main")) stock_main else colnames(stock_return)[1]
  stopifnot(bank %in% colnames(stock_return))
  r <- stock_return[, bank]
  i <- which.min(r)
  d <- dates[-1][i]
  list(date = d, yyyymmdd = as.integer(format(d, "%Y%m%d")), drop = r[i], bank = bank)
}

# ---------- user-facing wrappers ----------
# 1) Plot by exact date
plot_adj_for_date <- function(date_input, fixed = TRUE, nearest_ok = TRUE,
                              out_dir = file.path(output_path, "Network")){
  date_int <- to_int_yyyymmdd(date_input)
  f <- pick_file_for_date(date_int, fixed = fixed, nearest_ok = nearest_ok)
  if (is.na(f)) stop("No adjacency CSV found (fixed=", fixed, ") for/near date ", date_int)
  message("Using ", basename(f), " (fixed=", fixed, ")")
  plot_adj_matrix(f, out_dir = out_dir)
}

# 2) Auto: plot the FRM crash day (largest daily FRM drop)
plot_adj_for_crash_day <- function(fixed = TRUE, nearest_ok = TRUE,
                                   out_dir = file.path(output_path, "Network")){
  info <- find_crash_day_frm()
  message(sprintf("Crash day by FRM drop: %s (ΔFRM = %.6f)", format(info$date, "%Y-%m-%d"), info$drop))
  plot_adj_for_date(info$yyyymmdd, fixed = fixed, nearest_ok = nearest_ok, out_dir = out_dir)
}

# 3) Auto: plot the crash day for a specific bank (largest daily return drop)
plot_adj_for_bank_crash <- function(bank = NULL, fixed = TRUE, nearest_ok = TRUE,
                                    out_dir = file.path(output_path, "Network")){
  info <- find_crash_day_bank(bank)
  message(sprintf("Crash day for %s: %s (return = %.6f)",
                  info$bank, format(info$date, "%Y-%m-%d"), info$drop))
  plot_adj_for_date(info$yyyymmdd, fixed = fixed, nearest_ok = nearest_ok, out_dir = out_dir)
}

# ----------------------- examples -----------------------

# 1. Criza post-dot-com / instabilitate Romania
plot_adj_for_date("2002-07-06", fixed=TRUE)

# 2. Criza financiara globala
plot_adj_for_date("2009-03-04", fixed=TRUE)

# 3. Criza datorii suverane europene
plot_adj_for_date("2011-01-26", fixed=TRUE)

# 4. COVID — prima zi de panica
plot_adj_for_date("2020-02-21", fixed=TRUE)

# 5. COVID — crash maxim (ai deja aceasta dar o repeti pentru consistenta)
plot_adj_for_date("2020-03-16", fixed=TRUE)

# 6. Inflatie / razboi Ucraina
plot_adj_for_date("2022-07-27", fixed=TRUE)

# 7. Criza Silicon Valley Bank
plot_adj_for_date("2023-03-14", fixed=TRUE)
#
# 2) Auto FRM crash-day (Fixed):
# plot_adj_for_crash_day(fixed = TRUE)
#
# 3) Auto bank crash-day (e.g., BT) using stage_20 returns (Fixed):
# plot_adj_for_bank_crash("BT", fixed = TRUE)

# ===============================================================
# Build FIXED adjacency matrices directly from CSV inputs
# - Uses FRM_SC_config.R (paths, dates, s/tau/I/J, channel)
# - Ensures stage_20_loaded.RData exists (runs loader once if missing)
# - LOCKS Top-J banks at the START of the fixed window (by mktcap)
# - Writes per-day adjacency CSVs into: Output/<channel>/Adj_Matrices/Fixed/
# - Writes helper lambdas file:        Output/<channel>/Lambda/Fixed/lambdas_fixed_<start>_<end>.csv
# ===============================================================

rm(list = setdiff(ls(), character(0)))
suppressPackageStartupMessages({ library(zoo) })

# ---- config + utils + FRM algorithm ----
stopifnot(file.exists("FRM_SC_config.R"))
source("FRM_SC_config.R")
if (file.exists("FRM_SC_utils.R")) source("FRM_SC_utils.R")
if (!file.exists("FRM_Statistics_Algorithm.R"))
  stop("Missing FRM_Statistics_Algorithm.R")
source("FRM_Statistics_Algorithm.R")

# ---- ensure stage_20 exists (or build once) ----
st20 <- file.path(output_path, "stage_20_loaded.RData")
if (!file.exists(st20)) {
  message("stage_20_loaded.RData not found — running FRM_SC_load_data.R once...")
  stopifnot(file.exists("FRM_SC_load_data.R"))
  source("FRM_SC_load_data.R")
}
load(st20)
# expected objects now:
# dates, ticker, stock_return, macro_return, M_stock, M_macro, mktcap_index, s, tau, I, J

# ---- output folders ----
fixed_dir  <- file.path(output_path, "Adj_Matrices", "Fixed")
lambda_dir <- file.path(output_path, "Lambda", "Fixed")
dir.create(fixed_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(lambda_dir, recursive = TRUE, showWarnings = FALSE)

# ---- helpers for date indexing ----
idx_start_at_or_after <- function(key, target){ w <- which(key >= target); if (length(w)) w[1] else NA_integer_ }
idx_end_at_or_before  <- function(key, target){ w <- which(key <= target); if (length(w)) tail(w,1) else NA_integer_ }

# ---- locate fixed window ----
N0_fixed <- idx_start_at_or_after(ticker, date_start_fixed)
N1_fixed <- idx_end_at_or_before( ticker, date_end_fixed)
if (is.na(N0_fixed) || is.na(N1_fixed)) {
  stop("Fixed window not found in ticker range. Check date_start_fixed/date_end_fixed in FRM_SC_config.R")
}
if ((N1_fixed - N0_fixed + 1) < s) {
  stop(sprintf("Fixed window [%s, %s] too short for s=%d.", date_start_fixed, date_end_fixed, s))
}

# ---- lock Top-J once at the START of the fixed window ----
J_eff <- min(J, M_stock)
if (!is.matrix(mktcap_index) || ncol(mktcap_index) < (J_eff + 1))
  stop("mktcap_index shape unexpected; rebuild stage_20.")
biggest_index_fixed <- as.integer(mktcap_index[N0_fixed, 2:(J_eff + 1), drop = FALSE])
if (any(!is.finite(biggest_index_fixed))) stop("Invalid Top-J selection at fixed start.")

bank_names  <- colnames(stock_return)
macro_names <- if (M_macro > 0) colnames(macro_return) else character(0)
col_full    <- c(bank_names[biggest_index_fixed], macro_names)
M_J         <- length(col_full)

message("=== FIXED universe (Top-", J_eff, " @ start) ===")
message("Banks: ", paste(bank_names[biggest_index_fixed], collapse = ", "))
if (length(macro_names)) message("Macros: ", paste(macro_names, collapse = ", "))

# ---- iterate days in fixed window and create CSV adjacencies ----
N_fixed <- N1_fixed - N0_fixed + 1
lambdas_fixed <- matrix(NA_real_, N_fixed, J_eff + 1)  # first col = date (yyyymmdd)
lambdas_fixed[, 1] <- ticker[N0_fixed:N1_fixed]

n_written <- 0L
for (t in N0_fixed:N1_fixed) {
  rows_ret <- (t - s):(t - 1)
  if (min(rows_ret) < 1) next
  
  X <- cbind(
    stock_return[rows_ret, biggest_index_fixed, drop = FALSE],
    if (M_macro > 0) macro_return[rows_ret, , drop = FALSE] else NULL
  )
  mode(X) <- "numeric"; X[!is.finite(X)] <- 0
  
  # Drop only columns that are entirely zero — keep columns with low but nonzero variance
  keep <- colSums(X != 0) > 0
  X <- X[, keep, drop = FALSE]
  if (ncol(X) == 0) next
  
  A <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  lam_row <- rep(NA_real_, J_eff)
  
  for (k in 1:ncol(X)) {
    est <- tryCatch(FRM_Quantile_Regression(as.matrix(X), k, tau, I), error = function(e) NULL)
    if (is.null(est)) next
    gacv   <- est$Cgacv
    best   <- if (!any(is.finite(gacv))) length(gacv) else which.min(gacv[is.finite(gacv)])
    est_lambda <- abs(data.matrix(est$lambda[best]))
    est_beta   <- t(as.matrix(est$beta[best, ]))
    A[k, -k] <- est_beta
    if (k <= J_eff) lam_row[k] <- est_lambda
  }
  
  colnames(A) <- colnames(X)
  rownames(A) <- colnames(X)
  
  out_csv <- file.path(fixed_dir, paste0("adj_matrix_", ticker[t], ".csv"))
  write.csv(A, out_csv, quote = FALSE)
  n_written <- n_written + 1L
  
  lambdas_fixed[t - N0_fixed + 1, 2:(J_eff + 1)] <- lam_row
}

# ---- save helper lambdas CSV for the fixed universe ----
colnames(lambdas_fixed) <- c("date", bank_names[biggest_index_fixed])
lam_csv <- file.path(lambda_dir, paste0("lambdas_fixed_", date_start_fixed, "_", date_end_fixed, ".csv"))
write.csv(lambdas_fixed, lam_csv, row.names = FALSE, quote = FALSE)

message("✅ Fixed adjacency CSVs written: ", n_written)
message("📁 Fixed folder: ", normalizePath(fixed_dir))
message("🧾 Lambdas CSV:  ", normalizePath(lam_csv))

# ================================================================
# Centrality indicators (CSV-only) for FIXED universe snapshots
# NO NORMALIZATION in overlays: FRM raw (left), centralities raw (right).
# - Input : Output/<channel>/Adj_Matrices/Fixed/adj_matrix_YYYYMMDD.csv
# - FRM   : Output/<channel>/Lambda/FRM_<channel>_index.csv
# - Output: Output/<channel>/Network/Fixed/*.csv
#           Website/<channel>/<date_end>/Fixed/*.png + *.pdf (optional)
# ================================================================

suppressPackageStartupMessages({
  library(igraph); library(dplyr); library(ggplot2); library(readr); library(scales)
})

source("FRM_SC_config.R")
source("FRM_SC_utils.R")

# ---- Paths ----
adj_dir <- file.path(output_path, "Adj_Matrices", "Fixed")
net_dir <- file.path(output_path, "Network", "Fixed")
dir.create(net_dir, recursive = TRUE, showWarnings = FALSE)

web_dir <- file.path(website_path, date_end, "Fixed")
dir.create(web_dir, recursive = TRUE, showWarnings = FALSE)

# ---- FRM index (CSV) ----
frm_csv <- file.path(output_path, "Lambda", paste0("FRM_", channel, "_index.csv"))
if (!file.exists(frm_csv)) stop("Missing FRM index CSV at: ", frm_csv)
FRM_index <- read.csv(frm_csv, stringsAsFactors = FALSE)
FRM_index$date <- suppressWarnings(as.Date(FRM_index$date))
FRM_index <- FRM_index[is.finite(FRM_index$frm) & !is.na(FRM_index$date), ]

# ---- Find Fixed adjacency CSVs ----
if (!dir.exists(adj_dir)) stop("Missing folder: ", normalizePath(adj_dir))
files <- list.files(adj_dir, pattern = "^[Aa]dj_?matrix_[0-9]{8}\\.csv$", full.names = TRUE)
if (!length(files)) {
  files <- list.files(adj_dir, pattern = "\\.csv$", full.names = TRUE)
  files <- files[grepl("adj.*matr?ix", basename(files), ignore.case = TRUE)]
}
if (!length(files)) stop("No Fixed adjacency CSVs found in: ", normalizePath(adj_dir))

# ---- Helpers ----
date_from_name <- function(fp){
  b <- basename(fp)
  m <- regexpr("[0-9]{8}", b, perl = TRUE)
  if (m[1] == -1) return(as.Date(NA))
  ds <- substr(b, m[1], m[1] + attr(m, "match.length")[1] - 1)
  suppressWarnings(as.Date(ds, "%Y%m%d"))
}
numify_matrix <- function(df){
  for (j in seq_along(df)) {
    if (!is.numeric(df[[j]])) {
      v <- df[[j]]; if (is.factor(v)) v <- as.character(v)
      v <- gsub("\\s+", "", v)
      df[[j]] <- suppressWarnings(as.numeric(v))
    }
  }
  as.matrix(df)
}
read_adj_csv <- function(fp){
  df <- tryCatch(read.csv(fp, row.names = 1, check.names = FALSE), error = function(e) NULL)
  if (!is.null(df)) {
    M <- numify_matrix(df); mode(M) <- "numeric"; M[!is.finite(M)] <- 0
    if (any(M != 0)) return(M)
  }
  df <- tryCatch(read.csv2(fp, row.names = 1, check.names = FALSE), error = function(e) NULL)
  if (!is.null(df)) {
    M <- numify_matrix(df); mode(M) <- "numeric"; M[!is.finite(M)] <- 0
    if (any(M != 0)) return(M)
  }
  stop("Cannot parse numeric adjacency: ", basename(fp))
}
inv_distance <- function(g) 1 / pmax(abs(E(g)$weight), 1e-8)
abs_weight  <- function(g) abs(E(g)$weight)
safe_mean   <- function(x) if (length(x)) mean(x, na.rm = TRUE) else NA_real_

# ---- Order files by date ----
map <- data.frame(file = files, stringsAsFactors = FALSE)
map$date <- vapply(map$file, date_from_name, as.Date(NA))
map <- map[!is.na(map$date), ]
map <- map[order(map$date), ]
if (!nrow(map)) stop("No parsable dates in Fixed adjacency filenames.")

# ---- Compute centralities per snapshot ----
by_node <- list()
avg_rows <- vector("list", nrow(map))

for (i in seq_len(nrow(map))) {
  A <- tryCatch(read_adj_csv(map$file[i]), error = function(e) NULL)
  if (is.null(A)) next
  mode(A) <- "numeric"
  A[!is.finite(A)] <- 0
  diag(A) <- 0
  if (ncol(A) < 2) next
  
  g <- igraph::graph_from_adjacency_matrix(A, mode = "directed", weighted = TRUE, diag = FALSE)
  
  if (ecount(g) == 0) {
    n <- vcount(g)
    outdeg <- indeg <- closeo <- betw <- eigv <- ostr <- istr <- rep(0, n)
  } else {
    w_d   <- inv_distance(g)            # distances ~ 1/|weight|
    w_abs <- abs_weight(g)              # strengths/eigen on |weight|
    
    outdeg <- tryCatch(degree(g, mode = "out", normalized = TRUE), error = function(e) rep(0, vcount(g)))
    indeg  <- tryCatch(degree(g, mode = "in",  normalized = TRUE), error = function(e) rep(0, vcount(g)))
    closeo <- tryCatch(closeness(g, mode = "out", normalized = TRUE, weights = w_d),
                       error = function(e) tryCatch(closeness(g, mode = "out", normalized = TRUE),
                                                    error = function(e2) rep(0, vcount(g))))
    betw   <- tryCatch(betweenness(g, directed = TRUE, normalized = TRUE, weights = w_d),
                       error = function(e) tryCatch(betweenness(g, directed = TRUE, normalized = TRUE),
                                                    error = function(e2) rep(0, vcount(g))))
    eigv   <- tryCatch(eigen_centrality(g, directed = TRUE, weights = w_abs)$vector,
                       error = function(e) rep(0, vcount(g)))
    ostr   <- tryCatch(strength(g, mode = "out", weights = w_abs), error = function(e) rep(0, vcount(g)))
    istr   <- tryCatch(strength(g, mode = "in",  weights = w_abs), error = function(e) rep(0, vcount(g)))
    
    fixv <- function(v){ v[!is.finite(v)] <- 0; as.numeric(v) }
    outdeg <- fixv(outdeg); indeg <- fixv(indeg); closeo <- fixv(closeo)
    betw   <- fixv(betw);   eigv  <- fixv(eigv);  ostr   <- fixv(ostr);  istr <- fixv(istr)
  }
  
  nodes <- igraph::V(g)$name
  by_node[[length(by_node) + 1L]] <- data.frame(
    date = map$date[i], node = nodes,
    outdegree = outdeg, indegree = indeg,
    closeness = closeo, betweenness = betw,
    eigenvector = eigv, out_strength = ostr, in_strength = istr,
    stringsAsFactors = FALSE
  )
  
  avg_rows[[i]] <- data.frame(
    date = map$date[i],
    outdegree_avg    = safe_mean(outdeg),
    indegree_avg     = safe_mean(indeg),
    closeness_avg    = safe_mean(closeo),
    betweenness_avg  = safe_mean(betw),
    eigenvector_avg  = safe_mean(eigv),
    out_strength_avg = safe_mean(ostr),
    in_strength_avg  = safe_mean(istr)
  )
}

centrality_by_node <- dplyr::bind_rows(by_node)
centrality_avg     <- dplyr::bind_rows(avg_rows) %>% dplyr::arrange(date)

# ---- Save CSV outputs ----
write.csv(centrality_by_node, file.path(net_dir, "Centrality_ByNode_Fixed.csv"), row.names = FALSE)
write.csv(centrality_avg,     file.path(net_dir, "Centrality_Averages_Fixed.csv"), row.names = FALSE)

# --- Fix date types before joining ---
to_Date <- function(x){
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(x, origin = "1970-01-01"))
  d <- suppressWarnings(as.Date(x))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, "%Y-%m-%d"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, "%m/%d/%Y"))
  d
}

FRM_index$date          <- to_Date(FRM_index$date)
centrality_avg$date     <- to_Date(centrality_avg$date)
centrality_by_node$date <- to_Date(centrality_by_node$date)  # if you use it later

# ---- Correlation with FRM (time series averages vs FRM) ----
aligned <- dplyr::full_join(
  FRM_index %>% dplyr::arrange(date),
  centrality_avg %>% dplyr::arrange(date),
  by = "date"
)

vars <- c("frm","outdegree_avg","indegree_avg","closeness_avg","betweenness_avg",
          "eigenvector_avg","out_strength_avg","in_strength_avg")
C <- aligned[, vars, drop = FALSE]
C <- C[stats::complete.cases(C), , drop = FALSE]

corr_mat <- if (nrow(C) > 3) stats::cor(C, use = "pairwise.complete.obs", method = "pearson") else NA
p_mat <- matrix(NA_real_, ncol = ncol(C), nrow = ncol(C),
                dimnames = list(colnames(C), colnames(C)))
if (nrow(C) > 3) {
  for (i in 1:(ncol(C) - 1)) for (j in (i + 1):ncol(C)) {
    ct <- tryCatch(stats::cor.test(C[[i]], C[[j]], method = "pearson"), error = function(e) NULL)
    if (!is.null(ct)) { p_mat[i, j] <- p_mat[j, i] <- ct$p.value }
  }
  diag(p_mat) <- 0
}
write.csv(corr_mat, file.path(net_dir, "Centrality_FRM_Corr_Fixed.csv"),  row.names = TRUE)
write.csv(p_mat,    file.path(net_dir, "Centrality_FRM_CorrP_Fixed.csv"), row.names = TRUE)

# ===================== OVERLAYS (RAW AXES) =====================
# Left axis: FRM raw
# Right axis: centrality raw
# We compute an affine map so the red line can be drawn on the left axis,
# but the right axis shows the *raw* centrality values (no normalization).

PNG_DPI <- 120
save_png <- function(plot_obj, base_name, out_dir = web_dir) {
  ggsave(
    filename = file.path(out_dir, paste0(base_name, ".png")),
    plot = plot_obj, device = "png",
    width = 1600/PNG_DPI, height = 900/PNG_DPI, dpi = PNG_DPI, bg = "white"
  )
}

plot_overlay <- function(df, colname, right_label,
                         frm_col = "blue",    # blue
                         cent_col = "red") { # red
  D <- df[, c("date","frm", colname), drop = FALSE]
  names(D)[3] <- "cent"
  D <- D[stats::complete.cases(D), , drop = FALSE]
  if (!nrow(D)) return(invisible(NULL))
  
  # Left axis range = FRM RAW
  frm_min <- min(D$frm, na.rm = TRUE); frm_max <- max(D$frm, na.rm = TRUE)
  if (!is.finite(frm_min) || !is.finite(frm_max) || frm_max <= frm_min) {
    frm_min <- 0; frm_max <- 1
  }
  
  # Right axis range = CENTRALITY RAW
  cen_min <- min(D$cent, na.rm = TRUE); cen_max <- max(D$cent, na.rm = TRUE)
  if (!is.finite(cen_min) || !is.finite(cen_max) || cen_max <= cen_min) {
    cen_min <- 0; cen_max <- 1
  }
  
  # Affine map to *draw* red series on left axis (axes themselves remain raw)
  to_left   <- function(x) (x - cen_min) * (frm_max - frm_min) / (cen_max - cen_min) + frm_min
  from_left <- function(y) (y - frm_min) * (cen_max - cen_min) / (frm_max - frm_min) + cen_min
  
  # constant legend label for the red series
  D$series_cent <- right_label
  
  p <- ggplot(D, aes(x = date)) +
    geom_line(aes(y = frm,           color = "FRM@Banci"), linewidth = 1.3, lineend = "round") +
    geom_line(aes(y = to_left(cent), color = series_cent),  linewidth = 1.3, lineend = "round") +
    scale_color_manual(
      values = c("FRM@Banci" = frm_col, setNames(cent_col, right_label)),
      name = NULL
    ) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(
      limits = c(frm_min, frm_max),
      name   = paste0("FRM@", channel, " (raw)"),
      sec.axis = sec_axis(
        trans = ~ from_left(.),
        name  = right_label,
        breaks = scales::pretty_breaks(n = 6)
      )
    ) +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  base <- paste0("FRM_vs_", gsub("_avg$","", colname), "_Fixed_DualAxis")
  save_png(p, base)
}

plot_overlay(aligned, "outdegree_avg",    "Out-degree")
plot_overlay(aligned, "indegree_avg",     "In-degree")
plot_overlay(aligned, "closeness_avg",    "Closeness")
plot_overlay(aligned, "betweenness_avg",  "Betweenness")
plot_overlay(aligned, "eigenvector_avg",  "Eigenvector")
plot_overlay(aligned, "out_strength_avg", "Out-strength")
plot_overlay(aligned, "in_strength_avg",  "In-strength")

message("✅ Centrality CSVs: ", normalizePath(net_dir))
message("✅ Overlay PNGs:    ", normalizePath(web_dir))

# =============== Correlation table (HTML — FRM, closeness, eigen) ===============
# Uses: FRM_<channel>_index.csv and Network/Fixed/Centrality_Averages_Fixed.csv
# Output: Website/<channel>/<date_end>/Fixed/Correlation_FRM_Closeness_Eigen.html

# Paths (reuse previously defined web_dir, net_dir, output_path, channel)
frm_csv  <- file.path(output_path, "Lambda", paste0("FRM_", channel, "_index.csv"))
cent_csv <- file.path(net_dir, "Centrality_Averages_Fixed.csv")
stopifnot(file.exists(frm_csv), file.exists(cent_csv))

FRM_index <- read.csv(frm_csv, stringsAsFactors = FALSE, check.names = FALSE)
CENT      <- read.csv(cent_csv, stringsAsFactors = FALSE, check.names = FALSE)

names(FRM_index) <- tolower(names(FRM_index))
names(CENT)      <- tolower(names(CENT))

# robust date coercer (handles Date, POSIXct, "YYYY-mm-dd", "mm/dd/YYYY", and numeric days since 1970)
to_Date <- function(x){
  if (inherits(x, "Date"))   return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x))         return(as.Date(round(x), origin = "1970-01-01"))
  d <- suppressWarnings(as.Date(x))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, "%Y-%m-%d"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, "%m/%d/%Y"))
  d
}

FRM_index$date <- to_Date(FRM_index$date)
CENT$date      <- to_Date(CENT$date)

dat <- merge(
  FRM_index[, c("date","frm")],
  CENT[, c("date","closeness_avg","eigenvector_avg")],
  by = "date", all = FALSE
)
dat <- dat[is.finite(dat$frm) & is.finite(dat$closeness_avg) & is.finite(dat$eigenvector_avg), ]
if (!nrow(dat)) stop("No overlapping rows to compute correlations.")

# helper: r and p, then format with stars
star <- function(p){ if (is.na(p)) "" else if (p < 0.01) "***" else if (p < 0.05) "**" else if (p < 0.10) "*" else "" }
rp   <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(c(r = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  c(r = unname(ct$estimate), p = ct$p.value)
}
fmt <- function(r, p){ if (is.na(r)) "" else paste0(formatC(r, format = "f", digits = 3), star(p)) }

# compute pairs
r_FRM_clo <- rp(dat$frm, dat$closeness_avg)
r_FRM_eig <- rp(dat$frm, dat$eigenvector_avg)
r_clo_eig <- rp(dat$closeness_avg, dat$eigenvector_avg)

# build 3×3 text matrix (blank diagonal)
cell <- matrix("", 3, 3, dimnames = list(c("FRM","closeness","eigen"), c("FRM","closeness","eigen")))
cell["FRM","closeness"]      <- fmt(r_FRM_clo["r"], r_FRM_clo["p"])
cell["closeness","FRM"]      <- cell["FRM","closeness"]
cell["FRM","eigen"]          <- fmt(r_FRM_eig["r"], r_FRM_eig["p"])
cell["eigen","FRM"]          <- cell["FRM","eigen"]
cell["closeness","eigen"]    <- fmt(r_clo_eig["r"], r_clo_eig["p"])
cell["eigen","closeness"]    <- cell["closeness","eigen"]

# HTML
header <- "<tr><th></th><th>FRM</th><th>closeness</th><th>eigen</th></tr>"
row1   <- paste0("<tr><th>FRM</th><td></td><td>", cell["FRM","closeness"], "</td><td>", cell["FRM","eigen"], "</td></tr>")
row2   <- paste0("<tr><th>closeness</th><td>", cell["closeness","FRM"], "</td><td></td><td>", cell["closeness","eigen"], "</td></tr>")
row3   <- paste0("<tr><th>eigen</th><td>", cell["eigen","FRM"], "</td><td>", cell["eigen","closeness"], "</td><td></td></tr>")

q <- "\""
html <- paste0(
  "<!doctype html>\n<html lang=", q, "en", q, ">\n<head>\n<meta charset=", q, "utf-8", q, ">\n",
  "<title>Correlation: FRM, closeness, eigen</title>\n",
  "<style>",
  "body{font-family:-apple-system,system-ui,Arial,sans-serif;color:#000;background:#fff;margin:24px}",
  "table{border-collapse:collapse}",
  "th,td{border:1px solid #000;padding:6px 10px;text-align:right}",
  "th:first-child{text-align:left}",
  "</style>\n</head>\n<body>\n",
  "<h1>Correlation (Pearson)</h1>",
  "<p>Entries are <em>r</em>; significance: <strong>***</strong> p &lt; 0.01, ",
  "<strong>**</strong> p &lt; 0.05, <strong>*</strong> p &lt; 0.10.</p>",
  "<table>", header, row1, row2, row3, "</table>\n</body>\n</html>\n"
)

out_html <- file.path(web_dir, "Correlation_FRM_Closeness_Eigen.html")
writeLines(html, out_html)
message("🧾 HTML written: ", normalizePath(out_html))


# ==========================================================
# HHI charts: Lambda + MarketCap (Daily & Cumulative-to-date)
# Raw values (no normalization), transparent background,
# no titles, no x-axis label "Date".
# Outputs:
#   CSVs -> Output/<channel>/Lambda/
#   PNGs -> Website/<channel>/<date_end>/
# ==========================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# Expect: FRM_SC_utils.R defines theme_transparent_bottom
source("FRM_SC_config.R")
source("FRM_SC_utils.R")

# Data from earlier stages
load(file.path(output_path,"stage_20_loaded.RData"))   # mktcap, dates, etc.
load(file.path(output_path,"stage_30_varying.RData"))  # FRM_individ (list)

dir.create(file.path(output_path,"Lambda"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(website_path, date_end), showWarnings=FALSE, recursive=TRUE)

# ---------------- Helpers ----------------
hhi_from_values <- function(x){
  x <- as.numeric(x); x[!is.finite(x) | x < 0] <- 0
  s <- sum(x)
  if (s <= 0) return(NA_real_)
  p <- x / s
  sum(p^2)
}

# ---------------- 1) HHI from lambdas ----------------
# names(FRM_individ) expected to be date strings (YYYY-MM-DD or YYYYMMDD)
lambda_dates <- suppressWarnings(as.Date(names(FRM_individ)))
if (any(is.na(lambda_dates))) {
  lambda_dates2 <- suppressWarnings(as.Date(names(FRM_individ), "%Y%m%d"))
  lambda_dates[is.na(lambda_dates)] <- lambda_dates2[is.na(lambda_dates)]
}

hhi_lambda <- vapply(
  FRM_individ,
  function(m){
    if (is.null(m) || !ncol(m)) return(NA_real_)
    # use first row of lambda vector as in your pipeline
    hhi_from_values(m[1, ])
  },
  numeric(1)
)

HHI_lambda_df <- data.frame(date = lambda_dates, HHI_lambda = round(hhi_lambda, 6))
HHI_lambda_df <- HHI_lambda_df[order(HHI_lambda_df$date), ]
write.csv(HHI_lambda_df, file.path(output_path,"Lambda","HHI_lambda.csv"), row.names=FALSE)

png(file.path(website_path, date_end, "HHI_Lambda.png"), width=1600, height=900, res=120, bg="white")
print(
  ggplot(HHI_lambda_df, aes(x=date, y=HHI_lambda, color="Lambda HHI")) +
    geom_line(linewidth=1.2, lineend="round") +
    scale_color_manual(values=c("Lambda HHI"="#007AFF"), name=NULL) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(title=NULL, x=NULL, y="HHI") +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
)
dev.off()

# Încarcă datele
load(file.path(output_path, "stage_20_loaded.RData")) 

# REPARĂ EROAREA: Dacă obiectul s-a încărcat ca 'mktcap_index', îl redenumim în 'mktcap'
if (exists("mktcap_index") && !exists("mktcap")) {
  mktcap <- mktcap_index
}

# 2. Creăm obiectul mktcap_dates folosind direct vectorul 'ticker'
# (În stage_20_loaded.RData, ticker este de obicei un vector independent)
mktcap_dates <- as.Date(as.character(ticker), format = "%Y%m%d")

# 3. Verificăm dacă obiectul a fost creat cu succes
if(exists("mktcap_dates")) {
  print("Succes! Obiectul mktcap_dates a fost creat.")
  head(mktcap_dates)
}
mktcap_dates <- as.Date(as.character(mktcap[, 1]), format = "%Y%m%d")

# Verificăm dacă a funcționat
head(mktcap_dates)
if (any(is.na(mktcap_dates))) {
  mktcap_dates2 <- suppressWarnings(as.Date(mktcap$ticker, "%Y%m%d"))
  mktcap_dates[is.na(mktcap_dates)] <- mktcap_dates2[is.na(mktcap_dates)]
}

mcap_mat <- as.matrix(mktcap[, -1, drop=FALSE]); storage.mode(mcap_mat) <- "numeric"
mcap_mat[!is.finite(mcap_mat) | mcap_mat < 0] <- 0

hhi_mktcap_daily <- apply(mcap_mat, 1, hhi_from_values)
HHI_mktcap_daily_df <- data.frame(date = mktcap_dates, HHI_mktcap_daily = round(hhi_mktcap_daily, 6))
HHI_mktcap_daily_df <- HHI_mktcap_daily_df[order(HHI_mktcap_daily_df$date), ]

# Backward-compat file name expected elsewhere
write.csv(
  data.frame(date = HHI_mktcap_daily_df$date, HHI_mktcap = HHI_mktcap_daily_df$HHI_mktcap_daily),
  file.path(output_path,"Lambda","HHI_mktcap_all.csv"),
  row.names = FALSE
)
# Explicit daily file
write.csv(HHI_mktcap_daily_df, file.path(output_path,"Lambda","HHI_mktcap_daily.csv"), row.names=FALSE)

png(file.path(website_path, date_end, "HHI_MarketCap_Daily.png"), width=1600, height=900, res=120, bg="white")
print(
  ggplot(HHI_mktcap_daily_df, aes(x=date, y=HHI_mktcap_daily, color="Market-Cap HHI (Daily)")) +
    geom_line(linewidth=1.2, lineend="round") +
    scale_color_manual(values=c("Market-Cap HHI (Daily)"="#FF3B30"), name=NULL) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(title=NULL, x=NULL, y="HHI") +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
)
dev.off()


# ---------------- 3) HHI from Market Cap — CUMULATIVE-TO-DATE ----------------
mcap_cum_mat <- apply(mcap_mat, 2, function(col){
  col[!is.finite(col) | col < 0] <- 0
  cumsum(col)
})

hhi_mktcap_cum <- apply(mcap_cum_mat, 1, hhi_from_values)
HHI_mktcap_cum_df <- data.frame(date = mktcap_dates, HHI_mktcap_cum = round(hhi_mktcap_cum, 6))
HHI_mktcap_cum_df <- HHI_mktcap_cum_df[order(HHI_mktcap_cum_df$date), ]
write.csv(HHI_mktcap_cum_df, file.path(output_path,"Lambda","HHI_mktcap_cumulative.csv"), row.names=FALSE)

png(file.path(website_path, date_end, "HHI_MarketCap_Cumulative.png"), width=1600, height=900, res=120, bg="white")
print(
  ggplot(HHI_mktcap_cum_df, aes(x=date, y=HHI_mktcap_cum, color="Market-Cap HHI (Cumulative)")) +
    geom_line(linewidth=1.2, lineend="round") +
    scale_color_manual(values=c("Market-Cap HHI (Cumulative)"="#FF9500"), name=NULL) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(title=NULL, x=NULL, y="HHI") +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
)
dev.off()

# ---------------- 4) Combined CSV (optional) ----------------
HHI_combo <- merge(HHI_mktcap_daily_df, HHI_mktcap_cum_df, by="date", all=TRUE)
write.csv(HHI_combo, file.path(output_path,"Lambda","HHI_mktcap_variants.csv"), row.names=FALSE)

# ================================================================
# HHI (MarketCap) vs FRM@Banci — RAW axes (no normalization)
# Left axis: FRM raw  — blue "#0A84FF"
# Right axis: HHI in [0,1] — red "#FF3B30"
# Saves Daily, Cumulative, and Combined PNGs + aligned CSVs
# ================================================================

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

# --- config / paths ---
if (!exists("output_path") || !exists("website_path") || !exists("date_end") || !exists("channel")) {
  stopifnot(file.exists("FRM_SC_config.R")); source("FRM_SC_config.R")
}
dir.create(file.path(output_path, "Lambda"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(website_path, date_end), showWarnings = FALSE, recursive = TRUE)

# --- theme from utils (for transparent background) ---
if (!exists("theme_transparent_bottom")) source("FRM_SC_utils.R")

# --- helpers ---
parse_date_robust <- function(x){
  if (inherits(x,"Date"))   return(x)
  if (inherits(x,"POSIXt")) return(as.Date(x))
  if (is.numeric(x))        return(as.Date(round(x), origin = "1970-01-01"))
  d <- suppressWarnings(as.Date(x))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, "%Y-%m-%d"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, "%m/%d/%Y"))
  d
}
clamp01 <- function(v) pmin(pmax(v, 0), 1)

# --- FRM index (RAW; use column `frm`) ---
frm_rds  <- file.path(output_path, "stage_50_index.RData")
frm_csv  <- file.path(output_path, "Lambda", paste0("FRM_", channel, "_index.csv"))
if (file.exists(frm_rds)) {
  load(frm_rds)  # FRM_index (date, frm)
} else {
  stopifnot(file.exists(frm_csv))
  FRM_index <- read.csv(frm_csv, stringsAsFactors = FALSE, check.names = FALSE)
}
FRM_index <- as.data.frame(FRM_index)
names(FRM_index) <- tolower(names(FRM_index))
stopifnot(all(c("date","frm") %in% names(FRM_index)))
FRM_index$date <- parse_date_robust(FRM_index[["date"]])
FRM_index <- FRM_index[is.finite(FRM_index$frm) & !is.na(FRM_index$date), ] %>% arrange(date)

# --- HHI files (daily & cumulative) ---
daily_file <- file.path(output_path, "Lambda", "HHI_mktcap_all.csv")         # columns: date, HHI_mktcap
cum_file   <- file.path(output_path, "Lambda", "HHI_mktcap_cumulative.csv")  # columns: date, HHI_mktcap_cum
stopifnot(file.exists(daily_file), file.exists(cum_file))

HHI_d <- read.csv(daily_file, stringsAsFactors = FALSE, check.names = FALSE)
HHI_c <- read.csv(cum_file,   stringsAsFactors = FALSE, check.names = FALSE)
names(HHI_d) <- tolower(names(HHI_d))
names(HHI_c) <- tolower(names(HHI_c))

# tidy daily
stopifnot(all(c("date","hhi_mktcap") %in% names(HHI_d)))
HHI_d$date      <- parse_date_robust(HHI_d[["date"]])
HHI_d$hhi_daily <- clamp01(as.numeric(HHI_d$hhi_mktcap))
HHI_d <- HHI_d[, c("date","hhi_daily")] %>% arrange(date)

# tidy cumulative (accept _cum or _cumulative)
cum_col <- if ("hhi_mktcap_cum" %in% names(HHI_c)) "hhi_mktcap_cum" else "hhi_mktcap_cumulative"
stopifnot(cum_col %in% names(HHI_c))
HHI_c$date           <- parse_date_robust(HHI_c[["date"]])
HHI_c$hhi_mktcap_cum <- clamp01(as.numeric(HHI_c[[cum_col]]))
HHI_c <- HHI_c[, c("date","hhi_mktcap_cum")] %>% arrange(date)

# --- align with FRM (inner join) & save aligned CSVs ---
frm_df   <- FRM_index[, c("date","frm")] %>% arrange(date)
aligned_d <- full_join(frm_df, HHI_d, by = "date")
aligned_c <- full_join(frm_df, HHI_c, by = "date")
write.csv(aligned_d, file.path(output_path, "Lambda", "HHI_vs_FRM_Banci_DAILY_aligned_raw.csv"), row.names = FALSE)
write.csv(aligned_c, file.path(output_path, "Lambda", "HHI_vs_FRM_Banci_CUMULATIVE_aligned_raw.csv"), row.names = FALSE)

# --- generic dual-axis plotter (FRM raw left; HHI raw [0,1] right) ---
plot_dual_raw <- function(D, hhi_col, right_label, out_png,
                          frm_col = "blue",   # blue
                          hhi_colr = "red"   # red
){
  if (!nrow(D)) return(invisible(NULL))
  stopifnot(all(c("date","frm", hhi_col) %in% names(D)))
  d <- D[, c("date","frm", hhi_col)]
  names(d)[3] <- "hhi"
  d <- d[is.finite(d$frm) & is.finite(d$hhi) & !is.na(d$date), ]
  if (!nrow(d)) return(invisible(NULL))
  
  # left axis = FRM (raw range)
  fmin <- min(d$frm, na.rm = TRUE); fmax <- max(d$frm, na.rm = TRUE)
  if (!is.finite(fmin) || !is.finite(fmax) || fmax == fmin) { fmin <- 0; fmax <- 1 }
  
  # right axis fixed [0,1]; map HHI to left for drawing
  hmin <- 0; hmax <- 1
  to_left   <- function(x) (x - hmin) * (fmax - fmin) / (hmax - hmin) + fmin
  from_left <- function(y) (y - fmin) * (hmax - hmin) / (fmax - fmin) + hmin
  
  # ✅ correct color mapping keyed to actual legend labels
  pal <- setNames(c(frm_col, hhi_colr), c("FRM@Banci", right_label))
  
  png(out_png, width = 1600, height = 900, res = 120, bg = "white"); print(
    ggplot(d, aes(x = date)) +
      geom_line(aes(y = frm,          color = "FRM@Banci"), linewidth = 1.4, lineend = "round") +
      geom_line(aes(y = to_left(hhi), color = right_label),  linewidth = 1.4, lineend = "round") +
      scale_color_manual(values = pal, breaks = names(pal), limits = names(pal), name = NULL) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      scale_y_continuous(
        limits = c(fmin, fmax),
        name   = paste0("FRM@", channel, " (raw)"),
        sec.axis = sec_axis(~ from_left(.),
                            breaks = seq(0, 1, 0.25),
                            labels = sprintf("%.2f", seq(0, 1, 0.25)),
                            name   = "HHI (0–1)")
      ) +
      theme_transparent_bottom +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ); dev.off()
}

out_dir <- file.path(website_path, date_end)

# Daily (FRM blue, HHI red)
plot_dual_raw(
  aligned_d, "hhi_daily",
  right_label = "Market-Cap HHI (Daily)",
  out_png = file.path(out_dir, "HHI_vs_FRM_Banci_DAILY_Raw_DualAxis.png")
)

# Cumulative (FRM blue, HHI red)
plot_dual_raw(
  aligned_c, "hhi_mktcap_cum",
  right_label = "Market-Cap HHI (Cumulative)",
  out_png = file.path(out_dir, "HHI_vs_FRM_Banci_CUMULATIVE_Raw_DualAxis.png")
)

# Combined (FRM blue, HHI daily red, HHI cumulative orange)
combo <- full_join(aligned_d, aligned_c, by = c("date","frm"))
if (nrow(combo)) {
  fmin <- min(combo$frm, na.rm = TRUE); fmax <- max(combo$frm, na.rm = TRUE)
  if (!is.finite(fmin) || !is.finite(fmax) || fmax == fmin) { fmin <- 0; fmax <- 1 }
  hmin <- 0; hmax <- 1
  to_left   <- function(x) (x - hmin) * (fmax - fmin) / (hmax - hmin) + fmin
  from_left <- function(y) (y - fmin) * (hmax - hmin) / (fmax - fmin) + hmin
  
  pal_comb <- c("FRM@Banci"   = "#0A84FF",
                "HHI (Daily)"  = "#FF3B30",
                "HHI (Cumulative)" = "#FF9500")
  
  png(file.path(out_dir, "HHI_vs_FRM_Banci_BOTH_Raw_DualAxis.png"),
      width = 1600, height = 900, res = 120, bg = "white"); print(
        ggplot(combo, aes(x = date)) +
          geom_line(aes(y = frm,                    color = "FRM@Banci"),       linewidth = 1.4, lineend = "round") +
          geom_line(aes(y = to_left(hhi_daily),     color = "HHI (Daily)"),      linewidth = 1.4, lineend = "round") +
          geom_line(aes(y = to_left(hhi_mktcap_cum),color = "HHI (Cumulative)"), linewidth = 1.4, lineend = "round") +
          scale_color_manual(values = pal_comb, breaks = names(pal_comb), limits = names(pal_comb), name = NULL) +
          scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
          scale_y_continuous(
            limits = c(fmin, fmax),
            name   = paste0("FRM@", channel, " (raw)"),
            sec.axis = sec_axis(~ from_left(.),
                                breaks = seq(0, 1, 0.25),
                                labels = sprintf("%.2f", seq(0, 1, 0.25)),
                                name   = "HHI (0–1)")
          ) +
          theme_transparent_bottom +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
      ); dev.off()
}

message("✓ HHI vs FRM PNGs at: ", normalizePath(out_dir))



# ---- dynamic portfolio (GMV daily, 90-day windows; robust x-sec stats; all-bank table) ----

source("FRM_SC_utils.R"); load(file.path(output_path,"stage_20_loaded.RData"))

suppressPackageStartupMessages({
  library(quadprog); library(ggplot2)
  library(zoo); library(dplyr)
})

# ========== SETTINGS ==========
s_opt <- 90            # rolling window for GMV weights & rolling Sharpe
rebal_every <- 1       # rebalance every trading day
ANN_DAYS <- 365        # annualization days
K_VOL <- 90            # rolling volatility window
MIN_COV_FRAC <- 0.95   # require >=95% non-NA obs in window for a bank/date
KMIN_CROSS <- 5        # need at least 5 banks to compute x-sectional stat

# ========== HELPERS ==========
ann_mu   <- function(x) mean(x,na.rm=TRUE)*ANN_DAYS
ann_vol  <- function(x) sd(x,na.rm=TRUE)*sqrt(ANN_DAYS)
sharpe   <- function(x) if (sd(x,na.rm=TRUE)>0) ann_mu(x)/ann_vol(x) else NA_real_
to_wealth <- function(r) exp(cumsum(r))
roll_apply <- function(x,k,FUN) zoo::rollapply(x, k, FUN, by=1, align="right", fill=NA, na.rm=TRUE)
roll_vol    <- function(x,k) roll_apply(x,k,sd)
roll_mean   <- function(x,k) roll_apply(x,k,mean)
roll_sd     <- function(x,k) roll_apply(x,k,sd)
roll_sharpe <- function(x,k,ann=ANN_DAYS){
  m <- roll_mean(x,k); s <- roll_sd(x,k)
  out <- (m*ann)/(s*sqrt(ann)); out[!is.finite(out)] <- NA_real_; out
}
solve_gmv_longonly <- function(S){
  p <- ncol(S); Dmat <- 2 * (S + t(S))/2; dvec <- rep(0,p)
  Amat <- cbind(rep(1,p), diag(p)); bvec <- c(1, rep(0,p)); meq <- 1
  out <- tryCatch(quadprog::solve.QP(Dmat,dvec,Amat,bvec,meq=meq), error=function(e) NULL)
  if (is.null(out)) return(rep(1/p,p))
  w <- out$solution; w[w<0] <- 0; w/sum(w)
}
DD <- function(wealth){ peak <- cummax(wealth); draw <- (wealth/peak) - 1; min(draw, na.rm=TRUE) }
# finite-count in a rolling window
roll_count_finite <- function(x, k) rollapply(x, k, function(v) sum(is.finite(v)), by=1, align="right", fill=NA)
# robust row-wise stat requiring at least kmin values
row_stat <- function(M, FUN, kmin=KMIN_CROSS){
  apply(M, 1, function(row){
    v <- row[is.finite(row)]
    if (length(v) < kmin) NA_real_ else FUN(v)
  })
}
df_baza <- read.csv("Baza_FRM_Banci.csv")

# Setam data si tickerul
dates <- as.Date(as.character(df_baza$Date), format="%Y%m%d")
# ========== DATA ==========
dates_ret <- dates
R <- as.matrix(stock_return); mode(R) <- "numeric"
P <- ncol(R)

# ========== ROLLING GMV (DAILY) ==========
w_gmv <- matrix(NA_real_, nrow=nrow(R), ncol=P)
for (t in seq_len(nrow(R))) {
  if (t < s_opt) next
  if (((t - s_opt) %% rebal_every) != 0) { w_gmv[t,] <- w_gmv[t-1,]; next }
  Rw <- R[(t - s_opt + 1):t, , drop=FALSE]
  S  <- stats::cov(Rw, use="pairwise.complete.obs"); S[!is.finite(S)] <- 0
  w_gmv[t,] <- solve_gmv_longonly(S)
}
# forward-fill any gaps
w0 <- rep(1/P, P)
for (i in seq_len(nrow(w_gmv))) {
  if (i==1 && any(!is.finite(w_gmv[i,]))) w_gmv[i,] <- w0
  if (i>1  && any(!is.finite(w_gmv[i,])))  w_gmv[i,] <- w_gmv[i-1,]
}

# Portfolio returns & wealth
rp_gmv <- rowSums(R * w_gmv)
W_gmv  <- to_wealth(rp_gmv)

# Individual Banks wealth (for plots & MaxDD per bank)
W_single <- apply(R, 2, to_wealth)

# Cross-sectional wealth summaries each date (no coverage filter needed)
W_min <- apply(W_single, 1, function(x) min(x, na.rm=TRUE))
W_med <- apply(W_single, 1, function(x) stats::median(x, na.rm=TRUE))
W_max <- apply(W_single, 1, function(x) max(x, na.rm=TRUE))
W_avg <- rowMeans(W_single, na.rm=TRUE)

# Also track equal-weight average return across Individual Banks
rp_avg <- rowMeans(R, na.rm=TRUE)

# ========== ROLLING VOLATILITY (90d) with coverage filter ==========
vol_single_all <- apply(R, 2, function(col) roll_vol(col, K_VOL))
cover_mat <- sapply(seq_len(ncol(R)), function(j) roll_count_finite(R[, j], K_VOL))
MIN_COV <- ceiling(MIN_COV_FRAC * K_VOL)
vol_single_all[ cover_mat < MIN_COV ] <- NA

vol_min <- row_stat(vol_single_all, min)
vol_med <- row_stat(vol_single_all, median)
vol_max <- row_stat(vol_single_all, max)
vol_avg <- roll_vol(rp_avg, K_VOL)
vol_gmv <- roll_vol(rp_gmv, K_VOL)

# ========== ROLLING SHARPE (90d, ANN=365) with coverage filter ==========
# Compute per-bank rolling Sharpe first
sh_single_all <- apply(R, 2, function(col) roll_sharpe(col, s_opt, ann=ANN_DAYS))
# also enforce coverage for Sharpe (use the same cover_mat)
sh_single_all[ cover_mat < MIN_COV ] <- NA

sh_min <- row_stat(sh_single_all, min)
sh_med <- row_stat(sh_single_all, median)
sh_max <- row_stat(sh_single_all, max)
sh_avg <- roll_sharpe(rp_avg, s_opt, ann=ANN_DAYS)
sh_gmv <- roll_sharpe(rp_gmv, s_opt, ann=ANN_DAYS)

# ========== COLORS & STYLES ==========
col_GMV <- "dodgerblue"
col_min <- "firebrick2"
col_med <- "forestgreen"
col_max <- "purple"
col_avg <- "darkorange"
lw_gmv <- 2.2
lw_dash <- 1.3

# ========== PLOTS ==========
# 1) Cumulative wealth
png(file.path(output_path,"CAPM_Portfolio_vs_Individual Banks_Wealth.png"), width=1400, height=900, bg="white")
print(
  ggplot() +
    geom_line(aes(x=dates_ret, y=W_min, color="Individual Banks (min)"),     linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=W_med, color="Individual Banks (median)"),  linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=W_max, color="Individual Banks (max)"),     linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=W_avg, color="Individual Banks (average)"), linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=W_gmv, color="GMV"), linewidth=lw_gmv) +
    scale_color_manual(values=c("GMV"=col_GMV,
                                "Individual Banks (min)"=col_min,
                                "Individual Banks (median)"=col_med,
                                "Individual Banks (max)"=col_max,
                                "Individual Banks (average)"=col_avg)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(title=NULL, x=NULL, y="Wealth (start = 1)", color=NULL) +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
)
dev.off()

# 2) Rolling volatility (90-day)
png(file.path(output_path,"CAPM_Portfolio_vs_Individual_Banks_RollVol.png"), width=1400, height=800, bg="white")
print(
  ggplot() +
    geom_line(aes(x=dates_ret, y=vol_min, color="Individual Banks (min)"),     linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=vol_med, color="Individual Banks (median)"),  linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=vol_max, color="Individual Banks (max)"),     linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=vol_avg, color="Individual Banks (average)"), linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=vol_gmv, color="GMV"), linewidth=lw_gmv) +
    scale_color_manual(values=c("GMV"=col_GMV,
                                "Individual Banks (min)"=col_min,
                                "Individual Banks (median)"=col_med,
                                "Individual Banks (max)"=col_max,
                                "Individual Banks (average)"=col_avg)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(title=NULL, x=NULL, y="σ (90-day)", color=NULL) +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
)
dev.off()


# ==============================================================================
# BLOCUL: CALCULUL PORTOFOLIULUI GMV ȘI REZULTATE INVESTIȚIONALE
# ==============================================================================

# ---------- 1. Alinierea Datelor ----------
# CORECȚIE: Folosim tot vectorul 'dates' (5456 zile) pentru că stock_return are 5456 rânduri.
# Nu eliminăm prima zi (dates[-1]) deoarece randamentele sunt deja calculate în CSV.
dates_ret <- dates 
R <- as.matrix(stock_return)
mode(R) <- "numeric"
P <- ncol(R)

# ---------- 2. Optimizarea Portofoliului GMV (Global Minimum Variance) ----------
# Această secțiune caută să minimizeze riscul total prin rebalansare zilnică.

# Inițializăm matricea pentru ponderile investițiilor (weights)
w_gmv <- matrix(NA_real_, nrow=nrow(R), ncol=P)

for (t in seq_len(nrow(R))) {
  # Avem nevoie de un istoric minim de 90 de zile pentru a estima corelațiile
  if (t < s_opt) next
  
  # Rebalansăm conform frecvenței setate (zilnic în cazul nostru)
  if (((t - s_opt) %% rebal_every) != 0) { 
    w_gmv[t,] <- w_gmv[t-1,]
    next 
  }
  
  # Selectăm fereastra rulantă de date pentru calculul covarianței
  Rw <- R[(t - s_opt + 1):t, , drop=FALSE]
  
  # Calculăm matricea de covarianță (cum tind acțiunile să se miște împreună)
  S <- stats::cov(Rw, use="pairwise.complete.obs")
  S[!is.finite(S)] <- 0 
  
  # Rezolvăm problema de optimizare: Găsirea ponderilor care minimizează varianța portofoliului
  # sub restricția că suma ponderilor este 1 și nu avem vânzări în lipsă (long-only).
  w_gmv[t,] <- solve_gmv_longonly(S)
}

# Completăm ponderile pentru primele 90 de zile cu o strategie egal-ponderată (1/7)
w0 <- rep(1/P, P)
for (i in seq_len(nrow(w_gmv))) {
  if (i==1 && any(!is.finite(w_gmv[i,]))) w_gmv[i,] <- w0
  if (i>1  && any(!is.finite(w_gmv[i,])))  w_gmv[i,] <- w_gmv[i-1,]
}

# ---------- 3. Calculul Performanței (Randamente și Wealth) ----------

# Calculăm randamentul zilnic al strategiei noastre GMV
rp_gmv <- rowSums(R * w_gmv)

# Transformăm randamentele în bogăție cumulată (Wealth): pornind de la 1 unitate monetară
W_gmv  <- to_wealth(rp_gmv) 

# Calculăm evoluția bogăției pentru fiecare bancă individuală
W_single <- apply(R, 2, to_wealth)

# Statistici transversale: Minimul pieței, Mediana și Maximul fiecărei zile
W_min <- apply(W_single, 1, function(x) min(x, na.rm=TRUE))
W_med <- apply(W_single, 1, function(x) stats::median(x, na.rm=TRUE))
W_max <- apply(W_single, 1, function(x) max(x, na.rm=TRUE))
W_avg <- rowMeans(W_single, na.rm=TRUE)

# Randamentul mediu al pieței (fără optimizare)
rp_avg <- rowMeans(R, na.rm=TRUE)

# ---------- 4. Analiza Riscului (Volatilitatea și Sharpe Rulant) ----------

# Volatilitatea rulantă pe 90 de zile
vol_single_all <- apply(R, 2, function(col) roll_vol(col, K_VOL))
cover_mat <- sapply(seq_len(ncol(R)), function(j) roll_count_finite(R[, j], K_VOL))
MIN_COV <- ceiling(MIN_COV_FRAC * K_VOL)
vol_single_all[ cover_mat < MIN_COV ] <- NA

vol_min <- row_stat(vol_single_all, min)
vol_med <- row_stat(vol_single_all, median)
vol_max <- row_stat(vol_single_all, max)
vol_avg <- roll_vol(rp_avg, K_VOL)
vol_gmv <- roll_vol(rp_gmv, K_VOL) # Volatilitatea strategiei noastre

# Sharpe Ratio rulant (Randamentul ajustat la risc)
sh_single_all <- apply(R, 2, function(col) roll_sharpe(col, s_opt, ann=ANN_DAYS))
sh_single_all[ cover_mat < MIN_COV ] <- NA

sh_min <- row_stat(sh_single_all, min)
sh_med <- row_stat(sh_single_all, median)
sh_max <- row_stat(sh_single_all, max)
sh_avg <- roll_sharpe(rp_avg, s_opt, ann=ANN_DAYS)
sh_gmv <- roll_sharpe(rp_gmv, s_opt, ann=ANN_DAYS)

# ---------- 5. REPREZENTĂRI GRAFICE ----------

# Culori oficiale pentru lucrare
col_GMV <- "dodgerblue"   # Strategia GMV (Albastru intens)
col_min <- "firebrick2"   # Pragul minim de performanță (Roșu)
col_med <- "forestgreen"  # Evoluția mediană (Verde)
col_max <- "purple"       # Performanța maximă teoretică (Mov)
col_avg <- "darkorange"   # Media simplă (Portocaliu)

# Grafic 1: Bogăția Cumulată (Wealth)
png(file.path(output_path,"CAPM_Portfolio_vs_Individual_Banks_Wealth.png"), width=1400, height=900, bg="white")
print(
  ggplot() +
    # Benzile de variație ale acțiunilor individuale
    geom_line(aes(x=dates_ret, y=W_min, color="Individual Banks (min)"), linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=W_med, color="Individual Banks (median)"), linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=W_max, color="Individual Banks (max)"), linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=W_avg, color="Individual Banks (average)"), linewidth=lw_dash, linetype="dashed") +
    # Strategia noastră optimizată (Linie groasă continuă)
    geom_line(aes(x=dates_ret, y=W_gmv, color="GMV"), linewidth=lw_gmv) +
    
    scale_color_manual(values=c("GMV"=col_GMV,
                                "Individual Banks (min)"=col_min,
                                "Individual Banks (median)"=col_med,
                                "Individual Banks (max)"=col_max,
                                "Individual Banks (average)"=col_avg)) +
    scale_x_date(date_breaks="1 year", date_labels="%Y") + # Afișăm pe ani pentru claritate
    labs(title="Analiza Bogăției: Portofoliul GMV vs. Acțiuni Bancare Individuale", 
         x=NULL, y="Wealth (Index Start = 1)", color="Legendă Strategii") +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
)
dev.off()

# 2) Rolling volatility (90-day)
png(file.path(output_path,"CAPM_Portfolio_vs_Individual_Banks_RollVol.png"), width=1400, height=800, bg="white")
print(
  ggplot() +
    geom_line(aes(x=dates_ret, y=vol_min, color="Individual Banks (min)"),     linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=vol_med, color="Individual Banks (median)"),  linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=vol_max, color="Individual Banks (max)"),     linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=vol_avg, color="Individual Banks (average)"), linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=vol_gmv, color="GMV"), linewidth=lw_gmv) +
    scale_color_manual(values=c("GMV"=col_GMV,
                                "Individual Banks (min)"=col_min,
                                "Individual Banks (median)"=col_med,
                                "Individual Banks (max)"=col_max,
                                "Individual Banks (average)"=col_avg)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(title=NULL, x=NULL, y="σ (90-day)", color=NULL) +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
)
dev.off()


# 3) Rolling Sharpe (90-day, ANN=365)
png(file.path(output_path,"CAPM_Portfolio_vs_Individual Banks_RollSharpe.png"), width=1400, height=800, bg="white")
print(
  ggplot() +
    geom_line(aes(x=dates_ret, y=sh_min, color="Individual Banks (min)"),     linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=sh_med, color="Individual Banks (median)"),  linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=sh_max, color="Individual Banks (max)"),     linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=sh_avg, color="Individual Banks (average)"), linewidth=lw_dash, linetype="dashed") +
    geom_line(aes(x=dates_ret, y=sh_gmv, color="GMV"), linewidth=lw_gmv) +
    scale_color_manual(values=c("GMV"=col_GMV,
                                "Individual Banks (min)"=col_min,
                                "Individual Banks (median)"=col_med,
                                "Individual Banks (max)"=col_max,
                                "Individual Banks (average)"=col_avg)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(title=NULL, x=NULL, y="Rolling Sharpe (ANN=365, 90-day)", color=NULL) +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
)
dev.off()

# ========== STATS (ANN=365) — GMV + EVERY BANK ==========
stats_tbl <- data.frame(
  Strategy  = c("GMV", colnames(R)),
  AnnReturn = c(ann_mu(rp_gmv), apply(R, 2, ann_mu)),
  AnnVol    = c(ann_vol(rp_gmv), apply(R, 2, ann_vol)),
  Sharpe    = c(sharpe(rp_gmv),  apply(R, 2, sharpe)),
  MaxDD     = c(DD(W_gmv),       apply(W_single, 2, DD))
)
num_cols <- vapply(stats_tbl, is.numeric, logical(1))
stats_tbl[num_cols] <- lapply(stats_tbl[num_cols], function(x) round(x, 4))

# ====== SUMMARY TABLE (time-averages & medians across days) ======
avg <- function(x) mean(x, na.rm=TRUE)
med <- function(x) median(x, na.rm=TRUE)

summary_tbl <- data.frame(
  Indicator = rep(c("Volatility (90d)", "Sharpe (90d)"), each = 5),
  Series    = rep(c("GMV", "Singles (min)", "Singles (median)", "Singles (max)", "Singles (average)"), times = 2),
  Mean      = c(avg(vol_gmv), avg(vol_min), avg(vol_med), avg(vol_max), avg(vol_avg),
                avg(sh_gmv),  avg(sh_min),  avg(sh_med),  avg(sh_max),  avg(sh_avg)),
  Median    = c(med(vol_gmv), med(vol_min), med(vol_med), med(vol_max), med(vol_avg),
                med(sh_gmv),  med(sh_min),  med(sh_med),  med(sh_max),  med(sh_avg))
)

# round and save
summary_tbl$Mean   <- round(summary_tbl$Mean,   6)
summary_tbl$Median <- round(summary_tbl$Median, 6)
write.csv(summary_tbl, file.path(output_path, "CAPM_Rolling_Metrics_Summary.csv"), row.names = FALSE)

# ====== BOX PLOTS (full distribution across time) ======
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

# --- Use a strong, system-like blue for GMV (you can change if you prefer) ---
col_GMV <- "#007AFF"   # iOS System Blue (fully blue fill + outline for GMV)

# If you already defined col_min / col_med / col_max / col_avg above, this will use them.
# Otherwise, give them quick defaults:
if (!exists("col_min"))  col_min  <- "#8E8E93"  # muted grey
if (!exists("col_med"))  col_med  <- "#FF9F0A"  # orange
if (!exists("col_max"))  col_max  <- "#FF3B30"  # red
if (!exists("col_avg"))  col_avg  <- "#34C759"  # green

# ---------- Common high-res PNG helper (transparent) ----------
png_hr <- function(path, width_px = 2400, height_px = 1600, res = 300) {
  png(filename = path, width = width_px, height = height_px, res = res,
      units = "px", bg = "white")
}

# ---------- Aggregate (GMV vs Singles stats) ----------
vol_df <- data.frame(
  GMV                = vol_gmv,
  `Singles (min)`    = vol_min,
  `Singles (median)` = vol_med,
  `Singles (max)`    = vol_max,
  `Singles (average)`= vol_avg
) |>
  pivot_longer(cols = everything(), names_to = "Series", values_to = "Value") |>
  filter(is.finite(Value))

sh_df <- data.frame(
  GMV                = sh_gmv,
  `Singles (min)`    = sh_min,
  `Singles (median)` = sh_med,
  `Singles (max)`    = sh_max,
  `Singles (average)`= sh_avg
) |>
  pivot_longer(cols = everything(), names_to = "Series", values_to = "Value") |>
  filter(is.finite(Value))

pal_agg <- c(
  "GMV"               = col_GMV,
  "Singles (min)"     = col_min,
  "Singles (median)"  = col_med,
  "Singles (max)"     = col_max,
  "Singles (average)" = col_avg
)

# 1) 90-day Volatility (aggregate)
png_hr(file.path(output_path, "CAPM_Boxplot_RollVol_90d.png"))
print(
  ggplot(vol_df, aes(x = Series, y = Value, color = Series, fill = Series)) +
    geom_boxplot(outlier.alpha = 0.35, linewidth = 1) +
    scale_color_manual(values = pal_agg, guide = "none") +
    scale_fill_manual(values  = pal_agg, guide = "none") +
    labs(title = NULL, x = NULL, y = "σ (90-day)") +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust = 1))
)
dev.off()

# 2) 90-day Sharpe (aggregate)
png_hr(file.path(output_path, "CAPM_Boxplot_RollSharpe_90d.png"))
print(
  ggplot(sh_df, aes(x = Series, y = Value, color = Series, fill = Series)) +
    geom_boxplot(outlier.alpha = 0.35, linewidth = 1) +
    scale_color_manual(values = pal_agg, guide = "none") +
    scale_fill_manual(values  = pal_agg, guide = "none") +
    labs(title = NULL, x = NULL, y = "Rolling Sharpe (ANN=365, 90-day)") +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust = 1))
)
dev.off()

# ====== BOX PLOTS PER BANK (+ GMV) FOR 90d VOL & 90d SHARPE ======

# --- Ensure bank names from your return matrix R ---
bank_names <- colnames(R)

# Long data for per-bank volatility + GMV
vol_df_bank <- as.data.frame(vol_single_all)
colnames(vol_df_bank) <- bank_names
vol_df_bank$GMV <- vol_gmv
vol_long <- vol_df_bank |>
  pivot_longer(cols = everything(), names_to = "Series", values_to = "Value") |>
  filter(is.finite(Value))

# Long data for per-bank Sharpe + GMV
sh_df_bank <- as.data.frame(sh_single_all)
colnames(sh_df_bank) <- bank_names
sh_df_bank$GMV <- sh_gmv
sh_long <- sh_df_bank |>
  pivot_longer(cols = everything(), names_to = "Series", values_to = "Value") |>
  filter(is.finite(Value))

# Order with GMV first, then banks alphabetically
lvl <- c("GMV", sort(bank_names))
vol_long$Series <- factor(vol_long$Series, levels = lvl)
sh_long$Series  <- factor(sh_long$Series,  levels = lvl)

# ---------- Distinct colors for individual banks (GMV stays full blue) ----------
# A quick, bright, distinct palette generator
distinct_cols <- function(n, exclude = NULL) {
  # HCL-based distinct hues
  base <- grDevices::hcl(
    h = seq(15, 375, length.out = n + 1)[1:n],
    c = 100, l = 55
  )
  # Make sure none equals the excluded blue (rare anyway)
  if (!is.null(exclude)) base[base == exclude] <- "#5856D6"  # swap if collision (purple)
  base
}

bank_cols <- distinct_cols(length(bank_names), exclude = col_GMV)
names(bank_cols) <- sort(bank_names)

pal_by_bank <- c("GMV" = col_GMV, bank_cols)

# ---------- Boxplot: 90d Rolling Volatility (per bank + GMV) ----------
png_hr(file.path(output_path, "CAPM_Boxplot_RollVol_90d_ByBank.png"),
       width_px = 2800, height_px = 1600, res = 300)
print(
  ggplot(vol_long, aes(x = Series, y = Value, color = Series, fill = Series)) +
    geom_boxplot(outlier.alpha = 0.30, linewidth = 0.9) +
    scale_color_manual(values = pal_by_bank, guide = "none") +
    scale_fill_manual(values  = pal_by_bank, guide = "none") +
    labs(title = NULL, x = NULL, y = "σ (90-day)") +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
)
dev.off()

# ---------- Boxplot: 90d Rolling Sharpe (per bank + GMV) ----------
png_hr(file.path(output_path, "CAPM_Boxplot_RollSharpe_90d_ByBank.png"),
       width_px = 2800, height_px = 1600, res = 300)
print(
  ggplot(sh_long, aes(x = Series, y = Value, color = Series, fill = Series)) +
    geom_boxplot(outlier.alpha = 0.30, linewidth = 0.9) +
    scale_color_manual(values = pal_by_bank, guide = "none") +
    scale_fill_manual(values  = pal_by_bank, guide = "none") +
    labs(title = NULL, x = NULL, y = "Rolling Sharpe (ANN=365, 90-day)") +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
)
dev.off()

# ---------- Per-bank summary tables (time-mean & time-median) ----------
avg <- function(x) mean(x, na.rm = TRUE)
med <- function(x) median(x, na.rm = TRUE)

# Volatility summary
vol_summary <- data.frame(
  Series = lvl,
  Mean   = c(avg(vol_gmv), sapply(vol_df_bank[ , sort(bank_names), drop = FALSE], avg)),
  Median = c(med(vol_gmv), sapply(vol_df_bank[ , sort(bank_names), drop = FALSE], med))
)
vol_summary$Mean   <- round(vol_summary$Mean,   6)
vol_summary$Median <- round(vol_summary$Median, 6)
write.csv(vol_summary, file.path(output_path, "CAPM_RollVol_90d_ByBank_Summary.csv"), row.names = FALSE)

# Sharpe summary
sh_summary <- data.frame(
  Series = lvl,
  Mean   = c(avg(sh_gmv), sapply(sh_df_bank[ , sort(bank_names), drop = FALSE], avg)),
  Median = c(med(sh_gmv), sapply(sh_df_bank[ , sort(bank_names), drop = FALSE], med))
)
sh_summary$Mean   <- round(sh_summary$Mean,   6)
sh_summary$Median <- round(sh_summary$Median, 6)
write.csv(sh_summary, file.path(output_path, "CAPM_RollSharpe_90d_ByBank_Summary.csv"), row.names = FALSE)

# If you computed a stats_tbl elsewhere, this persists it:
if (exists("stats_tbl")) {
  write.csv(stats_tbl,
            file.path(output_path, "CAPM_Portfolio_vs_Individual Banks_Stats.csv"),
            row.names = FALSE)
}

# ===== GMV weights over time (exports + plots) =====

# ==== GMV Weights — All Coins, Stacked Area Over Time =========================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

# Fallback theme if theme_transparent_bottom isnt loaded
if (!exists("theme_transparent_bottom")) {
  theme_transparent_bottom <- theme_minimal(base_size = 12) +
    theme(
      panel.background      = element_rect(fill = "white", colour = NA),
      plot.background       = element_rect(fill = "white", colour = NA),
      legend.background     = element_rect(fill = "white"),
      legend.box.background = element_rect(fill = "white"),
      legend.position       = "bottom",
      axis.text.x           = element_text(angle = 30, vjust = 1, hjust = 1)
    )
}

# --- Ensure CAPS labels for banks ---
bank_names <- toupper(colnames(R))
colnames(w_gmv) <- bank_names

# --- Wide weights (Date first, then sorted banks) ---
weights_wide <- as.data.frame(w_gmv)
weights_wide$Date <- as.Date(dates_ret)
weights_wide <- weights_wide[, c("Date", sort(bank_names))]

# Export CSVs
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
write.csv(weights_wide, file.path(output_path, "GMV_Weights_Daily_wide.csv"), row.names = FALSE)

weights_long <- weights_wide |>
  tidyr::pivot_longer(-Date, names_to = "Bank", values_to = "Weight") |>
  dplyr::arrange(Date, Bank)

write.csv(weights_long, file.path(output_path, "GMV_Weights_Daily_long.csv"), row.names = FALSE)

# --- Color palette (distinct hues, stable mapping) ---
distinct_cols <- function(n){
  grDevices::hcl(h = seq(15, 375, length.out = n + 1)[1:n], c = 100, l = 55)
}
all_banks <- sort(unique(weights_long$Bank))
pal_banks <- distinct_cols(length(all_banks))
names(pal_banks) <- all_banks

# --- Order layers by average weight (largest first) ---
bank_order <- weights_long |>
  dplyr::group_by(Bank) |>
  dplyr::summarise(Avg = mean(Weight, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(desc(Avg)) |>
  dplyr::pull(Bank)

weights_all_long <- weights_long
weights_all_long$Bank <- factor(weights_all_long$Bank, levels = bank_order)

pal_all <- pal_banks[bank_order]
# start plotting at first "real" weight date (after the 90-day window)
start_plot_date <- as.Date(dates_ret[s_opt])   # s_opt = 90 in your script
weights_all_long_plot <- dplyr::filter(weights_all_long, Date >= start_plot_date)

png(file.path(output_path, "GMV_Weights_AllBanks_StackedArea.png"),
    width = 3000, height = 1700, res = 300, bg = "white")
print(
  ggplot(weights_all_long_plot, aes(Date, Weight, fill = Bank)) +
    geom_area(position = "stack", alpha = 0.98, linewidth = 0, na.rm = TRUE) +
    scale_fill_manual(values = pal_all, breaks = bank_order, limits = bank_order, name = NULL) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = "GMV Weights (stack = 1)") +
    theme_transparent_bottom +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
)
dev.off()

# ==============================================================================




# ===============================================================
# Build FIXED adjacency matrices from CSV inputs (no RData needed)
# - Uses FRM_SC_config.R (paths, dates, s/tau/I/J, channel)
# - Reads your input CSVs via FRM_SC_load_data.R (if needed)
# - Writes: Output/<channel>/Adj_Matrices/Fixed/adj_matrix_YYYYMMDD.csv
# - Writes: Output/<channel>/Lambda/Fixed/lambdas_fixed_<start>_<end>.csv
# ===============================================================

rm(list = setdiff(ls(), c()))
suppressPackageStartupMessages({ library(zoo) })

# ---- config + utils + FRM algorithm ----
source("FRM_SC_config.R")
source("FRM_SC_utils.R")
if (!file.exists("FRM_Statistics_Algorithm.R"))
  stop("Missing FRM_Statistics_Algorithm.R")
source("FRM_Statistics_Algorithm.R")

# ---- load stage_20 or build it from CSVs once ----
st20 <- file.path(output_path, "stage_20_loaded.RData")
if (!file.exists(st20)) {
  message("stage_20_loaded.RData not found — running FRM_SC_load_data.R once...")
  source("FRM_SC_load_data.R")
}
load(file.path(output_path, "stage_20_loaded.RData"))
# objects we expect now: dates, ticker, stock_return, macro_return, M_stock, M_macro, mktcap_index, s, tau, I, J

# ---- ensure output folders ----
fixed_dir  <- file.path(output_path, "Adj_Matrices", "Fixed")
lambda_dir <- file.path(output_path, "Lambda", "Fixed")
dir.create(fixed_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(lambda_dir, recursive = TRUE, showWarnings = FALSE)

# ---- helpers for date indexing ----
idx_start_at_or_after <- function(key, target){ w <- which(key >= target); if (length(w)) w[1] else NA_integer_ }
idx_end_at_or_before  <- function(key, target){ w <- which(key <= target); if (length(w)) tail(w,1) else NA_integer_ }

# ---- fixed window from config ----
N0_fixed <- idx_start_at_or_after(ticker, date_start_fixed)
N1_fixed <- idx_end_at_or_before( ticker, date_end_fixed)
if (is.na(N0_fixed) || is.na(N1_fixed) || (N1_fixed - N0_fixed + 1) < s) {
  stop(sprintf("Fixed window [%s, %s] too short for s=%d.", date_start_fixed, date_end_fixed, s))
}

# ---- lock Top-J once at the START of the fixed window ----
J_eff <- min(J, M_stock)
biggest_index_fixed <- as.integer(mktcap_index[N0_fixed, 2:(J_eff + 1), drop = FALSE])
M_J <- J_eff + M_macro

bank_names  <- colnames(stock_return)
macro_names <- if (M_macro > 0) colnames(macro_return) else character(0)
col_full    <- c(bank_names[biggest_index_fixed], macro_names)

# ---- iterate days in fixed window and create CSV adjacencies ----
N_fixed <- N1_fixed - N0_fixed + 1
lambdas_fixed <- matrix(NA_real_, N_fixed, J_eff + 1)
lambdas_fixed[, 1] <- ticker[N0_fixed:N1_fixed]

n_written <- 0L
for (t in N0_fixed:N1_fixed) {
  rows_ret <- (t - s):(t - 1)
  if (min(rows_ret) < 1) next
  
  X <- cbind(
    stock_return[rows_ret, biggest_index_fixed, drop = FALSE],
    if (M_macro > 0) macro_return[rows_ret, , drop = FALSE] else NULL
  )
  X[!is.finite(X)] <- 0
  if (!all(colSums(X != 0) > 0)) next
  
  A <- matrix(0, M_J, M_J)
  lam_row <- rep(NA_real_, J_eff)
  
  for (k in 1:M_J) {
    est <- tryCatch(FRM_Quantile_Regression(as.matrix(X), k, tau, I), error = function(e) NULL)
    if (is.null(est)) next
    gacv   <- est$Cgacv
    k_best <- if (!any(is.finite(gacv))) length(gacv) else which(gacv == min(gacv[is.finite(gacv)], na.rm = TRUE))[1]
    est_lambda <- abs(data.matrix(est$lambda[k_best]))
    est_beta   <- t(as.matrix(est$beta[k_best, ]))
    A[k, -k] <- est_beta
    if (k <= J_eff) lam_row[k] <- est_lambda
  }
  
  colnames(A) <- colnames(X)
  rownames(A) <- colnames(X)
  
  out_csv <- file.path(fixed_dir, paste0("adj_matrix_", ticker[t], ".csv"))
  write.csv(A, out_csv, quote = FALSE)
  n_written <- n_written + 1L
  
  lambdas_fixed[t - N0_fixed + 1, 2:(J_eff + 1)] <- lam_row
}

# ---- save helper lambdas CSV for the fixed universe ----
colnames(lambdas_fixed) <- c("date", bank_names[biggest_index_fixed])
lam_csv <- file.path(lambda_dir, paste0("lambdas_fixed_", date_start_fixed, "_", date_end_fixed, ".csv"))
write.csv(lambdas_fixed, lam_csv, row.names = FALSE, quote = FALSE)

message("✅ Fixed adjacency CSVs written: ", n_written)
message("📁 Fixed folder: ", normalizePath(fixed_dir))
message("🧾 Lambdas CSV:  ", normalizePath(lam_csv))


# ---- network GIF across time (Adj_Matrices -> GIF) ----

# ===============================================================
# Network GIF (FIXED) — circular layout, bottom DATE — FRM label
# Nodes sized by lambda (global scaling), UPPERCASE labels
# Edges: into main bank = blue; out of main bank = orange; others grey
# Transparent background
# ===============================================================

suppressPackageStartupMessages({
  library(igraph)
  library(magick)
  library(scales)
})

# ---------- 0) Config ----------
stopifnot(file.exists("FRM_SC_config.R"))
source("FRM_SC_config.R")     # output_path, website_path, channel, date_end, date_start_fixed, date_end_fixed, etc.
if (file.exists("FRM_SC_utils.R")) source("FRM_SC_utils.R")  

fixed_adj_dir <- file.path(output_path, "Adj_Matrices", "Fixed")
frm_csv       <- file.path(output_path, "Lambda", paste0("FRM_", channel, "_index.csv"))
out_dir       <- file.path(website_path, date_end, "Fixed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
gif_path      <- file.path(out_dir, paste0("Network_", date_start_fixed, "_", date_end_fixed, "_", channel, ".gif"))

# ---------- 1) FRM index & lookup ----------
stopifnot(file.exists(frm_csv))
FRM_index <- read.csv(frm_csv, stringsAsFactors = FALSE, check.names = FALSE)
names(FRM_index) <- tolower(names(FRM_index))
stopifnot(all(c("date","frm") %in% names(FRM_index)))
FRM_index$date <- suppressWarnings(as.Date(FRM_index$date))
if (any(is.na(FRM_index$date))) {
  d1 <- suppressWarnings(as.Date(FRM_index$date, "%Y-%m-%d"))
  d2 <- suppressWarnings(as.Date(FRM_index$date, "%Y/%m/%d"))
  FRM_index$date[is.na(FRM_index$date)] <- d1[is.na(FRM_index$date)]
  FRM_index$date[is.na(FRM_index$date)] <- d2[is.na(FRM_index$date)]
}
frm_lookup <- function(d){  # d is Date
  if (is.na(d)) return(NA_real_)
  w <- which(FRM_index$date == d)
  if (!length(w)) return(NA_real_)
  FRM_index$frm[w[1]]
}

# ---------- 2) Ensure ticker/N0/N1 ----------
if (!(exists("ticker") && length(ticker) > 0)) {
  f <- file.path(output_path, "stage_20_loaded.RData")
  if (file.exists(f)) load(f)
}
stopifnot(exists("ticker") && length(ticker) > 0)

idx_start_at_or_after <- function(key, target){ w <- which(key >= target); if (length(w)) w[1] else NA_integer_ }
idx_end_at_or_before  <- function(key, target){ w <- which(key <= target); if (length(w)) tail(w,1) else NA_integer_ }

if (!exists("N0_fixed") || !is.numeric(N0_fixed)) N0_fixed <- idx_start_at_or_after(ticker, date_start_fixed)
if (!exists("N1_fixed") || !is.numeric(N1_fixed)) N1_fixed <- idx_end_at_or_before( ticker, date_end_fixed)
stopifnot(is.finite(N0_fixed), is.finite(N1_fixed), N1_fixed >= N0_fixed)

# ---------- 3) Lambda table (prefer in-memory; else latest CSV) ----------
if (!exists("FRM_individ_fixed")) {
  lf_dir <- file.path(output_path, "Lambda", "Fixed")
  lf_all <- list.files(lf_dir, pattern = "^lambdas_fixed_.*\\.csv$", full.names = TRUE)
  if (!length(lf_all)) stop("No lambdas_fixed_*.csv found in: ", normalizePath(lf_dir))
  lf_path <- lf_all[which.max(file.info(lf_all)$mtime)]
  LFR     <- read.csv(lf_path, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot("date" %in% names(LFR))
  key <- gsub("\\D", "", as.character(LFR[[1]]))  # keep digits only (YYYYMMDD)
  FRM_individ_fixed <- cbind(as.integer(key), LFR[ , setdiff(names(LFR), "date"), drop = FALSE])
  colnames(FRM_individ_fixed)[1] <- "date"
}



if (!is.numeric(FRM_individ_fixed[,1])) {
  FRM_individ_fixed[,1] <- as.integer(gsub("\\\\D","", as.character(FRM_individ_fixed[,1])))
}

lambdas_fixed_vals <- as.matrix(FRM_individ_fixed[, -1, drop = FALSE])
storage.mode(lambdas_fixed_vals) <- "numeric"
lambdas_fixed_vals[!is.finite(lambdas_fixed_vals)] <- 0

if (!exists("J_eff_fixed") || !is.numeric(J_eff_fixed) || J_eff_fixed < 1) {
  J_eff_fixed <- ncol(lambdas_fixed_vals)
}

# Global scaling for node sizes (stable across frames)
rng <- range(lambdas_fixed_vals, na.rm = TRUE)
rescale_safe <- function(v) {
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) return(rep(10, length(v)))
  scales::rescale(v, to = c(6, 22), from = rng)
}

# ---------- 4) Labels & main bank ----------
if (!exists("stock_main") || !is.character(stock_main) || !length(stock_main)) {
  means <- colMeans(lambdas_fixed_vals, na.rm = TRUE)
  stock_main <- names(sort(means, decreasing = TRUE))[1]
  if (!length(stock_main)) stock_main <- colnames(lambdas_fixed_vals)[1]
}
label_vec <- function(cols){
  cols
}
gif_start <- idx_start_at_or_after(ticker, 20200101)
gif_end   <- idx_end_at_or_before(ticker,  20201231)
# ---------- 5) Render frames (bottom margin text) ----------
stopifnot(dir.exists(fixed_adj_dir))
fig <- image_graph(width = 1200, height = 1200, res = 120, bg = "white")
options(show.error.messages = FALSE)

for (t in gif_start:gif_end) {
  fcsv <- file.path(fixed_adj_dir, paste0("adj_matrix_", ticker[t], ".csv"))
  if (!file.exists(fcsv)) next

  adj0 <- tryCatch(read.csv(fcsv, header = TRUE, sep = ",", row.names = 1, check.names = FALSE),
                   error = function(e) NULL)
  if (is.null(adj0)) next
  adj0 <- as.matrix(adj0)
  storage.mode(adj0) <- "numeric"
  adj0[!is.finite(adj0)] <- 0
  J_here <- min(J_eff_fixed, nrow(adj0), ncol(adj0))
  adj0   <- adj0[1:J_here, 1:J_here, drop = FALSE]
  diag(adj0) <- 0

  netw1 <- graph_from_adjacency_matrix(adj0, mode = "directed", weighted = TRUE)

  row_id <- which(as.integer(FRM_individ_fixed[, 1]) == as.integer(ticker[t]))
  if (length(row_id)) row_id <- row_id[1]
  sizes <- if (length(row_id) && row_id >= 1 && row_id <= nrow(lambdas_fixed_vals)) {
    v <- as.numeric(lambdas_fixed_vals[row_id, 1:J_here]); v[!is.finite(v)] <- 0; rescale_safe(v)
  } else rep(10, J_here)

  cols_now <- colnames(adj0)
  lbl_now  <- toupper(label_vec(cols_now))  # UPPERCASE
  main_lbl <- toupper(stock_main)

  V(netw1)$label <- lbl_now
  V(netw1)$size  <- sizes
  V(netw1)$color <- ifelse(lbl_now == main_lbl, "orange", "lightgrey")
  V(netw1)$frame.color <- NA
  V(netw1)$label.cex <- 0.9

  ecols <- rep("gray", length(E(netw1)))
  if (ecount(netw1) > 0) {
    heads <- head_of(netw1, E(netw1))$name
    tails <- tail_of(netw1, E(netw1))$name
    ecols[heads == stock_main] <- "blue"    # into main bank
    ecols[tails == stock_main] <- "orange"  # out of main bank
  }
  E(netw1)$color <- ecols
  E(netw1)$width <- 0.8
  E(netw1)$arrow.size  <- 0.9
  E(netw1)$arrow.width <- 1

  # Bottom margin text (outside the circle)
  par(mar = c(7, 0, 0, 0), xpd = NA, bg = NA)
  plot(netw1,
       layout = layout_in_circle,
       edge.curved = 0.15,
       asp = 0)

  this_date <- as.Date(as.character(ticker[t]), "%Y%m%d")
  frm_val   <- frm_lookup(this_date)
  lab <- if (is.finite(frm_val)) {
    paste0(format(this_date, "%Y-%m-%d"), "  —  FRM: ", formatC(frm_val, format = "f", digits = 5))
  } else {
    format(this_date, "%Y-%m-%d")
  }
  mtext(lab, side = 1, line = 2.6, cex = 1.1, font = 2)
}

options(show.error.messages = TRUE)
dev.off()

# ---------- 6) Animate & save ----------
animation <- image_animate(fig, fps = 5, loop = 0)
image_write(animation, gif_path)
message("Network GIF saved to: ", normalizePath(gif_path))

