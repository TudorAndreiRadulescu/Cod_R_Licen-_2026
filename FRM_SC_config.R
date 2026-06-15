rm(list = ls(all = TRUE)); graphics.off(); options(digits=6, warn=-1)

# --- Edit this single line to your local project folder ---
wdir <- "C:/Users/troiu/Desktop/fotos/licenta/R/Claude Version"
setwd(wdir)

# --- Channel & dates ---
channel <- "Banci"
date_start_source <- 20010103; date_end_source <- 20251231
date_start        <- 20010103; date_end        <- 20251231
date_start_fixed  <- 20010103; date_end_fixed  <- 20251231

# --- Model knobs ---
s   <- 90        # rolling window (trading days)
tau <- 0.05      # quantile (5% tail risk)
I   <- 25        # max LASSO iterations
L   <- 5         # (reserved, not used in this pipeline)
J   <- 7         # number of banks to include
stock_main <- "BT"
quantiles  <- c(0.99,0.95,0.90,0.85,0.80,0.75,0.70,0.65,0.60,0.55,0.50,0.25)

# --- Paths ---
input_path  <- file.path("Input",  channel, paste0(date_start_source, "-", date_end_source))
output_path <- if (tau == 0.05 & s == 90) file.path("Output", channel) else
  file.path("Output", channel, paste0("Sensitivity/tau=", 100*tau, "/s=", s))
website_path <- file.path("Website", channel)

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
invisible(lapply(dirs, dir.create, recursive=TRUE, showWarnings=FALSE))

# --- Helper functions needed by all scripts ---
idx_start_at_or_after <- function(key, target) { w <- which(key >= target); if (length(w)) w[1] else NA_integer_ }
idx_end_at_or_before  <- function(key, target) { w <- which(key <= target); if (length(w)) tail(w,1) else NA_integer_ }

# --- Packages ---
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(lubridate); library(zoo)
  library(ggplot2); library(data.table); library(igraph); library(magick); library(scales)
  library(stringr); library(plotly); library(reshape2); library(quadprog); library(MASS)
})

# --- FRM Algorithm ---
if (!file.exists("FRM_Statistics_Algorithm.R"))
  stop("Missing FRM_Statistics_Algorithm.R in project root.")
source("FRM_Statistics_Algorithm.R")
