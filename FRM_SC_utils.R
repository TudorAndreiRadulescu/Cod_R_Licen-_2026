# ============================================================
# FRM_SC_utils.R — shared helpers for the Banci FRM pipeline
# ============================================================

# --- ggplot2 theme: white background, legend at bottom ---
theme_transparent_bottom <- theme(
  panel.background      = element_rect(fill="white", colour=NA),
  plot.background       = element_rect(fill="white", colour=NA),
  legend.box.background = element_rect(fill="white"),
  legend.background     = element_rect(fill="white"),
  legend.position       = "bottom",
  legend.direction      = "horizontal",
  legend.box.just       = "center",
  axis.line             = element_line(colour="black"),
  axis.text             = element_text(size=12, colour="black"),
  axis.title            = element_text(size=14, colour="black"),
  panel.grid            = element_blank()
)

# --- Risk color palette ---
risk_colors <- c(
  "1. Low risk"      = "green",
  "2. General risk"  = "blue",
  "3. Elevated risk" = "yellow",
  "4. High risk"     = "orange",
  "5. Severe risk"   = "red"
)

# --- Bank display name map (replaces old crypto COIN_NAME_MAP) ---
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

# --- Flexible file picker: tries strong regex first, then weak ---
pick_file_flexible <- function(path, strong_regex, weak_regex, role="file") {
  f <- list.files(path, pattern=strong_regex, full.names=TRUE, ignore.case=TRUE)
  if (!length(f)) f <- list.files(path, pattern=weak_regex, full.names=TRUE, ignore.case=TRUE)
  if (!length(f)) stop(sprintf("No %s found in %s", role, path))
  extract_max_date_from_name <- function(x) {
    b <- basename(x)
    ds <- gregexpr("[0-9]{8}", b, perl=TRUE)
    vals <- regmatches(b, ds)[[1]]
    if (length(vals)) as.integer(max(vals)) else NA_integer_
  }
  name_dates <- vapply(f, extract_max_date_from_name, integer(1))
  if (any(!is.na(name_dates))) return(f[which.max(name_dates)])
  f[1]
}

# --- Robust date parser (handles multiple formats) ---
parse_date_robust <- function(x) {
  if (inherits(x, "Date"))   return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x))         return(as.Date(round(x), origin="1970-01-01"))
  d <- suppressWarnings(as.Date(x))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, "%Y-%m-%d"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, "%m/%d/%Y"))
  d
}

# --- Date index helpers ---
idx_start_at_or_after <- function(key, target) { w <- which(key >= target); if (length(w)) w[1] else NA_integer_ }
idx_end_at_or_before  <- function(key, target) { w <- which(key <= target); if (length(w)) tail(w,1) else NA_integer_ }

# --- Clamp to [0,1] ---
clamp01 <- function(v) pmin(pmax(v, 0), 1)
