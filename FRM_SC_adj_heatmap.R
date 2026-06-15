# =====================================================================
# FRM_SC_adj_heatmap.R — Adjacency matrix heatmap for a chosen date
# =====================================================================
source("FRM_SC_config.R")
source("FRM_SC_utils.R")

suppressPackageStartupMessages({ library(ggplot2); library(reshape2) })

# ---- Core heatmap function ----
plot_adj_matrix <- function(
    adj_file,
    out_dir   = file.path(output_path, "Network"),
    digits    = 2,
    zero_eps  = 1e-3,   # |value| < zero_eps -> shown as blank
    text_size = 3
) {
  stopifnot(file.exists(adj_file))

  date_str <- gsub("adj_matrix_|\\.csv", "", basename(adj_file))
  mat      <- as.matrix(read.csv(adj_file, row.names=1, check.names=FALSE))
  mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0

  L <- max(abs(mat), na.rm=TRUE)
  if (!is.finite(L) || L == 0) L <- 1e-8

  # Use toupper of the column name (pretty_bank for human-readable labels)
  rn     <- rownames(mat); cn <- colnames(mat)
  rn_lbl <- toupper(rn);   cn_lbl <- toupper(cn)

  df_long <- reshape2::melt(mat, varnames=c("From","To"), value.name="Weight")
  df_long$From <- factor(df_long$From, levels=rn, labels=rn_lbl)
  df_long$To   <- factor(df_long$To,   levels=cn, labels=cn_lbl)

  # Blank near-zero cells
  df_long$Weight_plot <- df_long$Weight
  df_long$Weight_plot[abs(df_long$Weight_plot) < zero_eps] <- NA
  df_long$lab     <- ifelse(abs(df_long$Weight) < zero_eps, "",
                            sprintf(paste0("%.", digits, "f"), df_long$Weight))
  df_long$lab_col <- ifelse(abs(df_long$Weight) > 0.6 * L, "white", "black")

  p <- ggplot(df_long, aes(x=To, y=From, fill=Weight_plot)) +
    geom_tile(color=NA) +
    geom_text(aes(label=lab, colour=I(lab_col)), size=text_size) +
    scale_fill_gradient2(
      low="#3B82F6", mid="white", high="#EF4444",
      midpoint=0, limits=c(-L, L), oob=scales::squish,
      na.value="transparent", name=NULL
    ) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    labs(title=NULL, x="To", y="From") +
    coord_fixed() +
    theme_minimal(base_size=14) +
    theme(
      panel.grid.major      = element_blank(),
      panel.grid.minor      = element_blank(),
      axis.text.x           = element_text(angle=45, hjust=1),
      panel.background      = element_rect(fill="white", colour=NA),
      plot.background       = element_rect(fill="white", colour=NA),
      legend.box.background = element_rect(fill="white"),
      legend.background     = element_rect(fill="white")
    )

  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
  out_png <- file.path(out_dir, paste0("AdjMatrix_", date_str, ".png"))
  ggsave(out_png, p, width=8, height=6, bg="white", dpi=300)
  message("Saved: ", normalizePath(out_png, mustWork=FALSE))
  invisible(p)
}

# ---- Helper: convert any date format to integer YYYYMMDD ----
to_int_yyyymmdd <- function(x) {
  if (is.numeric(x)) return(as.integer(x))
  d <- suppressWarnings(as.Date(x))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, "%Y-%m-%d"))
  if (is.na(d)) d <- suppressWarnings(as.Date(x, "%m/%d/%Y"))
  if (is.na(d)) stop("Cannot parse date: ", x)
  as.integer(format(d, "%Y%m%d"))
}

# ---- List available adjacency CSVs ----
list_adj_files <- function(fixed=FALSE) {
  dir <- if (fixed) file.path(output_path, "Adj_Matrices", "Fixed") else
                    file.path(output_path, "Adj_Matrices")
  if (!dir.exists(dir)) return(character(0))
  list.files(dir, pattern="^adj_matrix_[0-9]{8}\\.csv$", full.names=TRUE)
}

extract_date_from_name <- function(fp) {
  b <- basename(fp)
  m <- regexpr("[0-9]{8}", b, perl=TRUE)
  if (m[1] == -1) return(NA_integer_)
  as.integer(substr(b, m[1], m[1] + attr(m, "match.length")[1] - 1))
}

# ---- Pick the file for a given date (exact or nearest) ----
pick_file_for_date <- function(date_int, fixed=FALSE, nearest_ok=TRUE) {
  files <- list_adj_files(fixed)
  if (!length(files)) return(NA_character_)
  d     <- vapply(files, extract_date_from_name, integer(1))
  exact <- which(d == date_int)
  if (length(exact)) return(files[exact[1]])
  if (!nearest_ok) return(NA_character_)
  ok <- is.finite(d)
  if (!any(ok)) return(NA_character_)
  ii <- which.min(abs(d[ok] - date_int))
  files[which(ok)[ii]]
}

# ---- Public wrapper: plot by date ----
plot_adj_for_date <- function(date_input, fixed=FALSE, nearest_ok=TRUE,
                              out_dir=file.path(output_path, "Network")) {
  date_int <- to_int_yyyymmdd(date_input)
  f <- pick_file_for_date(date_int, fixed=fixed, nearest_ok=nearest_ok)
  if (is.na(f)) stop("No adjacency CSV found for/near date ", date_int)
  message("Using: ", basename(f))
  plot_adj_matrix(f, out_dir=out_dir)
}

# ---- Auto: find the day with the largest FRM drop ----
find_crash_day_frm <- function() {
  frm_csv <- file.path(output_path, "Lambda", paste0("FRM_", channel, "_index.csv"))
  if (!file.exists(frm_csv)) stop("FRM index CSV not found: ", frm_csv)
  X <- read.csv(frm_csv, stringsAsFactors=FALSE)
  names(X) <- tolower(names(X))
  X$date <- suppressWarnings(as.Date(X$date))
  X <- X[is.finite(X$frm) & !is.na(X$date), ]
  if (nrow(X) < 2) stop("Not enough FRM observations.")
  ch <- diff(X$frm)
  i  <- which.min(ch)
  list(date=X$date[i+1], yyyymmdd=as.integer(format(X$date[i+1], "%Y%m%d")), drop=ch[i])
}

plot_adj_for_crash_day <- function(fixed=FALSE, nearest_ok=TRUE,
                                   out_dir=file.path(output_path, "Network")) {
  info <- find_crash_day_frm()
  message(sprintf("Crash day by FRM drop: %s (ΔFRM = %.6f)",
                  format(info$date, "%Y-%m-%d"), info$drop))
  plot_adj_for_date(info$yyyymmdd, fixed=fixed, nearest_ok=nearest_ok, out_dir=out_dir)
}

# =====================================================================
# --- Example calls (uncomment to use) ---
# Plot the adjacency heatmap for a specific date:
# plot_adj_for_date("2020-03-16")   # COVID crash day
# plot_adj_for_date("2022-02-24")   # Russia-Ukraine invasion
#
# Auto-find and plot the worst FRM crash day:
# plot_adj_for_crash_day()
# =====================================================================
