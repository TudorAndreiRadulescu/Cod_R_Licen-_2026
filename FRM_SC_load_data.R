# =====================================================================
# FRM_SC_load_data.R — Load bank return data and build stage_20
# =====================================================================
source("FRM_SC_config.R")
source("FRM_SC_utils.R")

cat("\n[1] Reading Baza_FRM_Banci.csv...\n")
df <- read.csv("Baza_FRM_Banci.csv")

# Date vector (integer YYYYMMDD, used as the 'ticker' key throughout)
ticker <- as.integer(df$Date)
N      <- length(ticker)

# Bank columns (must match CSV column names exactly)
banci_cols <- c("AlphaBank", "BNPPariBas", "BRD", "BT", "Erste", "ING", "Unicredit")

# Everything else (except Date and banks) is treated as a macro variable
macro_cols <- setdiff(names(df), c("Date", banci_cols))

cat("[2] Banks (", length(banci_cols), "):", paste(banci_cols, collapse=", "), "\n")
cat("    Macro vars (", length(macro_cols), "):", paste(macro_cols, collapse=", "), "\n")

# Build return matrices
stock_return <- as.matrix(df[, banci_cols])
macro_return <- if (length(macro_cols) > 0) as.matrix(df[, macro_cols]) else matrix(0, N, 0)

M_stock <- ncol(stock_return)
M_macro <- ncol(macro_return)
M       <- M_stock + M_macro

# Synthetic market-cap index: all 7 banks are always included (static ranking 1:7).
# In the original crypto pipeline this was a daily market-cap ranking; here we
# fix it so every bank enters the model on every day.
mktcap_index <- cbind(ticker, matrix(rep(1:M_stock, each=N), nrow=N, ncol=M_stock))

cat("[3] Saving stage_20_loaded.RData...\n")
dir.create(output_path, recursive=TRUE, showWarnings=FALSE)
save(ticker, stock_return, macro_return, mktcap_index, M_stock, M_macro, M,
     file=file.path(output_path, "stage_20_loaded.RData"))
cat("Done. Ready for LASSO estimation.\n\n")
