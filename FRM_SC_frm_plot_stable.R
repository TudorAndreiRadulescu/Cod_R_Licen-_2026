# =====================================================================
# FRM_SC_frm_plot_stable.R — Blue FRM index line plot
# =====================================================================
source("FRM_SC_config.R")
source("FRM_SC_utils.R")

load(file.path(output_path, "stage_50_index.RData"))
# Loaded: FRM_index (data.frame: date, frm)

dir.create(file.path(website_path, date_end), recursive=TRUE, showWarnings=FALSE)
out_png <- file.path(website_path, date_end, paste0("FRM_", channel, "_Index.png"))

png(out_png, width=1200, height=700, bg="white")
print(
  ggplot(FRM_index, aes(x=date, y=frm)) +
    geom_line(linewidth=0.9, color="blue") +
    scale_x_date(date_breaks="1 year", date_labels="%Y") +
    labs(title=NULL, x=NULL, y=paste0("FRM@", channel)) +
    theme_transparent_bottom +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
)
dev.off()
message("Done. Saved: ", normalizePath(out_png, mustWork=FALSE))
