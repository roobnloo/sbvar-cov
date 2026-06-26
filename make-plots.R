# make-plots.R
# Generate ggplot boxplots comparing SBVAR methods across simulation settings.
#
# Variants compared (Stage 2 / post-BEA):
#   sb_s2_rss, sb_s2_sig, sb_s2_pro  (SB variants)
#   nll_s2                             (NLL)
#   dt_s2                              (DT)
#
# Output: <sim_dir>/plots-<label>.pdf  (one file per setting)
#
# Usage:  Rscript make-plots.R <sim_dir>
#         (sim_dir defaults to "sim" if not provided)

library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================
# Parse arguments
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
sim_dir <- if (length(args) >= 1L) args[[1L]] else "sim"
if (!dir.exists(sim_dir)) stop(sprintf("Directory not found: %s", sim_dir))

# ============================================================
# Locate result files
# ============================================================

rds_files <- list.files(sim_dir,
  pattern = "^run-sim-sbvar-.*\\.rds$",
  full.names = TRUE
)
if (length(rds_files) == 0L) stop(sprintf("No results .rds files found in %s", sim_dir))

# ============================================================
# Flatten one results list into a long data frame
# ============================================================

VARIANTS <- c("sb_s2_rss", "sb_s2_sig", "sb_s2_pro", "nll_s2", "dt_s2")
METHOD_LABELS <- c(
  sb_s2_rss = "SB RSS",
  sb_s2_sig = "SB sigma",
  sb_s2_pro = "SB profile",
  nll_s2    = "NLL",
  dt_s2     = "DT"
)

SCALAR_METRICS <- c("hd", "B_err", "Sigma_err")
METRIC_LABELS <- c(
  hd        = "log(1 + Hausdorff distance)",
  B_err     = expression(B[err]),
  Sigma_err = expression(Sigma[err])
)

flatten_results <- function(results, setting_label) {
  rows <- vector("list", length(results) * length(VARIANTS))
  idx <- 1L
  for (i in seq_along(results)) {
    rep <- results[[i]]
    for (v in VARIANTS) {
      entry <- rep[[v]]
      if (is.null(entry)) {
        row <- as.list(rep(NA_real_, length(SCALAR_METRICS)))
        names(row) <- SCALAR_METRICS
      } else {
        row <- lapply(SCALAR_METRICS, function(m) as.numeric(entry[[m]]))
        names(row) <- SCALAR_METRICS
      }
      row$rep     <- i
      row$variant <- v
      row$setting <- setting_label
      rows[[idx]] <- row
      idx <- idx + 1L
    }
  }
  df <- bind_rows(rows)
  df$method <- METHOD_LABELS[df$variant]
  df
}

# ============================================================
# Parse setting label from filename
# ============================================================

parse_label <- function(path) {
  base <- sub("^run-sim-sbvar-(.+)\\.rds$", "\\1", basename(path))
  # e.g. "scenario1_sigscale1.5" -> "Scenario 1  sigscale=1.5"
  base <- gsub("_", "  ", base)
  base <- sub("scenario(\\d+)\\s+(sigscale|sigma)(.+)", "Scenario \\1  \\2=\\3", base, perl = TRUE)
  base <- sub("scenario(\\d+)$", "Scenario \\1", base, perl = TRUE)
  base <- sub("\\s+sigscale=1$", "", base)
  base <- sub("\\s+sigma=1$",    "", base)
  base
}

# ============================================================
# Build full data frame
# ============================================================

all_df <- bind_rows(lapply(rds_files, function(f) {
  label   <- parse_label(f)
  results <- readRDS(f)
  flatten_results(results, label)
}))

all_df$scenario    <- sub("^(Scenario \\d+).*$", "\\1", all_df$setting)
all_df$sigma_label <- trimws(sub("^Scenario \\d+\\s*", "", all_df$setting))
all_df$sigma_label[all_df$sigma_label == ""] <- "(base)"

# ============================================================
# Plotting helpers
# ============================================================

METHOD_COLORS <- c(
  "SB RSS"     = "#92b7dc",
  "SB sigma"   = "#4DAC26",
  "SB profile" = "#cb8ee5",
  "NLL"        = "#D6604D",
  "DT"         = "#F4A749"
)

all_df$method <- factor(all_df$method, levels = names(METHOD_COLORS))

make_boxplot <- function(df, metric, y_label, setting) {
  sub_df <- df[, c("method", metric)]
  colnames(sub_df)[2] <- "value"
  sub_df <- sub_df[!is.na(sub_df$value), ]
  if (metric == "hd") sub_df$value <- log1p(sub_df$value)

  ggplot(sub_df, aes(x = method, y = value, fill = method)) +
    geom_boxplot(outlier.size = 0.8, linewidth = 0.4, width = 0.55) +
    scale_fill_manual(values = METHOD_COLORS, guide = "none") +
    labs(x = NULL, y = y_label) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(size = 10, face = "bold"),
      axis.text.x   = element_text(angle = 20, hjust = 1, size = 9),
      panel.grid.minor = element_blank()
    )
}

make_combined_plot <- function(df, setting) {
  plot_list <- lapply(seq_along(SCALAR_METRICS), function(i) {
    m   <- SCALAR_METRICS[i]
    lbl <- METRIC_LABELS[[m]]
    make_boxplot(df, m, lbl, if (i == 1L) setting else "")
  })

  if (requireNamespace("gridExtra", quietly = TRUE)) {
    library(gridExtra)
    gridExtra::arrangeGrob(
      grobs = plot_list, nrow = 1L,
      top = grid::textGrob(setting, gp = grid::gpar(fontsize = 13, fontface = "bold"))
    )
  } else if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    wrap_plots(plot_list, nrow = 1L) +
      plot_annotation(title = setting)
  } else {
    stop("Install 'patchwork' or 'gridExtra' to combine panels.")
  }
}

# ============================================================
# Produce one PDF per setting
# ============================================================

settings <- unique(all_df$setting)
out_dir  <- sim_dir

for (s in settings) {
  sub <- all_df[all_df$setting == s, ]
  p   <- make_combined_plot(sub, s)

  fname <- file.path(
    out_dir,
    paste0("plots-", gsub("[[:space:]]+", "_", gsub("=", "", s)), ".pdf")
  )
  ggsave(fname, p, width = 14, height = 4, device = "pdf")
  cat(sprintf("Saved %s\n", fname))
}

cat("\nDone.\n")
