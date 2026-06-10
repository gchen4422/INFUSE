# Plotting utilities for INFUSE website
# Source this file at the top of analysis Rmd files:
#   source(here::here("code/plot_utils.R"))

library(ggplot2)
library(dplyr)
library(ggrepel)

# ── Shared theme ──────────────────────────────────────────────────────────────
custom_theme <- function() {
  theme(
    axis.text.x      = element_text(size = 7.5),
    axis.text.y      = element_text(size = 7.5),
    axis.title.x     = element_text(size = 10, face = "bold"),
    axis.title.y     = element_text(size = 10, face = "bold"),
    strip.text.x     = element_text(size = 12, face = "bold"),
    strip.text.y     = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 12, face = "bold"),
    plot.title       = element_text(size = 14, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    plot.tag         = element_text(size = 12),
    axis.line        = element_line(color = "black")
  )
}

drop_guides <- guides(color = "none", shape = "none", fill = "none", linetype = "none")

r2_scale <- scale_color_stepsn(
  colors      = c("navy", "lightskyblue", "green", "orange", "red"),
  breaks      = seq(0.2, 0.8, by = 0.2),
  limits      = c(0, 1),
  show.limits = TRUE,
  na.value    = "grey50",
  name        = expression(R^2)
)

# Theme for genomic annotation tracks (no x-axis, compact y-axis)
track_theme <- theme(
  axis.title.y       = element_text(angle = 0, hjust = 1, vjust = 0.5, margin = margin(r = 12)),
  axis.text.y        = element_blank(),
  axis.ticks.y       = element_blank(),
  axis.line.y        = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  axis.title.x       = element_blank(),
  axis.text.x        = element_blank(),
  axis.ticks.x       = element_blank(),
  axis.line.x        = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank()
)

# ── Plotting functions ────────────────────────────────────────────────────────
gwas_plot_fun <- function(data_plot, xlab_name, ylab_name, yintercept) {
  pos_rng <- range(data_plot$POS, na.rm = TRUE)
  pos_mb  <- max(data_plot$POS) > 1e6
  ggplot() +
    geom_point(data = filter(data_plot, Lead_SNP == 0), aes(POS, PIP, color = r2), size = 1) +
    geom_point(data = filter(data_plot, Lead_SNP == 1), aes(POS, PIP), size = 1.5, color = "red") +
    geom_text_repel(data = filter(data_plot, Lead_SNP == 1), aes(POS, PIP, label = SNP),
                    vjust = 1.2, size = 7/14*5, show.legend = FALSE) +
    r2_scale +
    geom_hline(yintercept = yintercept, linetype = "dashed", color = "grey50", size = 0.5) +
    geom_vline(xintercept = filter(data_plot, Lead_SNP == 1)$POS,
               linetype = "dashed", color = "grey50", size = 0.5) +
    scale_x_continuous(
      labels = if (pos_mb) function(x) paste0(x/1e6, " MB") else function(x) paste0(x/1e3, " KB"),
      expand = expansion(mult = c(0.05, 0.05)), limits = pos_rng
    ) +
    xlab(xlab_name) + ylab(ylab_name) + theme_bw() + custom_theme()
}

finemap_plot_fun <- function(data_plot, xlab_name, ylab_name, yintercept) {
  pos_rng <- range(data_plot$POS, na.rm = TRUE)
  pos_mb  <- max(data_plot$POS) > 1e6
  ggplot() +
    geom_point(data = data_plot, aes(POS, PIP, color = r2, shape = cat)) +
    scale_shape_manual(name = "Category", drop = FALSE, values = c(20, 24, 25, 23, 22)) +
    geom_text_repel(data = filter(data_plot, Lead_SNP == 1), aes(POS, PIP, label = SNP),
                    vjust = 1.2, size = 7/14*5, show.legend = FALSE) +
    r2_scale +
    geom_hline(yintercept = yintercept, linetype = "dashed", color = "grey50", size = 0.5) +
    geom_vline(xintercept = filter(data_plot, Lead_SNP == 1)$POS,
               linetype = "dashed", color = "grey50", size = 0.5) +
    scale_x_continuous(
      labels = if (pos_mb) function(x) paste0(x/1e6, " MB") else function(x) paste0(x/1e3, " KB"),
      expand = expansion(mult = c(0.05, 0.05)), limits = pos_rng
    ) +
    xlab(xlab_name) + ylab(ylab_name) + theme_bw() + custom_theme()
}

gene_range_plot_fun <- function(gene_list_data, plot.range) {
  ggplot(gene_list_data) +
    geom_linerange(aes(x = Gene, ymin = Start, ymax = End)) +
    geom_text(aes(x = Gene, y = Start, label = Gene), hjust = "right", size = 5/14*5) +
    scale_y_continuous(limits = plot.range, labels = function(y) paste0(y/1e6, " MB"),
                       expand = expansion(mult = c(0.05, 0.05))) +
    coord_flip() +
    ylab(paste0("chr", unique(gsub("chr", "", gene_list_data$Chrom)))) + xlab("Gene") +
    theme_bw() + theme(
      axis.text.x      = element_text(size = 7.5),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.title.x     = element_text(size = 10, face = "bold"),
      axis.title.y     = element_text(size = 10, face = "bold"),
      strip.text       = element_text(size = 12, face = "bold"),
      strip.background = element_blank(),
      panel.grid       = element_blank(),
      panel.border     = element_blank(),
      axis.line        = element_line(color = "black")
    )
}

#' Plot a single binary annotation track, highlighting a specific SNP in red
plot_annotation_highlight <- function(data, annotation_col, highlight_snp_id,
                                      track_label = annotation_col,
                                      track_color = "gray30",
                                      snp_col = "SNP", pos_col = "POS", chr_col = "CHR",
                                      plot_range = NULL) {
  plot_data     <- data[data[[annotation_col]] == 1, , drop = FALSE]
  highlight_snp <- data[data[[snp_col]] == highlight_snp_id, , drop = FALSE]
  pos_rng       <- if (is.null(plot_range)) range(data[[pos_col]], na.rm = TRUE) else plot_range
  chr_label     <- paste0("chr", unique(data[[chr_col]])[1])

  ggplot(plot_data, aes(x = .data[[pos_col]], y = 1)) +
    geom_point(color = track_color, shape = "|", size = 5, alpha = 0.7) +
    geom_point(data = highlight_snp, aes(x = .data[[pos_col]], y = 1),
               color = "red", size = 1.6, inherit.aes = FALSE) +
    scale_x_continuous(limits = pos_rng, labels = function(x) paste0(round(x/1e6, 2), " Mb")) +
    scale_y_continuous(limits = c(0.8, 1.2), breaks = NULL) +
    labs(x = chr_label, y = track_label) +
    theme_bw() + custom_theme()
}

# ── Data helpers ──────────────────────────────────────────────────────────────

#' Extract r² vector relative to a lead SNP from an LD matrix.
#' Returns zeros if the lead SNP is absent from the matrix.
get_ld_vec <- function(cov_mat, lead_snp) {
  idx <- which(colnames(cov_mat) == lead_snp)
  if (length(idx) == 0) setNames(rep(0, ncol(cov_mat)), colnames(cov_mat))
  else                  setNames(cov_mat[, idx]^2,       colnames(cov_mat))
}

#' Assign method-specific PIP categories for INFUSE vs. comparator methods
assign_cats <- function(df) {
  df %>% mutate(
    Paintor_cat      = case_when(Paintor_PIP > 0.5 ~ 4, TRUE ~ 0),
    Paintor_annot_cat = case_when(Paintor_PIP_all > 0.5 ~ 4, TRUE ~ 0),
    MESuSiE_cat      = case_when(MESuSiE_PIP_WB > 0.5 ~ 1, MESuSiE_PIP_BB > 0.5 ~ 2,
                                 MESuSiE_PIP_Shared > 0.5 ~ 3, TRUE ~ 0),
    INFUSE_cat       = case_when(MESuSiE_PIP_WB_all_mlk > 0.5 ~ 1,
                                 MESuSiE_PIP_BB_all_mlk > 0.5 ~ 2,
                                 MESuSiE_PIP_Shared_all_mlk > 0.5 ~ 3, TRUE ~ 0)
  )
}

#' Build a plot-ready data frame for a fine-mapping method
make_plot_data <- function(df, r2_col, pip_col, cat_col, lead_snp) {
  df %>% mutate(
    r2       = .data[[r2_col]],
    PIP      = .data[[pip_col]],
    Lead_SNP = as.integer(SNP == lead_snp),
    POS      = as.numeric(POS),
    cat      = factor(.data[[cat_col]],
                      levels = c("0", "1", "2", "3", "4"),
                      labels = c("Non", "EUR", "AFR", "Shared", "Paintor"))
  ) %>% select(SNP, POS, r2, PIP, Lead_SNP, cat)
}

#' Build a plot-ready data frame for a GWAS association panel
make_gwas_data <- function(df, r2_col, zscore_col, lead_snp) {
  df %>% mutate(
    r2       = .data[[r2_col]],
    PIP      = -log10(2 * pnorm(-abs(.data[[zscore_col]]))),
    Lead_SNP = as.integer(SNP == lead_snp),
    POS      = as.numeric(POS)
  ) %>% dplyr::select(SNP, POS, r2, PIP, Lead_SNP)
}
