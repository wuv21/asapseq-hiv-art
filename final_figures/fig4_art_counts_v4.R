source("defaultFigureSettings.R")
library(ggtext)

figA <- readRDS("../outs/rds/art_umap_labeledCluster_withPlotLabels.rds")
figB <- readRDS("../outs/rds/art_umap_haystackout.rds")
figC <- readRDS("../outs/rds/art_discreteAbsolute_matched_byDonor.rds")
figD <- readRDS("../outs/rds/art_cluster_proportions_by_donor.rds")
# figD <- readRDS("../outs/rds/art_discreteHivOnly_matched_byDonor.rds")

a_legend <- as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), nrow = 11,
      title = "Legend for (a)",
      title.position = "top",
      title.hjust = 0)) +
    theme(legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      legend.justification = "left",
      legend.text = element_text(margin = margin(r = 12, unit = "pt"), size = 6),
      legend.title = element_text(size = 7, color = "#000000", hjust = 0))))

b_legend <- as_ggplot(get_legend(figB + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1, title = "Legend for (b)",
      title.position = "top",
      title.hjust = 0)) +
    theme(
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      legend.justification = "left",
      legend.text = element_text(margin = margin(r = 12, unit = "pt"), size = 6),
      legend.title = element_text(size = 7, color = "#000000"))))


figABLegend <- a_legend / b_legend & subplotTheme 

figABLegend[[1]] <- figABLegend[[1]] +
  theme(plot.margin = margin(t = 20, l = 10))

figABLegend[[2]] <- figABLegend[[2]] +
  theme(plot.margin = margin(t = -40, b = 70, l = 10))

figA <- figA + umapPlotThemeNoLeg +
  subplotTheme

figB <- figB + umapPlotThemeNoLeg + 
  subplotTheme +
  theme(axis.title.y = element_text(color = "#FFFFFF00"),
    plot.margin = unit(c(0, 0, 0, 0), "pt"))

figAB <- figA + figB + plot_annotation(tag_levels = "a")

figD <- figD + subplotTheme +
  theme(
    strip.background = element_rect(color = NA, fill = "transparent"),
    legend.position = "bottom",
    legend.margin = margin(t = -10),
    plot.margin = margin(t = -20, l = 4, r = 20),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    panel.background = element_rect(colour = "black", fill = NA, size = 1)
  )

layout <- c(
  area(1, 1, 2, 4), #ab
  area(1, 5, 3, 9), #legend
  area(3, 1, 5, 9), #c
  area(6, 1, 7, 9)
)

discreteLollipopTheme <- subplotTheme + theme(
  axis.title.y = element_blank(),
  axis.title.x = element_markdown(margin = margin(t = 2.5)),
  axis.text.x = element_text(size = 6),
  plot.subtitle = element_text(size = 7, hjust = 0.5, margin = margin(b = -1)),
  plot.title.position = "panel",
  plot.margin = margin(l = -30, t = -20, r = 20)
)

decreasePointAndText <- function(x) {
  x$layers[[3]]$aes_params$size <- 5 / ggplot2:::.pt
  x$layers[[2]]$aes_params$size <- 0.75
  
  return(x)
}

figC <- lapply(figC, decreasePointAndText)
for (i in c(2:length(figC))) {
  figC[[i]] <- figC[[i]] + theme(axis.text.y = element_blank())
  
  if (i == 4) {
    figC[[i]] <- figC[[i]] +
      labs(x = "<span style='font-size:7pt'>Number of cells (<span style='color:#999999;'>HIV-</span> | <span style='color:#e63946;'>HIV+</span>)</span>")
  } else {
    figC[[i]] <- figC[[i]] + labs(x = "")
  }
}
figC[[1]] <- figC[[1]] + 
  labs(x = "") +
  theme(axis.text.y = element_text(angle = 45, hjust = 1, size = 5))

figC <- wrap_plots(figC, nrow = 1) &
  (list(
    scale_x_continuous(labels = scales::label_number(suffix = "K", scale = 1e-3, accuracy = 1),
    expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)),
    discreteLollipopTheme))


p <- wrap_elements(full = figAB, ignore_tag = TRUE) +
  wrap_elements(panel = figABLegend, ignore_tag = TRUE, clip = FALSE) +
  wrap_elements(figC) +
  wrap_elements(figD) +
  plot_annotation(tag_levels = list(c("c", "d"))) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig4_art_count", devices = c("png", "pdf"), gwidth = 8.2, gheight = 9)

