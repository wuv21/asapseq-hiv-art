source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/art_umap_labeledCluster_withPlotLabels.rds")
figB <- readRDS("../outs/rds/art_umap_haystackout.rds")
figC <- readRDS("../outs/rds/art_discreteAbsolute_matched_byDonor.rds")
figD <- readRDS("../outs/rds/art_discreteHivOnly_matched_byDonor.rds")

a_legend <- as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), nrow = 10,
      title = "Legend for (a)",
      title.position = "top",
      title.hjust = 0)) +
    theme(legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      legend.justification = "left",
      legend.text = element_text(margin = margin(r = 12, unit = "pt"), size = 6),
      legend.title = element_text(size = 8, color = "#000000", hjust = 0))))

b_legend <- as_ggplot(get_legend(figB + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1, title = "Legend for (b)",
      title.position = "top",
      title.hjust = 0)) +
    theme(
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      legend.justification = "left",
      legend.text = element_text(margin = margin(r = 12, unit = "pt"), size = 6),
      legend.title = element_text(size = 8, color = "#000000"))))


figABLegend <- a_legend / b_legend & subplotTheme 

figABLegend[[1]] <- figABLegend[[1]] +
  theme(plot.margin = margin(t = 20, l = 10))

figABLegend[[2]] <- figABLegend[[2]] +
  theme(plot.margin = margin(t = -50, b = 75, l = 10))

figA <- figA + umapPlotThemeNoLeg +
  subplotTheme

figB <- figB + umapPlotThemeNoLeg + 
  subplotTheme +
  theme(axis.title.y = element_text(color = "#FFFFFF00"),
    plot.margin = unit(c(0, 0, 0, 0), "pt"))

figAB <- figA + figB + plot_annotation(tag_levels = "a")

layout <- c(
  area(1, 1, 2, 4), #ab
  area(1, 5, 3, 9), #legend
  area(3, 1, 5, 9), #c
  area(6, 1, 7, 9)
  # area(6, 1, 7, 4), #d
  # area(6, 5, 7, 7), #e
  # area(6, 8, 7, 9) #f
)

discreteLollipopTheme <- subplotTheme + theme(
  axis.title.y = element_blank(),
  axis.title.x = element_blank(),
  axis.text.x = element_text(size = 5),
  plot.subtitle = element_text(size = 7, hjust = 0.5, margin = margin(b = -1)),
  plot.title.position = "panel",
  plot.margin = margin(l = -20, t = -20, r = 20)
)

decreasePointAndText <- function(x) {
  x$layers[[3]]$aes_params$size <- 5 / ggplot2:::.pt
  x$layers[[2]]$aes_params$size <- 0.75
  
  return(x)
}

figC <- lapply(figC, decreasePointAndText)
for (i in c(2:length(figC))) {
  figC[[i]] <- figC[[i]] + theme(axis.text.y = element_blank())
}
figC[[1]] <- figC[[1]] + theme(axis.text.y = element_text(angle = 45, hjust = 1, size = 5))

figC <- wrap_plots(figC, nrow = 1) &
  (list(
    scale_x_continuous(labels = scales::label_number(suffix = "K", scale = 1e-3, accuracy = 1),
    expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)),
    discreteLollipopTheme))

figD <- lapply(figD, decreasePointAndText)
for (i in c(2:length(figD))) {
  figD[[i]] <- figD[[i]] + theme(axis.text.y = element_blank())
}
figD[[1]] <- figD[[1]] + theme(axis.text.y = element_text(angle = 45, hjust = 1, size = 5))

figD <- wrap_plots(figD, nrow = 1) & discreteLollipopTheme
# figD[[4]] <- figD[[4]] + plot_layout(tag_level = "new")

p <- wrap_elements(full = figAB, ignore_tag = TRUE) +
  wrap_elements(panel = figABLegend, ignore_tag = TRUE, clip = FALSE) +
  figC +
  (figD + theme(plot.margin = margin(t = 5))) +
  # (wrap_plots(figD[c(1:3)]) & discreteLollipopTheme) +
  # (wrap_plots(figD[c(4:5)]) & discreteLollipopTheme) +
  # (wrap_plots(figD[c(6)]) & discreteLollipopTheme) +
  plot_annotation(tag_levels = list(c("c", "", "", "", "", "d", "", "", "e", "", "f"))) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig3_art", devices = c("png", "pdf"), gwidth = 8.2, gheight = 9)

