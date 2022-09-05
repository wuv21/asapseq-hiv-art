source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/art_upset_cbc.rds")
figB <- readRDS("../outs/rds/art_fragCoverage_byAnnot.rds")
figC <- readRDS("../outs/rds/art_umap_unannotCluster_withPlotLabels.rds")


figB <- figB +
  theme(
    panel.spacing = unit(1, "pt"),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(t = 0, b = -10),
    legend.position = "top",
    strip.text = element_blank())

figC <- figC + 
  guides(color = guide_legend(ncol = 5, direction = "horizontal", override.aes = list(size = 3))) +
  umapPlotThemeLeg

layout <- c(
  area(1, 1, 6, 6), #a
  area(7, 1, 8, 12), #b
  area(1, 7, 6, 12) #c
)

p <- wrap_elements(figA + subplotTheme) +
  wrap_elements(figB + subplotTheme) +
  wrap_elements(figC + subplotTheme) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = c("a")) &
  theme(plot.tag = element_text(size = 12),
    plot.tag.position = c(0, 1))

saveFinalFigure(fn = "supfig_art_frags", devices = c("png", "pdf"), plot = p, gwidth = 8.25, gheight = 7)

