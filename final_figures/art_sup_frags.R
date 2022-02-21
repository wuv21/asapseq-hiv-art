source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/art_upset_cbc.rds")
figB <- readRDS("../outs/rds/art_frags2.rds")

figC <- readRDS("../outs/rds/art_qc_haystack_frags.rds")
figD <- readRDS("../outs/rds/art_qc_haystack_TSS.rds")

figE <- readRDS("../outs/rds/art_umap_unlabeledCluster.rds")
figE <- figE + guides(color = guide_legend(ncol = 5, direction = "horizontal", override.aes = list(size = 3))) + umapPlotThemeLeg

layout <- c(
  area(1, 1, 4, 6), #a
  area(5, 1, 14, 12), #b
  area(1, 7, 4, 12)
  # area(1, 6, 2, 7), #c
  # area(3, 6, 4, 7), #d
  # area(1, 8, 4, 12) #e
)

p <- wrap_elements(figA + subplotTheme) +
  wrap_elements(figB + subplotTheme + theme(plot.margin = margin(0,0,0,-100))) + 
  # wrap_elements(figC + subplotTheme) +
  # wrap_elements(figD + subplotTheme) +
  wrap_elements(figE + subplotTheme) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = c("a")) &
  theme(plot.tag = element_text(size = 12),
    plot.tag.position = c(0, 1))

saveFinalFigure(fn = "supfig_art_frags", devices = c("png", "pdf"), plot = p, gwidth = 8.25, gheight = 10.5)

