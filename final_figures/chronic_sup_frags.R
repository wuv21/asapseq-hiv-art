source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/chronic_upset_cbc.rds")

figB <- readRDS("../outs/rds/chronic_frags.rds")
figB <- figB & subplotTheme
figB[[1]] <- figB[[1]] + theme(strip.background = element_rect(fill = "#DDDDDD60"))
figB[[2]] <- figB[[2]] + theme(plot.margin = margin(0, 0, 0.5, 0, unit = "lines"))

figC <- readRDS("../outs/rds/chronic_qc_haystack_frags.rds")
figD <- readRDS("../outs/rds/chronic_qc_haystack_TSS.rds")

figE <- readRDS("../outs/rds/chronic_umap_unlabeledCluster.rds")
figE <- figE + guides(color = guide_legend(ncol = 5, direction = "horizontal", override.aes = list(size = 3))) + umapPlotThemeLeg

layout <- c(
  area(1, 1, 3, 6), #a
  area(1, 7, 9, 12), #bcd
  area(4, 1, 5, 3), #e
  area(4, 4, 5, 6), #f
  area(6, 1, 9, 6) #g
)

p <- wrap_elements(figA + subplotTheme) +
  wrap_elements(figB) + 
  wrap_elements(figC + subplotTheme) +
  wrap_elements(figD + subplotTheme) +
  wrap_elements(figE + subplotTheme) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = c("a")) &
  theme(plot.tag = element_text(size = 12),
    plot.tag.position = c(0, 1))

saveFinalFigure(fn = "supfig_chronic_frags", devices = c("png", "pdf"), plot = p, gwidth = 8, gheight = 9)

