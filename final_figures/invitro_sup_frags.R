source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/invitro_upset_cbc.rds")
figB <- readPNG("../outs/png/p24pos.png", native = TRUE)

figC <- readRDS("../outs/rds/invitro_frags.rds")
figC <- figC & subplotTheme
figC[[2]] <- figC[[2]] + theme(plot.margin = margin(0, 0, 0.5, 0, unit = "lines"))

figD <- readRDS("../outs/rds/invitro_qc_haystack_frags.rds")
figE <- readRDS("../outs/rds/invitro_qc_haystack_TSS.rds")

figF <- readRDS("../outs/rds/InVitro_umap_unlabeledCluster.rds")
figF <- figF + guides(color = guide_legend(ncol = 5, direction = "horizontal", override.aes = list(size = 3))) + umapPlotThemeLeg

layout <- c(
  area(1, 1, 3, 4), #a
  area(1, 5, 3, 6), #flow
  area(1, 7, 9, 12), #bcd
  area(4, 1, 5, 3), #e
  area(4, 4, 5, 6), #f
  area(6, 1, 9, 6) #g
)

p <- wrap_elements(figA + subplotTheme) +
  (wrap_elements(full = figB) + subplotTheme) +
  wrap_elements(figC) + 
  wrap_elements(figD + subplotTheme) +
  wrap_elements(figE + subplotTheme) +
  wrap_elements(figF + subplotTheme) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = c("a")) &
  theme(plot.tag = element_text(size = 12),
    plot.tag.position = c(0, 1))

saveFinalFigure(plot = p, fn = "supfig_invitro_frags", devices = c("png", "pdf"), gwidth = 8, gheight = 9)

