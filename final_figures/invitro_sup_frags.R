source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/invitro_upset_cbc.rds")
figA[[1]] <- figA[[1]] + ggtitle("a") +
  theme(plot.title.position = "plot",
    plot.title = element_text(size = 12, margin = margin(0,0,0,0)))

figB <- readRDS("../outs/rds/invitro_frags.rds")
figB[[1]] <- figB[[1]] + ggtitle("b")
figB[[2]] <- figB[[2]] + ggtitle("c")
figB[[3]] <- figB[[3]] + ggtitle("d")

figE <- readRDS("../outs/rds/invitro_qc_haystack_frags.rds") + ggtitle("e")
figF <- readRDS("../outs/rds/invitro_qc_haystack_TSS.rds") + ggtitle("f")

figG <- readRDS("../outs/rds/InVitro_umap_unlabeledCluster.rds") + ggtitle("g")
figG <- figG + umapPlotThemeLeg

layout <- c(
  area(1, 1, 3, 4), #a
  area(1, 5, 9, 8), #bcd
  area(4, 1, 5, 2), #e
  area(4, 3, 5, 4), #f
  area(6, 1, 9, 4) #g
)

p <- wrap_elements(figA + subplotTheme) +
  wrap_elements(figB & subplotTheme) + 
  wrap_elements(figE + subplotTheme) +
  wrap_elements(figF + subplotTheme) +
  wrap_elements(figG + subplotTheme) +
  plot_layout(design = layout)

ggsave(filename = "../outs/pdf/supfig_invitro_frags.pdf", device = cairo_pdf, plot = p, width = 8, height = 8, units = "in")

