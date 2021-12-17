source("defaultFigureSettings.R")

figA <- readRDS("outs/rds/art_umap_labeledCluster_noPlotLabels.rds") + ggtitle("a")
figB <- readRDS("outs/rds/art_umap_haystackOut.rds") + ggtitle("b")
figC <- readRDS("outs/rds/art_discreteAbsolute_matched.rds") + ggtitle("c")
figD <- readRDS("outs/rds/art_discreteHivOnly_matched.rds") + ggtitle("d")
figE <- readRDS("outs/rds/art_lollipop.rds") + ggtitle("e")
figF <- readRDS("outs/rds/art_volcano_motifs_chromVAR.rds") + ggtitle("f")

a_legend <- as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1)) +
    theme(legend.margin=margin(0,0,0,0), legend.box.margin=margin(-10,-10,-10,-10))))

b_legend <- as_ggplot(get_legend(figB + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1)) +
    theme(
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10))))

figA <- figA + umapPlotThemeNoLeg
figB <- figB + umapPlotThemeNoLeg

layout <- c(
  area(1, 1, 3, 3), #a
  area(1, 4, 3, 6), #b
  area(1, 7, 3, 9), #legend
  area(4, 1, 6, 3), #c
  area(4, 4, 6, 6), #d
  area(4, 7, 8, 9), #e
  area(7, 1, 8, 2)  #f
  # area(7, 3, 9, 7) #g and h
)

p <- wrap_elements(figA + subplotTheme) +
  wrap_elements(figB + subplotTheme) + 
  wrap_elements((a_legend + subplotTheme + theme(plot.margin = unit(c(25, 0, 15, 0), "pt"))) / 
      (b_legend + subplotTheme + theme(plot.margin = unit(c(15, 0, 30, 0), "pt")))) +
  wrap_elements(figC + subplotTheme + theme(plot.margin = unit(c(0, 12, 0, 0), "pt"))) +
  wrap_elements(figD + subplotTheme) +
  wrap_elements(figE + subplotTheme) +
  wrap_elements(figF + subplotTheme + theme(legend.margin=margin(3,0,0,0), legend.box.margin=margin(-10,-10,-10,-10))) +
  plot_layout(design = layout)

ggsave(filename = "outs/pdf/fig4.pdf", device = cairo_pdf, plot = p, width = 8, height = 7.5, units = "in")

