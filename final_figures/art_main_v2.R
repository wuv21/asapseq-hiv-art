source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/art_umap_labeledCluster_noPlotLabels.rds") + ggtitle("a")
figB <- readRDS("../outs/rds/art_umap_haystackOut.rds") + ggtitle("b")
figC <- readRDS("../outs/rds/art_discreteAbsolute_matched.rds") + ggtitle("c")
figD <- readRDS("../outs/rds/art_discreteHivOnly_matched.rds") + ggtitle("d")
figE <- readRDS("../outs/rds/art_lollipop.rds") + ggtitle("e")
figF <- readRDS("../outs/rds/art_volcano_motifs_chromVAR.rds") + ggtitle("f")
legIndicator <- readPNG("asapseq_legend_photos/3 - art.png", native = TRUE)

a_legend <- as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 3)) +
    theme(
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      legend.justification = "left")))

b_legend <- as_ggplot(get_legend(figB + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1)) +
    theme(
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      legend.justification = "left")))

figA <- figA + umapPlotThemeNoLeg
figB <- figB + umapPlotThemeNoLeg

xInit <- -7.5
yInit <- -9
adj <- 2.5
arrow_1 <- aes(x = xInit, y = yInit, xend = xInit, yend = yInit + adj)
arrow_2 <- aes(x = xInit, y = yInit, xend = xInit + adj, yend = yInit)

figA <- figA + 
  umapCleanPlotTheme +
  geom_segment(arrow_1, arrow = umapArrowSettings, size = umapArrowSize) +
  geom_segment(arrow_2, arrow = umapArrowSettings, size = umapArrowSize)

figB <- figB +
  umapCleanPlotTheme +
  geom_segment(arrow_1, arrow = umapArrowSettings, size = umapArrowSize, color = "#FFFFFF") +
  geom_segment(arrow_2, arrow = umapArrowSettings, size = umapArrowSize, color = "#FFFFFF") +
  theme(axis.title = element_text(color = "#FFFFFF"))

figABLegend <- (a_legend + subplotTheme + theme(plot.margin = unit(c(25, 0, 15, 0), "pt"))) / 
  (b_legend + subplotTheme + theme(plot.margin = unit(c(15, 0, 0, 0), "pt")))

layout <- c(
  area(1, 1, 2, 4), #ab
  area(1, 5, 2, 9), #legend
  area(3, 1, 7, 3), #c
  area(3, 4, 5, 6), #d
  area(3, 7, 5, 9), #e
  area(6, 4, 7, 5)  #f
)

p <- wrap_elements(figA + figB & subplotTheme) +
  wrap_elements(figABLegend) +
  inset_element(p = legIndicator, left = 0.9, right = 1, top = 0.6, bottom = 0, clip = FALSE, align_to = "full") +
  wrap_elements(figC + subplotTheme + 
      scale_x_continuous(labels = scales::label_number(suffix = "K", scale = 1e-3, accuracy = 1),
        expand = expansion(mult = c(0, 0.05)), limits = c(0, NA))) +
  wrap_elements(figD + subplotTheme) +
  wrap_elements(figE + subplotTheme) +
  wrap_elements(figF + subplotTheme + theme(legend.margin=margin(3,0,0,0), legend.box.margin=margin(-10,-10,-10,-10))) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig4", devices = c("png", "pdf"), gwidth = 8, gheight = 7)

