source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/chronic_umap_labeledCluster_noPlotLabels.rds") + ggtitle("a")
figB <- readRDS("../outs/rds/chronic_umap_haystackOut.rds") + ggtitle("b")
figC <- readRDS("../outs/rds/chronic_discreteAbsolute_matched.rds") + ggtitle("c")
figD <- readRDS("../outs/rds/chronic_discreteHivOnly_matched.rds") + ggtitle("d")
figE <- readRDS("../outs/rds/chronic_lollipop.rds") + ggtitle("e")
figF <- readRDS("../outs/rds/chronic_lollipop_tfh.rds") + ggtitle("f")
figG <- readRDS("../outs/rds/chronic_volcano_motifs_chromVAR.rds") + ggtitle("g")
figH <- readRDS("../outs/rds/chronic_chromVAR_motifsUp_HIVneg.rds") + ggtitle("h")
figI <- readRDS("../outs/rds/chronic_chromVAR_motifsUp_HIVpos.rds") + ggtitle("i")
legIndicator <- readPNG("asapseq_legend_photos/2 - chronic.png", native = TRUE)

a_legend <- as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
    theme(legend.margin=margin(0,0,0,0),
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

arrow_1 <- aes(x = -9.25, y = -7, xend = -9.25, yend = -4.5)
arrow_2 <- aes(x = -9.25, y = -7, xend = -6.75, yend = -7)

figA <- figA + 
  umapCleanPlotTheme +
  geom_segment(arrow_1, arrow = umapArrowSettings, size = umapArrowSize) +
  geom_segment(arrow_2, arrow = umapArrowSettings, size = umapArrowSize)

figB <- figB +
  umapCleanPlotTheme +
  geom_segment(arrow_1, arrow = umapArrowSettings, size = umapArrowSize, color = "#FFFFFF") +
  geom_segment(arrow_2, arrow = umapArrowSettings, size = umapArrowSize, color = "#FFFFFF") +
  theme(axis.title = element_text(color = "#FFFFFF"))

figABLegend <- (a_legend + subplotTheme + theme(plot.margin = unit(c(10, 20, 0, 10), "pt"))) / 
  (b_legend + subplotTheme + theme(plot.margin = unit(c(10, 20, 0, 10), "pt")))
  
motifTheme <- subplotTheme + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(size = 5.5, angle = 45, vjust = 1, hjust = 1, margin = margin(0, 0, 0, 0)),
  plot.title = element_text(margin = margin(0, 0, 5, 0)))

figEF <- ((figE + subplotTheme + theme(plot.margin = margin(0, 0, 0, 0))) / 
  (figF + subplotTheme) + theme(plot.margin = margin(0,0,0,0))) + plot_layout(heights = c(3, 1)) + theme(plot.margin = margin(0,0,-40,0))
  

figHI <- ((figH + motifTheme) + plot_spacer() + theme(plot.background = element_rect(fill = "transparent"))) / 
  (figI + motifTheme) + plot_layout(widths = c(3, 1))

layout <- c(
  area(1, 1, 2, 4), #ab
  area(1, 5, 2, 9), #legend
  area(3, 1, 5, 3), #c
  area(3, 4, 5, 6), #d
  area(3, 7, 6, 9), #ef
  area(6, 1, 7, 2), #g
  area(6, 3, 7, 9) #hi
)

p <- wrap_elements(figA + figB & subplotTheme) +
  wrap_elements(figABLegend) +
  inset_element(p = legIndicator, left = 0.9, right = 1, top = 1, bottom = 0.5, clip = FALSE, align_to = "full") +
  wrap_elements(figC + subplotTheme) +
  wrap_elements(figD + subplotTheme) +
  wrap_elements(figEF & subplotTheme + theme(plot.margin = margin(0, 0, 10, 0))) +
  wrap_elements(figG + subplotTheme + theme(legend.margin=margin(10,0,0,0), legend.box.margin=margin(-20,-20,-5,-20))) +
  wrap_elements(figHI & theme(plot.background = element_rect(fill = "transparent", color = NA), plot.margin = margin(0, 0, -10, 0))) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig3", devices = c("png", "pdf"), gwidth = 8, gheight = 8)

