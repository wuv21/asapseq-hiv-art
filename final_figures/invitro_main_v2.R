source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/InVitro_umap_labeledCluster_noPlotLabels.rds") + ggtitle("a")
figB <- readRDS("../outs/rds/InVitro_umap_haystackOut.rds") + ggtitle("b")
figC <- readRDS("../outs/rds/inVitro_discreteAbsolute_matched.rds") + ggtitle("c")
figD <- readRDS("../outs/rds/inVitro_discreteHivOnly_matched.rds") + ggtitle("d")
figE <- readRDS("../outs/rds/invitro_lollipop_activeLate.rds") + ggtitle("e")
figF <- readRDS("../outs/rds/invitro_peaks_volcano_activatedHIV.rds") + ggtitle("f")

figTopPeaks <- readRDS("../outs/rds/invitro_markerTest_activated_topPeaks.rds") + ggtitle("g")
figCCR5 <- readRDS("../outs/rds/invitro_activated_ccr5.rds")

figMotifDown <- readRDS("../outs/rds/invitro_motifsUp_activatedHIVneg.rds") + ggtitle("i")
figMotifUp <- readRDS("../outs/rds/invitro_motifsUp_activatedHIVpos.rds") + ggtitle("j")

legIndicator <- readPNG("asapseq_legend_photos/1 - invitro.png", native = TRUE)

a_legend <- as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1)) +
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

xInit <- -5
yInit <- -8.5
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

figABLegend <- (a_legend + subplotTheme + theme(plot.margin = unit(c(50, 0, 0, 0), "pt"))) /
  (b_legend + subplotTheme + theme(plot.margin = unit(c(70, 0, -30, 0), "pt")))

figCCR5[[1]] <- figCCR5[[1]] + 
  subplotTheme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  ggtitle("h")
figCCR5[[2]] <- figCCR5[[2]] + theme(plot.title = element_blank()) + ggtitle("")
figCCR5 <- wrap_plots(figCCR5, ncol = 1, heights = c(3,2,1)) & subplotTheme

figTopPeaks <- figTopPeaks + 
  subplotTheme + 
  theme(legend.position = "right",
    legend.direction = "vertical",
    legend.text.align = 1,
    legend.box.margin=margin(0,0,0,0),
    plot.margin = margin(0,10,30,10))

figMotifs <- figMotifDown + figMotifUp & subplotTheme + theme(plot.margin = margin(-20, 0, 0, 0))

layout <- c(
  area(1, 1, 2, 4), #ab
  area(1, 5, 2, 6), #legend
  area(3, 1, 5, 3), #c
  area(3, 4, 5, 6), #d
  area(1, 7, 6, 9), #e
  area(8, 1, 10, 5), #track
  area(6, 1, 7, 2), #f
  area(6, 3, 8, 8), #topPeaks
  area(9, 6, 10, 9) #motifs
)

p <- wrap_elements(figA + figB & subplotTheme) +
  wrap_elements(figABLegend) + 
  inset_element(p = legIndicator, left = 0.8, right = 1, top = 1, bottom = 0.8, clip = FALSE, align_to = "full", ignore_tag = TRUE) + 
  wrap_elements(figC + subplotTheme) +
  wrap_elements(figD + subplotTheme + theme(plot.margin = unit(c(0, 0, 0, 10), "pt"))) +
  wrap_elements(figE + subplotTheme) +
  wrap_elements(figCCR5) +
  wrap_elements(figF + subplotTheme + theme(legend.margin=margin(0,0,5,0), legend.box.margin=margin(-5,-10,-10,-10))) +
  wrap_elements(figTopPeaks) +
  wrap_elements(figMotifs) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig2", devices = c("png", "pdf"), gwidth = 8, gheight = 9)

