source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/InVitro_umap_labeledCluster_noPlotLabels.rds")
figB <- readRDS("../outs/rds/InVitro_umap_haystackOut.rds")
figC <- readRDS("../outs/rds/inVitro_discreteAbsolute_matched.rds")
figD <- readRDS("../outs/rds/inVitro_discreteHivOnly_matched.rds")
figE <- readRDS("../outs/rds/invitro_lollipop_activeLate.rds")
figF <- readRDS("../outs/rds/invitro_peaks_volcano_activatedHIV.rds")

figTopPeaks <- readRDS("../outs/rds/invitro_markerTest_activated_topPeaks.rds")
figCCR5 <- readRDS("../outs/rds/invitro_activated_ccr5.rds")

figMotifDown <- readRDS("../outs/rds/invitro_motifsUp_activatedHIVneg.rds")
figMotifUp <- readRDS("../outs/rds/invitro_motifsUp_activatedHIVpos.rds")

a_legend <- as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1, title = "Legend for (A)")) +
    theme(legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      legend.justification = "left",
      legend.title = element_text(size = 8, color = "#000000")) +
    labs(color = "Legend for (A)")))

b_legend <- as_ggplot(get_legend(figB + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1, title = "Legend for (B)")) +
    theme(
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      legend.justification = "left",
      legend.title = element_text(size = 8, color = "#000000"))))

figA <- figA + umapPlotThemeNoLeg 
figB <- figB + umapPlotThemeNoLeg + theme(axis.title.y = element_blank())

figAB <- figA + figB & subplotTheme
figAB[[2]] <- figAB[[2]] + plot_layout(tag_level = "new")

figABLegend <- (a_legend + subplotTheme + theme(plot.margin = unit(c(50, 0, 0, 0), "pt"))) /
  (b_legend + subplotTheme + theme(plot.margin = unit(c(70, 0, -30, 0), "pt")))

figCCR5[[1]] <- figCCR5[[1]] + 
  subplotTheme +
  ggtitle(" ") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

figCCR5[[2]] <- figCCR5[[2]] + theme(plot.title = element_blank()) + ggtitle("")
figCCR5 <- wrap_plots(figCCR5, ncol = 1, heights = c(3,2,1)) & subplotTheme

figTopPeaks <- figTopPeaks + 
  subplotTheme + 
  theme(legend.position = "right",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    legend.direction = "vertical",
    legend.text.align = 1,
    legend.box.margin=margin(0,0,0,0),
    plot.margin = margin(0,10,30,10))

figMotifs <- figMotifDown + figMotifUp & subplotTheme + theme(plot.margin = margin(-20, 0, 0, 0))

layout <- c(
  area(1, 1, 2, 2), #a
  area(1, 3, 2, 4), #b
  area(1, 5, 2, 6), #legend
  area(3, 1, 5, 3), #c
  area(3, 4, 5, 6), #d
  area(1, 7, 6, 9), #e
  area(6, 1, 7, 2), #f
  area(6, 3, 8, 8), #topPeaks
  area(8, 1, 10, 5), #track
  area(9, 6, 10, 9) #motifs
)

p <- wrap_elements(figA + subplotTheme) +
  wrap_elements(figB + subplotTheme) +
  wrap_elements(figABLegend) + 
  wrap_elements(figC + subplotTheme) +
  wrap_elements(figD + subplotTheme + theme(plot.margin = unit(c(0, 0, 0, 10), "pt"))) +
  wrap_elements(figE + subplotTheme) +
  wrap_elements(figF + subplotTheme + theme(legend.margin=margin(0,0,5,0), legend.box.margin=margin(-5,-10,-10,-10))) +
  wrap_elements(figTopPeaks) +
  wrap_elements(figCCR5) +
  wrap_elements(figMotifs + plot_annotation(theme = theme(plot.title = element_text(margin = margin(-50,0,150,0))))) +
  plot_annotation(tag_levels = list(c(letters[1:2], "", letters[3:10]))) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig1", devices = c("png", "pdf"), gwidth = 8, gheight = 9)

