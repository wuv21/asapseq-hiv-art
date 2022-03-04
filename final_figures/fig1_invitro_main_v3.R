source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/InVitro_umap_labeledCluster_withPlotLabels.rds")
figB <- readRDS("../outs/rds/InVitro_umap_haystackOut.rds")
figC <- readRDS("../outs/rds/inVitro_discreteAbsolute_matched.rds")
figD <- readRDS("../outs/rds/inVitro_discreteHivOnly_matched.rds")
figE <- readRDS("../outs/rds/invitro_lollipop_activeLate.rds")

a_legend <- as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2, title = "Legend\nfor (a)")) +
    theme(legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      legend.justification = "left",
      legend.title = element_text(size = 8, color = "#000000"))))

b_legend <- as_ggplot(get_legend(figB + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1, title = "Legend\nfor (b)")) +
    theme(
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      legend.justification = "left",
      legend.title = element_text(size = 8, color = "#000000"))))

figA <- figA + umapPlotThemeNoLeg +
  subplotTheme

figB <- figB + umapPlotThemeNoLeg + 
  subplotTheme +
  theme(axis.title.y = element_text(color = "#FFFFFF00"),
    plot.margin = unit(c(0, 0, 0, 0), "pt"))

figAB <- figA + figB + plot_annotation(tag_levels = "a")
figAB[[2]] <- figAB[[2]] + plot_layout(tag_level = "keep")


figABLegend <- (a_legend + subplotTheme + theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))) /
  (b_legend + subplotTheme +
  theme(plot.margin = unit(c(-50, 0, 0, 10), "pt")))

layout <- c(
  area(1, 1, 2, 4), #ab
  area(1, 5, 2, 8), #legend
  area(3, 1, 5, 3), #c
  area(3, 4, 5, 6), #d
  area(2, 7, 5, 9) #e
)

p <- wrap_elements(full = figAB, ignore_tag = TRUE) +
  wrap_elements(figABLegend, ignore_tag = TRUE) + 
  wrap_elements(figC + subplotTheme) +
  wrap_elements(figD + subplotTheme + theme(plot.margin = unit(c(0, 0, 0, 10), "pt"))) +
  wrap_elements(figE + subplotTheme + theme(plot.margin = unit(c(-20, 0, 0, 0), "pt"))) +
  plot_annotation(tag_levels = list(c(letters[3:5]))) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig1", devices = c("png", "pdf"), gwidth = 8, gheight = 6)

