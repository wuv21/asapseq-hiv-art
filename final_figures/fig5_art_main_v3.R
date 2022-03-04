source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/art_umap_labeledCluster_withPlotLabels.rds")
figB <- readRDS("../outs/rds/art_umap_haystackOut.rds")
figC <- readRDS("../outs/rds/art_discreteAbsolute_matched.rds")
figD <- readRDS("../outs/rds/art_discreteHivOnly_matched.rds")
figE <- readRDS("../outs/rds/art_lollipop.rds")
figF <- readRDS("../outs/rds/art_volcano_motifs_chromVAR.rds")


a_legend <- as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), nrow = 9, title = "Legend\nfor (a)")) +
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


figABLegend <- a_legend / b_legend & subplotTheme 

figABLegend[[1]] <- figABLegend[[1]] +
  theme(plot.margin = margin(t = 20, l = 20))

figABLegend[[2]] <- figABLegend[[2]] +
  theme(plot.margin = margin(b = 50, l = 20))

figA <- figA + umapPlotThemeNoLeg +
  subplotTheme

figB <- figB + umapPlotThemeNoLeg + 
  subplotTheme +
  theme(axis.title.y = element_text(color = "#FFFFFF00"),
    plot.margin = unit(c(0, 0, 0, 0), "pt"))

figAB <- figA + figB + plot_annotation(tag_levels = "a")

layout <- c(
  area(1, 1, 2, 4), #ab
  area(1, 5, 3, 9), #legend
  area(3, 1, 7, 3), #c
  area(3, 4, 5, 6), #d
  area(3, 7, 5, 9), #e
  area(6, 4, 7, 5)  #f
)

p <- wrap_elements(full = figAB, ignore_tag = TRUE) +
  wrap_elements(panel = figABLegend, ignore_tag = TRUE, clip = FALSE) +
  wrap_elements(figC + subplotTheme + 
      scale_x_continuous(labels = scales::label_number(suffix = "K", scale = 1e-3, accuracy = 1),
        expand = expansion(mult = c(0, 0.05)), limits = c(0, NA))) +
  wrap_elements(figD + subplotTheme) +
  wrap_elements(figE + subplotTheme) +
  wrap_elements(figF + subplotTheme + theme(legend.margin=margin(3,0,0,0), legend.box.margin=margin(-10,-10,-10,-10))) +
  plot_annotation(tag_levels = list(c(letters[3:7]))) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig5", devices = c("png", "pdf"), gwidth = 8, gheight = 7)

