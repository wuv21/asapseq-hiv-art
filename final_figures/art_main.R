source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/art_umap_labeledCluster_noPlotLabels.rds") + ggtitle("a")
figB <- readRDS("../outs/rds/art_umap_haystackOut.rds") + ggtitle("b")
figC <- readRDS("../outs/rds/art_discreteAbsolute_matched.rds") + ggtitle("c")
figD <- readRDS("../outs/rds/art_discreteHivOnly_matched.rds") + ggtitle("d")
figE <- readRDS("../outs/rds/art_lollipop.rds") + ggtitle("e")
figF <- readRDS("../outs/rds/art_volcano_motifs_chromVAR.rds") + ggtitle("f")
legIndicator <- readPNG("asapseq_legend_photos/3 - art.png", native = TRUE)

a_legend <- as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
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

layout <- c(
  area(1, 1, 3, 3), #a
  area(1, 4, 3, 6), #b
  area(1, 7, 3, 9), #legend
  area(4, 1, 8, 3), #c
  area(4, 4, 6, 6), #d
  area(4, 7, 6, 9), #e
  area(7, 4, 8, 5)  #f
)

p <- wrap_elements(figA + subplotTheme) +
  wrap_elements(figB + subplotTheme) + 
  wrap_elements((a_legend + subplotTheme + theme(plot.margin = unit(c(25, 0, 15, 0), "pt"))) / 
      (b_legend + subplotTheme + theme(plot.margin = unit(c(15, 0, 30, 0), "pt")))) +
  inset_element(p = legIndicator, left = 0.7, right = 1, top = 0.6, bottom = 0, clip = FALSE, align_to = "full") +
  wrap_elements(figC + subplotTheme + 
      scale_x_continuous(labels = scales::label_number(suffix = "K", scale = 1e-3, accuracy = 1),
        expand = expansion(mult = c(0, 0.05)), limits = c(0, NA))) +
  wrap_elements(figD + subplotTheme) +
  wrap_elements(figE + subplotTheme) +
  wrap_elements(figF + subplotTheme + theme(legend.margin=margin(3,0,0,0), legend.box.margin=margin(-10,-10,-10,-10))) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig4", devices = c("png", "pdf"), gwidth = 8, gheight = 8.5)

