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

layout <- c(
  area(1, 1, 3, 3), #a
  area(1, 4, 3, 6), #b
  area(1, 7, 3, 9), #legend
  area(4, 1, 6, 3), #c
  area(4, 4, 6, 6), #d
  area(4, 7, 8, 9), #ef
  area(7, 1, 8, 2), #g
  area(7, 3, 9, 9) #hi
)

motifTheme <- subplotTheme + 
    theme(axis.title.x = element_blank(),
      axis.text.x = element_text(size = 5.5, angle = 90, vjust = 0.5, hjust = 1, margin = margin(0, 0, 0, 0)),
      plot.title = element_text(margin = margin(0, 0, 3, 0)))

figEF <- ((figE + subplotTheme + theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))) / 
  (figF + subplotTheme)) + plot_layout(heights = c(3, 1))
  

figHI <- ((figH + motifTheme) + plot_spacer() + theme(plot.background = element_rect(fill = "transparent"))) / 
  (figI + motifTheme) + plot_layout(widths = c(3, 1))


p <- wrap_elements(figA + subplotTheme) +
  wrap_elements(figB + subplotTheme) + 
  wrap_elements((a_legend + subplotTheme + theme(plot.margin = unit(c(80, 0, 30, 50), "pt"))) / 
      (b_legend + subplotTheme + theme(plot.margin = unit(c(60, 0, 30, 50), "pt")))) +
  inset_element(p = legIndicator, left = 0.7, right = 1, top = 1, bottom = 0.8, clip = FALSE, align_to = "full") +
  wrap_elements(figC + subplotTheme) +
  wrap_elements(figD + subplotTheme) +
  wrap_elements(figEF & subplotTheme + theme(plot.margin = margin(0, 0, 10, 0))) +
  wrap_elements(figG + subplotTheme + theme(legend.margin=margin(10,0,0,0), legend.box.margin=margin(-20,-20,-20,-20))) +
  wrap_elements(figHI & theme(plot.background = element_rect(fill = "transparent", color = NA)) & subplotTheme) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig3", devices = c("png", "pdf"), gwidth = 8, gheight = 8.5)

