source("defaultFigureSettings.R")

figA <- readRDS("outs/rds/chronic_umap_labeledCluster_noPlotLabels.rds") + ggtitle("a")
figB <- readRDS("outs/rds/chronic_umap_haystackOut.rds") + ggtitle("b")
figC <- readRDS("outs/rds/chronic_discreteAbsolute_matched.rds") + ggtitle("c")
figD <- readRDS("outs/rds/chronic_discreteHivOnly_matched.rds") + ggtitle("d")
figE <- readRDS("outs/rds/chronic_lollipop_tfh.rds") + ggtitle("e")
figF <- readRDS("outs/rds/chronic_volcano_motifs_chromVAR.rds") + ggtitle("f")
figG <- readRDS("outs/rds/chronic_chromVAR_motifsUp_HIVneg.rds") + ggtitle("g")
figH <- readRDS("outs/rds/chronic_chromVAR_motifsUp_HIVpos.rds") + ggtitle("h")
legIndicator <- readPNG("final_figures/asapseq_legend_photos/2 - chronic.png", native = TRUE)

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
  area(3, 7, 8, 9), #e
  area(7, 1, 8, 2), #f
  area(7, 3, 9, 7) #g and h
)

figGH <- (figG + labs(x = "")) / figH & 
  (subplotTheme + 
      theme(plot.margin = unit(c(0, 7, 0, 0), "pt"),
        plot.title = element_text(margin = margin(0, 0, 10, 0)),
        axis.title.x = element_text(margin = margin(0,0,0,0)),
        axis.text.x = element_text(size = 5.5, angle = 90, vjust = 0.5, hjust = 1, margin = margin(0, 0, 0, 0))))

p <- wrap_elements(figA + subplotTheme) +
  wrap_elements(figB + subplotTheme) + 
  wrap_elements((a_legend + subplotTheme + theme(plot.margin = unit(c(60, 0, 30, 50), "pt"))) / 
      (b_legend + subplotTheme + theme(plot.margin = unit(c(60, 0, 40, 50), "pt")))) +
  inset_element(p = legIndicator, left = 0.7, right = 1, top = 1, bottom = 0.8, clip = FALSE, align_to = "full") +
  wrap_elements(figC + subplotTheme) +
  wrap_elements(figD + subplotTheme) +
  wrap_elements(figE + subplotTheme + theme(plot.margin = unit(c(40, 0, 0, 0), "pt"))) +
  wrap_elements(figF + subplotTheme + theme(legend.margin=margin(3,0,0,0), legend.box.margin=margin(-10,-10,-10,-10))) +
  wrap_elements(figGH) +
  plot_layout(design = layout)

ggsave(filename = "outs/pdf/fig3.pdf", device = cairo_pdf, plot = p, width = 8, height = 8.5, units = "in")

