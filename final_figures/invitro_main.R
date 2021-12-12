# library(ggplot2)
# library(patchwork)

figA <- readRDS("outs/rds/InVitro_umap_labeledCluster_noPlotLabels.rds") + ggtitle("a")
figB <- readRDS("outs/rds/InVitro_umap_haystackOut.rds") + ggtitle("b")
figC <- readRDS("outs/rds/inVitro_discreteAbsolute_matched.rds") + ggtitle("c")
figD <- readRDS("outs/rds/inVitro_discreteHivOnly_matched.rds") + ggtitle("d")
figE <- readRDS("outs/rds/invitro_lollipop.rds") + ggtitle("e")
figF <- readRDS("outs/rds/invitro_peaks_volcano_activatedHIV.rds") + ggtitle("f")
figG <- readRDS("outs/rds/invitro_motifsUp_activatedHIVneg.rds") + ggtitle("g")
figH <- readRDS("outs/rds/invitro_motifsUp_activatedHIVpos.rds") + ggtitle("h")

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

umapPlotTheme <- theme(
  legend.position = "none",
  axis.ticks = element_blank(),
  axis.text = element_blank()
)

figA <- figA + umapPlotTheme
figB <- figB + umapPlotTheme

layout <- c(
  area(1, 1, 3, 3), #a
  area(1, 4, 3, 6), #b
  area(1, 7, 3, 9), #legend
  area(4, 1, 6, 3), #c
  area(4, 4, 6, 6), #d
  area(3, 7, 8, 9), #e
  area(7, 1, 8, 2), #f
  area(7, 3, 8, 4), #g
  area(7, 5, 8, 6) #h
)

subplotTheme <- theme(
  plot.title.position = "plot",
  plot.title = element_text(size = 12, margin = margin(0,0,0,0)),
  plot.margin = unit(c(0,0,0,0), "pt"),
  text = element_text(family = "Arial"),
  rect = element_rect(fill = "transparent", colour = NULL)
)

p <- wrap_elements(figA + subplotTheme) +
  wrap_elements(figB + subplotTheme) + 
  wrap_elements((a_legend + subplotTheme + theme(plot.margin = unit(c(10, 0, 0, 30), "pt"))) / 
      (b_legend + subplotTheme + theme(plot.margin = unit(c(0, 0, 40, 30), "pt")))) +
  wrap_elements(figC + subplotTheme) +
  wrap_elements(figD + subplotTheme + theme(plot.margin = unit(c(0, 0, 0, 10), "pt"))) +
  wrap_elements(figE + subplotTheme) +
  wrap_elements(figF + subplotTheme + theme(legend.margin=margin(0,0,5,0), legend.box.margin=margin(-10,-10,-10,-10))) +
  wrap_elements(figG + subplotTheme) +
  wrap_elements(figH + subplotTheme) +
  plot_layout(design = layout)

ggsave(filename = "outs/pdf/fig2.pdf", device = cairo_pdf, plot = p, width = 8, height = 8, units = "in")

