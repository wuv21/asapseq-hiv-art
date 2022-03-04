source("defaultFigureSettings.R")

figG <- readRDS("../outs/rds/chronic_volcano_motifs_chromVAR.rds")
figH <- readRDS("../outs/rds/chronic_chromVAR_motifsUp_HIVneg.rds")
figI <- readRDS("../outs/rds/chronic_chromVAR_motifsUp_HIVpos.rds")
figJ <- readRDS("../outs/rds/chronic_models.rds")

figG <- figG + 
  subplotTheme + 
  theme(legend.margin=margin(10,0,0,0),
    legend.box.margin=margin(-20,-20,-5,-20))

motifTheme <- subplotTheme + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(size = 5.5, angle = 45, vjust = 1, hjust = 1, margin = margin(0, 0, 0, 0)),
  plot.title = element_text(margin = margin(0, 0, 5, 0)))


figHI <- ((figH + motifTheme) + plot_spacer() + theme(plot.background = element_rect(fill = "transparent"))) / 
  (figI + motifTheme) + plot_layout(widths = c(3, 1))

figHI <- figHI +
  plot_annotation(tag_levels = list(c("b", "c")),
    theme = theme(
      plot.margin = unit(c(0, -25, -5, 7.5), "pt")
    ))

layout <- c(
  area(1, 1, 2, 2), #g
  area(1, 3, 2, 6), #hi
  area(1, 7, 2, 9) #j
)

p <- wrap_elements(figG) +
  wrap_elements(full = figHI, ignore_tag = TRUE) +
  wrap_elements(panel = figJ + subplotTheme) +
  plot_annotation(tag_levels = list(c("a", "d"))) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig4", devices = c("png", "pdf"), gwidth = 8, gheight = 2.5)

