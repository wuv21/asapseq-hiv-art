source("defaultFigureSettings.R")
library(patchwork)

basePanelPlots <- c(
  "../outs/rds/invitro_basepanel_separate.rds",
  "../outs/rds/chronic_basepanel_separate.rds",
  "../outs/rds/art_basepanel_separate.rds"
)

imputePanelPlots <- c(
  "../outs/rds/invitro_imputePanel.rds",
  "../outs/rds/chronic_imputePanel.rds",
  "../outs/rds/art_imputePanel.rds"
)

finalFns <- c(
  "supfig_invitro_annotPanels",
  "supfig_chronic_annotPanels",
  "supfig_art_annotPanels"
)


for (i in c(1:length(basePanelPlots))) {
  baseFig <- readRDS(basePanelPlots[i])
  
  baseFig <- wrap_plots(baseFig)
  baseFig <- baseFig &
    theme(
      plot.margin = unit(c(0,0,2.5,0), "pt"),
      plot.title = element_text(size = 6, margin = margin(0, 0, 2, 0)),
      axis.text.y = element_text(size = 4.75),
      axis.text.x = element_text(size = 4.75)
    )
  
  # saveFinalFigure(plot = p, fn = basePanelFns[i], devices = c("png", "pdf"), gwidth = 8, gheight = ifelse(i == 3, 7, 6))

  impFig <- readRDS(imputePanelPlots[i]) + patchwork::plot_layout(ncol = 6)

  impFig <- impFig &
    theme(
      axis.title = element_blank(),
      plot.title = element_text(size = 6, margin = margin(2, 0, 0, 0)),
      legend.margin = margin(0,0,0,0),
      legend.box.margin = margin(-5,2.5,0,-7.5),
      legend.box.spacing = unit(0, "pt"),
      legend.text = element_text(size = 3, margin = margin(t = -3)),
      legend.title = element_blank(),
      legend.key.size = unit(0.3, "lines"),
      legend.key.width = unit(0.5, "lines"),
      legend.position = "top",
      legend.justification = "right",
      plot.margin = unit(c(0,5,5,5), "pt"))
  
  basePanelHeight <- ifelse(i == 3, 7, 6)
  tmp <- (wrap_elements(baseFig) / wrap_elements(impFig)) +
    plot_layout(heights = c(basePanelHeight, 11 - basePanelHeight)) +
    plot_annotation(tag_levels = "A") 
  
  saveFinalFigure(
    plot = tmp,
    fn = finalFns[i],
    devices = c("png", "pdf"),
    gwidth = 8,
    gheight = 11)

  # p <- wrap_elements(fig)
  # saveFinalFigure(plot = p, fn = basePanelFns[i], devices = c("png", "pdf"), gwidth = 8, gheight = 4)
}



