source("defaultFigureSettings.R")

basePanelPlots <- c(
  "../outs/rds/invitro_basepanel.rds",
  "../outs/rds/chronic_basepanel.rds", 
  "../outs/rds/art_basepanel.rds"
)

basePanelFns <- c(
  "supfig_invitro_basepanel",
  "supfig_chronic_basepanel",
  "supfig_art_basepanel"
)

for (i in seq_along(basePanelPlots)) {
  fig <- readRDS(basePanelPlots[i])
  p <- wrap_elements(fig + subplotTheme)
  saveFinalFigure(plot = p, fn = basePanelFns[i], devices = c("png", "pdf"), gwidth = 8, gheight = 8)
}


basePanelPlots <- c(
  "../outs/rds/invitro_imputePanel.rds",
  "../outs/rds/chronic_imputePanel.rds", 
  "../outs/rds/art_imputePanel.rds"
)

basePanelFns <- c(
  "supfig_invitro_imputepanel",
  "supfig_chronic_imputepanel",
  "supfig_art_imputepanel"
)


for (i in seq_along(basePanelPlots)) {
  fig <- readRDS(basePanelPlots[i])
  p <- wrap_elements(fig)
  saveFinalFigure(plot = p, fn = basePanelFns[i], devices = c("png", "pdf"), gwidth = 8, gheight = 10)
}