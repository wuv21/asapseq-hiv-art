source("defaultFigureSettings.R")

fig <- readRDS("../outs/rds/InVitro_basepanel.rds")
p <- wrap_elements(fig + subplotTheme)

ggsave(filename = "../outs/pdf/supfig_invitro_basepanel.pdf", device = cairo_pdf, plot = p, width = 8, height = 8, units = "in")


fig <- readRDS("../outs/rds/invitro_imputePanel.rds")
p <- wrap_elements(fig + subplotTheme)

ggsave(filename = "../outs/pdf/supfig_invitro_imputepanel.pdf", device = cairo_pdf, plot = p, width = 8, height = 10, units = "in")

