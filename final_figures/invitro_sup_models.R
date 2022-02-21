source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/invitro_models.rds")
figB <- readRDS("../outs/rds/invitro_logitRegression_volcano.rds")
figC <- readRDS("../outs/rds/invitro_models_activated.rds")

figA <- figA + guides(color = guide_legend(nrow = 3)) + labs(subtitle = "Total CD4+ T-cells")
figC <- figC + guides(color = guide_legend(nrow = 3)) + labs(subtitle = "Activated CD4+ T-cells")

rocTheme <- theme(plot.subtitle = element_text(hjust = 0.5), plot.title.position = "panel")

layout <- c(
  area(1, 1, 4, 3), #a
  area(1, 4, 6, 8),#b
  area(5, 1, 8, 3) #c
)

p <- wrap_elements(figA + subplotTheme + rocTheme) +
  wrap_elements(figB + subplotTheme) +
  wrap_elements(figC + subplotTheme + rocTheme) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = c("a")) &
  theme(plot.tag = element_text(size = 12),
    plot.tag.position = c(0, 1))

saveFinalFigure(fn = "supfig_invitro_model", devices = c("pdf", "png"), plot = p, gwidth = 8, gheight = 5)

