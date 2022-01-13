source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/invitro_lollipop.rds")
figB <- readRDS("../outs/rds/invitro_lollipop_early.rds")

layout <- c(
  area(1, 1, 5, 5), #a
  area(1, 6, 3, 10) #b
)

p <- wrap_elements(figA + labs(subtitle = "All CD4+ T-cells") + supplementalLollipopTheme) +
  wrap_elements(figB + labs(subtitle = "Early differentiated/naive CD4+ T-cells") + supplementalLollipopTheme) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = c("a")) &
  theme(plot.tag = element_text(size = 12),
    plot.tag.position = c(0, 1))

saveFinalFigure(fn = "supfig_invitro_adt", devices = c("pdf", "png"), plot = p, gwidth = 6, gheight = 7)

