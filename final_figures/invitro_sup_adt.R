source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/invitro_lollipop2.rds")
figB <- readRDS("../outs/rds/invitro_lollipop_early2.rds")

layout <- c(
  area(1, 1, 5, 5), #a
  area(1, 6, 3, 10) #b
)

titleTheme <- subplotTheme + theme(
  plot.subtitle = element_text(size = 8, hjust = 0.5, margin = margin(b = 4)),
  plot.title.position = "panel",
  axis.text = element_text(size = 6),
  strip.background = element_rect(fill = "transparent", colour = NA_character_),
  strip.text = element_text(size = 8, color = "#000000"),
  panel.background = element_rect(fill = "transparent", colour = NA_character_),
  plot.background = element_rect(fill = "transparent", colour = NA_character_),
  axis.title.y = element_blank(),
  legend.text = element_text(margin = margin(-4), hjust = 1),
  plot.margin = margin(b = 10))


p <- wrap_elements(figA + labs(subtitle = "All CD4+ T-cells") + titleTheme) +
  wrap_elements(figB + labs(subtitle = "Early differentiated/naive CD4+ T-cells") + titleTheme) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = c("a")) &
  theme(plot.tag = element_text(size = 12),
    plot.tag.position = c(0, 1))

saveFinalFigure(fn = "supfig_invitro_adt", devices = c("pdf", "png"), plot = p, gwidth = 6, gheight = 7)

