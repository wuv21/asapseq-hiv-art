library(ggplot2)
library(patchwork)
library(png)
library(ggpubr)

subplotTheme <- theme(
  plot.title.position = "plot",
  plot.title = element_text(size = 12, margin = margin(0,0,0,0)),
  plot.margin = unit(c(0,0,0,0), "pt"),
  text = element_text(family = "Arial"),
  rect = element_rect(fill = "transparent", colour = NULL)
)

umapPlotThemeNoLeg <- theme(
  legend.position = "none",
  axis.ticks = element_blank(),
  axis.text = element_blank()
)

umapPlotThemeLeg <- theme(
  legend.position = "none",
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(-10,-10,-10,-10),
)