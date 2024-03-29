library(ggplot2)
library(patchwork)
library(png)
library(ggpubr)
library(ggrastr)
library(glue)
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

subplotTheme <- theme(
  plot.title.position = "plot",
  plot.title = element_text(size = 12, margin = margin(0,0,0,0), hjust = 0),
  plot.margin = unit(c(0,0,0,0), "pt"),
  panel.background = element_rect(fill = "transparent", colour = NA_character_),
  plot.background = element_rect(fill = "transparent", colour = NA_character_),
  text = element_text(family = "Arial"),
  rect = element_rect(fill = "transparent", colour = NULL)
)

supplementalLollipopTheme <- theme(
  plot.title.position = "plot",
  plot.subtitle = element_text(size = 7, margin = margin(4,0,0,0), hjust = 0.5),
  plot.margin = unit(c(0,0,0,0), "pt"),
  panel.background = element_rect(fill = "transparent", colour = NA),
  plot.background = element_rect(fill = "transparent", colour = NA),
  text = element_text(family = "Arial"),
  rect = element_rect(fill = "transparent", colour = NULL)
)

umapPlotThemeNoLeg <- list(
  # coord_fixed(),
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank())
)

umapPlotThemeLeg <- theme(
  legend.position = "bottom",
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  legend.margin = margin(0,0,0,0)
)


wrapNewAnnotLevel <- function(p) {
  return(p + plot_layout(tag_level = "new"))
}

unclipStripText <- function(g) {
  g2 <- patchworkGrob(g)
  for (i in which(grepl("strip-r", g2$layout$name))){

    g2$grobs[[i]]$grobs[[1]]$layout$clip <- "off"
  }
  
  return(g2)
}

saveFinalFigure <- function(
  plot,
  prefixDir = "../outs",
  fn,
  devices = c("png", "pdf"),
  gheight,
  gwidth) {
  
  if (!is.vector(devices)) {
    devices <- c(devices)
  }
  
  for (d in devices) {
    gfn <- glue("{prefixDir}/{d}/{fn}.{d}")
    
    if (d == "rds") {
      saveRDS(plot, gfn)
    } else if (d == "pdf") {
      ggsave(gfn, plot = plot, dpi = "retina", device = cairo_pdf, width = gwidth, height = gheight, units = "in")
    } else {
      ggsave(gfn, plot = plot, dpi = "retina", device = d, width = gwidth, height = gheight, units = "in")  
    }
  }
}