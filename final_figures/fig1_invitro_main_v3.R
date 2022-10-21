source("defaultFigureSettings.R")
library(ggtext)

figA <- readRDS("../outs/rds/InVitro_umap_labeledCluster_withPlotLabels.rds")
figB <- readRDS("../outs/rds/invitro_umap_haystackout.rds")
figC <- readRDS("../outs/rds/inVitro_discreteAbsolute_matched.rds")
figD <- readRDS("../outs/rds/inVitro_discreteHivOnly_matched.rds")
figE <- readRDS("../outs/rds/invitro_lollipop_activeLate2.rds")

a_legend <- as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2,
      title.position = "top",
      title.hjust = 0,
      title = "Legend for (a)")) +
    theme(legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      legend.justification = "left",
      legend.text = element_text(margin = margin(r = 12, unit = "pt")),
      legend.title = element_text(size = 7, color = "#000000"))))

b_legend <- as_ggplot(get_legend(figB + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1,
      title = "Legend for (b)",
      title.position = "top",
      title.hjust = 0)) +
    theme(
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      legend.justification = "left",
      legend.text = element_text(margin = margin(r = 12, unit = "pt")),
      legend.title = element_text(size = 7, color = "#000000"))))

figA <- figA + umapPlotThemeNoLeg +
  subplotTheme +
  coord_fixed(ratio = 1)

figB <- figB + umapPlotThemeNoLeg + 
  subplotTheme +
  theme(axis.title.y = element_text(color = "#FFFFFF00"),
    plot.margin = unit(c(0, 0, 0, 0), "pt")) +
  coord_fixed(ratio = 1)

figAB <- figA + figB + plot_annotation(tag_levels = "a")
figAB[[2]] <- figAB[[2]] + plot_layout(tag_level = "keep")

figC <- figC +
  labs(x = "<span style='font-size:7pt'>Number of cells (<span style='color:#999999;'>HIV-</span> | <span style='color:#e63946;'>HIV+</span>)</span>") +
  subplotTheme +
  theme(axis.title.x = element_markdown())

figD <- figD +
  labs(x = "<span style='font-size:7pt'>Proportion of <span style='color:#e63946;'>HIV+</span> cells</span>") +
  subplotTheme +
  theme(axis.title.x = element_markdown())

figE <- figE +
  subplotTheme + theme(
    plot.title = element_text(size = 7, hjust = 0.5, margin = margin(b = 4)),
    axis.text = element_text(size = 6),
    strip.background = element_rect(fill = "transparent", colour = NA_character_),
    strip.text = element_text(size = 7, color = "#000000"),
    panel.background = element_rect(fill = "transparent", colour = NA_character_),
    plot.background = element_rect(fill = "transparent", colour = NA_character_),
    axis.title.y = element_blank(),
    legend.text = element_text(margin = margin(-8)),
    plot.margin = margin(b = 10))

figABLegend <- (
  a_legend + 
    subplotTheme +
    theme(
      plot.margin = unit(c(0, 0, 0, -15), "pt"))) /
  (b_legend +
    subplotTheme +
    theme(
      plot.margin = unit(c(-50, 0, 0, -15), "pt")))

layout <- c(
  area(1, 1, 2, 4), #ab
  area(1, 5, 2, 8), #legend
  area(3, 1, 5, 3), #c
  area(3, 4, 5, 6), #d
  area(2, 7, 5, 9) #e
)

p <- wrap_elements(full = figAB, ignore_tag = TRUE) +
  wrap_elements(figABLegend, ignore_tag = TRUE) + 
  wrap_elements(figC) +
  wrap_elements(figD + theme(plot.margin = unit(c(0, 0, 0, 10), "pt"))) +
  wrap_elements(figE) +
  plot_annotation(tag_levels = list(c(letters[3:5]))) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig1", devices = c("png", "pdf"), gwidth = 8, gheight = 6)

