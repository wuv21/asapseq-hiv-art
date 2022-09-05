source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/chronic_umap_labeledCluster_withPlotLabels.rds")
figB <- readRDS("../outs/rds/chronic_umap_haystackout.rds")
figC <- readRDS("../outs/rds/chronic_discreteAbsolute_matched.rds")
figD <- readRDS("../outs/rds/chronic_discreteHivOnly_matched.rds")
figE <- readRDS("../outs/rds/chronic_lollipop_deseq2.rds")
figF <- readRDS("../outs/rds/chronic_lollipop_tfh.rds")
figG <- readRDS("../outs/rds/chronic_volcano_motifs_chromVAR.rds")
figH <- readRDS("../outs/rds/chronic_chromVAR_motifsUp_HIVneg_enhanced.rds")
figI <- readRDS("../outs/rds/chronic_chromVAR_motifsUp_HIVpos_enhanced.rds")
figJ <- readRDS("../outs/rds/chronic_models.rds")


a_legend <- as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2, title = "Legend\nfor (a)")) +
    theme(legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      legend.justification = "left",
      legend.title = element_text(size = 8, color = "#000000"))))

b_legend <- as_ggplot(get_legend(figB + 
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1, title = "Legend\nfor (b)")) +
    theme(
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      legend.justification = "left",
      legend.title = element_text(size = 8, color = "#000000"))))

figA <- figA + umapPlotThemeNoLeg +
  subplotTheme

figB <- figB + umapPlotThemeNoLeg + 
  subplotTheme +
  theme(axis.title.y = element_text(color = "#FFFFFF00"),
    plot.margin = unit(c(0, 0, 0, 0), "pt"))

figAB <- figA + figB + plot_annotation(tag_levels = "a")
figAB[[2]] <- figAB[[2]] + plot_layout(tag_level = "keep")

figABLegend <- a_legend / b_legend & subplotTheme &
  theme(plot.margin = unit(c(0, 10, 0, 0), "pt"))

figD <- figD + 
  subplotTheme + 
  theme(plot.margin = margin(0,10,0,0, "pt"),
    axis.title.y = element_blank())

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

figEF <- (figE + labs(subtitle = "All T cells") + titleTheme) /
  (figF + labs(subtitle = "Tfh cells") + titleTheme) +
  plot_layout(heights = c(2.6, 1)) +
  plot_annotation(tag_levels = list(c("e", "f"))) &
  theme(
    plot.margin = margin(5,0,0,0),
    plot.tag = element_text(size = 12, margin = margin(l = 4))
  )

figG <- figG +
  subplotTheme +
  theme(legend.margin=margin(10,0,0,0),
    legend.box.margin=margin(-20,-20,-5,-20))

motifTheme <- subplotTheme + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(size = 5.5, angle = 45, vjust = 1, hjust = 1, margin = margin(0, 0, 0, 0)),
  legend.key.size = unit(0.5, 'line'),
  legend.key.height = unit(0.4, 'line'),
  legend.box.margin = margin(r = -5),
  plot.subtitle = element_text(size = 8, hjust = 0.5, margin = margin(b = 4)),
  plot.title = element_text(size = 12, margin = margin(0, 0, 5, 0)))


figH <- figH +
  labs(
    x = "chromVAR motif",
    color = "Mean\ndifference",
    y = "-log10(FDR)") +
  motifTheme +
  theme(plot.margin = unit(c(1, 0, -5, 7.5), "pt"))


figI <- figI + 
  labs(
    x = "chromVAR motif",
    color = "Mean\ndifference",
    y = "-log10(FDR)") +
  motifTheme +
  theme(plot.margin = unit(c(0, 0, -5, 7.5), "pt"))


figHI <- (figH / figI) +
  plot_annotation(tag_levels = list(c("h", "i"))) &
  theme(
    plot.tag = element_text(size = 12),
  )


layout <- c(
  area(1, 1, 2, 4), #ab
  area(1, 5, 2, 8), #legend
  area(3, 1, 5, 3), #c
  area(3, 4, 5, 6), #d
  area(2, 7, 5, 9), #ef
  area(6, 1, 7, 2), #g
  area(6, 3, 7, 6), #hi
  area(6, 7, 7, 9) #j
)

p <- wrap_elements(full = figAB, ignore_tag = TRUE) +
  wrap_elements(figABLegend, ignore_tag = TRUE) +
  wrap_elements(figC + subplotTheme) +
  wrap_elements(figD) +
  wrap_elements(full = figEF, ignore_tag = TRUE) +
  wrap_elements(figG) +
  wrap_elements(full = figHI, ignore_tag = TRUE) +
  wrap_elements(panel = figJ + subplotTheme) +
  plot_annotation(tag_levels = list(c(letters[3:4], "g", "j"))) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "supfig_chronic_main", devices = c("png", "pdf"), gwidth = 8, gheight = 8)