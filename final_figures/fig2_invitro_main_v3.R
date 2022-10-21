source("defaultFigureSettings.R")
library(ggtext)

figVolcano <- readRDS("../outs/rds/invitro_peaks_volcano_activatedHIV.rds")

figTopPeaks <- readRDS("../outs/rds/invitro_markerTest_activated_topPeaks.rds")
figCCR5 <- readRDS("../outs/rds/invitro_activated_ccr5.rds")

figMotifDown <- readRDS("../outs/rds/invitro_motifsUp_activatedHIVneg.rds")
figMotifUp <- readRDS("../outs/rds/invitro_motifsUp_activatedHIVpos.rds")

figModel <- readRDS("../outs/rds/invitro_models.rds")
figROCVolcano <- readRDS("../outs/rds/invitro_logitRegression_volcano.rds")
figActModel <- readRDS("../outs/rds/invitro_models_activated.rds")

figModel <- figModel + guides(color = guide_legend(ncol = 2)) + labs(subtitle = "Total CD4+ T-cells")
figActModel <- figActModel + guides(color = guide_legend(ncol = 2)) + labs(subtitle = "Activated CD4+ T-cells")

figFancyUpset <- readRDS("../outs/rds/invitro_fancyUpset_activatedLaterCells.rds")

figROCVolcano <- figROCVolcano + 
  subplotTheme +
  guides(color = guide_legend(ncol = 3)) +
  theme(
    legend.title = element_blank(),
    plot.margin = unit(c(0, 5, -20, -30), "pt"),
    legend.position = "bottom",
    # legend.key.height = unit(10, "points"),
    # legend.spacing.y = unit(2, "points"),
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-5,-10,20,-10))

rocTheme <- theme(
  plot.subtitle = element_text(hjust = 0.5),
  legend.text = element_text(size = 6),
  plot.title.position = "panel")

figCCR5[[1]] <- figCCR5[[1]] + 
  subplotTheme +
  labs(title = "<span style='font-size:7pt'> 
    <span style='color:#999999;'>HIV-</span> | 
    <span style='color:#e63946;'>HIV+</span> | Gene (<span style='color:#FFA500;'>positive</span> vs <span style='color:#0000FF;'>negative</span> orientation)
    </span>") +
  theme(
    plot.title = element_markdown(),
    plot.title.position = "panel",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

figCCR5[[2]] <- figCCR5[[2]] + subplotTheme + theme(plot.title = element_blank())
figCCR5 <- wrap_plots(figCCR5, ncol = 1, heights = c(3,2,1))

figTopPeaks <- figTopPeaks + 
  subplotTheme + 
  theme(legend.position = "right",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    legend.direction = "vertical",
    legend.text.align = 1,
    legend.box.margin=margin(0,0,0,0),
    plot.margin = margin(0,0,10,0))

motifTheme <- theme(
  plot.subtitle = element_text(size = 7, hjust = 0.5, margin = margin(0,0,1,0)),
  plot.title.position = "panel"
)

figMotifDown <- figMotifDown + labs(subtitle = "HIV-") + motifTheme
figMotifUp <- figMotifUp + labs(subtitle = "HIV+") + motifTheme
  
figMotifs <- figMotifDown + (figMotifUp + theme(axis.title.y = element_text(color = "#FFFFFF00"))) &
    subplotTheme + theme(plot.title.position = "panel", plot.margin = margin(0, 0, 0, 0))

layout <- c(
  area(1, 1, 3, 2), #volcano
  area(1, 3, 3, 9), #topPeaks
  area(4, 1, 6, 5), #track
  area(4, 6, 6, 9), #motifs
  area(7, 1, 9, 3), #model
  area(7, 4, 9, 6), #volcano
  area(7, 7, 9, 9), #act Model
  area(10, 1, 12, 9) # fancy upset
)

p <- wrap_elements(panel = figVolcano + subplotTheme +
    theme(legend.margin=margin(0,0,20,0), legend.box.margin=margin(-5,-10,0,-10))) +
  wrap_elements(figTopPeaks + theme(plot.margin = margin(t = -10), strip.text = element_text(size = 8))) +
  wrap_elements(figCCR5 & theme(plot.margin = margin(t = -20))) +
  wrap_elements(figMotifs) +
  wrap_elements(full = figModel + subplotTheme + rocTheme + theme(plot.margin = unit(c(0,0,0,-10), "pt"))) +
  wrap_elements(panel = figROCVolcano) +
  wrap_elements(full = figActModel + subplotTheme + rocTheme) +
  wrap_elements(full = figFancyUpset + subplotTheme) +
  plot_annotation(tag_levels = c("a")) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig2", devices = c("png", "pdf"), gwidth = 8.1, gheight = 9)

