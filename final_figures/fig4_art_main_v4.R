source("defaultFigureSettings.R")
library(ComplexHeatmap)

figA <- readRDS("../outs/rds/art_lollipop_deseq2.rds")
figB <- readRDS("../outs/rds/art_lollipop_wilcox.rds")
figC <- readRDS("../outs/rds/art_tcmttm_markers.rds")
figD <- readRDS("../outs/rds/art_tem_markers.rds")
figERds <- readRDS("../outs/rds/artLess_markerTest_deviations_heatmap.rds")
figETxtRds <- readRDS("../outs/rds/artLess_markerTest_deviations_heatmapTxt.rds")
figF <- readRDS("../outs/rds/art_chromVAR_tcmttm_motifsUp_HIVpos_enhanced.rds")

titleTheme <- subplotTheme + theme(
  plot.title = element_text(size = 8, hjust = 0.5, margin = margin(b = 4)),
  axis.text = element_text(size = 6))

figA <- figA + labs(title = "All T cells (DESeq2)") + titleTheme
figB <- figB + labs(title = "All T cells (Wilcoxon)") + titleTheme
figC <- figC + labs(title = "Tcm/Ttm (DESeq2)") + titleTheme
figD <- figD + labs(title = "Tem (DESeq2)") + titleTheme

figF <- figF + subplotTheme +
  labs(
    title = "Tcm/Ttm HIV+ Enriched Motifs",
    color = "Mean\ndifference",
    y = "FDR") +
  titleTheme +
  theme(plot.subtitle = element_blank(),
    axis.title.y = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.5, 'line'),
    plot.margin = margin(r = 5, unit = "line"))

# figD <- ggplotGrob(figD)
# for(i in which(grepl("strip-r", figD$layout$name))){
#   figD$grobs[[i]]$layout$clip <- "off"
# }

set.seed(21)
figE <- ComplexHeatmap::Heatmap(
  matrix = figERds,
  name = "Mean z deviation",
  border = FALSE,
  cluster_rows = TRUE,
  show_row_names = FALSE,
  row_title_gp = gpar(fontsize = 6),
  row_dend_side = "left",
  show_row_dend = TRUE,
  cluster_columns = FALSE,
  column_split = c(1,1,2,2,3,3),
  show_column_names = FALSE,
  show_heatmap_legend = TRUE,
  column_dend_side = "top",
  column_title = NULL,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 6),
    title_position = "topleft",
    legend_direction = "horizontal",
    legend_width = unit(1, "cm")),
  show_column_dend = TRUE,
  column_km = 3,
  row_km = 10,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    grid.text(figETxtRds[i, j], x, y, gp = gpar(fontsize = 3, hjust = 0.5, vjust = 0.5))
  },
  bottom_annotation = ComplexHeatmap::HeatmapAnnotation(
    `Phenotype` = str_split_fixed(colnames(figERds), "_", n = 2)[, 1],
    `HIV` = ifelse(str_split_fixed(colnames(figERds), "_", n = 2)[, 2] == "TRUE", "HIV+", "HIV-"),
    col = list(
      `Phenotype` = c("Tem/effector" = "#1b9e77", "MAIT" = "#d95f02", "Tcm/Ttm" = "#7570b3"),
      `HIV` = c("HIV+" = "#e63946", "HIV-" = "#999999")
    ),
    simple_anno_size = convertHeight(unit(1/nrow(figERds) * 2, "npc"), "mm"),
    annotation_name_gp = gpar(fontsize = 6),
    annotation_legend_param = list(
      `Phenotype` = list( 
        title_gp = gpar(fontsize = 6, 
          fontface = "bold"),
        legend_direction = "horizontal", 
        labels_gp = gpar(fontsize = 6)),
      `HIV` = list( 
        title_gp = gpar(fontsize = 6,
          fontface = "bold"), 
        legend_direction = "horizontal",
        labels_gp = gpar(fontsize = 6)))
  ),
  right_annotation = rowAnnotation(
    foo = anno_block(gp = gpar(fill = "#000000"), width = unit(1, "mm")),
    bar = anno_block(
      graphics = function(index, levels) {
        lbls <- rownames(figERds)[index]
        lbls <- gsub("_\\d+", "", lbls)
        lbls <- str_wrap(paste(lbls, collapse = " "), width = 30)
        
        grid.rect(gp = gpar(fill = NA, col = NA))
        txt = paste(lbls, collapse = ",")
        grid.text(txt, 0.01, 0.5, rot = 0, hjust = 0, gp = gpar(fontsize = 6))
      },
      width = unit(4, "cm"))
  )
)

layout <- c(
  area(1, 1, 4, 2), #a
  area(1, 3, 3, 4), #b
  area(5, 1, 6, 2), #c
  area(4, 3, 7, 4),  #d
  area(1, 5, 8, 12),  #e
  area(7, 1, 8, 3)  #f
)

p <- figA + figB + figC + figD +
  # wrap_elements(figB) +
  # wrap_elements(figC) +
  # wrap_elements(figD) +
  wrap_elements(grid.grabExpr(draw(figE, heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legend = TRUE))) +
  wrap_elements(plot = figF) +
  plot_annotation(tag_levels = list(c(letters[1:7]))) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig4_art", devices = c("png", "pdf"), gwidth = 8.25, gheight = 7)

warnings()

