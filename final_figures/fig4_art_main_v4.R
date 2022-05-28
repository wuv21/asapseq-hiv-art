source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/art_lollipop_deseq2.rds")
figB <- readRDS("../outs/rds/art_lollipop_negbinom.rds")
figC <- readRDS("../outs/rds/art_volcano_motifs_chromVAR.rds")
figD <- readRDS("../outs/rds/art_volcano_motifs_chromVAR_lola_encode.rds")

titleTheme <- subplotTheme + theme(plot.title = element_text(size = 8, hjust = 0.5, margin = margin(t = -10, b = 10)))
figA <- figA + labs(title = "DESeq2") + titleTheme
figB <- figB + labs(title = "Neg. binomial") + titleTheme

figC <- figC + 
  labs(title = "CISBP") + 
  titleTheme + 
  theme(legend.margin=margin(3,0,0,0), legend.box.margin=margin(-10,-10,-10,-10))
  
figD <- figD + 
  labs(title = "ENCODE (from LOLA)") + 
  titleTheme + theme(legend.margin=margin(3,0,0,0), legend.box.margin=margin(-10,-10,-10,-10))

decreaseText <- function(x) {
  x$layers[[3]]$aes_params$size <- 4 / ggplot2:::.pt
  
  return(x)
}

figC <- decreaseText(figC)
figD <- decreaseText(figD)



layout <- c(
  area(1, 1, 4, 3), #a
  area(1, 4, 4, 6), #b
  area(5, 1, 7, 3), #c
  area(5, 4, 7, 6)  #d
)

p <- wrap_elements(figA) +
  wrap_elements(figB) +
  wrap_elements(figC) +
  wrap_elements(figD) +
  plot_annotation(tag_levels = list(c(letters[1:7]))) +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig4_art", devices = c("png", "pdf"), gwidth = 8, gheight = 7)

