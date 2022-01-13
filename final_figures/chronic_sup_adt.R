source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/chronic_tcells_vln.rds")
figB <- readRDS("../outs/rds/chronic_tfh_vln.rds")

ncols <- 5

figAB <- c(figA, figB)
figAB <- lapply(seq_along(figAB), function(i) {
  if (i == 1 | i == length(figA) + 1) {
    figAB[[i]] <- figAB[[i]] + plot_layout(tag_level = "keep")
  } else {
    figAB[[i]] <- figAB[[i]] + plot_layout(tag_level = "new")
  }
  
  return(figAB[[i]])
})

pLayout <- c()
for (i in c(0:(length(figAB) - 1))) {
  row <- i %/% ncols * 2 + 1 
  
  if (length(figA) %% ncols != 0) {
    row <- row + (i %/% length(figA))
  }

  col <- i %% ncols * 2 + 1
  
  pLayout <- append(pLayout, patchwork::area(row, col, row + 1, col + 1))
}

p <- wrap_plots(figAB) +
  plot_annotation(tag_levels = c("a")) +
  plot_layout(design = pLayout) +
  theme(plot.tag = element_text(size = 12),
    plot.tag.position = c(0, 1))


saveFinalFigure(fn = "supfig_chronic_adt", devices = c("pdf", "png"), plot = p, gwidth = 8, gheight = 6)

