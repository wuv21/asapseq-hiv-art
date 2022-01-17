source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/invitro_activated_vln_top10.rds")
ncols <- 5

nPerGroup <- 10

figA <- lapply(seq_along(figA), function(i) {
  if (i == 1 | i == nPerGroup + 1) {
    figA[[i]] <- figA[[i]] + plot_layout(tag_level = "keep")
  } else {
    figA[[i]] <- figA[[i]] + plot_layout(tag_level = "new")
  }
  
  return(figA[[i]])
})

pLayout <- c()
for (i in c(0:(length(figA) - 1))) {
  row <- i %/% ncols * 2 + 1 
  
  if (nPerGroup %% ncols != 0) {
    row <- row + (i %/% length(nPerGroup))
  }
  
  col <- i %% ncols * 2 + 1
  
  pLayout <- append(pLayout, patchwork::area(row, col, row + 1, col + 1))
}

p <- wrap_plots(figA) +
  plot_annotation(tag_levels = c("a")) +
  plot_layout(design = pLayout) +
  theme(plot.tag = element_text(size = 12),
    plot.tag.position = c(0, 1))


saveFinalFigure(fn = "supfig_invitro_adt_vln", plot = p, gwidth = 8, gheight = 8)
