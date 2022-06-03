source("defaultFigureSettings.R")

# figA <- readRDS("../outs/rds/art_frags2_A01.rds")
# figB <- readRDS("../outs/rds/art_frags2_A08.rds")
# figC <- readRDS("../outs/rds/art_frags2_B45.rds")

figA <- readRDS("../outs/rds/art_frags2.rds")

figA[[1]] <- figA[[1]] +
  guides(colour = guide_legend(nrow = 2))

figA[[2]] <- figA[[2]] +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), nrow = 2))

figA <- figA + subplotTheme + theme(plot.margin = margin(0,0,0,-100))

# figB <- figB + subplotTheme + theme(plot.margin = margin(0,0,0,-100))
# figC <- figC + subplotTheme + theme(plot.margin = margin(0,0,0,-100))

# layout <- c(
#   area(1, 1, 4, 6), #a
#   area(5, 1, 14, 12), #b
#   area(1, 7, 4, 12)
# )


#   # plot_layout(design = layout) +
#   
p <- (figA)
  # plot_annotation(tag_levels = c("a")) &
  # theme(plot.tag = element_text(size = 12),
    # plot.tag.position = c(0, 1))

saveFinalFigure(fn = "supfig_art_frags_all", devices = c("png", "pdf"), plot = p, gwidth = 8, gheight = 11)

