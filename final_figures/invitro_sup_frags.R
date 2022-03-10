source("defaultFigureSettings.R")

figA <- readRDS("../outs/rds/invitro_upset_cbc.rds")
figB <- readPNG("../outs/png/p24pos.png", native = TRUE)

figC <- readRDS("../outs/rds/invitro_frags.rds")
figC <- figC & subplotTheme
figC[[2]] <- figC[[2]] + theme(plot.margin = margin(0, 0, 0.5, 0, unit = "lines"))


figF <- readRDS("../outs/rds/InVitro_umap_unannotCluster_withPlotLabels.rds")
figF <- figF + guides(color = guide_legend(ncol = 5, direction = "horizontal", override.aes = list(size = 3))) + umapPlotThemeLeg

idxDir <- "../data/idxstats/"
fns <- list.files(idxDir, pattern = ".tsv", full.names = FALSE)
mimitouDataset <- lapply(paste0(idxDir, fns), read.table, sep = "\t", header = FALSE)
names(mimitouDataset) <- gsub("_stats.tsv", "", fns)

figIdxStats <- bind_rows(mimitouDataset, .id = "genome") %>%
  mutate(chrClean = ifelse(grepl("^(GL|KI)", V1), "Unlocalized", V1)) %>% 
  mutate(chrClean = gsub("^chr", "", chrClean)) %>%
  mutate(chrClean = ifelse(chrClean == "*", "Unaligned", chrClean)) %>%
  mutate(chrClean = factor(chrClean, levels = str_sort(unique(chrClean), numeric = TRUE))) %>%
  group_by(genome, chrClean) %>% 
  summarize(Mapped = sum(V3), Unmapped = sum(V4)) %>% 
  pivot_longer(cols = all_of(c("Mapped", "Unmapped")), names_to = "metric", values_to = "reads") %>%
  ggplot(aes(y = chrClean, color = metric, x = reads)) +
  geom_point(stat = "identity", position = position_dodge(0.8)) +
  geom_linerange(aes(xmin = 0, xmax = reads, y = chrClean),
    linetype = "dashed", size = 0.5, position = position_dodge(0.8)) +
  theme_classic() +
  labs(y = "Chromosome/region", x = "Reads", color = "Metric") +
  theme(legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.title = element_text(size = 6),
    axis.text.y = element_text(size = 6, color = "#000000"),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1, color = "#000000"),
    panel.grid.major.y = element_line(color = "#BCBCBC30", size = 5),
    legend.margin = margin(-10, 0, -10, 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 6)) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~ genome, nrow = 1, scales = "free")


layout <- c(
  area(1, 1, 3, 4), #a
  area(1, 5, 3, 6), #flow/b
  area(4, 1, 8, 6),
  area(1, 7, 11, 12), #cde
  area(9, 1, 11, 6) #f
)

p <- wrap_elements(figA + subplotTheme) +
  (wrap_elements(full = figB) + subplotTheme) +
  wrap_elements(figIdxStats + subplotTheme) + 
  wrap_elements(figC) + 
  wrap_elements(figF + subplotTheme) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = c("a")) &
  theme(plot.tag = element_text(size = 12),
    plot.tag.position = c(0, 1))

saveFinalFigure(plot = p, fn = "supfig_invitro_frags", devices = c("png", "pdf"), gwidth = 8, gheight = 10.5)



