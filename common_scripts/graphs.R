######
# Generate UMAP df from ArchRProject
######
generateUmapDfFromArchR <- function(
  proj,
  secondGraphColumn = "haystackOut",
  cluster = "Clusters",
  embedding = "UMAP") {
  
  umapFromArchr <- getEmbedding(proj, embedding = embedding)
  
  df <- data.frame(x = umapFromArchr[, 1],
    y = umapFromArchr[, 2],
    secondMetadata = getCellColData(proj, select = secondGraphColumn)[, 1],
    sample = proj$Sample,
    cluster = getCellColData(proj, select = cluster)[, 1])
  
  return(df)
}

######
# UMAP plot themes
######
umapTheme <- theme(
  legend.position = "bottom",
  plot.background = element_rect(fill = "transparent", colour = NA))

######
# Plot dual umap
######
plotDualUmap <- function(proj,
  fn,
  ggtheme = umapTheme,
  secondGraphColumn = "haystackOut",
  cluster = "Clusters",
  embedding = "UMAP") {
  
  df <- generateUmapDfFromArchR(proj,
    secondGraphColumn = secondGraphColumn, cluster = cluster, embedding = embedding)
  
  clusterLabelUmapPos <- df %>% 
    group_by(cluster) %>% 
    summarize(medianX = median(x),
      medianY = median(y),
      iqrXFactor = IQR(x) / 1.25,
      iqrYFactor = IQR(y) / 1.25,
      x = mean(x[x > (medianX - iqrXFactor) & x < (iqrXFactor + medianX)]),
      y = mean(y[y > (medianY - iqrYFactor) & y < (iqrYFactor + medianY)]))
  
  p1 <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.8, aes(color = cluster), size = 0.5) +
    labs(x = "UMAP 1",
      y = "UMAP 2",
      color = "Cluster") +
    theme_classic() +
    ggtheme +
    guides(colour = guide_legend(override.aes = list(size=6), nrow = 3)) +
    ggrepel::geom_label_repel(data = clusterLabelUmapPos,
      aes(x = x, y = y, label = cluster),
      label.size = 0.1,
      size = 2.5,
      segment.color = NA)
  
  p2 <- df %>%
    arrange(secondMetadata) %>%
    ggplot(., aes(x = x, y = y)) +
    geom_point(alpha = 0.7, aes(color = secondMetadata), size = 0.5) +
    labs(x = "UMAP 1",
      y = "UMAP 2",
      color = "HIV+") +
    theme_classic() +
    ggtheme +
    scale_color_manual(values = c("#cccccc30", "#ff000060")) +
    guides(colour = guide_legend(override.aes = list(size=6)))
  
  ggsave(fn, plot = p1 + p2, dpi = "retina", width = 10, height = 6)  
}

######
# Discrete bar plot theme
######
ggDiscreteBarTheme <- list(
  theme_classic(),
  theme(legend.position = "none",
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 5),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.margin = unit(c(0.1, 0.3, 0.1, 0.1), "in")
  ),
  coord_flip(clip = "off"),
  scale_fill_manual(values = c("#999999", "#e63946")),
  scale_color_manual(values = c("#444444", "#e63946"))
)


plotDiscreteBar <- function(proj, fn,
  cluster,
  secondGraphColumn = "haystackOut",
  graphType = "proportion") {
  df <- generateUmapDfFromArchR(proj, cluster = cluster, secondGraphColumn = secondGraphColumn)
  
  stopifnot(graphType %in% c("proportion", "absolute", "hivOnly"))
  
  if (graphType == "proportion") {
    p <- df %>%
      group_by(cluster) %>%
      summarize(propPos = sum(secondMetadata) / n(),
        propNeg = sum(!secondMetadata) / n()) %>%
      pivot_longer(cols = starts_with("prop"),
        names_to = "proviralStatus",
        values_to = "proportion") %>%
      mutate(proviralStatus = gsub("prop", "", proviralStatus)) %>%
      ggplot(., aes(x = cluster, y = proportion, fill = proviralStatus)) +
      geom_bar(stat = "identity", width = 0.7, position = position_dodge(0.7)) +
      labs(x = "Cluster",
        y = "Proportion of cluster") +
      geom_text(aes(y = ifelse(proportion > 0.85, proportion, proportion + 0.02),
        label = scales::label_percent(accuracy = 0.1)(proportion),
        color = proviralStatus,
        hjust = ifelse(proportion > 0.85, 1, 0)),
        position = position_dodge(0.7),
        vjust = 0.5,
        size = 1.5) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
      ggDiscreteBarTheme    
  } else if (graphType == "absolute") {
    p <- df %>%
      mutate(cluster = factor(cluster),
        hivPos = factor(ifelse(secondMetadata, "Pos", "Neg"))) %>%
      dplyr::count(cluster, hivPos, .drop = FALSE) %>%
      {ggplot(., aes(x = cluster, y = n, fill = hivPos)) +
          geom_bar(stat = "identity", width = 0.7, position = position_dodge(0.7)) +
          labs(x = "Cluster",
            y = "Number of cells") +
          geom_text(aes(y = n + max(.$n) * 0.02, label = n, color = hivPos),
            position = position_dodge(0.7),
            vjust = 0.5,
            hjust = 0,
            size = 1.5) +
          scale_y_continuous(expand = c(0, 0), limits = c(0, max(.$n + 500))) +
          scale_x_discrete(drop = FALSE) +
          ggDiscreteBarTheme}
    
  } else {
    p <- df %>%
      mutate(cluster = factor(cluster),
        hivPos = ifelse(secondMetadata, "Pos", "Neg")) %>%
      dplyr::filter(hivPos != "Neg") %>%
      dplyr::count(cluster) %>%
      mutate(proportion = n / sum(n)) %>%
      {ggplot(., aes(x = cluster, y = proportion)) +
          geom_bar(stat = "identity", width = 0.7, fill = "#e63946") +
          labs(x = "Cluster",
            y = "Proportion of HIV+ cells") +
          geom_text(aes(y = ifelse(proportion > 0.85, proportion, proportion + 0.02),
            label = scales::label_percent(accuracy = 0.1)(proportion),
            hjust = ifelse(proportion > 0.85, 1, 0)),
            vjust = 0.5,
            color = "#e63946",
            size = 1.5) +
          scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
          ggDiscreteBarTheme +
          scale_fill_manual(values = "#e63946")}
  }

  ggsave(fn, plot = p, width = 2.5, height = 2.5, dpi = "retina")
}

# violin plot cleaner
clean_VlnPlot <- function(gList, finalplotCol = 5, newTitle = c()) {
  g_grob <- lapply(seq_along(gList), function(i) {
    x <- gList[[i]]
    x <- x +
      theme(axis.title = element_blank(),
        legend.position = "none")
    
    if (length(newTitle) > 0) {
      x <- x +
        labs(title = newTitle[i]) +
        theme(plot.title = element_text(size = 10))
    }
    
    x <- ggplot_build(x)
    x$data[[2]]$colour <- "#00000020"
    x$data[[2]]$size <- 0.2
    
    return(ggplot_gtable(x))
  })
  
  return(gridExtra::grid.arrange(grobs = g_grob, ncol = finalplotCol))
}



# cell plot
makeCellPlot <- function(seu, tsa_catalog, cellID,
  outputDir, outputFn = glue("cellplot_{cellID}.png")) {
  
  cellIDData <- seu@assays$tsa@data[, cellID]
  
  singleAdtDf <- data.frame(marker = names(cellIDData),
    normExp = cellIDData) %>%
    left_join(tsa_catalog, by = c("marker" = "DNA_ID")) %>%
    filter(!isCtrl) %>%
    dplyr::select(marker, normExp, cleanName) %>%
    mutate(cleanName = factor(cleanName, levels = cleanName),
      id = row_number(),
      angle = 90 - 360 * (id - 0) / n(),
      hjust = ifelse(angle < -90, 1, 0),
      angle = ifelse(angle < -90, angle + 180, angle),
      label = ifelse(normExp > 2, as.character(cleanName), ""),
      barFill = ifelse(label != "", "#e63946", "#999999"))
  
  
  cellPlt <- ggplot(singleAdtDf, aes(x = cleanName, y = normExp)) +
    geom_bar(stat = "identity", alpha = 0.8, aes(fill = barFill), width = 0.4) +
    # scale_fill_distiller(palette = "Accent") +
    scale_fill_identity() +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-0.5, 4), "in")
    ) +
    # geom_segment(x = 0, y = 0, xend = 360, yend = 0, colour = "#e6394610", size=0.1, inherit.aes = FALSE) +
    ylim(-1, max(singleAdtDf$normExp, na.rm = T) + 0.5) +
    geom_text(aes(x = cleanName, y = normExp + 0.1, label = label, hjust = hjust, angle = angle),
      color="black", alpha = 0.6, size = 2.25, inherit.aes = FALSE) +
    coord_polar()
  
  ggsave(filename = paste0(outputDir, "/", outputFn),
    plot = cellPlt, width = 8, height = 8,
    dpi = "retina", bg = "#FFFFFF")
  
  return(cellPlt)
}

makeDifferentialLollipop <- function(markers, fn) {
  pSapPi <- markers %>%
    arrange(desc(abs(piScore)), .by_group = TRUE) %>%
    mutate(cleanName = factor(cleanName, levels = rev(cleanName))) %>%
    mutate(markerColor = ifelse(piScore > 0, "#e63946", "#999999")) %>%
    {ggplot(., aes(x = abs(piScore), y = cleanName, color = markerColor)) +
        geom_point(size = 1.5) + 
        geom_segment(aes(x = 0, xend = abs(piScore), y = cleanName, yend = cleanName), linetype = "dotted") +
        labs(y = "Surface marker",
          x = "π score (log2FC * -log10(adj pval))") +
        scale_x_continuous(expand = c(0, 0), limits = c(0, max(abs(.$piScore) + 0.3))) +
        theme_classic() +
        theme(axis.text.y = element_text(size = 6),
          legend.position = "none",
          panel.background = element_rect(fill = NULL, colour = NULL),
          strip.placement = "outside") +
        scale_color_identity() +
        facet_grid(Status ~ ., scales = "free", space = "free")} 
  
  ggsave(filename = fn, plot = pSapPi, width = 3, height = 2 + (nrow(markers) * 0.08), dpi = "retina")
  
  return(pSapPi)
}



