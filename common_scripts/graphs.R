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
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 6),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA)
  ),
  coord_flip(),
  scale_fill_manual(values = c("#999999", "#e63946")),
  scale_color_manual(values = c("#444444", "#e63946"))
)


plotDiscreteBar <- function(proj, fn,
  cluster,
  secondGraphColumn = "haystackOut",
  graphType = "proportion") {
  df <- generateUmapDfFromArchR(proj, cluster = cluster, secondGraphColumn = secondGraphColumn)
  
  stopifnot(graphType %in% c("proportion", "absolute"))
  
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
  } else {
    p <- df %>%
      mutate(cluster = factor(cluster),
        hivPos = factor(ifelse(secondMetadata, "Pos", "Neg"))) %>%
      dplyr::count(cluster, hivPos, .drop = FALSE) %>%
      {ggplot(., aes(x = cluster, y = n, fill = hivPos)) +
          geom_bar(stat = "identity", width = 0.7, position = position_dodge(0.7)) +
          labs(x = "Cluster",
            y = "Number of cells") +
          geom_text(aes(y = n + 30, label = n, color = hivPos),
            position = position_dodge(0.7),
            vjust = 0.5,
            hjust = 0,
            size = 1.5) +
          scale_y_continuous(expand = c(0, 0), limits = c(0, max(.$n + 500))) +
          scale_x_discrete(drop = FALSE) +
          ggDiscreteBarTheme}
    
  }

  ggsave(fn, width = 2.5, height = 2.5, dpi = "retina")
}

