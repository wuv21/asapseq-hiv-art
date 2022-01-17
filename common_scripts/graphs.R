library(gghalves)

######
# COLOR SCHEMES
######
HIVPOSCOLOR <- "#e63946"
HIVNEGCOLOR <- "#aaaaaa"
BASEPTAXISFONTSIZE <- 7
BASEPTFONTSIZE <- 7
BASEFONTSIZE <- BASEPTFONTSIZE / ggplot2:::.pt


savePlot <- function(plot, fn, devices, gheight, gwidth, rdsPlot = NULL, scale = 1, customSavePlot = NULL) {
  if (!is.vector(devices)) {
    devices <- c(devices)
  }
  
  for (d in devices) {
    gfn <- glue("outs/{d}/{fn}.{d}")
    
    if (d == "rds" & !is.null(rdsPlot)) {
      saveRDS(rdsPlot, gfn)
    } else if (d == "rds") {
      saveRDS(plot, gfn)
    } else if (!is.null(customSavePlot)) {
      ggsave(gfn, plot = customSavePlot, dpi = "retina", device = d, width = gwidth, height = gheight, scale = scale)
    } else {
      ggsave(gfn, plot = plot, dpi = "retina", device = d, width = gwidth, height = gheight, scale = scale)  
    }
  }
}


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
  legend.text = element_text(size = BASEPTFONTSIZE),
  legend.title = element_blank(),
  axis.title = element_text(size = BASEPTFONTSIZE),
  axis.text = element_text(size = BASEPTAXISFONTSIZE),
  legend.spacing.x = unit(BASEPTFONTSIZE / 2, 'points'),
  # legend.spacing.y = unit(BASEPTFONTSIZE / 4, 'points'),
  legend.key.size = unit(BASEPTFONTSIZE * 1.1, 'points'),
  legend.background = element_rect(fill = "transparent", colour = NA),
  panel.background = element_rect(fill = "transparent", colour = NA),
  plot.background = element_rect(fill = "transparent", colour = NA))


######
# Plot umap
######
plotUmap <- function(
  proj,
  fn,
  devices = c("png", "rds"),
  ggtheme = umapTheme,
  colorBy = "Clusters",
  colorByLabel = colorBy,
  colorScheme = NULL,
  labelColors = TRUE,
  rasterize = TRUE,
  bringToTop = FALSE,
  propInLegend = FALSE,
  embedding = "UMAP") {
  
  df <- generateUmapDfFromArchR(proj, cluster = colorBy, embedding = embedding)
  
  if (bringToTop) {
    df <- df %>%
      arrange(cluster)
  }
  
  if (colorBy == "haystackOut") {
    df <- df %>%
      mutate(cluster = ifelse(cluster, "HIV+", "HIV-"))
  }
  
  if (propInLegend) {
    df <- df %>%
      group_by(cluster) %>%
      mutate(cluster = glue("{cluster} ({round(n() / nrow(.) * 100, digits = 1)}%)"))
  }
  
  p1 <- ggplot(df, aes(x = x, y = y))
  
  if (rasterize) {
    p1 <- p1 + rasterize(geom_point(alpha = 0.6, aes(color = cluster), size = 0.25), dpi = 300)
  } else {
    p1 <- p1 + geom_point(alpha = 0.8, aes(color = cluster), size = 0.25)
  }
  
  p1 <- p1 +
  labs(x = "UMAP 1",
    y = "UMAP 2",
    color = colorByLabel) +
  theme_classic() +
  ggtheme +
  guides(colour = guide_legend(override.aes = list(size = 4), nrow = 5))
    
  if (is.null(colorScheme) & colorBy == "haystackOut") {
    p1 <- p1 + scale_color_manual(values = c(HIVNEGCOLOR, HIVPOSCOLOR))
  } else if (!is.null(colorScheme)) {
    p1 <- p1 + colorScheme
  }
  
  if (labelColors) {
    clusterLabelUmapPos <- df %>% 
      group_by(cluster) %>% 
      summarize(medianX = median(x),
        medianY = median(y),
        iqrXFactor = IQR(x) / 1.5,
        iqrYFactor = IQR(y) / 1.5,
        x = mean(x[x > (medianX - iqrXFactor) & x < (iqrXFactor + medianX)]),
        y = mean(y[y > (medianY - iqrYFactor) & y < (iqrYFactor + medianY)]))
    
    p1 <- p1 + 
      ggrepel::geom_label_repel(data = clusterLabelUmapPos,
        aes(x = x, y = y, label = cluster),
        label.size = 0.05,
        size = BASEFONTSIZE,
        segment.color = NA)
  }
  
  savePlot(plot = p1, fn = fn, devices = devices, gheight = 4, gwidth = 3.5)
}


ggDiscreteLollipopTheme <- list(
  theme_classic(),
  coord_cartesian(clip = "off"),
  theme(legend.position = "none",
    axis.title = element_text(size = BASEPTFONTSIZE),
    axis.text = element_text(size = BASEPTAXISFONTSIZE, color = "#000000"),
    # plot.background = element_rect(fill = "transparent", colour = NA),
    # panel.background = element_rect(fill = "transparent", colour = NA),
    plot.margin = unit(c(0.1, 0.3, 0.1, 0.1), "in")
  )
)


plotDiscreteLollipop <- function(proj,
  fn,
  devices = c("png", "rds"),
  cluster,
  gheight = 3,
  gwidth = 3,
  secondGraphColumn = "haystackOut",
  graphType = "proportion") {
  df <- generateUmapDfFromArchR(proj, cluster = cluster, secondGraphColumn = secondGraphColumn)
  
  stopifnot(graphType %in% c("absolute", "hivOnly"))
  
  if (graphType == "absolute") {
    p <- df %>%
      mutate(cluster = factor(cluster),
        hivPos = factor(ifelse(secondMetadata, "Pos", "Neg"))) %>%
      dplyr::count(cluster, hivPos, .drop = FALSE) %>%
      {ggplot(., aes(y = cluster, fill = hivPos)) +
          geom_linerange(aes(xmin = 0, xmax = n, y = cluster, color = hivPos),
            linetype = "dotted", position = position_dodge(0.8)) +
          geom_point(aes(x = n, color = hivPos), size = 1, position = position_dodge(0.8)) + 
          labs(y = "Cluster",
            x = "Number of cells") +
          geom_text(aes(x = n + max(n) * 0.05, label = n, color = hivPos),
            position = position_dodge(0.8),
            vjust = 0.5,
            hjust = 0,
            size = BASEFONTSIZE) +
          scale_x_continuous(expand = c(0, 0), limits = c(0, max(.$n * 1.05))) +
          scale_y_discrete(drop = FALSE) +
          scale_fill_manual(values = c(HIVNEGCOLOR, HIVPOSCOLOR)) +
          scale_color_manual(values = c(HIVNEGCOLOR, HIVPOSCOLOR)) +
          ggDiscreteLollipopTheme}
    
  } else {
    p <- df %>%
      mutate(cluster = factor(cluster),
        hivPos = ifelse(secondMetadata, "Pos", "Neg")) %>%
      dplyr::filter(hivPos != "Neg") %>%
      dplyr::count(cluster) %>%
      mutate(proportion = n / sum(n)) %>%
      {ggplot(., aes(y = cluster)) +
          geom_segment(aes(x = 0, xend = proportion, y = cluster, yend = cluster),
            color = HIVPOSCOLOR, linetype = "dotted") +
          geom_point(aes(x = proportion), color = HIVPOSCOLOR, size = 1) + 
          labs(y = "Cluster",
            x = "Proportion of HIV+ cells") +
          geom_text(aes(x = proportion + 0.05,
              label = scales::label_percent(accuracy = 0.1)(proportion)),
            hjust = 0,
            vjust = 0.5,
            color = HIVPOSCOLOR,
            size = BASEFONTSIZE) +
          scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
          ggDiscreteLollipopTheme}
  }
  
  savePlot(plot = p, fn = fn, devices = devices, gheight = gheight, gwidth = gwidth)
}

plotVlnEnhanced <- function(
  seu,
  cells,
  feats,
  titles,
  fn,
  colorVar = "haystackOut",
  xVar = "individual",
  separator = NULL,
  showAggregate = TRUE,
  ncols = 5) {
  
  df <- data.frame(
    x = FetchData(seu, cells = cells, vars = xVar)[, 1],
    colorData = FetchData(seu, cells = cells, vars = colorVar)[, 1]
  )
  
  titleDict <- titles
  names(titleDict) <- feats
  
  df[, feats] <- t(seu@assays$tsa@data[feats, cells])
  df <- df %>%
    pivot_longer(cols = all_of(feats), names_to = "DNA_ID", values_to = "y") %>%
    mutate(cleanName = titleDict[DNA_ID]) %>%
    mutate(cleanName = factor(cleanName, levels = titles))
  
  if (!is.null(separator)) {
    sepDict <- separator
    names(sepDict) <- feats
    
    df <- df %>%
      mutate(sepData = sepDict[DNA_ID])
  }
  
  if (showAggregate) {
    origX <- unique(df$x)
    origNRow <- nrow(df)
    
    df <- df %>%
      slice(c(rep(1:n()), rep(1:n())))
    
    df[c((origNRow + 1) : nrow(df)), "x"] <- "Aggr"
    
    df$x <- factor(df$x, levels = c(origX, "Aggr"))
  }
  
  manualPalette <- c(HIVNEGCOLOR, HIVPOSCOLOR)
  halfWidths <- 0.9
  gSkeleton <- function(d, title, titleColor) {
    g <- ggplot(d, aes(x = x, y = y)) +
      geom_half_violin(aes(fill = colorData), color = "#000000", side = "l", width = halfWidths, lwd = 0.3) +
      geom_half_point(aes(fill = colorData, color = colorData), side = "r", width = halfWidths, size = 0.1, alpha = 0.3) +
      theme_classic() +
      theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = BASEPTFONTSIZE),
        plot.title.position = "plot",
        plot.title = element_text(size = BASEPTFONTSIZE, hjust = 0.5, color = titleColor, margin = margin(0,0,0,0)),
        strip.text = element_text(size = BASEPTFONTSIZE)) +
      scale_color_manual(values = manualPalette) +
      scale_fill_manual(values = manualPalette) +
      # scale_y_continuous(limits = c(-0.25, NA), expand = expansion(mult = c(0, 0.05))) +
      labs(y = "Expression level",
        title = title)
      
    
    # if (length(unique(d$cleanName)) > 1) {
    #   g <- g + facet_wrap(~ cleanName, ncol = ncols, scales = "free")
    # }
    
    return(g)
  }
  
  separatorDict <- separator
  names(separatorDict) <- titles
  gList <- lapply(titles, function(z) {
    df2 <- df %>%
      dplyr::filter(cleanName == z)
    
    return(gSkeleton(df2, z, ifelse(separatorDict[z] == "HIV-", HIVNEGCOLOR, HIVPOSCOLOR)))
  })
  
  # g <- gSkeleton(df)
  # g <- ggplot_gtable(ggplot_build(g))
  # stript <- which(grepl("strip-t", g$layout$name))
  # 
  # separatorDict <- separator
  # names(separatorDict) <- titles
  # 
  # for (s in stript) {
  #   lbl <- g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$label
  #   stripColor <- ifelse(separatorDict[lbl] == "HIV-", HIVNEGCOLOR, HIVPOSCOLOR)
  #   g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- stripColor
  # }
  
  g <- wrap_plots(gList, ncol = ncols)
  
  savePlot(plot = g, rdsPlot = gList, fn = fn, devices = c("rds", "png"),
    gheight = ceiling(length(feats) / ncols) * (8/ncols),
    gwidth = ifelse(length(feats) > ncols, 8, 8 / ncols * length(feats)))
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
    ylim(-1, max(singleAdtDf$normExp, na.rm = T) + 0.5) +
    geom_text(aes(x = cleanName, y = normExp + 0.1, label = label, hjust = hjust, angle = angle),
      color="black", alpha = 0.6, size = 2.25, inherit.aes = FALSE) +
    coord_polar()
  
  ggsave(filename = paste0(outputDir, "/", outputFn),
    plot = cellPlt, width = 8, height = 8,
    dpi = "retina", bg = "#FFFFFF")
  
  return(cellPlt)
}

makeDifferentialLollipop <- function(
  markers,
  fn,
  devices = c("png", "svg", "rds")) {
  
  pSapPi <- markers %>%
    arrange(desc(abs(piScore)), .by_group = TRUE) %>%
    mutate(cleanName = factor(cleanName, levels = rev(cleanName))) %>%
    mutate(markerColor = ifelse(piScore > 0, HIVPOSCOLOR, HIVNEGCOLOR)) %>%
    {ggplot(., aes(x = abs(piScore), y = cleanName, color = markerColor)) +
        geom_segment(aes(x = 0, xend = abs(piScore), y = cleanName, yend = cleanName), linetype = "dotted") +
        geom_point(size = 1.5) + 
        labs(y = "Surface marker",
          x = "π score") +
        scale_x_continuous(expand = c(0, 0), limits = c(0, max(abs(.$piScore) + 0.3))) +
        theme_classic() +
        theme(axis.text = element_text(size = BASEPTFONTSIZE),
          axis.title = element_text(size = BASEPTFONTSIZE),
          legend.position = "none",
          panel.background = element_rect(fill = NULL, colour = NULL),
          strip.background = element_blank(),
          strip.placement = "outside") +
        scale_color_identity() +
        coord_cartesian(clip = "off") +
        facet_grid(Status ~ ., scales = "free", space = "free")} 
  
  savePlot(plot = pSapPi, fn = fn, devices = devices, gheight = 2 + (nrow(markers) * 0.08), gwidth = 3)
}

plotMotifRank <- function(
  motifChanged,
  fn,
  nLabel = 10,
  color = NULL,
  devices = c("png", "rds"),
  direction = NULL,
  chromVARmode = FALSE) {
  
  
  if (chromVARmode) {
    df <- data.frame(
      y = log10(assays(motifChanged)$FDR[, 1]) * -1,
      meandiff = assays(motifChanged)$MeanDiff[, 1],
      motif = elementMetadata(motifChanged)$name)
    
    if (direction == "negative") {
      df <- filter(df, meandiff < 0)
    } else {
      df <- filter(df, meandiff > 0)
    }
    
    df <- df[order(df$y, decreasing = TRUE), ]
    df$x <- seq_len(nrow(df))
    
    yLbl <- "-log10(FDR)"
  } else {
    df <- data.frame(motif = rownames(motifChanged), y = assay(motifChanged)[,1])
    df <- df[order(df$y, decreasing = TRUE), ]
    df$x <- seq_len(nrow(df))
    
    yLbl <- "-log10(p adj)"
  }

  df <- df %>%
    mutate(color = ifelse(y > -log10(0.05), color, "#DDDDDD"))
  
  g <- ggplot(df, aes(x, y)) + 
    geom_hline(yintercept = -log10(0.05), color = "#333333", linetype = "dotted") +
    geom_point(size = 1, shape = 16, aes(color = color)) +
    ggrepel::geom_text_repel(
      data = df[rev(seq_len(nLabel)), ], aes(x = x, y = y, label = motif), 
      size = 5 / ggplot2:::.pt,
      max.overlaps = 30,
      force = 9,
      family = "Arial",
      direction = "y",
      nudge_x = 500,
      hjust = 0,
      ylim = c(-log10(0.05), max(df$y) * 1.02),
      segment.color = color,
      color = color,
      min.segment.length = 0.1
    ) +
    scale_color_identity() +
    theme_classic() +
    theme(
      axis.title = element_text(size = BASEPTFONTSIZE),
      axis.text = element_text(size = BASEPTFONTSIZE)) +
    labs(x = "Rank Sorted Motifs Enriched",
      y = yLbl)
  
  savePlot(plot = g, fn = fn, devices = devices, gheight = 2.5, gwidth = 2.5)
}

plotMotifDot <- function(
  markerFeatures,
  fn,
  direction = "positive",
  devices = c("png", "rds")
) {
  df <- data.frame(
    y = log10(assays(markerFeatures)$FDR[, 1]) * -1,
    meandiff = assays(markerFeatures)$MeanDiff[, 1],
    motif = elementMetadata(markerFeatures)$name)
  
  if (direction == "negative") {
    df <- filter(df, meandiff < 0)
    pointColor  <- HIVNEGCOLOR
    
  } else {
    df <- filter(df, meandiff > 0)
    pointColor <- HIVPOSCOLOR
  }
  
  df <- filter(df, y > -log10(0.05))
  df <- df[order(df$y, decreasing = TRUE), ]
  df$motif <- factor(df$motif, levels = df$motif)
  
  g <- ggplot(df, aes(x = motif, y = y)) +
    geom_point(color = pointColor) +
    theme_classic() +
    coord_cartesian(clip = "off") +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.text = element_text(size = BASEPTFONTSIZE),
      axis.title = element_text(size = BASEPTFONTSIZE)
    ) +
    labs(x = "Motif",
      y = "-log10(FDR)")
  
  savePlot(plot = g, fn = fn, devices = devices, gheight = 2, gwidth = nrow(df) / 5)
}


plotVolcanoFromGetMarkerFeatures <- function(
  markerFeatures,
  fn,
  posFoldChangeName = "Upregulated in HIV+",
  negFoldChageName = "Downregulated in HIV+",
  devices = c("png", "rds"),
  chromVARmode = FALSE
) {
  if (chromVARmode) {
      df <- data.frame(
        y = log10(assays(markerFeatures)$FDR[, 1]) * -1,
        x = assays(markerFeatures)$MeanDiff[, 1],
        motif = elementMetadata(markerFeatures)$name)
        
    xLbl <- "Mean difference in chromVAR\nbias-corrected deviations"
    yLbl <- "-log10(FDR)"
    
  } else {
    df <- data.frame(x = markerFeatures@assays@data@listData$Log2FC$x,
      y = log10(markerFeatures@assays@data@listData$FDR$x) * -1)
    
    xLbl <- "log2(Fold change)"
    yLbl <- "-log10(FDR)"
  }
  
  df <- df %>%
    mutate(pointColor = case_when(
      y > -1 * log10(0.05) & x > 0 ~ posFoldChangeName,
      y > -1 * log10(0.05) & x < 0 ~ negFoldChageName,
      TRUE ~ "Not significant"
    ))
  
  g <- df %>%  
    {ggplot(data = ., aes(x = x, y = y, color = pointColor)) +
        rasterize(geom_point(size = 0.3, alpha = 0.8), dpi = 300) +
        theme_classic() +
        labs(x = xLbl,
          y = yLbl,
          color = "") +
        scale_y_continuous(expand = expansion(mult = c(0,0.05)), limits = c(0, NA)) +
        scale_color_manual(values = c(HIVNEGCOLOR, "#dddddd", HIVPOSCOLOR))
    } +
    theme(legend.position = "bottom",
      axis.text = element_text(size = BASEPTFONTSIZE),
      axis.title = element_text(size = BASEPTFONTSIZE),
      legend.key.size = unit(BASEPTFONTSIZE, 'points'),
      legend.text = element_text(size = BASEPTFONTSIZE)) +
    guides(colour = guide_legend(override.aes = list(size = 1), ncol = 1))
  
  savePlot(plot = g, fn = fn, devices = devices, gheight = 3, gwidth = 3)
}

plotGeneAnnot <- function(geneAnnots) {
  g <- ggplot(geneAnnots, aes(x = startPos, xend = endPos, y = annotation, yend = annotation)) +
    geom_line(aes(group = annotation2), color = "#CCCCCC") +
    geom_segment(size = 1.5) +
    # geom_text(aes(x = (startPos + endPos) / 2, label = annotation), hjust = 0.5, position = position_nudge(y = 0.6), size = 2) +
    scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    labs(y = "Annotation") +
    theme(
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 6),
      axis.title.x = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(size = 6))
  
  return(g)
}

greaterThanZero <- function(d1, d2) {
  d2Check <- sapply(d2, function(x) {
    return(x > 0)
  })
  
  return(d1 + d2Check)
}

calculateCoverage <- function(frags, coverage, startCol = "startBp", endCol = "endBp") {
  uniqCbc <- unique(frags$newCbc)
  cbcCoverage <- vector("list", length(uniqCbc))
  names(cbcCoverage) <- uniqCbc
  
  for (i in seq_along(frags[, 1])) {
    cbc <- frags$newCbc[i]
    start <- frags[i, startCol] + 1
    end <- frags[i, endCol] + 1
    
    if (is.null(cbcCoverage[[cbc]])) {
      cbcCoverage[[cbc]] <- coverage
    }
    
    cbcCoverage[[cbc]][c(start:end)] <- cbcCoverage[[cbc]][c(start:end)] + 1
  }
  
  condensedCoverage <- Reduce(greaterThanZero, cbcCoverage, init = coverage)
  df <- data.frame(
    y = condensedCoverage / length(uniqCbc),
    x = seq(0, length(coverage) - 1)
  )
  
  return(df)
}

plotFragDistribution <- function(frags, coverage, inferredCoverage) {
  seqCoverage <- calculateCoverage(frags, coverage)
  inferredCoverage <- calculateCoverage(as.data.frame(inferredCoverage), coverage, startCol = "inferredStart", endCol = "inferredEnd")
  
  df <- bind_rows(list("sequenced" = seqCoverage, "inferred" = inferredCoverage), .id = "coverageType")
  
  g <- df %>%
    ggplot(aes(x = x, y = y, color = coverageType)) +
    geom_line(size = 0.3) +
    theme_classic() +
    coord_cartesian(clip = "off") +
    scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,NA)) +
    scale_color_manual(values = c("#1b9e7770", "#d95f0290")) +
    labs(y = "Prop. of cells") +
    theme(axis.line.y = element_blank(),
      axis.text = element_text(size = 6),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      legend.key = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = "#ffffff00", color = "#ffffff00"),
      legend.margin=margin(-1, 0, -1, 0, unit='lines'),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 6))
  
  return(g)
}

viralFragGraphTheme <- list(
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA)),
  theme_classic(),
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    axis.title.x = element_text(size = 6),
    panel.grid.major.y = element_line(color = "#EEEEEE80"),
    strip.background = element_rect(color = "#ffffff", fill = "#ffffff"),
    strip.text = element_text(size = BASEPTFONTSIZE)),
  labs(x = "Read fragment aligned to provirus (bp)")
)

plotFragCoverage <- function(frags, inferredCoverage, separateByIndividual) {
  nCells <- length(unique(frags$newCbc))
  if (nCells > 100) {
    lineSize <- 0.1
  } else if (nCells > 50) {
    lineSize <- 0.5
  } else {
    lineSize <- 1
  }
  
  g <- ggplot(frags, aes(x = startBp, xend = endBp, y = newCbc, yend = newCbc)) +
    geom_segment(color = "#d95f0290", size = lineSize) +
    geom_segment(data = inferredCoverage, aes(x = inferredStart, xend = inferredEnd, y = newCbc, yend = newCbc),
      color = "#1b9e7770", size = lineSize) +
    labs(y = "Cell")
  
  if (separateByIndividual) {
    g <- g +
      facet_grid(individual ~ ., scales = "free_y", space = "free_y")
  }
  
  g <- g + viralFragGraphTheme
  
  return(g)
}

plotFragMultiGraph <- function(
  frags,
  coverage,
  geneAnnots,
  fn,
  gwidth = 4,
  gheight = 7.5,
  separateByIndividual = FALSE,
  devices = c("png", "rds")
) {
  
  if (separateByIndividual) {
    inferredCoverage <- frags %>% 
      group_by(individual, readname, newCbc)
  } else {
    inferredCoverage <- frags %>% 
      group_by(readname, newCbc)    
  }
  
  inferredCoverage <- inferredCoverage %>%
    filter(n() == 2) %>%
    arrange(min(startBp, endBp), .by_group = TRUE) %>%
    filter(first(endBp) < last(startBp)) %>%
    summarise(
      inferredStart = first(endBp) + 1,
      inferredEnd = last(startBp) - 1)
  
  annotG <- plotGeneAnnot(geneAnnots)
  distG <- plotFragDistribution(frags, coverage, inferredCoverage)
  fragG <- plotFragCoverage(frags, inferredCoverage, separateByIndividual)
  
  p <- fragG / distG / annotG + 
    plot_layout(heights = c(6, 0.5, 0.75)) & 
    theme(text = element_text(family = "Arial"))
  
  savePlot(plot = p, fn = fn, devices = devices, gheight = gheight, gwidth = gwidth)
}

plotArchRQCData <- function(
  proj,
  fn,
  xcol = "as.character(haystackOut)",
  ycol,
  xlbl = NULL,
  ylbl = ycol,
  gwidth = 3,
  gheight = 3,
  devices = c("png", "rds")
) {
  df <- data.frame(yaxis = getCellColData(proj, ycol)[,1],
    xaxis = getCellColData(proj, xcol)[,1],
    cell = getCellNames(proj))
  
  if (xcol == "as.character(haystackOut)") {
    df <- df %>%
      mutate(xaxis = ifelse(xaxis, "HIV+", "HIV-"))
  }
  
  g <- ggplot(df, aes(x = xaxis, y = yaxis, color = xaxis)) +
    geom_violin(aes(fill = xaxis), alpha = 0.4) +
    geom_boxplot(fill = "#FFFFFF00", outlier.size = 0.8) +
    theme_classic() +
    labs(y = ylbl) +
    theme(legend.position = "none",
      axis.text = element_text(size = BASEPTFONTSIZE, family = "Arial"),
      axis.title = element_text(size = BASEPTFONTSIZE, family = "Arial"),
      axis.title.x = element_blank())
  
  if (xcol == "as.character(haystackOut)") {
    g <- g + 
      stat_compare_means(method = "t.test", label.x.npc = 0.5, label.y.npc = 0.85, hjust = 0.5, size = BASEFONTSIZE, label = "p.format") +
      scale_color_manual(values = c(HIVNEGCOLOR, HIVPOSCOLOR)) +
      scale_fill_manual(values = c(HIVNEGCOLOR, HIVPOSCOLOR))
  }
  
  savePlot(plot = g, fn = fn, devices = devices, gheight = gheight, gwidth = gwidth)
  
  return(df)
}

plotUpsetCBC <- function(
  atacCBC,
  hivCBC,
  adtCBC,
  fn,
  gheight = 2.5,
  gwidth = 2.5,
  devices = c("png", "rds")
) {
  df <- data.frame(
    cbc = c(atacCBC, hivCBC, adtCBC),
    modal = c(rep("ATAC", length(atacCBC)), rep("HIV", length(hivCBC)), rep("ADT", length(adtCBC))),
    value = TRUE
  )
  
  df <- df %>%
    pivot_wider(id_cols = cbc, names_from = modal, values_from = value) %>%
    group_by(ATAC, HIV, ADT) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    mutate(rowID = factor(row_number()))
  
  barPlot <- df %>%
    ggplot(aes(x = rowID, y = n)) +
    geom_bar(stat = "identity", width = 0.5, fill = "#000000") +
    geom_text(aes(label = n), nudge_y = 250, size = BASEFONTSIZE) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(clip = "off") +
    labs(y = "Count") +
    theme_classic() +
    theme(
      axis.text = element_text(size = BASEPTFONTSIZE),
      axis.title.y = element_text(size = BASEPTFONTSIZE),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank())
  
  axisPlot <- df %>%
    dplyr::select(-n) %>%
    pivot_longer(cols = c(ATAC, HIV, ADT), names_to = "modals", values_to = "presence") %>%
    mutate(lineGroup = ifelse(is.na(presence), paste(rowID, modals, NA), rowID)) %>%
    mutate(colorGroup = ifelse(is.na(presence), "#CCCCCC", "#000000")) %>%
    ggplot(aes(x = rowID, y = modals)) +
    geom_point(aes(color = colorGroup), size = 3) +
    geom_line(aes(group = lineGroup)) +
    theme_classic() +
    theme(axis.line = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = BASEPTFONTSIZE),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major.y = element_line(size = 7, color = "#EFEFEF80")) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_color_identity() +
    coord_cartesian(clip = "off")
  
  p <- barPlot / axisPlot + plot_layout(heights = c(5, 1)) & theme(axis.text = element_text(family = "Arial"))
  
  savePlot(plot = p, fn = fn, devices = devices, gheight = gheight, gwidth = gwidth)
}

cleanUpTrackAndSave <- function(
  archrTrack,
  fn,
  gheight = 2.25,
  gwidth = 4,
  devices = c("png", "rds")
) {
  # plot the raw figure
  p <- wrap_plots(archrTrack, ncol = 1, heights = c(3,2,1))
  savePlot(p, fn = paste0(fn, "_orig"), devices = "png", gheight = gheight, gwidth = gwidth)
  
  # clean up theme
  cleanUpTheme <- theme(
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    text = element_text(family = "Arial", size = 6),
    axis.text = element_text(family = "Arial", size = 6),
    panel.border = element_blank())
  
  # clean up
  archrTrack$bulktrack <- archrTrack$bulktrack +
    labs(y = "Normalized signal") +
    cleanUpTheme
  
  archrTrack$sctrack <- archrTrack$sctrack +
    labs(y = "Binarized signal") +
    cleanUpTheme +
    scale_color_manual(values = c(HIVNEGCOLOR, HIVPOSCOLOR))
  
  archrTrack$genetrack <- archrTrack$genetrack +
    cleanUpTheme +
    labs(y = "Gene") +
    scale_color_manual(values = c("blue", "orange"))
  
  savePlot(archrTrack, customSavePlot = wrap_plots(archrTrack, ncol = 1, heights = c(2,1.5,1.5)),
    fn = fn, devices = devices, gheight = gheight, gwidth = gwidth)
}