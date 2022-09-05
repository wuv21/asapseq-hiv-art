library(gghalves)

######
# COLOR SCHEMES
######
HIVPOSCOLOR <- "#e63946"
HIVNEGCOLOR <- "#999999"
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


customSortAnnotation <- function(x, reverse = FALSE, ignoreClusterPrefix = TRUE) {
  priority <- c("MAIT", "Tfh", "Naive", "Tcm", "Tem", "Activated", "Treg")
  
  origNames <- x
  if (ignoreClusterPrefix) {
    x <- gsub("C\\d+: ", "", x)
    
    mapNames <- origNames
    names(mapNames) <- x
  }
  
  x <- sort(unique(x))
  
  cd4Indices <- grepl("^CD4", x)
  
  oldCd4 <- x[cd4Indices]
  newCd4 <- c()
  
  for (p in priority) {
    pIndices <- grepl(p, oldCd4)
    newCd4 <- append(newCd4, oldCd4[pIndices])
    oldCd4 <- oldCd4[!pIndices]
  }
  
  if (length(oldCd4) > 0) {
    newCd4 <- append(newCd4, oldCd4)
  }

  finalSort <- c(newCd4, x[!cd4Indices])
  
  if (ignoreClusterPrefix) {
    finalSort <- mapNames[finalSort]
  }
  
  return(finalSort)
}

######
# Generate UMAP df from ArchRProject
######
generateUmapDfFromArchR <- function(
  proj,
  secondGraphColumn = "haystackOut",
  cluster = "Clusters",
  customSort = FALSE,
  donorColumn = NULL,
  colorLabelCluster = NULL,
  embedding = "UMAP") {
  
  umapFromArchr <- getEmbedding(proj, embedding = embedding)
  
  # check to make sure df is the same size as expected
  stopifnot(nrow(umapFromArchr) == length(proj$cellNames))
  
  # check to make sure embedding order is same as project order
  if (sum(rownames(umapFromArchr) == proj$cellNames) != length(proj$cellNames)) {
    message("Embedding rownames from getEmbedding() not in the same order as ArchR project's cell names. Will fix now.")
    
    # sort embedding df by proj$cellNames
    umapFromArchr <- umapFromArchr[proj$cellNames, ]
  }
  
  # check to make sure that all of getCellColData is in the same order too
  stopifnot(
    sum(rownames(getCellColData(proj, select = secondGraphColumn)) == proj$cellNames) == length(proj$cellNames),
    sum(rownames(getCellColData(proj, select = cluster)) == proj$cellNames) == length(proj$cellNames)
  )
  
  df <- data.frame(x = umapFromArchr[, 1],
    y = umapFromArchr[, 2],
    secondMetadata = getCellColData(proj, select = secondGraphColumn)[, 1],
    sample = proj$Sample)
    
  if (is.vector(cluster) && length(cluster) > 1) {
    tmp <- getCellColData(proj, select = cluster)
    tmp <- tidyr::unite(as.data.frame(tmp), col = "cluster", sep = ": ")
    
    df$cluster <- tmp$cluster
    
  } else {
    tmp <- data.frame(cluster = getCellColData(proj, select = cluster)[, 1])
    
    if (customSort) {
      tmp$cluster <- factor(tmp$cluster, levels = customSortAnnotation(tmp$cluster))
    }
    
    df$cluster <- tmp$cluster
  }
  
  if (!is.null(colorLabelCluster)) {
    tmp <- getCellColData(proj, select = colorLabelCluster)[, 1]
  
    stopifnot(rownames(tmp) != proj$cellNames)
    df$colorLabelCluster <- tmp
  }
  
  if (!is.null(donorColumn)) {
    tmp <- getCellColData(proj, select = donorColumn)[, 1]
    
    stopifnot(rownames(tmp) != proj$cellNames)
    df$donor <- tmp
  }
  
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
  colorLabelBy = colorBy,
  colorLabel = colorBy,
  colorScheme = NULL,
  rasterize = TRUE,
  bringToTop = FALSE,
  customClusterSort = FALSE,
  propInLegend = FALSE,
  propDigits = 1,
  embedding = "UMAP") {
  
  df <- generateUmapDfFromArchR(proj,
    cluster = colorBy,
    embedding = embedding,
    colorLabelCluster = colorLabelBy,
    customSort = customClusterSort)
  
  if (bringToTop) {
    df <- df %>%
      dplyr::arrange(cluster)
  }
  
  if (length(colorBy) == 1 && colorBy == "haystackOut") {
    df <- df %>%
      mutate(cluster = ifelse(cluster, "HIV+", "HIV-"))
  }
  
  if (propInLegend) {
    df <- df %>%
      group_by(cluster) %>%
      mutate(cluster = glue("{cluster} ({round(n() / nrow(.) * 100, digits = propDigits)}%)"))
  }
  
  if (customClusterSort) {
    # TODO fix because cluster has more info...
    
    df$cluster <- factor(df$cluster, levels = customSortAnnotation(df$cluster))
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
      color = colorLabel) +
    theme_classic() +
    ggtheme +
    guides(colour = guide_legend(override.aes = list(size = 4), nrow = 5))
  
  if (is.null(colorScheme) && colorBy == "haystackOut") {
    p1 <- p1 + scale_color_manual(values = c(HIVNEGCOLOR, HIVPOSCOLOR))
  } else if (!is.null(colorScheme)) {
    p1 <- p1 + colorScheme
  }
  
  if (!is.null(colorBy) & !is.null(colorLabelBy)) {
    clusterLabelUmapPos <- df %>% 
      group_by(colorLabelCluster) %>% 
      dplyr::summarize(
        topX = quantile(x, c(.6)),
        topY = quantile(y, c(.6)),
        bottomX = quantile(x, c(.4)),
        bottomY = quantile(y, c(.4)),
        x = (topX + bottomX) / 2,
        y = (topY + bottomY) / 2)
    
    p1 <- p1 + 
      ggrepel::geom_label_repel(data = clusterLabelUmapPos,
        aes(x = x, y = y, label = colorLabelCluster),
        label.size = 0.05,
        force = 300,
        max.time = 5,
        max.iter = 25000,
        size = BASEFONTSIZE,
        seed = 21,
        min.segment.length = 0.25,
        segment.color = "#BBBBBB")
  }
  
  savePlot(plot = p1, fn = fn, devices = devices, gheight = 4, gwidth = 3.5)
}


ggDiscreteLollipopTheme <- list(
  theme_classic(),
  coord_cartesian(clip = "off"),
  theme(legend.position = "none",
    axis.title = element_text(size = BASEPTFONTSIZE),
    axis.text = element_text(size = BASEPTAXISFONTSIZE, color = "#000000"),
    strip.background = element_rect(fill = "transparent", colour = "#000000", size = 1),
    panel.spacing = unit(1.5, "line"),
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
  donorColumn = NULL,
  donorColumnOrder = NULL,
  showAggregate = FALSE,
  graphType = "proportion") {
  
  if (is.null(donorColumn)) {
    df <- generateUmapDfFromArchR(proj,
                                  cluster = cluster,
                                  secondGraphColumn = secondGraphColumn,
                                  customSort = TRUE)
  } else {
    df <- generateUmapDfFromArchR(proj,
                                  cluster = cluster,
                                  secondGraphColumn = secondGraphColumn,
                                  donorColumn = donorColumn,
                                  customSort = TRUE)
  }
  
  stopifnot(graphType %in% c("absolute", "hivOnly"))
  
  if (!is.null(donorColumn) & !is.null(donorColumnOrder)) {
    df <- df %>%
      mutate(donor = factor(donor, levels = donorColumnOrder))
  } else if(!is.null(donorColumn)) {
    uniqueDonors <- sort(unique(df$donor))
    df <- df %>%
      mutate(donor = factor(donor, levels = uniqueDonors))
  }

  if (graphType == "absolute") {
    dfg <- df %>%
      mutate(hivPos = factor(ifelse(secondMetadata, "Pos", "Neg")))
    
    absGraphSettings <- list(
      geom_linerange(aes(xmin = 0, xmax = n, y = cluster, color = hivPos),
        linetype = "dotted", position = position_dodge(0.8)),
      geom_point(aes(x = n, color = hivPos), size = 1, position = position_dodge(0.8)),
      labs(y = "Cluster",
        x = "Number of cells"),
      geom_text(aes(x = n + max(n) * 0.05, label = n, color = hivPos),
        position = position_dodge(0.8),
        vjust = 0.5,
        hjust = 0,
        size = BASEFONTSIZE),
      scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)),
      scale_y_discrete(drop = FALSE, limits = rev),
      scale_fill_manual(values = c(HIVNEGCOLOR, HIVPOSCOLOR)),
      scale_color_manual(values = c(HIVNEGCOLOR, HIVPOSCOLOR)),
      ggDiscreteLollipopTheme
    )
    
    if (is.null(donorColumn)) {
      dfg <- dfg %>% dplyr::count(cluster, hivPos, .drop = FALSE)
      p <- ggplot(dfg, aes(y = cluster, fill = hivPos)) +
        absGraphSettings
      
    } else {
      dfg <- dfg %>% dplyr::count(donor, cluster, hivPos, .drop = FALSE)
      
      pList <- lapply(seq_along(levels(dfg$donor)), function(i) {
        d <- levels(dfg$donor)[i]
        
        dfgDonor <- dfg %>%
          dplyr::filter(donor == d)
        
        g <- ggplot(dfgDonor, aes(y = cluster)) +
          absGraphSettings +
          labs(subtitle = d) +
          theme(
            plot.subtitle = element_text(size = BASEPTFONTSIZE, hjust = 0.5),
            axis.title.x = element_blank())
        
        if (i != 1) {
          g <- g +
            theme(
              axis.text.y = element_blank(),
              axis.title.y = element_blank())
        }
        
        return(g)
      })
    }
    
  } else {
    dfg <- df %>%
      mutate(hivPos = ifelse(secondMetadata, "Pos", "Neg")) %>%
      dplyr::filter(hivPos != "Neg")
    
    if (is.null(donorColumn)) {
      dfg <- dfg %>% dplyr::count(cluster)
    } else {
      dfg <- dfg %>% group_by(donor) %>% dplyr::count(cluster)
    }
    
    if (showAggregate) {
      dfAggr <- dfg %>%
        ungroup() %>%
        group_by(cluster) %>%
        dplyr::summarize(n = sum(n)) %>%
        mutate(donor = "Aggregate")
      
      donorLvls <- c(levels(dfg$donor), "Aggregate")
      
      dfg <- dfg %>% 
        bind_rows(dfAggr) %>%
        mutate(donor = factor(donor, levels = donorLvls))
    }
    
    dfg <- dfg %>%
      mutate(proportion = n / sum(n))
    
    propGraphSettings <- list(
      geom_segment(aes(x = 0, xend = proportion, y = cluster, yend = cluster),
        color = HIVPOSCOLOR, linetype = "dotted"),
      geom_point(aes(x = proportion), color = HIVPOSCOLOR, size = 1),
      labs(y = "Cluster",
        x = "Proportion of HIV+ cells"),
      geom_text(aes(x = proportion + 0.05,
        label = scales::label_percent(accuracy = 0.1)(proportion)),
        hjust = 0,
        vjust = 0.5,
        color = HIVPOSCOLOR,
        size = BASEFONTSIZE),
      ggDiscreteLollipopTheme
    )
    
    if (!is.null(donorColumn)) {
      pList <- lapply(seq_along(levels(dfg$donor)), function(i) {
        d <- levels(dfg$donor)[i]
        
        dfgDonor <- dfg %>%
          dplyr::filter(donor == d)
        
        g <- ggplot(dfgDonor, aes(y = cluster)) +
          propGraphSettings +
          labs(subtitle = d) +
          scale_y_discrete(limits = rev, drop = FALSE) +
          scale_x_continuous(limits = c(0, 1.05), labels = scales::label_percent(suffix = "")) +
          theme(
            plot.subtitle = element_text(size = BASEPTFONTSIZE, hjust = 0.5),
            axis.title.x = element_blank())
        
        if (i != 1) {
          g <- g +
            theme(
              axis.text.y = element_blank(),
              axis.title.y = element_blank())
        }
        
        return(g)
      })
      
    } else {
      p <- ggplot(dfg, aes(y = cluster)) +
        propGraphSettings +
        scale_y_discrete(limits = rev)
    }
  }
  
  if (!is.null(donorColumn)) {
    savePlot(plot = wrap_plots(pList, nrow = 1),
      rdsPlot = pList, fn = fn, devices = devices, gheight = gheight, gwidth = gwidth)
  } else {
    savePlot(plot = p, fn = fn, devices = devices, gheight = gheight, gwidth = gwidth)
  }
  
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
    tidyr::pivot_longer(cols = all_of(feats), names_to = "DNA_ID", values_to = "y") %>%
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
      rasterize(geom_half_point(aes(fill = colorData, color = colorData), side = "r", width = halfWidths, size = 0.1, alpha = 0.3), dpi = 300) +
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
makeCellPlot <- function(seu,
  tsa_catalog,
  cellID,
  outputDir,
  labelFeatures = FALSE,
  colorByHaystack = NULL,
  plotMargin = -0.5,
  yLimMin = -1,
  outputFn = glue("cellplot_{cellID}.png")) {
  
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
  
  
  if (!is.null(colorByHaystack)) {
    singleAdtDf$barFill <- ifelse(colorByHaystack, HIVPOSCOLOR, HIVNEGCOLOR)
  }
  
  cellPlt <- ggplot(singleAdtDf, aes(x = cleanName, y = normExp)) +
    geom_bar(stat = "identity", alpha = 0.8, aes(fill = barFill), width = 0.4) +
    scale_fill_identity() +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(plotMargin, 4), "in")
    ) +
    ylim(yLimMin, max(singleAdtDf$normExp, na.rm = T) + 0.5) +
    coord_polar()
    
  if (labelFeatures) {
    cellPlt <- cellPlt +
      geom_text(aes(x = cleanName, y = normExp + 0.1, label = label, hjust = hjust, angle = angle),
        color="black", alpha = 0.6, size = 2.25, inherit.aes = FALSE)
  }
  
  if (!is.null(outputFn)) {
    ggsave(filename = paste0(outputDir, "/", outputFn),
      plot = cellPlt, width = 5, height = 5,
      dpi = "retina", bg = "transparent")
  }
  
  return(cellPlt)
}

makeDifferentialLollipop <- function(
  markers,
  fn,
  devices = c("png", "rds")) {
  
  pSapPi <- markers %>%
    dplyr::arrange(dplyr::desc(abs(piScore)), .by_group = TRUE) %>%
    mutate(cleanName = factor(cleanName, levels = rev(cleanName))) %>%
    mutate(markerColor = ifelse(piScore > 0, HIVPOSCOLOR, HIVNEGCOLOR)) %>%
    {ggplot(., aes(x = abs(piScore), y = cleanName, color = markerColor)) +
        geom_segment(aes(x = 0, xend = abs(piScore), y = cleanName, yend = cleanName), linetype = "dotted") +
        geom_point(size = 1.5) + 
        labs(y = "Surface marker",
          x = "π score") +
        scale_x_continuous(expand = c(0, 0), limits = c(0, max(abs(.$piScore) + 0.3))) +
        theme_classic() +
        theme(axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
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

makeDifferentialLollipop2 <- function(
    markers,
  fn,
  devices = c("png", "rds")) {

  pSapPi <- markers %>%
    dplyr::arrange(dplyr::desc(abs(piScore)), .by_group = TRUE) %>%
    mutate(cleanName = factor(cleanName, levels = rev(cleanName))) %>%
    mutate(markerColor = ifelse(piScore > 0, HIVPOSCOLOR, HIVNEGCOLOR)) %>%
    mutate(graphingP = -log10(p_val_adj)) %>%
    {ggplot(., aes(x = graphingP, y = cleanName, fill = avg_log2FC)) +
        geom_segment(aes(x = 0, xend = graphingP, y = cleanName, yend = cleanName), linetype = "dotted", color = "#555555") +
        geom_point(size = 1.75, shape = 21, color = "#000000") + 
        labs(y = "Surface marker",
          x = "-log10(p adj)",
          fill = "log2(fold change)") +
        guides(fill = guide_colorbar(
          title.position = "left",
          title.hjust = 1,
          title.vjust = 1,
          barwidth = 4,
          barheight = 0.36)) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, max(abs(.$graphingP) + 0.3))) +
        theme_classic() +
        theme(
          axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
          axis.title = element_text(size = BASEPTFONTSIZE),
          legend.position = "bottom",
          legend.margin = margin(-5, 10, -5, -30, "pt"),
          plot.margin = unit(c(-20, 0, 0, 0), "pt"),
          legend.title = element_text(size = BASEPTFONTSIZE, margin = margin(0, 0, 0, 0)),
          legend.text = element_text(size = BASEPTFONTSIZE - 1, angle = 90, hjust = 1, vjust = 0.5),
          panel.background = element_rect(fill = NULL, colour = NULL),
          strip.background = element_blank(),
          strip.placement = "outside") +
        scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") +
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
  
  yMax <- max(df$y)
  yMin <- -log10(0.05) * 2
  yDiff <- (yMax - yMin) / (nLabel - 1)
  xMin <- 500
  
  dfLbl <- df[rev(seq_len(nLabel)), ] %>%
    mutate(yLbl = yDiff * (nLabel - x) + yMin,
      xLbl = xMin)
  
  g <- ggplot(df, aes(x, y)) + 
    geom_hline(yintercept = -log10(0.05), color = "#333333", linetype = "dotted") +
    geom_segment(data = dfLbl, aes(y = y, yend = yLbl, x = x, xend = xLbl - 25), color = color, size = 0.5) +
    geom_point(size = 1, shape = 16, aes(color = color)) +
    geom_text(data = dfLbl, aes(label = motif, y = yLbl, x = xLbl), hjust = 0, color = color, size = 2) +
    scale_color_identity() +
    theme_classic() +
    theme(
      axis.title = element_text(size = BASEPTFONTSIZE),
      axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000")) +
    labs(x = "Rank Sorted Motifs Enriched",
      y = yLbl)
  
  savePlot(plot = g, fn = fn, devices = devices, gheight = 2.5, gwidth = 2.5)
}

plotEnhancedMotifDot <- function(
  markerFeatures,
  fn,
  xVar,
  yVar,
  xLbl = rlang::as_name(xVar),
  yLbl = rlang::as_name(yVar), 
  colorVar,
  subtitle = paste0(xLbl, " vs ", yLbl),
  preFilter = NULL,
  devices = c("png", "rds"),
  alphaByPVal = NULL,
  pValMax = 0.05,
  showTopNByPVal = NULL,
  palette = "YlOrRd",
  paletteDirection = 1
) {
  
  .xVar <- rlang::parse_expr(xVar)
  .yVar <- rlang::parse_expr(yVar)
  .colorVar <- rlang::parse_expr(colorVar)
  pValMax <- -1 * log10(pValMax)
  
  df <- data.frame(
    fdr = -1 * log10(assays(markerFeatures)$FDR[, 1]),
    meandiff = assays(markerFeatures)$MeanDiff[, 1],
    motif = elementMetadata(markerFeatures)$name)
  
  if (!is.null(preFilter)) {
    preFilter <- rlang::parse_expr(preFilter)
    df <- dplyr::filter(df, !!preFilter)
  }
  
  if (!is.null(showTopNByPVal)) {
    df <- df[order(df$fdr, decreasing = TRUE), ]
    df <- df[c(1:showTopNByPVal), ]
  } else {
    df <- filter(df, fdr > pValMax)
    df <- df[order(df$fdr, decreasing = TRUE), ]
  }
  df$motif <- factor(df$motif, levels = df$motif)
  
  g <- ggplot(df, aes(x = !!.xVar, y = !!.yVar))
  
  if (!is.null(showTopNByPVal) & yLbl == "fdr") {
    g <- g + geom_hline(yintercept = pValMax, linetype = "dotted", color = "#00000075")
  }
  
  g <- g +
    geom_point(aes(color = !!.colorVar)) +
    theme_classic() +
    coord_cartesian(clip = "off") +
    scale_color_distiller(palette = palette, direction = paletteDirection) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
      axis.title = element_text(size = BASEPTFONTSIZE),
      text = element_text(size = BASEPTFONTSIZE),
      legend.title = element_text(size = BASEPTFONTSIZE),
      legend.text = element_text(size = BASEPTFONTSIZE),
      legend.key.width = unit(BASEPTFONTSIZE, 'pt'),
      legend.key.height = unit(BASEPTFONTSIZE * 2, 'pt'),
      legend.position = "right"
    ) +
    labs(
      x = xLbl,
      y = yLbl,
      subtitle = subtitle)
  
  gwidth <- ifelse(nrow(df) >= 5, nrow(df) / 4.5, 3)
  savePlot(plot = g, fn = fn, devices = devices, gheight = 2, gwidth = gwidth)
}



plotMotifDot <- function(
  markerFeatures,
  fn,
  direction = "positive",
  devices = c("png", "rds"),
  alphaBySig = FALSE,
  colorBySig = FALSE,
  pValMax = 0.05,
  showTopNByPVal = NULL
) {
  df <- data.frame(
    y = log10(assays(markerFeatures)$FDR[, 1]) * -1,
    meandiff = assays(markerFeatures)$MeanDiff[, 1],
    motif = elementMetadata(markerFeatures)$name)
  
  # TODO change to colorByMeanDiff...
  
  # add fill....
  if (direction == "negative") {
    df <- filter(df, meandiff < 0)
    pointColor <- HIVNEGCOLOR
    
  } else {
    df <- filter(df, meandiff > 0)
    pointColor <- HIVPOSCOLOR
  }
  
  if (alphaBySig) {
    df <- df %>%
      mutate(pointColor = ifelse(y > -log10(pValMax), pointColor, paste0(pointColor, "50")))
  }
  
  if (!is.null(showTopNByPVal)) {
    df <- df[order(df$y, decreasing = TRUE), ]
    df <- df[c(1:showTopNByPVal), ]
  } else {
    df <- filter(df, y > -log10(pValMax))
    df <- df[order(df$y, decreasing = TRUE), ]
  }
  df$motif <- factor(df$motif, levels = df$motif)
  
  g <- ggplot(df, aes(x = motif, y = y))
  
  if (!is.null(showTopNByPVal)) {
    g <- g + geom_hline(yintercept = -log10(pValMax), linetype = "dotted", color = "#00000075")
  }
  
  g <- g +
    geom_point(aes(color = pointColor)) +
    theme_classic() +
    coord_cartesian(clip = "off") +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
      axis.title = element_text(size = BASEPTFONTSIZE)
    ) +
    scale_color_identity() +
    labs(x = "Motif",
      y = "-log10(FDR)")
  
  gwidth <- ifelse(nrow(df) >= 5, nrow(df) / 5, 3)
  savePlot(plot = g, fn = fn, devices = devices, gheight = 2, gwidth = gwidth)
}


plotVolcanoFromGetMarkerFeatures <- function(
  markerFeatures,
  fn,
  posFoldChangeName = "Upregulated in HIV+",
  negFoldChageName = "Downregulated in HIV+",
  devices = c("png", "rds"),
  chromVARmode = FALSE,
  returnDf = FALSE,
  labelTopN = NULL,
  filterOutlier = FALSE,
  isEncode = FALSE
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
  
  if (isEncode) {
    motifSep <- str_match(df$motif, "(\\d+).(.*)-(.*)...")
    df$motif <- paste0(motifSep[, 3], " (", motifSep[, 4], ")")
  }
  
  df <- df %>%
    mutate(pointColor = case_when(
      y > -1 * log10(0.05) & x > 0 ~ posFoldChangeName,
      y > -1 * log10(0.05) & x < 0 ~ negFoldChageName,
      TRUE ~ "Not significant"
    )) %>%
    mutate(pointColor = factor(pointColor, levels = c(negFoldChageName, "Not significant", posFoldChangeName))) %>%
    filter(!is.na(y) & !is.na(x) & !is.nan(x))
  
  if (returnDf) {
    return(df)
  }
  
  if (filterOutlier) {
    origLength <- nrow(df)
    df <- df %>%
      dplyr::filter(y < 10)
    
    print(paste0("Note that ", origLength - nrow(df), " outliers were filtered out for graphical purposes."))
  }
  
  hivColors <- c(HIVNEGCOLOR, "#dddddd", HIVPOSCOLOR)
  names(hivColors) <- levels(df$pointColor)
  hivColorScale <- scale_colour_manual(values = hivColors)
  
  g <- df %>%  
    {ggplot(data = ., aes(x = x, y = y)) +
        geom_hline(yintercept = -log10(0.05), color = "#00000070", linetype = "dotted") +
        rasterize(geom_point(aes(color = pointColor), size = 0.3, alpha = 0.8), dpi = 300) +
        theme_classic() +
        labs(x = xLbl,
          y = yLbl,
          color = "") +
        scale_y_continuous(expand = expansion(mult = c(0,0.05)), limits = c(0, NA)) +
        hivColorScale
    } +
    theme(legend.position = "bottom",
      axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
      axis.title = element_text(size = BASEPTFONTSIZE),
      legend.key.size = unit(BASEPTFONTSIZE, 'points'),
      legend.text = element_text(size = BASEPTFONTSIZE)) +
    guides(colour = guide_legend(override.aes = list(size = 1), ncol = 1))
  
  if (!is.null(labelTopN)) {
    dfTopN <- df %>%
      mutate(direction = x > 0) %>%
      group_by(direction) %>%
      dplyr::top_n(labelTopN, wt = y)
    
    yMax <- max(df$y) * 1.05
    yMin <- 0.1
    
    dfLbl <- dfTopN %>%
      mutate(xLbl = ifelse(x > 0, max(df$x) + 0.1, min(df$x) - 0.1)) %>%
      group_by(direction) %>%
      dplyr::arrange(desc(y)) %>%
      mutate(gid = seq_along(motif)) %>%
      mutate(yLbl = (yMax - yMin) / (max(gid) - 1) * (max(gid) - gid) + yMin)
    
    g <- g +
      scale_x_continuous(expand = expansion(mult = c(0.4,0.4))) +
      geom_text(data = dfLbl,
        aes(x = xLbl, y = yLbl, label = motif, hjust = ifelse(x > 0, 0, 1)),
        color = "#000000",
        size = 1.75) +
      geom_segment(data = dfLbl,
        aes(x = xLbl * 0.98, xend = x, y = yLbl, yend = y),
        color = "#00000030", size = 0.25)
  }
  
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

calculateCoverage <- function(frags, coverage, totalUniqCells = NULL, startCol = "startBp", endCol = "endBp") {
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
  
  denominator <- ifelse(is.null(totalUniqCells), length(uniqCbc), totalUniqCells)
  
  condensedCoverage <- Reduce(greaterThanZero, cbcCoverage, init = coverage)
  df <- data.frame(
    y = condensedCoverage / denominator,
    x = seq(0, length(coverage) - 1)
  )
  
  return(df)
}

plotFragDistribution <- function(frags, coverage, inferredCoverage) {
  seqCoverage <- calculateCoverage(frags, coverage)
  inferredCoverage <- calculateCoverage(
    as.data.frame(inferredCoverage),
    coverage,
    totalUniqCells = length(unique(frags$newCbc)),
    startCol = "inferredStart",
    endCol = "inferredEnd")
  
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
  gwidth = 3.5,
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
    dplyr::arrange(min(startBp, endBp), .by_group = TRUE) %>%
    filter(dplyr::first(endBp) < dplyr::last(startBp)) %>%
    dplyr::summarise(
      inferredStart = dplyr::first(endBp) + 1,
      inferredEnd = dplyr::last(startBp) - 1)
  
  annotG <- plotGeneAnnot(geneAnnots)
  distG <- plotFragDistribution(frags, coverage, inferredCoverage)
  fragG <- plotFragCoverage(frags, inferredCoverage, separateByIndividual)
  
  p <- fragG / distG / annotG + 
    plot_layout(heights = c(6, 0.5, 0.75)) & 
    theme(text = element_text(family = "Arial"))
  
  savePlot(plot = p, fn = fn, devices = devices, gheight = gheight, gwidth = gwidth)
}


.spaceFacetRelative <- function(g, relHeights) {
  tmpG <- ggplotGrob(g)
  ind <- tmpG$layout$t[grepl("panel", tmpG$layout$name)]
  tmpG$heights[ind] <- unit(relHeights, "null")
  
  for(i in which(grepl("strip-l", tmpG$layout$name))){
    tmpG$grobs[[i]]$layout$clip <- "off"
  }
  
  return(tmpG)
}

plotFragByAnnot <- function(
  fragsAnnot,
  fn,
  devices = c("png", "rds"),
  gwidth = 8,
  gheight = 1) {
  
  g <- ggplot(fragsAnnot, aes(x = xVar, y = annot, fill = individual)) +
    geom_tile() +
    scale_y_discrete(limits = rev, drop = FALSE) +
    theme_bw() +
    labs(
      x = "Cell",
      y = "Annotation"
    ) +
    # coord_equal() +
    facet_grid(~ individual, scales = "free_x", space = "free") +
    theme(
      strip.background = element_rect(fill = "transparent", color = NA),
      legend.position = "none",
      axis.title = element_text(family = "Arial", size = BASEPTFONTSIZE),
      text = element_text(family = "Arial", size = BASEPTFONTSIZE),
      legend.key.size = unit(5, "pt"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  savePlot(g, fn = fn, devices = devices, gheight = gheight, gwidth = gwidth)
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
  devices = c("png", "rds"),
  separateByIndividual = FALSE
) {
  df <- data.frame(
    cbc = c(atacCBC, hivCBC, adtCBC),
    modal = c(rep("ATAC", length(atacCBC)), rep("HIV", length(hivCBC)), rep("ADT", length(adtCBC))),
    value = TRUE
  )
  
  if (separateByIndividual) {
    df <- df %>%
      separate(cbc, sep = "#", remove = FALSE, into = c("sample", "oligoBC")) %>%
      mutate(individual = str_match(sample, "(.*)(?=_.*$)")[, 1])
    
    df <- df %>%
      pivot_wider(id_cols = all_of(c("individual", "cbc")), names_from = modal, values_from = value) %>%
      group_by(individual, ATAC, HIV, ADT)

  } else {
    df <- df %>%
      pivot_wider(id_cols = cbc, names_from = modal, values_from = value) %>%
      group_by(ATAC, HIV, ADT)
  }
  
  df <- df %>%
    dplyr::summarize(n = n()) %>%
    ungroup()
  
  nudgeBuffer <- 0.05
  if (separateByIndividual) {
    df <- df %>%
      group_by(individual) %>%
      mutate(rowID = paste0(ATAC, HIV, ADT)) %>%
      mutate(rowID = factor(rowID, levels = unique(rowID))) %>%
      mutate(nudge = n + max(n) * (n_groups(.) * nudgeBuffer))
    
  } else {
    df <- df %>%
      mutate(rowID = factor(row_number())) %>%
      mutate(nudge = n + max(n) * nudgeBuffer)
  }
  
  barPlot <- df %>%
    ggplot(aes(x = rowID, y = n)) +
    geom_bar(stat = "identity", width = 0.5, fill = "#000000") +
    geom_text(aes(label = n, y = nudge), size = BASEFONTSIZE) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(clip = "off") +
    labs(y = "Count") +
    theme_classic() +
    theme(
      axis.text = element_text(size = BASEPTFONTSIZE),
      axis.title.y = element_text(size = BASEPTFONTSIZE),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = BASEPTFONTSIZE))
  
  if (separateByIndividual) {
    barPlot <- barPlot + facet_wrap(~ individual, scales = "free", strip.position = "right", ncol = 1)
  }
  
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

plotGgRoc <- function(
  perf,
  auc,
  fn,
  returnDf = FALSE,
  devices = c("rds", "png"),
  gheight = 4,
  gwidth = 4
) {
  if (typeof(perf) == "list") {
    dfs <- lapply(perf, function(pe) {
      return(data.frame(
        x = pe@x.values[[1]],
        y = pe@y.values[[1]],
        alpha = pe@alpha.values[[1]]
      ))
    })
    
    names(dfs) <- names(perf)
    df <- bind_rows(dfs, .id = "colorVar")
  } else {
    df <- data.frame(
      x = perf@x.values[[1]],
      y = perf@y.values[[1]],
      alpha = perf@alpha.values[[1]],
      colorVar = ""
    )
  }
  
  
  df$colorVar <- paste0(df$colorVar, " (AUC=", round(unlist(auc[df$colorVar]), digits = 2), ")")
  
  if (returnDf) {
    return(df)
  }
  
  g <- ggplot(df, aes(x = x, y = y, color = colorVar)) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linetype = "dashed", color = "#cccccc") +
    geom_path() +
    labs(x = "False positive rate", y = "True positive rate") +
    theme_classic() +
    scale_x_continuous(limits = c(0,1), expand = c(0, 0.01)) +
    scale_y_continuous(limits = c(0,1), expand = c(0, 0.01)) +
    guides(color = guide_legend(nrow = 3)) +
    coord_fixed() +
    theme(
      text = element_text(family = "Arial", size = BASEPTFONTSIZE, color = "#000000"),
      axis.text = element_text(family = "Arial", size = BASEPTFONTSIZE, color = "#000000"),
      legend.title = element_blank(),
      legend.key.size = unit(0.5, "lines"),
      legend.margin = margin(-5,0,0,0),
      legend.position = "bottom"
    )
  
  savePlot(g, fn = fn, devices = devices, gheight = gheight, gwidth = gwidth)
}

plotLogitRegressionVolcano <- function(
  logModel,
  tsaCatalog,
  adtToUse,
  fn,
  devices = c("rds", "png"),
  gheight = 3.5,
  gwidth = 3.5
) {
  
  df <- data.frame(cleanName = tsaCatalog[tsaCatalog$DNA_ID %in% adtToUse, "cleanName"],
    coeff = logModel$coefficients[-1],
    p = coef(summary(logModel))[-1, 4]) %>% 
    mutate(sig = case_when(
      p < 0.05 & coeff > 0 ~ "Up in HIV+",
      p < 0.05 & coeff < 0 ~ "Up in HIV-",
      TRUE ~ "Not significant"
    )) %>%
    mutate(direction = coeff > 0)

  
  yMin <- -log10(0.05) * 1.3
  yMax <- -log10(min(df$p)) * 1.5
  
  dfLbl <- df %>%
    filter(p < 0.05) %>%
    mutate(x = ifelse(coeff > 0, 1.25, -1.25)) %>%
    group_by(direction) %>%
    dplyr::arrange(p, .by_group = TRUE) %>%
    mutate(gid = seq_along(cleanName)) %>%
    mutate(y = (yMax - yMin) / (max(gid) - 1) * (max(gid) - gid) + yMin)
  
  p <- ggplot(df, aes(x = coeff, y = -log10(p), color = sig)) +
    geom_point() +
    theme_classic() +
    scale_color_manual(values = c("#dddddd", HIVNEGCOLOR, HIVPOSCOLOR)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#cccccc") +
    geom_text(data = dfLbl, aes(x = x, y = y, label = cleanName, hjust = ifelse(coeff > 0, 0, 1)), size = 1.5) +
    geom_segment(data = dfLbl, aes(x = x * 0.98, xend = coeff, y = y, yend = -log10(p)), size = 0.25, alpha = 0.5) +
    xlim(-2.25, 2.25) +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "none",
      axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000", family = "Arial"),
      text = element_text(size = BASEPTFONTSIZE, color = "#000000", family = "Arial")) +
    labs(x = "log(odds ratio)",
      y = "log(p)")
  
  savePlot(p, fn = fn, devices = devices, gheight = gheight, gwidth = gwidth)
}

makeFancyUpsetPlotHIV <- function(
  seu,
  cellsOfInterest,
  featuresOfInterest,
  featuresOfInterestNames = featuresOfInterest,
  metadata,
  thresholds,
  fn
) {
  df <- t(seu@assays$tsa@data[featuresOfInterest, ]) %>%
    as.data.frame(.)
  
  colnames(df) <- featuresOfInterestNames
  df$meta <- "All cells"
  df[cellsOfInterest, "meta"] <- metadata
  
  dfInterest <- df[cellsOfInterest, ]
  
  for (i in seq(1:length(featuresOfInterest))) {
    dfInterest[, i] <- dfInterest[, i] > thresholds[i]
  }
  dfInterest <- dfInterest %>%
    unite(col = "celltype", all_of(featuresOfInterestNames), sep = "___") %>%
    mutate(
      celltype = factor(celltype),
      meta = factor(meta)) %>%
    group_by(meta) %>%
    mutate(totalGroup = n()) %>%
    group_by(meta, celltype, .drop = FALSE) %>%
    dplyr::summarize(totalN = n()) %>%
    group_by(meta) %>%
    mutate(totalMeta = sum(totalN)) %>%
    mutate(prop = totalN / totalMeta) %>%
    mutate(comboID = as.character(row_number()))
  
  lolliPlot <- dfInterest %>%
    ggplot(., aes(x = comboID, color = meta, y = prop)) +
    geom_linerange(aes(ymin = 0, ymax = prop),
      linetype = "dashed", position = position_dodge(0.7)) +
    geom_hline(yintercept = 0) +
    geom_point(position = position_dodge(0.7), size = 1.25) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    scale_color_manual(values = c(HIVNEGCOLOR, HIVPOSCOLOR)) +
    labs(y = "Proportion of\nHIV- or HIV+ cells") +
    theme(
      legend.position = "none",
      axis.text = element_text(family = "Arial", size = BASEPTFONTSIZE),
      axis.title.y = element_text(family = "Arial", size = BASEPTFONTSIZE),
      legend.title = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank())
  
  dotPlot <- dfInterest %>%
    separate(celltype, into = featuresOfInterestNames, sep = "___") %>%
    pivot_longer(cols = all_of(featuresOfInterestNames), names_to = "metric", values_to = "bin") %>%
    mutate(metric = factor(metric, levels = featuresOfInterestNames)) %>%
    ggplot(aes(x = comboID, y = metric)) +
    geom_point(size = 1.5, aes(color = bin)) +
    scale_color_manual(values = c("#EFEFEF", "#000000")) +
    labs(y = "Marker") +
    theme(
      panel.background = element_rect(fill = "#FFFFFF"),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text = element_text(family = "Arial", size = BASEPTFONTSIZE),
      axis.title.y = element_text(family = "Arial", size = BASEPTFONTSIZE),
      panel.grid.major.y = element_blank(),
      axis.text.x = element_blank(),
      axis.line.y = element_line(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank())
  
  ridgePlot <- df %>%
    pivot_longer(cols = all_of(featuresOfInterestNames), names_to = "metric", values_to = "bin") %>%
    mutate(metric = factor(metric, levels = featuresOfInterestNames)) %>%
    mutate(meta = factor(meta, levels = c("All cells", "HIV-", "HIV+"))) %>% #need to refactor this...
    ggplot(., aes(y = metric, x = bin, fill = meta)) +
    ggridges::geom_density_ridges(alpha = 0.7)
  
  for (i in seq(1:length(thresholds))) {
    ridgePlot <- ridgePlot +
      annotate("segment", x = thresholds[i], xend = thresholds[i], y = i, yend = i + 1, linetype = "dotted")
  }
  
  ridgePlot <- ridgePlot +
    scale_fill_manual(values = c("#1111ff", HIVNEGCOLOR, HIVPOSCOLOR)) +
    labs(x = "Expression") +
    theme_classic() +
    guides(fill = guide_legend(override.aes = list(size = 0.5))) +
    theme(axis.title = element_blank(),
      legend.title = element_blank(),
      legend.key.size = unit(0.8, "line"),
      axis.line.y = element_blank(),
      axis.text = element_text(family = "Arial", size = BASEPTFONTSIZE),
      legend.text = element_text(family = "Arial", size = BASEPTFONTSIZE),
      axis.title.x = element_text(family = "Arial", size = BASEPTFONTSIZE),
      axis.text.x = element_text(family = "Arial", size = BASEPTFONTSIZE),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank())
  
  
  finalPlot <- (((lolliPlot) / dotPlot + plot_layout(heights = c(5, 3))) |
      (as_ggplot(get_legend(ridgePlot)) / (ridgePlot + theme(legend.position = "none")) + plot_layout(heights = c(5, 3)))) + plot_layout(widths = c(4, 1))
  
  savePlot(finalPlot, fn = fn, gheight = 3, gwidth = 6, devices = c("rds", "png"))
}
