## Importer function (min 200 genes per cell)
```{r}
#' Import Single cell sequencing experiments into Seurat3and perform normalisation and scale Data 
#' @author David John
#' @param pathways A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids Vector of strings that are assigned to the concordant cells
#' @return Merged seurat object
Importer <- function(pathway,id, TenX=TRUE, performNormalisation=TRUE, performScaling = FALSE,performVariableGeneDetection=TRUE) {
  if (TenX) {
    Matrix <- Read10X(pathway)
  }  else{
    Matrix <- read.table(pathway,header = TRUE,sep = ",", dec = ".", row.names = 1)
  }
  seuratObject =CreateSeuratObject(counts = Matrix, project = id, min.cells = 5)
  seuratObject$sample <- id
  tmp<-unlist(strsplit(id,split = "-"))
  seuratObject$condition <- paste0(tmp[1:length(tmp)-1],collapse = "-")
  seuratObject <- subset(x = seuratObject, subset = nFeature_RNA > 200)
  if (performNormalisation==TRUE) {
    seuratObject<-NormalizeData(object = seuratObject,verbose = FALSE)
  }
  if(performVariableGeneDetection){
    seuratObject<-FindVariableFeatures(object = seuratObject, do.plot = FALSE, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  if (performScaling==TRUE) {
    seuratObject<-ScaleData(object = seuratObject)
  }
  cat("Imported ", length(seuratObject@meta.data$orig.ident), " cells from ", pathway, "with ID ", id, "\n")
  return(seuratObject)
}




## prepare data for cluster t.test from the deg list and do a cluster t-test
do_cluster_t_test <- function(seurat_subset, degs, m_group="dataset", m_cluster="seurat_clusters"){
  gene_names<- names(table(degs$gene))
  #print(head(gene_names))
  p_values <- vector("list",length(gene_names))
  names(p_values) <- gene_names
  #gene_names <- row.names(cluster_subset)
  #if (celltype=="Adipocytes"){
  #  seurat_subset <- seurat_subset[,!seurat_subset$orig.ident=="D7"]
  #}
  group <- as.numeric(seurat_subset[[m_group]][,1])
  #print(head(group))
  cluster <- seurat_subset[[m_cluster]][,1]
  #print(head(cluster))
  for (gene in gene_names){
    #print(gene)
    y <- c(t(as.matrix(seurat_subset@assays$RNA[gene,])))
    test_info <- t.test.cluster(y, cluster = cluster, group = group)
    p_values[[gene]] <- test_info[nrow(test_info)]
  }
  p_values
}

## added line 54-56 so that each group is tested if 
## only one ovbservation is present and throw an error
my.t.test.cluster <- function (y, cluster, group, conf.int = 0.95) 
{
  group <- as.factor(group)
  cluster <- as.factor(cluster)
  s <- !(is.na(y) | is.na(cluster) | is.na(group))
  y <- y[s]
  cluster <- cluster[s]
  group <- group[s]
  n <- length(y)
  if (n < 2) 
    stop("n<2")
  gr <- levels(group)
  if (length(gr) != 2) 
    stop("must have exactly two treatment groups")
  n <- table(group)
  nc <- tapply(cluster, group, function(x) length(unique(x)))
  bar <- tapply(y, group, mean)
  u <- unclass(group)
  y1 <- y[u == 1]
  y2 <- y[u == 2]
  c1 <- factor(cluster[u == 1])
  c2 <- factor(cluster[u == 2])
  b1 <- tapply(y1, c1, mean)
  b2 <- tapply(y2, c2, mean)
  m1 <- table(c1)
  m2 <- table(c2)
  if (any(names(m1) != names(b1)))
    stop("logic error 1")
  if (any(names(m2) != names(b2)))
    stop("logic error 2")
  if (any(m2 < 2))
    stop(paste("The following clusters contain only one observation:",
               paste(names(m2[m2 < 2]), collapse = " ")))
  if (any(m1 < 2))
    stop(paste("The following clusters contain only one observation:",
               paste(names(m1[m1 < 2]), collapse = " ")))
  M1 <- mean(y1)
  M2 <- mean(y2)
  ssc1 <- sum(m1 * ((b1 - M1)^2))
  ssc2 <- sum(m2 * ((b2 - M2)^2))
  if (nc[1] != length(m1))
    stop("logic error 3")
  if (nc[2] != length(m2))
    stop("logic error 4")
  df.msc <- sum(nc) - 2
  msc <- (ssc1 + ssc2)/df.msc
  v1 <- tapply(y1, c1, var)
  v2 <- tapply(y2, c2, var)
  ssw1 <- sum((m1 - 1) * v1)
  ssw2 <- sum((m2 - 1) * v2)
  df.mse <- sum(n) - sum(nc)
  mse <- (ssw1 + ssw2)/df.mse
  na <- (sum(n) - (sum(m1^2)/n[1] + sum(m2^2)/n[2]))/(sum(nc) - 
                                                        1)
  rho <- (msc - mse)/(msc + (na - 1) * mse)
  r <- max(rho, 0)
  C1 <- sum(m1 * (1 + (m1 - 1) * r))/n[1]
  C2 <- sum(m2 * (1 + (m2 - 1) * r))/n[2]
  v <- mse * (C1/n[1] + C2/n[2])
  v.unadj <- mse * (1/n[1] + 1/n[2])
  de <- v/v.unadj
  dif <- diff(bar)
  se <- sqrt(v)
  zcrit <- qnorm((1 + conf.int)/2)
  cl <- c(dif - zcrit * se, dif + zcrit * se)
  z <- dif/se
  P <- 2 * pnorm(-abs(z))
  stats <- matrix(NA, nrow = 20, ncol = 2, dimnames = list(c("N", 
                                                             "Clusters", "Mean", "SS among clusters within groups", 
                                                             "SS within clusters within groups", "MS among clusters within groups", 
                                                             "d.f.", "MS within clusters within groups", "d.f.", "Na", 
                                                             "Intracluster correlation", "Variance Correction Factor", 
                                                             "Variance of effect", "Variance without cluster adjustment", 
                                                             "Design Effect", "Effect (Difference in Means)", "S.E. of Effect", 
                                                             paste(format(conf.int), "Confidence limits"), "Z Statistic", 
                                                             "2-sided P Value"), gr))
  stats[1, ] <- n
  stats[2, ] <- nc
  stats[3, ] <- bar
  stats[4, ] <- c(ssc1, ssc2)
  stats[5, ] <- c(ssw1, ssw2)
  stats[6, 1] <- msc
  stats[7, 1] <- df.msc
  stats[8, 1] <- mse
  stats[9, 1] <- df.mse
  stats[10, 1] <- na
  stats[11, 1] <- rho
  stats[12, ] <- c(C1, C2)
  stats[13, 1] <- v
  stats[14, 1] <- v.unadj
  stats[15, 1] <- de
  stats[16, 1] <- dif
  stats[17, 1] <- se
  stats[18, ] <- cl
  stats[19, 1] <- z
  stats[20, 1] <- P
  attr(stats, "class") <- "t.test.cluster"
  stats
}

library(cowplot)
PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

DotPlot_costum <- function (object, assay = NULL, features, cols = c("lightgrey", 
                                                                     "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                            idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
                            scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA, cluster_order = NULL) 
{
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                                                                                          paste(sort(x = levels(x = object), decreasing = TRUE), 
                                                                                                collapse = "|"), ")_)"), replacement = "", 
                         USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  if (!is.null(cluster_order)){
    for (i in names(cluster_order)){
      data.plot[as.numeric(as.character(data.plot$id)) %in% cluster_order[[i]],"celltype"] <- i
    }
    data.plot$celltype <- factor(data.plot$celltype, levels = names(cluster_order))
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(y = "features.plot", x = "id")) + 
    geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) + 
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
    theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + 
    guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Cluster", y = ifelse(test = is.null(x = split.by), 
                                   yes = "Marker genes", no = "Split Identity")) + theme_cowplot()
  if ( (!is.null(x = feature.groups)) & (!is.null(cluster_order)) ) {
    plot <- plot + facet_grid(vars(feature.groups), vars(celltype), scales = "free", 
                              space = "free") + theme(panel.spacing = unit(x = 1, 
                                                                           units = "lines"), strip.background = element_blank())
  }
  if ( (!is.null(x = feature.groups)) & (is.null(cluster_order)) ) {
    plot <- plot + facet_grid(rows = vars(feature.groups), scales = "free", 
                              space = "free") + theme(panel.spacing = unit(x = 1, 
                                                                           units = "lines"), strip.background = element_blank())
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient2(low = "blue", mid = "grey", high = "red")
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}
