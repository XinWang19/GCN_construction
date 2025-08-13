library(Seurat)
library(igraph)
# Filter genes by correlation and expression
MyFilterGenebyCor <- function(exp.data,
                              var.gene = NULL,
                              exp.cutoff = 1,
                              exp.prop.whole.max = 0.8,
                              exp.prop.whole.min = 0,
                              vector.group = NULL,
                              exp.prop.group.min = 0,
                              cor.method = c("pearson", "rho", "phs"),
                              cor.cutoff = 0.5,
                              partner.cutoff = 5,
                              iteration = FALSE) {
    cor.method <- cor.method[1]
    if (!is.null(var.gene)) {
        exp.data <- exp.data[rownames(exp.data) %in% var.gene, ]
    }
    if (cor.method == "pearson") {
        cor.data <- WGCNA::cor(t(exp.data), method = "pearson")
    } else if (cor.method == "rho") {
        cor.data <- propr:::lr2rho(t(exp.data))
        rownames(cor.data) <- rownames(exp.data)
        colnames(cor.data) <- rownames(exp.data)
    } else if (cor.method == "phs") {
        cor.data <- propr:::lr2phs(t(exp.data))
        rownames(cor.data) <- rownames(exp.data)
        colnames(cor.data) <- rownames(exp.data)
    }

    gene.cor <- rownames(cor.data)[rowSums(cor.data > cor.cutoff) >= partner.cutoff]
    if (iteration) {
        while (nrow(cor.data) > length(gene.cor)) {
            cor.data <- cor.data[gene.cor, gene.cor]
            gene.cor <- rownames(cor.data)[rowSums(cor.data > cor.cutoff) >= partner.cutoff]
        }
    }
    exp.data <- exp.data[rowSums(exp.data > exp.cutoff) <= exp.prop.whole.max * ncol(exp.data) &
        rowSums(exp.data > exp.cutoff) >= exp.prop.whole.min * ncol(exp.data), ]
    if (!is.null(vector.group)) {
        gene.group <- c()
        for (group.tmp in unique(vector.group)) {
            exp.data.tmp <- exp.data[, vector.group %in% exp.data.tmp]
            gene.group.tmp <- rownames(exp.data.tmp)[rowSums(exp.data.tmp >= exp.cutoff) >= exp.prop.group.min * ncol(exp.data.tmp)]
            gene.group <- c(
                gene.group,
                gene.group.tmp
            )
        }
        exp.data <- exp.data[rownames(exp.data) %in% gene.group, ]
    }
    return(cor.data)
}

# downsample
downsample <- function(x, n) {
  if (length(x) <= n) {
    x
  } else {
    sample(x, n)
  }
}

# convert name to color
MyName2Col <- function(name,
                       palette = WGCNA::labels2colors(1:20),
                       is.row = FALSE,
                       na.value = NULL) {
  name.factor <- as.factor(name)
  if (!is.null(names(palette))) {
    palette <- palette[levels(name.factor)]
  } else {
    print("your palette order must adapt to you names level order")
  }
  name.color <- palette[name.factor]
  name.color <- as.matrix(name.color)
  if (is.row) {
    name.color <- t(name.color)
  }
  if (!is.null(na.value)) {
    name.color[is.na(name.color)] <- na.value
  }
  return(name.color)
}