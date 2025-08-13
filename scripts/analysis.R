source("configs/configs.R")
source("src/utils.R")
seurat <- readRDS(Seurat_object)
gene_use <- intersect(rownames(seurat), gene_use)

# Create a list to store the downsampled Seurat objects
seurat_sub <- list()
for (i in 1:concensus_sample_time) {
    set.seed(i)
    sample_cell <- c()
    for (cell_type in seurat@meta.data[, group_name]) {
        sample_cell <- c(sample_cell, downsample(colnames(seurat)[seurat@meta.data[, group_name] == cell_type], sample_count_per_group))
    }
    seurat_sub[[i]] <- seurat[, sample_cell]
}
saveRDS(seurat_sub, file = "output/seurat_sub.rds")

# Filter TF list
TF_use <- list()
if (use_super_group) {
    for (i in 1:concensus_sample_time) {
        TF_use_sub <- list()
        for (j in seq_along(super_group)) {
            TF_use_sub[[j]] <- MyFilterGenebyCor(as.matrix(seurat_sub[[i]]@assays$RNA@data[, seurat_sub[[i]]@meta.data[, group_name] %in% super_group[[j]]]),
                var.gene = gene_use,
                exp.cutoff = exp.cutoff,
                exp.prop.whole.max = exp.prop.whole.max,
                exp.prop.whole.min = exp.prop.whole.min,
                cor.cutoff = cor.cutoff,
                partner.cutoff = partner.cutoff
            )
        }
        TF_use[[i]] <- unique(unlist(TF_use_sub))
    }
} else {
    for (i in 1:concensus_sample_time) {
        TF_use[[i]] <- TF_use_sub[[j]] <- MyFilterGenebyCor(as.matrix(seurat_sub[[i]]@assays$RNA@data),
            var.gene = gene_use,
            exp.cutoff = exp.cutoff,
            exp.prop.whole.max = exp.prop.whole.max,
            exp.prop.whole.min = exp.prop.whole.min,
            cor.cutoff = cor.cutoff,
            partner.cutoff = partner.cutoff
        )
    }
}

# Create GCN
genegraph <- list()
for (i in 1:concensus_sample_time) {
    cor_matrix <- WGCNA::cor(t(as.matrix(seurat_sub[[i]]@assays$RNA@data[TF_use[[i]], ])))
    cor_matrix[cor_matrix < cut_off_GCN] <- 0
    genegraph[[i]] <- graph.adjacency(cor_matrix, mode = "undirected", weighted = TRUE)
    genegraph[[i]] <- simplify(genegraph[[i]], remove.multiple = TRUE, remove.loops = TRUE)
    genegraph[[i]] <- delete.vertices(genegraph[[i]], V(genegraph[[i]])[igraph::clusters(genegraph[[i]])$membership %in% which(igraph::clusters(genegraph[[i]])$csize < 10)])
    if (length(exclude_gene) > 0) {
        genegraph[[i]] <- delete.vertices(genegraph[[i]], V(genegraph[[i]])[names(V(genegraph[[i]])) %in% exclude_gene])
    }
}

if (use_concensus) {
    genegraph <- Reduce("%s%", genegraph)
    genegraph <- delete.vertices(genegraph, V(genegraph)[igraph::clusters(genegraph)$membership %in% which(igraph::clusters(genegraph)$csize < 10)])
    genegraph <- igraph::simplify(genegraph)
}else {
   genegraph <- genegraph[[1]]
}

# plotting
layout <- layout_with_fr(genegraph)
cluster <- cluster_louvain(genegraph, resolution = louvain_resolution)
cluster_merge <- cluster$membership
# you can manually merge clusters if needed like cluster_merge[cluster_merge %in% c(1, 2, 3)] <- "c1"

pdf("output/gene_graph.pdf")
plot.igraph(genegraph,
              layout = layout,
              edge.color = "#EFEFEF88",
              edge.width = 2,
              vertex.label.cex = 1,
              vertex.frame.color = "#FFFFFF00",
              vertex.frame.width = 0.5,
              vertex.size = 6,
              vertex.color = MyName2Col(cluster_merge))
dev.off()

save(genegraph, layout, cluster_merge, file = "output/gene_graph.RData")
