# working directory
working_dir <- "../"

# path to Seurat object
Seurat_object <- paste0(working_dir, "data/src.0_32SS.use.rds")

# number of downsamples for each cell type
sample_count_per_group <- 100

# cutoff for sparsifying the GCN
cut_off_GCN <- 0.18

# name of Seurat metatable column of cell type
group_name <- "cluster"

# whether to use super groups (separately filter genes for each super group)
# for in vivo use True, because of the high batch effect; for in vitro use False
use_super_group <- TRUE

# super groups
super_group <- list(
  g1 = c("DE", "AL", "MG"),
  g2 = c("DPE", "VPE", "PPe", "PPl"),
  g3 = c("Tip", "Acinar", "Trunk"),
  g4 = c("Trunk", "Duct", "EP1", "EP2", "EP3"),
  g5 = c("EP3", "EP4", "alphaPP2"),
  g6 = c("EP4", "beta", "delta", "alphaPP1"),
  g7 = c("alpha", "alphaPP1", "gamma"),
  g8 = c("alpha", "alphaPP2", "gamma"),
)

# MyFilterGenebyCor parameters
exp.cutoff <- 0
exp.prop.whole.max <- 0.8
exp.prop.whole.min <- 0.1
cor.cutoff <- 0.25
partner.cutoff <- 5

# whether to sample multiple times and use consensus GCN
# for in vivo use True; for in vitro use False, we want higher sensitivity
use_concensus <- TRUE

# number of sampling times
concensus_sample_time <- 3

# louvain resolution
louvain_resolution <- 1

# gene use (We use TF)
gene_use <- read.delim("data/TFlist_human.tsv", header = FALSE)$V1

# exclude gene list due to batch effects, applied in human in vivo
exclude_gene <- c("ZNF124", "TCF3", "SMARCC1", "KDM5B", "SOX4", "HES6", "POU2F1", "MYCL", "ZFHX3", "TEAD1", "PBX1", "DMTF1", "ZNF91", "KDM5A", "ZNF652", "MLXIP",
                  "ARID1B", "FOXP1", "ZBTB18", "ID4", "RCOR2", "JARID2", "FOXN3", "ARID1A", "NFE2L1", "PBRM1", "UBTF", "ZNF426", "ZNF532", "TCF4", "SMAD4", "ZNF292",
                  "MGA", "CERS6", "SSRP1", "ZNF668", "PIAS3", "ZNF775", "ZNF384", "ZNF358", "ZNF579", "ZNF22", "ZNF205", "SALL2", "GATAD1", "HBP1", "FOXA3", "NR2F6",
                  "FOXJ3", "ZNF337", "KLF11", "ZNF513", "ESRRA", "ZNF32", "USF2", "SIM1", "BAZ2A", "YY1", "MTA1", "ZBED1", "YBX1", "CERS2", "MAZ", "FOXA2",
                  "MXD4", "JUND", "NFATC3", "E2F5", "NR4A2", "MAFF", "EGR2", "EGR3", "EGR4", "NR4A1", "ATF3", "SRF", "JUN", "JUNB", "FOS", "EGR1",
                  "FOSB", "NR0B2", "NCOR1", "SMARCE1", "MEIS1", "ZNF24", "ZNF252P", "FOXK1", "ZZZ3", "WHSC1")