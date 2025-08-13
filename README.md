# GCN Construction & Comparison

This script aims to construct gene coexpressed network (GCN) to identify coexpressed modules during in vivo and in vitro islet differentiation. Each coexpressed module corresponds to a unique regulatory programme driving one particular stage. Thereby, it enables to determine: 1) how many stable stages a developmental trajectory comprises, and 2) the extent of regulatory changes between two stable stages (quantified by connectivity between two modules). 

## Usage
1. Modify the parameters in `configs/configs.R`
2. Place Seurat object (v4) in the `data/` folder
3. Run `scripts/analysis.R`
4. Check results in the `outputs/` folder

## Dependencies

The project requires the following R packages:
- Seurat
- ggplot2
- igraph
- propr
- WGCNA