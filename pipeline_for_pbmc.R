#This code uses Seurat Version 2.4#


library(Matrix)
C1MF.data <- Read10X(data.dir = "/media/minion/Seagate Expansion Drive/SHAOBO scRNA/mutiple anaylsis/C1multi/outs/filtered_gene_bc_matrices/hg19")
C1MF <- CreateSeuratObject(raw.data = C1MF.data, min.cells = 3, min.genes = 200, project = "C1MultiFiltered")
mito.genes <- grep(pattern = "^MT-", x = rownames(x = C1MF@data), value = TRUE)
> percent.mito <- Matrix::colSums(C1MF@raw.data[mito.genes, ])/Matrix::colSums(C1MF@raw.data)
> C1MF <- AddMetaData(object = C1MF, metadata = percent.mito, col.name = "percent.mito")
> VlnPlot(object = C1MF, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
> par(mfrow = c(1, 2))
> GenePlot(object = C1MF, gene1 = "nUMI", gene2 = "percent.mito")
> GenePlot(object = C1MF, gene1 = "nUMI", gene2 = "nGene")
> C1MF <- FilterCells(object = C1MF, subset.names = c("nGene", "percent.mito"), 
                      +                     low.thresholds = c(200, -Inf), high.thresholds = c(3000, 0.15))
> C1MF <- NormalizeData(object = C1MF, normalization.method = "LogNormalize", 
                        +                       scale.factor = 10000)
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
  |----|----|----|----|----|----|----|----|----|----|
  **************************************************|
  > C1MF <- FindVariableGenes(object = C1MF, mean.function = ExpMean, dispersion.function = LogVMR, 
                              +                           x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
Calculating gene means
0%   10   20   30   40   50   60   70   80   90   100%
  |----|----|----|----|----|----|----|----|----|----|
  **************************************************|
  Calculating gene variance to mean ratios
0%   10   20   30   40   50   60   70   80   90   100%
  |----|----|----|----|----|----|----|----|----|----|
  **************************************************|
  > length(x = C1MF@var.genes)
[1] 1525

1  0.000000e+00 2.7623298 0.983 0.426  0.000000e+00       0        S100A9
2  0.000000e+00 2.6797697 0.981 0.316  0.000000e+00       0        S100A8
3 1.120898e-145 0.7701065 0.480 0.129 1.906760e-141       1          CCR7
4 5.675599e-142 0.6490500 0.828 0.559 9.654762e-138       1          LDHB
5 7.065132e-259 1.5147326 0.525 0.040 1.201850e-254       2          GZMK
6 2.444925e-174 1.3619179 0.703 0.210 4.159062e-170       2         KLRB1
7 1.191736e-198 1.0031776 0.693 0.151 2.027262e-194       3          CD8B
8 1.042083e-184 1.0070699 0.472 0.050 1.772687e-180       3 RP11-291B21.2
9  0.000000e+00 2.0378188 0.960 0.158  0.000000e+00       4          GZMH
10  0.000000e+00 1.8819221 0.993 0.311  0.000000e+00       4          CCL5

C1MF <- ScaleData(object = C1MF, vars.to.regress = c("nUMI", "percent.mito"))
[1] "Regressing out nUMI"         "Regressing out percent.mito"
|===========================================================| 100%
[1] "Scaling data matrix"
|===========================================================| 100%
> C1MF <- RunPCA(object = C1MF, pc.genes = C1MF@var.genes, do.print = TRUE, pcs.print = 1:5, 
                 +                genes.print = 5)
[1] "PC1"
[1] "IL32" "NKG7" "CCL5" "CST7" "GZMA"
[1] ""
[1] "CSTA"          "RP11-1143G9.4" "FCN1"          "CST3"         
[5] "LST1"         
[1] ""
[1] ""
[1] "PC2"
[1] "NKG7"   "TYROBP" "GNLY"   "CST7"   "GZMA"  
[1] ""
[1] "CD79A"     "MS4A1"     "LTB"       "LINC00926" "BANK1"    
[1] ""
[1] ""
[1] "PC3"
[1] "IL7R"      "CD3D"      "PRKCQ-AS1" "MAL"       "CCR7"     
[1] ""
[1] "CD74"     "HLA-DQA2" "HLA-DPB1" "HLA-DPA1" "HLA-DRB1"
[1] ""
[1] ""
[1] "PC4"
[1] "CD79B"     "CD79A"     "MS4A1"     "LINC00926" "BANK1"    
[1] ""
[1] "LILRA4"        "LRRC26"        "SERPINF1"      "SCT"          
[5] "RP11-117D22.2"
[1] ""
[1] ""
[1] "PC5"
[1] "SDPR"    "PF4"     "S100A12" "PPBP"    "VCAN"   
[1] ""
[1] "CDKN1C" "HES4"   "CKB"    "TCF7L2" "LYPD2" 
[1] ""
[1] ""
> VizPCA(object = C1MF, pcs.use = 1:2)
> PCAPlot(object = C1MF, dim.1 = 1, dim.2 = 2)
> C1MF <- ProjectPCA(object = C1MF, do.print = FALSE)
> PCHeatmap(object = C1MF, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
> C1MF <- JackStraw(object = C1MF, num.replicate = 100, do.print = FALSE)
> JackStrawPlot(object = C1MF, PCs = 1:12)
Warning message:
  Removed 12804 rows containing missing values (geom_point). 
> PCElbowPlot(object = C1MF)
> save.SNN = TRUE
> C1MF <- FindClusters(object = C1MF, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
> PrintFindClustersParams(object = C1MF)
Parameters used in latest FindClusters calculation run on: 2018-02-27 16:22:57
=============================================================================
  Resolution: 0.6
-----------------------------------------------------------------------------
  Modularity Function    Algorithm         n.start         n.iter
1                   1                 100             10
-----------------------------------------------------------------------------
  Reduction used          k.param          k.scale          prune.SNN
pca                 30                25              0.0667
-----------------------------------------------------------------------------
  Dims used in calculation
=============================================================================
  1 2 3 4 5 6 7 8 9 10

> C1MF <- RunTSNE(object = C1MF, dims.use = 1:10, do.fast = TRUE)
> TSNEPlot(object = pbmc)
Error in match(x, table, nomatch = 0L) : object 'pbmc' not found
> TSNEPlot(object = C1MF)
> C1MF.markers <- FindAllMarkers(object = C1MF, only.pos = TRUE, min.pct = 0.25, 
                                 +                                thresh.use = 0.25)
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 25s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 12s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 06s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 13s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 11s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 16s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 13s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 06s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 25s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 19s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 15s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 25s
> C1MF.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
# A tibble: 24 x 7
# Groups:   cluster [12]
p_val avg_logFC pct.1 pct.2     p_val_adj cluster
<dbl>     <dbl> <dbl> <dbl>         <dbl>  <fctr>
  1  0.000000e+00 2.7623298 0.983 0.426  0.000000e+00       0
2  0.000000e+00 2.6797697 0.981 0.316  0.000000e+00       0
3 1.120898e-145 0.7701065 0.480 0.129 1.906760e-141       1
4 5.675599e-142 0.6490500 0.828 0.559 9.654762e-138       1
5 7.065132e-259 1.5147326 0.525 0.040 1.201850e-254       2
6 2.444925e-174 1.3619179 0.703 0.210 4.159062e-170       2
7 1.191736e-198 1.0031776 0.693 0.151 2.027262e-194       3
8 1.042083e-184 1.0070699 0.472 0.050 1.772687e-180       3
9  0.000000e+00 2.0378188 0.960 0.158  0.000000e+00       4
10  0.000000e+00 1.8819221 0.993 0.311  0.000000e+00       4
# ... with 14 more rows, and 1 more variables: gene <chr>
> FeaturePlot(object = C1MF, features.plot = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols.use = c("grey", "blue"), 
              reduction.use = "tsne")
> top10 <- C1MF.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
> DoHeatmap(object = C1MF, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)




#T cell pipeline
Tcells.subset <- SubsetData(object = C1MF, ident.use = c("CD8 T cells", "CD4 T cells"))
mito.genes <- grep(pattern = "^MT-", x = rownames(x = Tcells.subset@data), value = TRUE)
percent.mito <- Matrix::colSums(Tcells.subset@raw.data[mito.genes, ])/Matrix::colSums(Tcells.subset@raw.data)
Tcells.subset <- AddMetaData(object = Tcells.subset, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = Tcells.subset, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
Tcells.subset <- AddMetaData(object = Tcells.subset, metadata = percent.mito, col.name = "percent.mito")
Tcells.subset <- FindVariableGenes(object = Tcells.subset, mean.function = ExpMean, dispersion.function = LogVMR,
                                   x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
Tcells.subset <- ScaleData(object = Tcells.subset, vars.to.regress = c("nUMI", "percent.mito"))
Tcells.subset <- RunPCA(object = Tcells.subset, pc.genes = Tcells.subset@var.genes, do.print = TRUE, pcs.print = 1:5,
                        genes.print = 5)
PrintPCA(object = Tcells.subset, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = Tcells.subset, pcs.use = 1:2)
PCAPlot(object = Tcells.subset, dim.1 = 1, dim.2 = 2)
Tcells.subset <- ProjectPCA(object = Tcells.subset, do.print = FALSE)
PCHeatmap(object = Tcells.subset, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = Tcells.subset, pc.use = 1:12, cells.use = 500, do.balanced = TRUE,
          label.columns = FALSE, use.full = FALSE)
Tcells.subset <- JackStraw(object = Tcells.subset, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = Tcells.subset, PCs = 1:12)
Tcells.subset <- FindClusters(object = Tcells.subset, reduction.type = "pca", dims.use = 1:10,
                              resolution = 0.8, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = Tcells.subset)
Tcells.subset <- RunTSNE(object = Tcells.subset, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = Tcells.subset, do.label = TRUE)
VlnPlot(object = Tcells.subset, point.size.use = 0, nCol = 3, features.plot= c("CD3D", "CD8A", "NKG7", "GZMA", "CD160", "TIGIT", "OASL", "IL7R", "CCR7", "LTB"))
Tcells.subset <- FindClusters(object = Tcells.subset, reduction.type = "pca", dims.use = 1:10, resolution = 1.0, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
TSNEPlot(object = Tcells.subset, do.label = TRUE)
Tcells.subset <- FindClusters(object = Tcells.subset, reduction.type = "pca", dims.use = 1:10, resolution = 0.8, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
TSNEPlot(object = Tcells.subset, do.label = TRUE)
table(Tcells.subset@ident)
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6)
new.cluster.ids <- c("CD4-Tn", "CD4-Tem", "CD8-Tn", "CD8-Tem","CD4-Tpm", "CD4-Tn", "CD8-Tn")