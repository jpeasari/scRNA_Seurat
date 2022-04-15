# script to integrate scRNA-Seq datasets to correct for batch effects
# setwd("~/Desktop/demo/single_cell_integrate")

dir <- getwd()
dir
# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# get data location
dirs <- list.dirs(path = 'all_data_new/', recursive = F, full.names = F)
dirs



for(x in dirs) {
  name <- gsub(x,x, x)
  
  url <- paste0('D:/Dave/snRNA/all_data_new/',x)
  cts <- Read10X(data.dir=url)
  
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts)) 
  
}


names <- ls()
names <- c('AD_1','AD_10', 'AD_13', 'AD_19', 'AD_2', 'AD_20','AD_21' ,'AD_4' , 'AD_6', 'AD_8','AD_9', 'NC_11','NC_12' ,'NC_14' ,'NC_15', 'NC_16', 'NC_17','NC_18' ,'NC_3' ,
          'NC_7')

# merge datasets

merged_seurat <- merge(AD_1, y = c(AD_10, AD_13, AD_19, AD_2, AD_20,AD_21 ,AD_4 , AD_6, AD_8,AD_9, NC_11,NC_12 ,NC_14 ,NC_15, NC_16, NC_17,NC_18 ,NC_3 ,
                                   NC_7),
                       add.cell.ids = names,
                       project = 'HB')


merged_seurat

# QC & filtering -----------------------

View(merged_seurat@meta.data)
# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('sample', 'No', 'Barcode'), 
                                    sep = '_')

# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')

# explore QC


# filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 500 &
                                   mitoPercent < 20)

merged_seurat_filtered

merged_seurat




# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)

tiff('plots/elbow_plot.tiff',units="in", width=5, height=6, res=300)
ElbowPlot(merged_seurat_filtered)
dev.off()

merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)


# plot
tiff('plots/dimplot_1.tiff',units="in", width=5, height=6, res=300)
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'sample', label = TRUE, repel = TRUE,raster=FALSE)
p1
dev.off()

#p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'No',cols = c('red','green','blue'),raster=FALSE)

#grid.arrange(p1, p2, ncol = 2, nrow = 2)


# perform integration to correct for batch effects ------
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'sample')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 1000)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features,dims = 1:30,reduction = "rpca",reference = c(1, length(obj.list)),k.anchor = 5)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors, dims = 1:20)
seurat.integrated_copy <- seurat.integrated
seurat.integrated <- seurat.integrated_copy

DefaultAssay(seurat.integrated) <- "integrated"


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated,npcs = 50)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:20)
seurat.integrated <- FindNeighbors(seurat.integrated, reduction = "pca", dims = 1:20)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 1)

seurat.integrated_1 <- seurat.integrated
DefaultAssay(seurat.integrated_1) <- "RNA"


tiff('plots/dimplot_2.tiff',units="in", width=6, height=6, res=300)
p3 <- DimPlot(seurat.integrated_1, reduction = 'umap', group.by = 'sample',raster=FALSE, label = TRUE, repel = TRUE)
p3
dev.off()

tiff('plots/dimplot_3.tiff',units="in", width=6, height=6, res=300)
p4 <- DimPlot(seurat.integrated_1, reduction = 'umap', label = TRUE, repel = TRUE,raster=FALSE)
p4
dev.off()


tiff('plots/dimplot_2_3.tiff',units="in", width=10, height=6, res=300)
combined <- grid.arrange(p3,p4, ncol = 2, nrow = 1)
combined
dev.off()





tiff('plots/dimplot_4.tiff',units="in", width=8, height=6, res=300)
p5 <- DimPlot(seurat.integrated_1, reduction = 'umap', split.by = 'sample',raster=FALSE, repel = TRUE)
p5
dev.off() 





nk.markers <- FindConservedMarkers(seurat.integrated_1, ident.1 = 6, grouping.var = "sample", verbose = FALSE)
View(nk.markers)


tiff('plots/feature_plot.tiff',units="in", width=12, height=8, res=300)
feature_plot <- FeaturePlot(seurat.integrated_1, features = c("ADGRV1", "GPC5", "RYR3","CLDN5","ABCB1","EBF1","CBLN2","LDB2","CAMK2A","LHFPL3", "GAD1","LRMDA","DOCK8","PLP1","ST18","MBP"), min.cutoff = "q9",raster=FALSE)
feature_plot
dev.off()

#immune.combined <- RenameIdents(seurat.integrated, `0` = "oligodendrocytes", `3` = "Astrocytes", `7` = "Microglia", `4`= "Inhibitory neurons" , `8` = "excitatory neurons")

tiff('plots/clusters_annotated.tiff',units="in", width=6, height=6, res=300)
seurat.integrated_1 <- RenameIdents(seurat.integrated_1,`3` = "Astro",`10` = "Astro",`35` = "Astro",`26` = "Astro",`31` = "Astro",`36` = "Astro",
                                  `20` = "Endo",`33` = "Endo",
                                  `2` = "Inhibit",`32` = "Inhibit",`37` = "Inhibit",`38` = "Inhibit", `39` = "Inhibit",`24` = "Inhibit",`11` = "Inhibit",`17` = "Inhibit",`6` = "Inhibit",`15` = "Inhibit" ,`29` = "Inhibit" ,
                                  `27` = "Micro",`9` = "Micro",
                                  `0` = "Oligo",`8` = "Oligo",`14` = "Oligo" ,`30` = "Oligo",
                                  `12` = "Excit",`13` = "Excit",`24` = "Excit", `23` = "Excit",`4` = "Excit",`16` = "Excit",`1` = "Excit",`7` = "Excit",`22` = "Excit",`5` = "Excit",`34` = "Excit",`28` = "Excit",`19` = "Excit",`18` = "Excit",`21` = "Excit",`25` = "Excit")
p6 <- DimPlot(seurat.integrated_1, label = TRUE,raster=FALSE)
p6
dev.off()

seurat.integrated_3 <- seurat.integrated_1

tiff('plots/conserved_markers.tiff',units="in", width=6, height=6, res=300)
Idents(seurat.integrated_3) <- factor(Idents(seurat.integrated_3), levels = c("Astro","Endo","Inhibit","Micro","Oligo", "Excit"))
markers.to.plot <- c("ELOVL1","ELOVL2","ELOVL4","ELOVL5","ELOVL6","ELOVL7")
p7 <- DotPlot(seurat.integrated_3, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "sample")  + RotatedAxis()
p7
dev.off()


tiff('plots/Elovl_genes.tiff',units="in", width=6, height=10, res=300)
p8 <- FeaturePlot(seurat.integrated_3, features = c("ELOVL1","ELOVL2","ELOVL4","ELOVL5","ELOVL6","ELOVL7"), split.by = "sample", max.cutoff = 3,
            cols = c("grey", "red"))
p8
dev.off()

seurat.integrated_4 <- seurat.integrated_3

#Identify differential expressed genes across conditions
tiff('plots/DEGs_clusters.tiff',units="in", width=6, height=6, res=300)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
t.cells <- subset(seurat.integrated_4, idents = "Astro")
Idents(t.cells) <- "sample"
avg.t.cells <- as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))
avg.t.cells$gene <- rownames(avg.t.cells)

cd14.mono <- subset(seurat.integrated_4, idents = "Oligo")
Idents(cd14.mono) <- "sample"
avg.cd14.mono <- as.data.frame(log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA))
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

genes.to.label = c("ELOVL1","ELOVL2","ELOVL4","ELOVL5","ELOVL6","ELOVL7")
p1 <- ggplot(avg.t.cells, aes(NC, AD)) + geom_point() + ggtitle("Astro")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE,xnudge  = 0,ynudge =0)
p2 <- ggplot(avg.cd14.mono, aes(NC, AD)) + geom_point() + ggtitle("Oligo")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE,xnudge  = 0,ynudge =0)
p9 <- p1 + p2
p9
dev.off()





tiff('plots/Violin_plots.tiff',units="in", width=10, height=15, res=300)
plots <- VlnPlot(seurat.integrated_4, features = c("ELOVL1","ELOVL2","ELOVL4","ELOVL5","ELOVL6","ELOVL7"), split.by = "sample",
                 pt.size = 0, combine = FALSE)

p10 <- wrap_plots(plots = plots, ncol = 1)
p10
dev.off()


##### End of the analysis 
