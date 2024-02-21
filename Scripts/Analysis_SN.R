library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(openxlsx)
library(tidyverse)
library(here)
library(naturalsort)
library(cowplot)
library(ggplot2)
library(rafalib)
library(DoubletFinder)
library(kableExtra)
library(ggplot2)

#instseurat.packages("Matrix", type = "source")
#install.packages("irlba", type = "source")

analysisfolder <- here::here("RNA")
setwd(analysisfolder)

data.dir <- "/Users/eliascrapa/Library/CloudStorage/Box-Box/Rat SAN different region snRNA seq/SN/Input"

# Load the library
library(Seurat)

# Define the file paths and Read data
data1 <- Read10X_h5("/Users/eliascrapa/Library/CloudStorage/Box-Box/Rat SAN different region snRNA seq/SN/Input/RATSANcluster3/filtered_feature_bc_matrix.h5")
data2 <- Read10X_h5("/Users/eliascrapa/Library/CloudStorage/Box-Box/Rat SAN different region snRNA seq/SN/Input/RATSANcluster4/filtered_feature_bc_matrix.h5")
# Create Seurat objects
seurat_object1 <- CreateSeuratObject(counts = data1, project = "RATSANcluster3")
seurat_object2 <- CreateSeuratObject(counts = data2, project = "RATSANcluster4")



seurat_object1 <- CreateSeuratObject(counts = data1, project = "RATSANcluster3")
seurat_object2 <- CreateSeuratObject(counts = data2, project = "RATSANcluster4")






seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat, features = rownames(seurat))
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 0.5)
seurat <- RunUMAP(seurat, dims = 1:10)
ElbowPlot(seurat)

#######DOUBLETREMOVAL    
sweep.res <- paramSweep(seurat, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats) # should return NULL??
pK.set <- unique(sweep.stats$pK)[2]
nExp_poi <- round(0.08*nrow(seurat@meta.data))
seurat <- doubletFinder(seurat, PCs = 1:20, pN = 0.25, pK = as.numeric(as.character(pK.set)), nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seurat <- subset(seurat,  DF.classifications_0.25_0.02_142 == "Singlet")

# Merge the two Seurat objects into one
seurat <- merge(seurat_object1, y = seurat_object2)
rm(data1,data2,seurat_object1,seurat_object2)



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RidgePlot(seurat, features=c("nFeature_RNA","nCount_RNA", "percent.mt"), log=T, ncol = 2)
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2

seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    seurat <- NormalizeData(seurat)
    seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2

#Scaling Data for DGE later on. 
all.genes <- rownames(seurat)
     seurat <- ScaleData(seurat, features = all.genes)


      #seurat <- SCTransform(seurat)

#Here we do Principal Component Analysis - You can see how your cells cluster as an overview - not really interesting for us right now.
      seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
# Examine and visualize PCA results a few different ways
print(seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
DimPlot(seurat, reduction = "pca") + NoLegend()
DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(seurat)

#Now the interesting part starts - we will do our own clustering first
    seurat <- FindNeighbors(seurat, dims = 1:10)
    seurat <- FindClusters(seurat, resolution = 0.5)
    seurat <- RunUMAP(seurat, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
    
    
    
    
librarDimPlot(seurat, reduction = "umap")

dir.create(paste0(here::here("RNA")))
dir.create(paste0(here::here("RNA","Rdata")))
saveRDS(seurat, file = paste0(here::here("RNA","Rdata"),"ouranalysis.rds"))





#ANALYSIS CLUSTERMARKER - THIS IS WHERE YOU_ FIND YOUR Cluster - We are using their Clustering
Idents(seurat) <- seurat@meta.data[["seurat_clusters"]]
my_levels <- unique(naturalsort(seurat@meta.data[["seurat_clusters"]]))
Idents(seurat) <- factor(Idents(seurat), levels= my_levels)
DimPlot(seurat, reduction = "umap") #+ NoLegend()

whichcluster <- "seurat_clusters"
whichreduction <- "umap"

dir.create(here::here("RNA","Plots"))
setwd(here::here("RNA","Plots"))


# plot this clustering
png("ClusterPlots.png", width = 1920, height = 960)
plot_grid(ncol = 3, DimPlot(seurat, label = T, reduction= whichreduction ) + NoAxes(), DimPlot(seurat, group.by = "orig.ident",  reduction= whichreduction) + 
            NoAxes(), DimPlot(seurat, group.by = whichcluster,  reduction= whichreduction) + NoAxes())
dev.off()


DefaultAssay(seurat) <- "RNA"

# Compute differentiseurat expression
markers_genes <- FindAllMarkers(seurat, logfc.threshold = 0.2, test.use = "negbinom", 
                                min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 1000, 
                                assay = "RNA", slot = "counts")

marker_genes_bp <- markers_genes
markers_genes <- markers_genes[!grepl("^mt-",markers_genes$gene),]

#We can now select the top 25 up regulated genes for plotting.

top25 <- markers_genes %>% group_by(cluster) %>% top_n(-25, p_val_adj)
factor(top25, levels= my_levels)
top25

write.xlsx(markers_genes,"ClusterMarkerGenes_negbinom.xlsx")
write.csv(top25,"ClusterMarkerGenes_negbinom_Top25.csv")
#We can now select the top 25 up regulated genes for plotting.

mypar(2, 5, mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)) {
  BarPlot <- paste0("BarPLot_ClusterGenes_C_",i,".png")
  
  png(BarPlot)
  barplot(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F), horiz = T,
          las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i",
          col = rainbow(length(top25$gene)),
          #col = colorRampPalette(c("red", "blue"))(length(top25$gene)),
          xlab = "Average Log2 Fold Change",
          ylab = "Gene")
  abline(v = c(0, 0.25), lty = c(1, 2))
  dev.off()
}


#barplot(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F), horiz = T, 
#        las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")


# We can visualize them as a heatmap. Here we are selecting the top 5.

top5 <- markers_genes %>% group_by(cluster) %>% slice_min(p_val_adj, n = 5, with_ties = FALSE) 

levels(seurat) <- my_levels
# create a scale.data slot for the selected genes
seurat <- ScaleData(seurat, features = as.character(unique(top5$gene)), assay = "RNA")
png("HearMapClusterGenes_Top3.png", width = 2400, height = 1440)

gg <- DoHeatmap(seurat, features = as.character(unique(top5$gene)), 
                assay = "RNA", group.by = naturalsort(`whichcluster`)) + scale_fill_gradientn(colors = c("blue", "black", "yellow")) 
gg
dev.off()

ggsave(plot=gg, filename ="HeatMapClusterGenes_Top5_300.pdf", device="pdf", dpi="print", width = 30, height = 25, unit="cm")

#Another way is by representing the overseurat group expression and detection rates in a dot-plot.

png("DotPlotClusters.png",width = 960, height = 1440)
DotPlot(seurat, features = rev(as.character(unique(top5$gene))), whichcluster, 
        assay = "RNA") + coord_flip() 
dev.off()



#We can also plot a violin plot for each gene.

# take top 3 genes per cluster/
top3 <- top5 %>% group_by(cluster) %>%  slice_min(p_val_adj, n = 3, with_ties = FALSE)

png("DotPlotClusters_Top3.png",width = 960, height = 1440)
DotPlot(seurat, features = rev(as.character(unique(top3$gene))), whichcluster, 
        assay = "RNA",  scale = FALSE, col.max = 1) + coord_flip() 
dev.off()

# take top 3 genes per cluster/
top1 <- top5 %>% group_by(cluster) %>%  slice_min(p_val_adj, n = 1, with_ties = FALSE)

png("DotPlotClusters_Top1.png",width = 960, height = 1440)
DotPlot(seurat, features = rev(as.character(unique(top1$gene))), whichcluster, 
        assay = "RNA",  scale = FALSE, col.max = 1) + coord_flip() 
dev.off()


# set pt.size to zero if you do not want seurat the points to hide the violin
# shapes, or to a smseurat value like 0.1

png("VlnPlots.png", width = 2880, height = 2160)
VlnPlot(seurat, features = as.character(unique(top3$gene)), ncol = 5, 
        assay = "RNA", pt.size = 0)#, group.by = whichcluster
dev.off()

