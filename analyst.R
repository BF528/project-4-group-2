library(dplyr)
library(Seurat)
cells <- readRDS("GSM2230758_seurat.rda")
#cells <- read.csv("combined_count_mtx.csv",header = T)
#1
cells.markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25)
#2
VlnPlot(cells, features = c("GCG")) #alpha 1 2 9
VlnPlot(cells, features = c("INS")) #beta 3 6
VlnPlot(cells, features = c("SST")) #delta 0
VlnPlot(cells, features = c("PPY")) #gamma 8
VlnPlot(cells, features = c("GHRL")) #Epsilon 0
VlnPlot(cells, features = c("KRT19")) #Ductal 4
VlnPlot(cells, features = c("CPA1")) #Acinar 5
VlnPlot(cells, features = c("PDGFRA", "PDGFRB")) #Activated stellate 7
VlnPlot(cells, features = c("RGS5","PDGFRB")) #Quiescent stellate 7
VlnPlot(cells, features = c("VWF", "PECAM1","CD34")) #Endothelial 12
VlnPlot(cells, features = c("SDS", "CD163","CD68","IgG")) #Macrophage 11
VlnPlot(cells, features = c("CD3", "CD8","TRAC")) #Cytotoxic T 
VlnPlot(cells, features = c("TPSAB1", "KIT","CPA3")) #Mast 11
VlnPlot(cells, features = c("SOX10")) #Schwann too few cells 7
##                  
VlnPlot(cells, features = c("KRT19", "CPA1", "VWF"))
FeaturePlot(cells, features = c("KRT19", "CPA1", "VWF"))
#unknown 9:Oligodendrocyte./ Cholangiocytes
#10:Basal cells
top5 <- cells.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
cluster0.marker <- top5[top5$cluster==0,]
cluster1.marker <- top5[top5$cluster==1,]
cluster2.marker <- top5[top5$cluster==2,]
cluster3.marker <- top5[top5$cluster==3,]
cluster4.marker <- top5[top5$cluster==4,]
cluster5.marker <- top5[top5$cluster==5,]
cluster6.marker <- top5[top5$cluster==6,]
cluster7.marker <- top5[top5$cluster==7,]
cluster8.marker <- top5[top5$cluster==8,]
cluster9.marker <- top5[top5$cluster==9,]
cluster10.marker <- top5[top5$cluster==10,]
cluster11.marker <- top5[top5$cluster==11,]
cluster12.marker <- top5[top5$cluster==12,]
new.cluster.ids <- c("Delta/Epsilon", "Alpha", "Alpha", "Beta", "Ductal","Acinar","Beta", 
                     "Activated/Quiescent stellate/Schwann", "Gamma","Alpha",
                     "Possible Basal cells","Mast/Macrophage", "Endothelial")
names(new.cluster.ids) <- levels(cells)
cells <- RenameIdents(cells, new.cluster.ids)
Umap <- DimPlot(cells, reduction = "umap", label = TRUE, pt.size = 1,
                label.size = 6) + NoLegend()
Umap
#ggsave(Umap,"Umap_clustered_cells.png")
#top5 <- cells.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(cells, features = top5$gene,angle = 90) + NoLegend()
#write.csv(top5,"top5marker.csv")
#write.csv(top10,"top10marker.csv")
top200 <- cells.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)
#write.csv(top200,"top200marker.csv")
