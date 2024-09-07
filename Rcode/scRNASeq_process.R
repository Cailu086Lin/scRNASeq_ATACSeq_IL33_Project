library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)
library(scales)
library(ggplot2)


# Load the CD11C_Control_IL33KO dataset
CD11C_Control_IL33KO.data <- Read10X(data.dir = "/home/cailu/data/DBH/CD11C_aggr/outs/count/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
CD11C_Control_IL33KO <- CreateSeuratObject(counts = CD11C_Control_IL33KO.data, project = "CD11C_Control_IL33KO", min.cells = 3, min.features = 200)
CD11C_Control_IL33KO


CD11C_Control_IL33KO[["percent.mt"]] <- PercentageFeatureSet(CD11C_Control_IL33KO, pattern = "^mt-")

plot1 <- FeatureScatter(CD11C_Control_IL33KO, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by="type", jitter = TRUE)+theme(legend.title=element_blank())+geom_smooth(method="loess")

plot2 <- FeatureScatter(CD11C_Control_IL33KO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by="type")+NoLegend()+geom_smooth(method="loess")


# We filter cells that have unique feature counts over 4000 or less than 500
#We filter cells that have >5% mitochondrial counts

CD11C_Control_IL33KO <- subset(CD11C_Control_IL33KO, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
CD11C_Control_IL33KO <- NormalizeData(CD11C_Control_IL33KO, normalization.method = "LogNormalize", scale.factor = 10000)

#CD11C_Control_IL33KO <- NormalizeData(CD11C_Control_IL33KO)

CD11C_Control_IL33KO <- FindVariableFeatures(CD11C_Control_IL33KO, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(CD11C_Control_IL33KO), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(CD11C_Control_IL33KO)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data
all.genes <- rownames(CD11C_Control_IL33KO)
CD11C_Control_IL33KO <- ScaleData(CD11C_Control_IL33KO, features = all.genes)

#Perform linear dimensional reduction

CD11C_Control_IL33KO <- RunPCA(CD11C_Control_IL33KO, features = VariableFeatures(object = CD11C_Control_IL33KO))


# Examine and visualize PCA results a few different ways
print(CD11C_Control_IL33KO[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(CD11C_Control_IL33KO, dims = 1:2, reduction = "pca")

DimPlot(CD11C_Control_IL33KO, reduction = "pca")

DimHeatmap(CD11C_Control_IL33KO, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(CD11C_Control_IL33KO, dims = 1:15, cells = 500, balanced = TRUE)

#Determine the ‘dimensionality’ of the dataset

CD11C_Control_IL33KO <- JackStraw(CD11C_Control_IL33KO, num.replicate = 100)
CD11C_Control_IL33KO <- ScoreJackStraw(CD11C_Control_IL33KO, dims = 1:20)


JackStrawPlot(CD11C_Control_IL33KO, dims = 1:15)

ElbowPlot(CD11C_Control_IL33KO)

CD11C_Control_IL33KO <- FindNeighbors(CD11C_Control_IL33KO, dims = 1:15)
CD11C_Control_IL33KO <- FindClusters(CD11C_Control_IL33KO, resolution = 0.2)

CD11C_Control_IL33KO <- RunUMAP(CD11C_Control_IL33KO, dims = 1:15)

p1<-DimPlot(CD11C_Control_IL33KO, reduction = "umap", label=T)

CD11C_Control_IL33KO@meta.data$type <-ifelse(gsub(".*-1","Control", colnames(CD11C_Control_IL33KO))=="Control", "Control", "IL33KO")
p2<-DimPlot(CD11C_Control_IL33KO, reduction = "umap", group="type")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(CD11C_Control_IL33KO, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)


# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(CD11C_Control_IL33KO, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)


# find markers for every cluster compared to all remaining cells, report only the positive
# ones
CD11C_Control_IL33KO.markers <- FindAllMarkers(CD11C_Control_IL33KO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
x<-CD11C_Control_IL33KO.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC)


DimPlot(ctr, reduction = "umap", label=T)    
CD11C_Control_IL33KO@meta.data$type <-ifelse(gsub(".*-1","Control", colnames(CD11C_Control_IL33KO))=="Control", "Control", "IL33KO")


#cell annotation
new.cluster.ids <- c("B cells","Immature cDC2/Monocyte like","cDC2","cDC1","Quiescent macrophages","Skin resident T cells","Langerhans cells", "Activated macrophages", "Progenitor cells","Stormal cells/Fibroblasts", "Plasmacytoid dendritic cells")
names(new.cluster.ids) <- levels(CD11C_Control_IL33KO)
CD11C_Control_IL33KO<- RenameIdents(CD11C_Control_IL33KO, new.cluster.ids)


saveRDS(CD11C_Control_IL33KO, file = "skin.rds")
saveRDS(CD11C_Control_IL33KO, file = "CD11_Control_IL33KO.rds")

p1<-DimPlot(CD11C_Control_IL33KO, reduction = "umap", label=T) +ggtitle("Control+IL33KO") +theme(plot.title=element_text(hjust=0.5))
p2<-DimPlot(CD11C_Control_IL33KO, reduction = "umap", group.by="type")+ggtitle("Control+IL33KO")  
p3=DimPlot(subset(CD11C_Control_IL33KO, type=="Control"), reduction = "umap", label=F)+ggtitle("Control")+theme(plot.title=element_text(hjust=0.5))    
p4=DimPlot(subset(CD11C_Control_IL33KO, type=="IL33KO"), reduction = "umap", label=F)+ggtitle("IL33KO")  +theme(plot.title=element_text(hjust=0.5))   



dev.new()
ggsave(filename ='Fig_IL33KO.tiff',p1+p2, width =32, height =10, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()  

#20221006 Figure in subset data (cDC2, cDC2/Monocyte like, Quiescent macrophages and Activated macrophages)
#Fig1

sub <- subset(CD11C_Control_IL33KO,idents=c("Immature cDC2/Monocyte like","cDC2","Quiescent macrophages","Activated macrophages"))
DimPlot(sub, reduction = "umap", label=T)


new.cluster.ids <- c("Immature cDC2\nMonocyte like","cDC2","Quiescent\nmacrophages", "Activated\nmacrophages")
names(new.cluster.ids) <- levels(sub)
sub<- RenameIdents(sub, new.cluster.ids)


######v1
px<-VlnPlot(sub, features = "Il1b",split.by = "type", split.plot=TRUE,y.max=5.5)+
  theme(axis.text.x=  element_blank(),axis.title.x = element_blank(), legend.position = "top", legend.justification='center',
        axis.title.y = element_text(size=20), plot.title = element_blank())+labs(y="Il1b") +scale_fill_manual(values=c("#a6cee3", "#b15928"))

px2<-VlnPlot(sub, features = "Il18", split.by = "type",split.plot=TRUE,y.max=3.5)+
  theme( axis.text.x=  element_blank(), axis.title.x = element_blank(),
         axis.title.y = element_text(size=20), plot.title = element_blank())+NoLegend()+labs(y="Il18")+scale_fill_manual(values=c("#a6cee3", "#b15928"))
 

px3<-VlnPlot(sub, features = "Il6", split.by = "type",split.plot=TRUE,y.max=4.1)+
  theme( axis.text.x=  element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_text(size=20), plot.title = element_blank())+NoLegend()+labs(y="Il6")+scale_fill_manual(values=c("#a6cee3", "#b15928"))



px4<-VlnPlot(sub, features = "Tnf", split.by = "type",split.plot=TRUE,y.max=4.5)+
  theme(axis.title.x=  element_blank(), axis.text.x = element_text(angle=0, hjust = 0.5),
        axis.title.y = element_text(size=20), plot.title = element_blank())+NoLegend()+labs(y="Tnf")+scale_fill_manual(values=c("#a6cee3", "#b15928"))


##
px<-VlnPlot(sub, features = "Il1b",split.by = "type", split.plot=TRUE,y.max=5.5)+
  theme(axis.text.x=  element_blank(),axis.title.x = element_blank(), legend.position = "top", legend.justification='center',
        axis.title.y = element_text(size=20), plot.title = element_blank())+labs(y="Il1b") +scale_fill_manual(values=c("#f9766e", "#068991"))

px2<-VlnPlot(sub, features = "Il18", split.by = "type",split.plot=TRUE,y.max=3.5)+
  theme( axis.text.x=  element_blank(), axis.title.x = element_blank(),
         axis.title.y = element_text(size=20), plot.title = element_blank())+NoLegend()+labs(y="Il18")+scale_fill_manual(values=c("#f9766e", "#068991"))


px3<-VlnPlot(sub, features = "Il6", split.by = "type",split.plot=TRUE,y.max=4.1)+
  theme( axis.text.x=  element_blank(), axis.title.x = element_blank(),
         axis.title.y = element_text(size=20), plot.title = element_blank())+NoLegend()+labs(y="Il6")+scale_fill_manual(values=c("#f9766e", "#068991"))



px4<-VlnPlot(sub, features = "Tnf", split.by = "type",split.plot=TRUE,y.max=4.5)+
  theme(axis.title.x=  element_blank(), axis.text.x = element_text(angle=0, hjust = 0.5),
        axis.title.y = element_text(size=20), plot.title = element_blank())+NoLegend()+labs(y="Tnf")+scale_fill_manual(values=c("#f9766e", "#068991"))



p<-cowplot::plot_grid(px, px2, px3, px4, ncol=1, rel_heights  = c(1.18,1,1,1.25))

p<-cowplot::plot_grid(px, px2,  px4, ncol=1, rel_heights  = c(1.18,1,1.25))

dev.new()
ggsave(filename ='Fig_SubsetVln_v2_coloeMatched.tiff',p, width =14, height =24, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()

##v220221031
px<-VlnPlot(sub, features = "Tgfb1",split.by = "type", split.plot=TRUE,y.max=4.1)+
  theme(axis.text.x=  element_blank(),axis.title.x = element_blank(), legend.position = "top", legend.justification='center',
        axis.title.y = element_text(size=20), plot.title = element_blank())+labs(y="Tgfb1") +scale_fill_manual(values=c("#f9766e", "#068991"))



#####
Tgfb1_data <- FetchData(sub, vars = c("Tgfb1", "type"))
Tgfb1_data$cluster<-sub@active.ident
stat_results <- compare_means(formula = Tgfb1 ~ type, data = Tgfb1_data, group.by="cluster", method = "wilcox.test",paired = FALSE)


# Save results to a CSV file
write.csv(stat_results, "suppl.Table_ms/stat_results_Tgfb1.csv", row.names = FALSE)

##stats differed between ggbur and stats, refer to https://github.com/satijalab/seurat/discussions/6811
####

px2<-VlnPlot(sub, features = "Ccl2", split.by = "type",split.plot=TRUE,y.max=4.3)+
  theme( axis.text.x=  element_blank(), axis.title.x = element_blank(),
         axis.title.y = element_text(size=20), plot.title = element_blank())+NoLegend()+labs(y="Ccl2")+scale_fill_manual(values=c("#f9766e", "#068991"))

Ccl2_data <- FetchData(sub, vars = c("Ccl2", "type"))
Ccl2_data$cluster<-sub@active.ident
stat_results <- compare_means(formula = Ccl2 ~ type, data = Ccl2_data, group.by="cluster", method = "wilcox.test",paired = FALSE)


# Save results to a CSV file
write.csv(stat_results, "suppl.Table_ms/stat_results_Ccl2.csv", row.names = FALSE)


px4<-VlnPlot(sub, features = "Ccl9", split.by = "type",split.plot=TRUE,y.max=5)+
  theme(axis.title.x=  element_blank(), axis.text.x = element_text(angle=0, hjust = 0.5),
        axis.title.y = element_text(size=20), plot.title = element_blank())+NoLegend()+labs(y="Ccl9")+scale_fill_manual(values=c("#f9766e", "#068991"))

####
Ccl9_data <- FetchData(sub, vars = c("Ccl9", "type"))
Ccl9_data$cluster<-sub@active.ident
stat_results <- compare_means(formula = Ccl9 ~ type, data = Ccl9_data, group.by="cluster", method = "wilcox.test",paired = FALSE)


# Save results to a CSV file
write.csv(stat_results, "suppl.Table_ms/stat_results_Ccl9.csv", row.names = FALSE)
####

p<-cowplot::plot_grid(px, px2,  px4, ncol=1, rel_heights  = c(1.18,1,1.25))
dev.new()
ggsave(filename ='/home/cailu/data/DBH/CD11C_aggr/Fig_SubsetVln_v3_coloeMatched_20221031.tiff',p, width =14, height =23, units ="cm",dpi = 600, bg="white",device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()

d <- FetchData(sub, vars = c("type", "Tgfb1", "Ccl2","Ccl9","Il1b" ,"Il18", "Il23a", "Tnf"))
d$cluster<-sub@active.ident

write.csv(d[c(9,1:8)], "suppl.Table_ms/NormalizedGeneExpressiondata.csv")
#Fig2 UMAP

p2<-DimPlot(CD11C_Control_IL33KO, reduction = "umap", group.by="type",raster=TRUE,cols=alpha(c("#f9766e","#068991"),0.5), pt.size = 0.1)+
  theme(legend.position = "bottom",legend.justification='center', plot.title = element_blank(), legend.text = element_text(size=12))+
  guides(color = guide_legend(override.aes = list(size=4)))
  
  
dev.new()
ggsave(filename ='Fig_UMAP_red_green.tiff',p2, width =8, height =8, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()

#Fig3 Cell Protion
d<-data.frame(table(CD11C_Control_IL33KO@meta.data$type, CD11C_Control_IL33KO@active.ident))%>%pivot_wider(names_from=Var1, values_from=Freq)
d$ctr<-d$Control/rowSums(d[,2:3])
d$IL33<-d$IL33KO/rowSums(d[,2:3])
d<-d%>%gather(key, value, -Var2, -Control,-IL33KO)%>%
  mutate(key=recode(key,"ctr"="Control", "IL33"="IL33KO"))

p5<-ggplot(d, aes(x = Var2, y =  value*100, fill = key)) +geom_bar(stat = "identity",color = "black") +
  scale_fill_manual(values=c("#f9766e", "#068991"))+theme_classic2()+
  ylab("Cell Proportion (%) ") +theme(axis.title.y=element_text(size=16), axis.title.x=element_blank(), 
                                      axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1),
                                      legend.text = element_text(size=14),legend.title=element_blank())


dev.new()
ggsave(filename ='Fig_CellProportion_red_green.tiff',p5, width =22, height =12, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()


###20220620
#Fig1
p3=DimPlot(subset(CD11C_Control_IL33KO, type=="Control"), reduction = "umap", label=F)+ggtitle("Control")+theme(plot.title=element_text(hjust=0.5)) +NoLegend()  
p4=DimPlot(subset(CD11C_Control_IL33KO, type=="IL33KO"), reduction = "umap", label=F)+ggtitle("IL33KO")  +theme(plot.title=element_text(hjust=0.5))   

p<-cowplot::plot_grid(p3, p4, ncol=2, rel_widths = c(1,1.8))

dev.new()
ggsave(filename ='Fig_UMAP_LINGOKO_Control_20220609.tiff',p1.0, width =20, height =14, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()


##volin
gv<-c('Il1a',
      'Il1b',
      'Il1r2',
      'Il6',
      'Il6st',
      'Il10',
      'Il11',
      'Il12a',
      'Il12b',
      'Il13',
      'Il15',
      'Il16',
      'Il17a',
      'Il17f',
      'Il18',
      'Il22',
      'Il23a',
      'Il27',
      'Tnf',
      'Ilf2',
      'Ilf3',
      'Ilrun')


d<-FetchData(CD11C_Control_IL33KO, vars=gv)
d$cluster<-CD11C_Control_IL33KO@active.ident
d$type<-CD11C_Control_IL33KO@meta.data[["type"]]

d2<-d%>%
  tidyr::gather(key, value, -cluster, -type)

d2$key<-factor(d2$key, levels=gv) 
d2$cluster<-factor(d2$cluster, levels=new.cluster.ids)
p<-ggplot(subset(d2, type=="Control"), aes(cluster, value, color=cluster))+geom_violin()+facet_grid(key~.,switch = "y")+
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(), axis.ticks.y=element_blank(), strip.text = element_text(size=16))+
  theme(legend.position="none", axis.text.x=element_text(angle=75, vjust=1, hjust=1, size=16), panel.background = element_blank(),strip.background = element_blank())+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black")+ggtitle("Control")+
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5))
  

p2<-ggplot(subset(d2, type=="IL33KO"), aes(cluster, value, color=cluster))+geom_violin()+facet_grid(key~.,switch = "y")+
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(), axis.ticks.y=element_blank(), strip.text = element_blank())+
  theme(legend.position="none", axis.text.x=element_text(angle=75, vjust=1, hjust=1, size=16), panel.background = element_blank(),strip.background = element_blank())+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black")+ggtitle("IL33KO")+
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5))


p3<-cowplot::plot_grid(p, p2, ncol=2, rel_widths=c(1.08,1))


p1<-VlnPlot(subset(CD11C_Control_IL33KO,type=="Control"), features = gv, ncol=1, pt.size = 0)
for (i in 1:21){
  p1[[i]]<-p1[[i]]+scale_y_continuous(expand = c(0, 0))+ theme(plot.title = element_blank(),axis.title.x = element_blank(), axis.text.x= element_blank(),axis.text.y= element_blank())+ylab(gv[i])
}
p1[[22]]<- p1[[22]]+scale_y_continuous(expand = c(0, 0))+ theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.y= element_blank())+ylab(gv[22])
p1a<-cowplot::plot_grid(p1[[1]],p1[[2]],p1[[3]],p1[[4]],p1[[5]],p1[[6]],p1[[7]],p1[[8]],p1[[9]],p1[[10]],
            p1[[11]],p1[[12]],p1[[13]],p1[[14]],p1[[15]],p1[1[6]],p1[[17]],p1[[18]],p1[[19]],p1[[20]],p1[[21]], p1[[22]], ncol=1,
            rel_heights =c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3))
  


p2<-VlnPlot(subset(CD11C_Control_IL33KO,type=="IL33KO"), features = gv, ncol=1, pt.size = 0)
for (i in 1:21){
  p2[[i]]<-p2[[i]]+scale_y_continuous(expand = c(0, 0))+ theme(plot.title = element_blank(),axis.title.x = element_blank(), axis.text.x = element_blank(),axis.text.y= element_blank())+ylab("")
}

p2[[22]]<- p2[[22]]+scale_y_continuous(expand = c(0, 0))+ theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.y= element_blank())+ylab("")
p2a<-cowplot::plot_grid(p2[[1]],p2[[2]],p2[[3]],p2[[4]],p2[[5]],p2[[6]],p2[[7]],p2[[8]],p2[[9]],p2[[10]],
                        p2[[11]],p2[[12]],p2[[13]],p2[[14]],p2[[15]],p2[1[6]],p2[[17]],p2[[18]],p2[[19]],p2[[20]],p2[[21]], p2[[22]], ncol=1,
                        rel_heights =c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3))

p3<-cowplot::plot_grid(p1, p2, ncol=2, label=c("Control", "IL33KO"), label_size = 16)


dev.new()
ggsave(filename ='Vlnplot_20220620.tiff',p3, width =25, height =55, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()


##### Neuropeptide receptor gene expression 

CD11C_Control_IL33KO<-readRDS("/home/cailu/data/DBH/CD11C_aggr/CD11_Control_IL33KO.rds")

g=c('Dbi', 'Calca', 'Scg5','Calcb', 'Gal', 'Nmb', 'Penk',  'Oxt', 'Galr3','Oprd1', 'Avpr2', 'Ramp3','Npffr1', 'Bdkrb1',  'Npr1','Crcp', 'Calcrl', 'Ramp1', 'Uts2r', 'Galr2','Gcgr',
'Bdkrb2', 'Vipr1', 'Tacr1', 'Ramp2', 'Npr2','Npr3', 'Vipr2',  'Npy1r','Oprm1')

p<-DoHeatmap(subset(CD11C_Control_IL33KO,idents = c("Immature cDC2/Monocyte like","Quiescent macrophages","Activated macrophages")), features = g)



cols<-c("black","black","black","green","green","red","red","black","black","black",
       "black","black","black","black","black","black","black","black","black","green",
       "black","black","red","black","black","black","black","green","green","black") ##3 clusters

cols<-c("black","green","red","black","black","black","black","black","black","black",
        "black","black","red","black","red","green","green","black","black","black",
        "black","black","green","black","green","black","black","green","black","black") ##all clusters

cols<-c("black","green","red","black","black","black","black","black","black","black",
        "black","black","red","black","red","green","green","black","black","black",
        "black","black","green","black","green","black","black","green","black","black") ##all clusters without Dbi



subset<-subset(CD11C_Control_IL33KO,idents = c("Immature cDC2/Monocyte like","Quiescent macrophages","Activated macrophages"))
subset<-CD11C_Control_IL33KO

d<-FetchData(subset, vars=g)
d$cluster<-subset@active.ident
d$type<-subset@meta.data[["type"]]

d2<-d%>%
  gather(key, value, -cluster, -type)%>%
  group_by(cluster,type, key)%>%
  summarize(exp=mean(value))

d2$exp<-rescale(log2(d2$exp+1))

library("RColorBrewer")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_fill_gradientn(colours = myPalette(100), limits=c(-0.1,1),breaks=c(-0.1, 1), labels=c("Min", "Max"))


p<-ggplot(d2%>%mutate(key = forcats::fct_reorder(key, desc(-exp))), aes(cluster, key, fill=exp)) + geom_tile(color="black")+facet_grid(.~type,switch = "x")+theme_classic()+sc+scale_x_discrete(position = "top") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle=60, size=14, vjust = 1, hjust=0), axis.title = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size=16), axis.text.y =  element_text(color = cols,size=14),
       legend.key.height = unit(2, "cm"),
        legend.text=element_text(size=14, angle=0))+expand_limits( y= 0)+
  guides(fill=guide_colourbar(title.theme = element_blank(), label.hjust = 0)) 


p1<-p+labs(title = "Gal and receptors",
            subtitle = "Cgrp and receptors")+theme(plot.title = element_text(size = 16, color="red",hjust=-0.1, vjust=-9),
                                                   plot.subtitle = element_text(size = 16, color="green", hjust=-0.1, vjust=-10))
  
  
  ###
g<-neuropeptide.receptors <- read.csv("~/data/DBH/Joy_pub/neuropeptide receptors", sep="", header=F)
d<-FetchData(CD11C_Control_IL33KO, vars=c("type",g$V1))
d$Cluster<-CD11C_Control_IL33KO@active.ident
write.csv(d,"~/data/DBH/CD11C_aggr/suppl.Table_ms/ExpressionNeuropeptideReceptors.csv") 
###v2
ctr<-readRDS("/home/cailu/data/DBH/datset/CD11CCRE/outs/CD11_Control.rds")
p1<-DimPlot(ctr, reduction = "umap", label=F)
new.cluster.ids <- c("B cells","immature cDC2/Monocyte like","Skin resident T cells","Quiescent macrophages","Plasmacytoid dendritic cells","Activated macrophages","cDC","Langerhans cells",  "not defined","Stormal cells/Fibroblasts" )
names(new.cluster.ids) <- levels(ctr)
ctr<- RenameIdents(ctr, new.cluster.ids)

subset<-subset(ctr,idents = c("Immature cDC2/Monocyte like","Quiescent macrophages","Activated macrophages"))

d<-FetchData(subset, vars=g)
d$cluster<-subset@active.ident
d$type<-"Control"

d2c<-d%>%
  gather(key, value, -cluster, -type)%>%
  group_by(cluster,type, key)%>%
  summarize(exp=mean(value))

###IL33KO

ko<-readRDS("/home/cailu/data/DBH/datset/CD11CCREIL33KO/outs/CD11_IL33KO.rds")
p1<-DimPlot(ko, reduction = "umap", label=F)
new.cluster.ids <- c("immature cDC2/Monocyte like","cDC1","cDC2","Quiescent macrophages","Progenitor cells","B cells","Langerhans cells",  "Stormal cells/Fibroblasts", "Skin resident T cells","Activated macrophages")
names(new.cluster.ids) <- levels(ko)
ko<- RenameIdents(ko, new.cluster.ids)

subset<-subset(ko,idents = c("Immature cDC2/Monocyte like","Quiescent macrophages","Activated macrophages"))

d<-FetchData(subset, vars=g)
d$cluster<-subset@active.ident
d$type<-"IL33KO"

d2k<-d%>%
  gather(key, value, -cluster, -type)%>%
  group_by(cluster,type, key)%>%
  summarize(exp=mean(value))

d2<-bind_rows(d2c,d2k)
d2$exp<-rescale(log10(d2$exp+1))


p<-ggplot(d2%>%mutate(key = forcats::fct_reorder(key, desc(-exp))), aes(cluster, key, fill=exp)) + geom_tile(color="black")+facet_grid(.~type,switch = "x")+theme_classic()+sc+scale_x_discrete(position = "top") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle=40, vjust = 1, hjust=0), axis.title = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size=16), axis.text = element_text(size=14), legend.key.height = unit(2.5, "cm"),
        legend.text=element_text(size=14, angle=0))+expand_limits( y= 0)+
  guides(fill=guide_colourbar(title.theme = element_blank(), label.hjust = 0))


##simpsin index plot
d<-data.frame(table(CD11C_Control_IL33KO@meta.data$type, CD11C_Control_IL33KO@active.ident))%>%pivot_wider(names_from=Var1, values_from=Freq)
d$ctr<-d$Control/sum(d$Control)
d$IL33<-d$IL33KO/sum(d$IL33KO) 
d<-d%>%gather(key, value, -Var2, -Control,-IL33KO)

p5<-ggplot(d, aes(x = key, y =  value, fill = Var2)) +geom_bar(stat = "identity",color = "black") +
  scale_fill_brewer(palette="Set3",direction=-1)+scale_y_continuous(labels = scales::percent)+theme_classic2()+
  ylab("Percentage of cell types in skin") +scale_x_discrete(labels=c("Control", "IL33KO"))+theme(axis.title.y=element_text(size=16), axis.title.x=element_blank(), axis.text=element_text(size=12),legend.title=element_blank())


#
d2<-data.frame(table(CD11C_Control_IL33KO@meta.data$type, CD11C_Control_IL33KO@active.ident))%>%pivot_wider(names_from=Var1, values_from=Freq)
d2$Fract1<-d2$Control/(d2$Control+d2$IL33KO) 
d2$Fract2<-1-d2$Fract1
d2<-d2%>%gather(key, value, -Var2, -Control,-IL33KO)%>%mutate(key=recode(key,"Fract1"="Control", "Fract2"="IL33KO"))

d2$key<-factor(d2$key, levels=c("IL33KO", "Control"))

p5<-ggplot(d2, aes(x = value, y =  Var2, fill = key)) +geom_bar(stat = "identity",color = "black") +
  scale_fill_brewer("blue",direction=-1, guide = guide_legend(reverse = TRUE)
  )+scale_x_continuous(labels = scales::percent)+theme_classic2()+xlab("Fraction of cells")+theme(axis.title.x=element_text(size=16), axis.title.y=element_blank(), axis.text=element_text(size=12),legend.title=element_blank(), legend.position="top", legend.text=element_text(size=14))

#
d3<-data.frame(table(CD11C_Control_IL33KO@active.ident))

p6<-ggplot(d3, aes(x = Freq, y =  Var1)) +geom_bar(stat = "identity",color = "black")+theme_classic2()+xlab("Number of cells") +theme(axis.title.x=element_text(size=16), axis.title.y=element_blank(), axis.text.y=element_blank(),axis.text=element_text(size=12))



dev.new()
ggsave(filename ='Fig_cellFract.tiff',p5+p6, width =20, height =13, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()  

##cell type diversity.

simpson<-function(data,type="complement"){
  simpson.diversity<-numeric(ncol(data))
  for(j in 1:ncol(data)){
    soma<-sum(data[,j])
    prop<-data[,j]/soma
    prop2<-prop^2
    D<-sum(prop2)
    if(type=="inverse") (simp<-1/D)
    if(type=="complement") (simp<-1-D)
    simpson.diversity[j]<-simp
  }
  plot.number<-1:ncol(data)
  return(rbind(plot.number,simpson.diversity))
}
plots<-table(CD11C_Control_IL33KO@active.ident,CD11C_Control_IL33KO@meta.data[["type"]])
d<-simpson(plots)
d2<-data.frame(key=c("Control", "IL33KO"), index=c(d[2,1], d[2,2]))

p<-ggplot(d2,aes(key, index, fill=key))+geom_bar(stat="identity", fill=c("black", "black")) + theme_classic2()+ylab("Simpson index")+
  theme(axis.text=element_text(size=12, face="bold"), axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold"),legend.title=element_blank(), legend.text=element_text(size=12))+ylim(0,1)+scale_x_discrete(limits =c("Control", "IL33KO"))

##heatmap

library(ComplexHeatmap)
# https://github.com/immunogenomics/presto
library(presto)
library(tictoc)

IL<-rownames(CD11C_Control_IL33KO)[str_detect(rownames(CD11C_Control_IL33KO), pattern ="Il|Tnf|Tlsp|Gmcsf")]

p7<-DoHeatmap(CD11C_Control_IL33KO, features=IL, label=F, group.by='ident')

dev.new()
ggsave(filename ='Fig_heatmap.tiff',p7, width =30, height =28, units ="cm",dpi = 300, device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()

mat<- CD11C_Control_IL33KO[["RNA"]]@data[IL, ] %>% as.matrix()

cluster_anno<- data.frame(CD11C_Control_IL33KO@meta.data$type,CD11C_Control_IL33KO@meta.data$seurat_clusters)
colnames(cluster_anno) <- c( 'type','ident')

d<-data.frame(t(mat))

d<-cbind(d, cluster_anno)

d2<-d%>%gather(gene, value, -type,-ident)


d3<-d2%>%
  filter(value !=0)%>%
  group_by(gene,type,ident)%>%
  summarise(av.ex= median(value), n=n())

d4<-d3%>%pivot_wider(names_from="type", values_from=c("av.ex", "n"))


tn<-data.frame(table(CD11C_Control_IL33KO@meta.data$type,CD11C_Control_IL33KO@meta.data$seurat_clusters))%>%
  pivot_wider(names_from="Var1", values_from="Freq")

d5<-merge(d4, tn, by.x="ident", by.y="Var2")
d5$prt_Cell_Control<-100*d5$n_Control/d5$Control 
d5$prt_Cell_IL33KO<-100*d5$n_IL33KO/d5$IL33KO 

d5[is.na(d5)]<-0

#wilcox.test(x, y, alternative = "two.sided")

d6<-d2%>%
  filter(value !=0)

results1<-data.frame()
for (i in unique(d6$gene)){
  for (j in unique(d6$ident)){
    result<-data.frame(statistic=NA,pvale=NA,ident=NA,gene=NA)
    d7<-subset(d6, gene==i&ident==j)
    if(length(unique(d7$type))>1){
      t<-wilcox.test(value~type,data=d7,alternative = "two.sided")
      result$statistic <-t$statistic[[1]] 
      result$pvale<-t$p.value[[1]] 
      result$ident<-j
      result$gene<-i
    }
    else{
      result$statistic<-"not_test"
      result$pvale<-"NA"
      result$ident<-j
      result$gene<-i
    }
    results1<-rbind(results1, result)
  }
}



Results_final<-merge(d5, results1, by=c("gene", "ident"))

##plot gene
tg<-c(
  'Il13',
  'Il16',
  'Il18',
  'Il1b',
  'Il1r2',
  'Il2rb',
  'Il3ra',
  'Il7r',
  'Ilf2',
  'Ilf3',
  'Tnfrsf18',
  'Tnfsf13')




d8<-d2[d2$gene %in%tg,]
d8$ident<-as.factor(d8$ident)
p10<- ggplot(d8, aes(x=ident, y=value, color=ident)) + 
  geom_violin(trim=FALSE) +facet_grid(gene~type,scales="free_y")+theme_classic2()+theme(legend.position="none")



## scale the rows
mat<- t(scale(t(mat)))

cluster_anno<- data.frame(CD11C_Control_IL33KO@meta.data$type,CD11C_Control_IL33KO@meta.data$seurat_clusters)
colnames(cluster_anno) <- c( 'type','ident')
colours <- list(
  'ident' = c('0'="#FF00FF",'1'= "#CA00CA",'2'= "#950095",'3'= "#6A006A",'4'= "#200020",'5'= "#202000",'6'="#6A6A00" , '7'="#B5B500" ,'8'= "#AA00AA",'9'="#FFFF00", '10'= "#9F009F"),
  'type' = c('Control' = 'red2', 'IL33KO' = 'royalblue'))
colAnn <- HeatmapAnnotation(df = cluster_anno,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))


#We can explicitly map the colors to the scaled expression values

# what's the value range in the matrix
quantile(mat, c(0.1, 0.95))

##        10%        95% 
## -0.5277426  2.3470090

Seurat::PurpleAndYellow()


###

##  [1] "#FF00FF" "#F400F4" "#EA00EA" "#DF00DF" "#D500D5" "#CA00CA" "#BF00BF"
##  [8] "#B500B5" "#AA00AA" "#9F009F" "#950095" "#8A008A" "#800080" "#750075"
## [15] "#6A006A" "#600060" "#550055" "#4A004A" "#400040" "#350035" "#2B002B"
## [22] "#200020" "#150015" "#0B000B" "#000000" "#000000" "#0B0B00" "#151500"
## [29] "#202000" "#2B2B00" "#353500" "#404000" "#4A4A00" "#555500" "#606000"
## [36] "#6A6A00" "#757500" "#808000" "#8A8A00" "#959500" "#9F9F00" "#AAAA00"
## [43] "#B5B500" "#BFBF00" "#CACA00" "#D4D400" "#DFDF00" "#EAEA00" "#F4F400"
## [50] "#FFFF00"

## make the black color map to 0. the yellow map to highest and the purle map to the lowest
col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))

#plot the heatmap

p7<-Heatmap(mat, name = "Expression",  
            column_split = factor(cluster_anno$ident),
            cluster_columns = TRUE,
            show_column_dend = FALSE,
            cluster_column_slices = TRUE,
            column_title_gp = gpar(fontsize = 8),
            column_gap = unit(0.5, "mm"),
            cluster_rows = TRUE,
            show_row_dend = FALSE,
            col = col_fun,
            row_names_gp = gpar(fontsize = 4),
            column_title_rot = 90,
            top_annotation=colAnn,
            show_column_names = FALSE,
            use_raster = TRUE,
            raster_quality = 4),
top_annotation_height=unit(1.0,"cm"))



p7 <- Heatmap(
  mat,
  name = "expression",show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  show_row_dend = TRUE,
  row_dend_reorder = TRUE,
  column_dend_reorder = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  width = unit(100, "mm"),
  
  top_annotation_height=unit(1.0,"cm"), top_annotation=colAnn)

#top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),


##Featureplor

g<-c("Cd79a","Cd79b","Ms4a1","Ear2","Cd209a","Lyz2","Apod","Xcr1","Cd207","Cxcl2","Csf1r","Nkg7","Trbc2","Cd207",
     "Epcam","Ccl5","Ube2c","Birc5", "Igfbp7","Siglech")

p<- FeaturePlot(CD11C_Control_IL33KO, c("Il12b",
                        "Il4i1",
                        "Tnfrsf4",
                        "Ilk",
                        "Il1r2"),
                split.by="type")


p1<- RidgePlot(subset(CD11C_Control_IL33KO,type=="Control"), features="Igkc")+NoLegend()+theme(axis.title.y = element_blank())+ggtitle("Control")
p2<- RidgePlot(subset(CD11C_Control_IL33KO,type=="IL33KO"), features="Igkc")+NoLegend()+theme(axis.title.y = element_blank(), axis.text.y=element_blank())+ggtitle("IL33KO")
p<-cowplot::plot_grid(p1,p2, ncol=2, rel_widths = c(2,1.1))


for (i in 1:10){
  p[[i]]<-p[[i]]+ NoLegend()+NoAxes()# +coord_fixed()#+theme(panel.border = element_rect(color="black"))#,plot.title = element_text(vjust = - 5, hjust=0.1), plot.margin=margin(0,0,0,0))
}

dev.new()
ggsave(filename ='Fig_featureplottiff',p, width =15, height =20, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()  


p<- VlnPlot(subset(CD11C_Control_IL33KO, type=="Control"), features = tg[c(2,4:7)], pt.size=0, combine = FALSE )
for(i in 1:(length(p)-1)){
  p[[i]] <- p[[i]] + NoLegend()+NoAxes()+coord_flip()
}
p[[5]] <- p[[5]]+ NoLegend()+theme(axis.title=element_blank(),axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+coord_flip()


p1<- VlnPlot(subset(CD11C_Control_IL33KO, type=="IL33KO"), features = tg[c(2,4:7)], pt.size=0, combine = FALSE)
for(i in 1:(length(p)-1)){
  p1[[i]] <- p1[[i]] + NoLegend()+NoAxes()
}

p1[[5]] <- p1[[5]]+ NoLegend()+theme(axis.title=element_blank(),axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

px<- VlnPlot(subset(CD11C_Control_IL33KO, type=="Control"), features = tg[c(8:12)], pt.size=0, combine = FALSE)
for(i in 1:(length(p)-1)){
  px[[i]] <- px[[i]] + NoLegend()+NoAxes()
}

px[[5]] <- px[[5]]+ NoLegend()+theme(axis.title=element_blank(),axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

p1x<- VlnPlot(subset(CD11C_Control_IL33KO, type=="IL33KO"), features = tg[c(8:12)], pt.size=0, combine = FALSE)
for(i in 1:(length(p)-1)){
  p1x[[i]] <- p1x[[i]] + NoLegend()+NoAxes()
}

p1x[[5]] <- p1x[[5]]+ NoLegend()+theme(axis.title=element_blank(),axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

title_gg <- ggplot() + 
  labs(title = "Control") + 
  theme(plot.title=element_text(size=20, face="bold", color="blue",hjust=0.5))

title_gg2 <- ggplot() + 
  labs(title = "IL33KO") + 
  theme(plot.title=element_text(size=20,face="bold",color="red", hjust=0.5))
py1<-cowplot::plot_grid(plotlist = c(p), ncol=1, rel_heights=c(1,1,1,1,2), scale=1.0)
py2<-cowplot::plot_grid(plotlist = c(p1), ncol=1,rel_heights=c(1,1,1,1,2), scale=1.0)
py3<-cowplot::plot_grid(plotlist = c(px), ncol=1,rel_heights=c(1,1,1,1,2), scale=1.0)
py4<-cowplot::plot_grid(plotlist = c(p1x), ncol=1, rel_heights=c(1,1,1,1,2),scale=1.0)

pz1<-cowplot::plot_grid(title_gg, py1, ncol = 1, rel_heights = c(0.028, 1))
pz2<-cowplot::plot_grid(title_gg2, py2, ncol = 1, rel_heights = c(0.028, 1))
pz3<-cowplot::plot_grid(title_gg, py3, ncol = 1, rel_heights = c(0.028, 1))
pz4<-cowplot::plot_grid(title_gg2, py4, ncol = 1, rel_heights = c(0.028, 1))


pz<- cowplot::plot_grid(pz1,pz2, pz3,pz4, ncol=4)  


dev.new()
ggsave(filename ='IL_top10.tiff',pz, width =30, height =28, units ="cm",dpi = 300, device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()

pdf( file ='IL_top10.pdf', width =24, height =16)
pz
dev.off() 


p<-c()
for (i in 1:length(IL)){  
  p[[i]]<-VlnPlot(CD11C_Control_IL33KO, features = IL[i], pt.size=0)+coord_flip()+theme(axis.title=element_blank(),
                                                                        axis.text=element_text(angle=0))
}  


#Seurat has several tests for differential expression which can be set with the test.use parameter (see our #DE vignette for details). For example, the ROC test returns the ‘classification power’ for any individual #marker (ranging from 0 - random, to 1 - perfect).

cluster0.markers <- FindMarkers(CD11C_Control_IL33KO, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(CD11C_Control_IL33KO, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(CD11C_Control_IL33KO, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(CD11C_Control_IL33KO, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

CD11C_Control_IL33KO.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(CD11C_Control_IL33KO, features = top10$gene) + NoLegend()

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(CD11C_Control_IL33KO)
CD11C_Control_IL33KO <- RenameIdents(CD11C_Control_IL33KO, new.cluster.ids)
DimPlot(CD11C_Control_IL33KO, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


saveRDS(CD11C_Control_IL33KO, file = "DBM_aggr.rds")

dev.new()
ggsave(filename ='DBM_aggr_UAMPPlot.tiff',p, width =8, height =8, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()


###Zoom in B cells

bc<-subset(CD11C_Control_IL33KO, ident="0")

bc <- RunPCA(bc, features = VariableFeatures(object = bc))
bc <- JackStraw(bc, num.replicate = 100)
bc <- ScoreJackStraw(bc, dims = 1:20)
JackStrawPlot(bc, dims = 1:15)
ElbowPlot(bc)
bc <- FindNeighbors(bc, dims = 1:15)
bc <- FindClusters(bc, resolution = 0.2)
bc <- RunUMAP(bc, dims = 1:15)

p1<-DimPlot(bc, reduction = "umap", label=T)


saveRDS(bc, file = "B_Cell_SubClusters.rds")

bc.markers <- FindAllMarkers(bc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
x<-bc.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC)

write.csv(x, "Markers_B_cell_subcluster.csv", quote = F, row.names = F)


##Cell anaotation
library(scCATCH)
bc <- findcelltype(bc)
mouse_skin <- rev_gene(data = mouse_skin, data_type = "data", species = "Mouse", geneinfo = geneinfo)



#
p1<-DimPlot(bc, reduction = "umap", label=T) 
p2<-DimPlot(bc, reduction = "umap", group.by="type") +theme(plot.title=element_blank())  

d2<-data.frame(table(bc@meta.data$type, bc@active.ident))%>%pivot_wider(names_from=Var1, values_from=Freq)
d2$Fract1<-d2$Control/(d2$Control+d2$IL33KO) 
d2$Fract2<-1-d2$Fract1
d2<-d2%>%gather(key, value, -Var2, -Control,-IL33KO)%>%mutate(key=recode(key,"Fract1"="Control", "Fract2"="IL33KO"))

d2$key<-factor(d2$key, levels=c("IL33KO", "Control"))

p5<-ggplot(d2, aes(x = value, y =  Var2, fill = key)) +geom_bar(stat = "identity",color = "black") +
  scale_fill_brewer("blue",direction=-1, guide = guide_legend(reverse = TRUE)
  )+scale_x_continuous(labels = scales::percent)+theme_classic2()+xlab("Fraction of cells")+theme(axis.title.x=element_text(size=16), axis.title.y=element_blank(), axis.text=element_text(size=12),legend.title=element_blank(), legend.position="top", legend.text=element_text(size=14))

#
d3<-data.frame(table(bc@active.ident))

p6<-ggplot(d3, aes(x = Freq, y =  Var1)) +geom_bar(stat = "identity",color = "black")+theme_classic2()+xlab("Number of cells") +theme(axis.title.x=element_text(size=16), axis.title.y=element_blank(),axis.text=element_text(size=12))




p<-cowplot::plot_grid(p1, p2,p5,p6, labels = c('A', 'B','C','D'), label_size=16,align = "h")


p<-ggarrange(ggarrange(p1,p2,labels = c("A","B"),ncol=2),
             ggarrange(p5,p6, labels = c("C", ""),ncol=2), ncol=1, align="h"                                 
)  


pdf(file = "IL33QC_S1.pdf", width=10, height = 12) 

dev.new()
ggsave(filename ='bc_subtype_UAMPPlot.tiff',p1+p2+p5+p6, width =18, height =15, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()

##
g<-c('Ly6e',
     'Cd83',
     'Ccr7',
     'Il2rg',
     'Cd86',
     'Ly6e',
     'Cd83',
     'Cd9',
     'Cd44',
     'Ahr',
     'Cd72',
     'Iglc2',
     'Iglc3',
     'Iglv3',
     'Klf2',
     'Ikzf3',
     'Cd55',
     'Ifit3',
     'Ifit2',
     'Irf7',
     'Tnfrsf17',
     'Igkv4-55',
     'Trf',
     'Ly6c2',
     'Iglv1',
     'Ighv2-3',
     'Iglv3')

p9<- DoHeatmap(bc, features=g, label=F, group.by='ident')

dev.new()
ggsave(filename ='bc_heatmap.tiff',p9, width =20, height =18, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()  




#featurePlot
g<-c('Cd79a',
     'Lyz2',
     'Tcf4',
     'Flt3',
     'Irf8',
     'Cxcl2',
     'Trbc2',
     'Rora',
     'Ctla2a',
     'Ccl5',
     'Cd207',
     'Pclaf',
     'Birc5',
     'Siglech',
     'Col1a2')

p<-FeaturePlot(CD11C_Control_IL33KO, features=g, ncol=5)
for (i in 1:15){
  p[[i]]<-p[[i]]+NoLegend()+NoAxes()+theme(panel.border = element_rect(color = "black", fill = NA))
}


###
g<-c('Hbb-bs',
     'Hba-a1',
     'Hbb-bt',
     'Rnaset2b',
     'Hba-a2',
     'Klk10',
     'Mid1',
     'Prss53',
     'Padi4',
     'Tubb3',
     'Hunk',
     'Cpm',
     'Cdh6',
     'Hoxc13',
     'Dlx3',
     'Gabrp',
     'Bambi',
     'Dct',
     'Krtap3-2',
     'Gm8276',
     'Krt36',
     'Tyrp1',
     'Car6',
     'Krt32',
     '2210017I01Rik',
     'Msx2',
     'Lef1',
     'Slc24a5',
     'S100a7a',
     'Dlx2',
     'Gjb2',
     'Dlx4',
     'Psors1c2',
     'Krt84',
     'Lrrc15',
     'Dlx1',
     'Pmel',
     'Otop2',
     'Adh6b',
     'Ly6g6g',
     'Ly6g6d',
     'Tagln3',
     'Gnmt',
     'A030005L19Rik',
     'Krtap17-1',
     'Pinlyp',
     'Crym',
     'Krt74',
     'Shh',
     'Tyr',
     'Foxe1',
     'Krt82',
     'Padi1',
     'Mlana',
     'Krtap3-3',
     'Krt28',
     'S100a3',
     'Capn8',
     'Krt83',
     'Krt87',
     'Fbp1',
     'Tchh',
     'Vsig8',
     'Krt31',
     'Prr9',
     'Krt73',
     'Krt35',
     'Krt27',
     'Krt71',
     'Krt25','Glo1-ps',
     'Rnps1-ps',
     'Klk1b22',
     'Pgk1-rs7',
     'Hmga1b',
     'Adgrg7',
     'Rps2-ps10',
     'Krt6b',
     'Aadacl3',
     'Gm45855',
     'Clca3a2',
     'Mup11',
     'Scd3',
     'Il33',
     'Clca3a1',
     'Gbp9',
     'Inka2',
     'Gbp7',
     'Tmem181b-ps',
     'Ikzf2',
     'H2-Q5',
     'Ppp2r2c',
     'Ano1')
  

scRNAP<-CD11C_Control_IL33KO
scRNAP$celltype<- paste0(Idents(scRNAP),scRNAP$type, "_") #dents sample subgruops to the cell type idents
new.idents <- scRNAP@meta.data$celltype
Idents(object = scRNAP) <- new.idents 
p<-DotPlot(scRNAP, cols = c("red"), dot.scale = 8, features =g) + RotatedAxis() +
  theme(axis.text = element_text(size = 14))+scale_colour_gradient(low =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Avg.Exp'))+ guides(size=guide_legend(title="Pert.Exp"))+
  scale_y_discrete(limits=c("B cellsIL33KO_", "B cellsControl_","Immature cDC2/Monocyte likeIL33KO_", "Immature cDC2/Monocyte likeControl_",	
                            "cDC2IL33KO_", "cDC2Control_",	"cDC1IL33KO_","cDC1Control_",	
                            "Quiescent macrophagesIL33KO_", "Quiescent macrophagesControl_",	"Skin resident T cellsIL33KO_","Skin resident T cellsControl_",
                            "Langerhans cellsIL33KO_","Langerhans cellsControl_","Activated macrophagesIL33KO_", "Activated macrophagesControl_",
                            "Progenitor cellsIL33KO_", "Progenitor cellsControl_", "Stormal cells/FibroblastsIL33KO_", "Stormal cells/FibroblastsControl_",
                            "Plasmacytoid dendritic cellsIL33KO_","Plasmacytoid dendritic cellsControl_"),
                   labels=c("B cells  IL33KO", "Control","Immature cDC2/Monocyte like  IL33KO", "Control",	
                            "cDC2  IL33KO", "Control",	"cDC1  IL33KO","Control",	
                            "Quiescent macrophages  IL33KO", "Control",	"Skin resident T cells  IL33KO","Control",
                            "Langerhans cells  IL33KO","Control","Activated macrophages  IL33KO", "Control",
                            "Progenitor cells  IL33KO", "Control", "Stormal cells/Fibroblasts  IL33KO", "Control",
                            "Plasmacytoid dendritic cells  IL33KO","Control"))+
  theme(axis.title = element_blank()) +
  geom_hline(yintercept =2.5, linetype=2, color='grey') +
  geom_hline(yintercept =4.5, linetype=2, color='grey')+
  geom_hline(yintercept =6.5, linetype=2, color='grey')+
  geom_hline(yintercept =8.5, linetype=2, color='grey')+
  geom_hline(yintercept =10.5, linetype=2, color='grey')+
  geom_hline(yintercept =12.5, linetype=2, color='grey')+
  geom_hline(yintercept =14.5, linetype=2, color='grey')+
  geom_hline(yintercept =16.5, linetype=2, color='grey')+
  geom_hline(yintercept =18.5, linetype=2, color='grey')+
  geom_hline(yintercept =20.5, linetype=2, color='grey')



######Done####
