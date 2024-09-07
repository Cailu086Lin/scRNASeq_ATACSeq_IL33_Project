library(ChIPQC)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(cinaR)
library(limma)
library(ggvenn)  
library(ggplot2)
library(eulerr)
library(fgsea)
library(data.table)
library(lattice)
library(enrichplot)
library(DOSE)
library(Rsubread)
library(GenomicAlignments)
library(tidyr)
library(DESeq2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tracktables)
library(clusterProfiler)
library(ChIPseeker)
library(DT)
library(patchwork)


peaks <- dir("/home/cailu/data/DBH/ATAC_novogene_240418Project2/cleanedData/usftp21.novogene.com/PeakCalling", pattern = "*._peaks.narrowPeak",
             full.names = TRUE)

myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)
allPeaksSet_nR <- reduce(unlist(GRangesList(myPeaks)))
overlap <- list()
for (i in 1:length(myPeaks)) {
  overlap[[i]] <- allPeaksSet_nR %over% myPeaks[[i]]
}
overlapMatrix <- do.call(cbind, overlap)
colnames(overlapMatrix) <- basename(peaks)
mcols(allPeaksSet_nR) <- overlapMatrix
allPeaksSet_nR[1:2, ]

names(myPeaks) <- c( "A1", "A2",  "A3", "B1","B2","B3", "C1", "C2", "C3")
Group <- factor(c( "A","A","A", "B","B","B", "C","C", "C"))
                  
myGRangesList<-GRangesList(myPeaks)  
reduced <- reduce(unlist(myGRangesList))
consensusIDs <- paste0("consensus_", seq(1, length(reduced)))
mcols(reduced) <- do.call(cbind, lapply(myGRangesList, function(x)(reduced %over% x) + 0))
reducedConsensus <- reduced
mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)), consensusIDs)

consensusToCount<-reducedConsensus

blklist <- import.bed("/home/cailu/data/DBH/ATAC_novogene_240418Project2/cleanedData/usftp21.novogene.com/ENCFF543DDX.bed.gz")
consensusToCount <- consensusToCount[!consensusToCount %over% blklist & !seqnames(consensusToCount) %in% 
                                       "chrM"]

consensusToCount
##

library(limma)
library(ggvenn)  
library(ggplot2)

x<-as.data.frame(elementMetadata(consensusToCount)) %>% 
  dplyr::select(A, B, C)

x[]<-lapply(x,as.logical)


p<-ggplot(x, aes(A=A,B =B, C=C))+    
  geom_venn(fill_color="white",text_size = 6)+
  ggtitle("Overlap for open regions")+theme(plot.title=element_text(size=18, vjust=0, hjust=0.5))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks= element_blank(),
        panel.border = element_blank()
  )
#####
fit2 <- euler(x)
dotplot(resid(fit2), xlab = "",
        panel = function(...) {
          panel.abline(v = 0, lty = 2)
          panel.dotplot(...)
        })

error_plot(fit2)

coef(fit2)
plot(fit2,quantities = TRUE,
fill = "transparent",
lty = 1:3,
labels = list(font = 4), main="Overlap for open regions of BMDC Samples")

postscript(colormodel="cmyk")
ggsave(filename ='/home/cailu/data/DBH/ATAC_novogene/DA/BMDC_venn_plot.tiff',p, width =16, height =14, units ="cm",dpi = 600, bg="white", device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()    


# PCA of overlaps (occupancy analysis)

myPlot <- as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(-consensusIDs) %>% 
  as.matrix %>% t %>% prcomp %>% .$x %>% data.frame %>% mutate(Samples = rownames(.)) %>% 
  mutate(Group = gsub("[0-9]", "", Samples)) %>% ggplot(aes(x = PC1, y = PC2, 
                                                           colour = Group)) + geom_point(size = 5)

myPlot

## Counting for differential ATACseq, echo=F

occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>% 
  rowSums

table(occurrences) %>% rev %>% cumsum

consensusToCount <- consensusToCount[occurrences >= 2, ]

consensusToCount


bamsToCount <- dir("/home/cailu/data/DBH/ATAC_novogene_240418Project2/cleanedData/usftp21.novogene.com/filteredBam/01shifted", full.names = TRUE, pattern = "*.\\.bam$")[c(1:9)]

regionsToCount <- data.frame(GeneID = paste("ID", seqnames(consensusToCount), 
                                            start(consensusToCount), end(consensusToCount), sep = "_"), Chr = seqnames(consensusToCount), 
                             Start = start(consensusToCount), End = end(consensusToCount), Strand = strand(consensusToCount))
fcResults <- featureCounts(bamsToCount, annot.ext = regionsToCount, isPairedEnd = TRUE, 
                           countMultiMappingReads = FALSE, maxFragLength = 100, annot.inbuilt = "mm10",GTF.attrType = "gene_id" )
myCounts <- fcResults$counts

colnames(myCounts) <-c("A1", "A2",  "A3", "B1","B2","B3", "D1", "D2",  "D3", "E1","E2","E3")
  
save(myCounts, file = "/home/cailu/data/DBH/ATAC_novogene_240418Project2/cleanedData/usftp21.novogene.com/ComparisonResults/countsFromATAC_ABDE.RData")

# DESeq2 for differential ATACseq

load( "/home/cailu/data/DBH/ATAC_novogene_240418Project2/cleanedData/usftp21.novogene.com/ComparisonResults/countsFromATAC_ABC.RData")
metaData <- data.frame(Group = c("A","A", "A", "B", "B", "B", "C", "C", "C"))
                         
atacDDS <- DESeqDataSetFromMatrix(myCounts, metaData, ~Group, rowRanges = consensusToCount)
atacDDS <- DESeq(atacDDS)
atac_Rlog <- rlog(atacDDS)
plotPCA(atac_Rlog, intgroup = "Group", ntop = nrow(atac_Rlog))


# DESeq2 for differential ATACseq

##A vs B
wt1 <- results(atacDDS, c("Group","A","B"), format = "GRanges")
wt1<-wt1[order(wt1$pvalue)]
wt1


library(TxDb.Mmusculus.UCSC.mm10.knownGene)
toOverLap <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene, 500, 500)
wt1 <- wt1[(!is.na(wt1$padj) & 
              wt1$padj < .2) & wt1 %over% toOverLap, 
]
wt1 <- annotatePeak(wt1, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
makebedtable(wt1, "BvsD.html", "/home/cailu/data/DBH/ATAC_novogene_240418Project2/cleanedData/usftp21.novogene.com/ComparisonResults")

# Annotation for differential ATACseq

go1 <- enrichGO(as.data.frame(as.GRanges(wt1)[as.GRanges(wt1)$log2FoldChange > 0])$geneId, OrgDb = "org.Mm.eg.db", ont = "BP", maxGSSize = 5000,pvalueCutoff = 0.05,  qvalueCutoff = 0.2)
go2 <- enrichGO(as.data.frame(as.GRanges(wt1)[as.GRanges(wt1)$log2FoldChange < 0])$geneId, OrgDb = "org.Mm.eg.db", ont = "BP", maxGSSize = 5000, pvalueCutoff = 0.05,  qvalueCutoff = 0.2)

head(go1, 10) %>% dplyr::select(ID, Description, pvalue, p.adjust) %>% datatable(elementId = "goEle1")
head(go2, 10) %>% dplyr::select(ID, Description, pvalue, p.adjust) %>% datatable(elementId = "goEle2")

up<-dotplot(go1, showCategory=10) + ggtitle("up_BvsD")
dw<-dotplot(go2, showCategory=10) + ggtitle("down_BvsD")
up+dw


###############
##############

geneList<- as.GRanges(wt1)$log2FoldChange
names(geneList)<-as.GRanges(wt1)$geneId
edo<-gseDO(geneList)
edo <- enrichDGN(de)
go <- enrichGO(as.data.frame(as.GRanges(wt1))$geneId, OrgDb = "org.Mm.eg.db", ont = "ALL", maxGSSize = 5000, pvalueCutoff = 0.01,  qvalueCutoff = .05)
de<-go@gene

edox <- setReadable(go, 'org.Mm.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
p2<-cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 

mutate(go1, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")+ggtitle("GO:up_A_B_BMDC")->up

mutate(go2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")+ggtitle("GO:down_A_B_BMDC")->dw

up1<-dotplot(go1, showCategory=10) + ggtitle("up_A_B_BMDC")
dw1<-dotplot(go2, showCategory=10) + ggtitle("dowp_A_B_BMDC")

up<-enrichWP(wt1$geneId, organism = "Mus musculus")
dotplot(up, showCategory=10) + ggtitle("BMDC_A_B")

#pathway: comaprsion_A vs B
go1 <- enrichGO(as.data.frame(as.GRanges(wt)[as.GRanges(wt)$log2FoldChange > 
                                               0])$geneId, OrgDb = "org.Mm.eg.db", ont = "ALL", maxGSSize = 5000, pvalueCutoff = 0.2,  qvalueCutoff = 0.2)
go2 <- enrichGO(as.data.frame(as.GRanges(wt)[as.GRanges(wt)$log2FoldChange < 
                                               0])$geneId, OrgDb = "org.Mm.eg.db", ont = "ALL", maxGSSize = 5000, pvalueCutoff = 0.1,  qvalueCutoff = 0.2)

head(go1, 10) %>% dplyr::select(ID, Description, pvalue, p.adjust) %>% datatable(elementId = "goEle1")
head(go2, 10) %>% dplyr::select(ID, Description, pvalue, p.adjust) %>% datatable(elementId = "goEle2")

up2<-dotplot(go1, showCategory=10) + ggtitle("up_wt_cKO2")
dw2<-dotplot(go2, showCategory=10) + ggtitle("dowp_wt_cKO2")



##############
df<-data.frame(myCounts)
df$Row.name<-rownames(df)
bed1<-df%>%
  tidyr::separate(Row.name, c("ID", "Chr", "Start","End"), sep = "_")

bed2<-bed1
bed2[,1:9]<-normalizeConsensus(bed2[1:9], norm.method = "cpm", log.option = FALSE)
colnames(bed2)[1:9]<-c("A1", "A2",  "A3", "B1","B2","B3", "C1", "C2", "C3")    
contrasts<- c("A", "A", "A","B", "B", "B","C", "C", "C")
results2 <- cinaR(bed2[c(11:13, 1:9)], contrasts=contrasts, reference.genome = "mm10", enrichment.method="GSEA") 

#A vs B"
p<-dot_plot(results2, fdr.cutoff =0.2, filter.pathways = F)+
  theme(plot.caption =  element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face = "bold"))

postscript(colormodel="cmyk")
ggsave(filename ='/home/cailu/data/DBH/ATAC_novogene_240418Project2/cleanedData/usftp21.novogene.com/ComparisonResults/pathway_ABDE.tiff',p, width =16, height =18, units ="cm",dpi = 600, bg="white",device='tiff', limitsize = TRUE, compression = "lzw")
dev.off() 

##
# Export enrichment results for A_B
write.csv(results2$Enrichment.Results$A_B, file = "~/data/DBH/CD11C_aggr/suppl.Table_ms/enrichment_A_B.csv")

# Export enrichment results for A_C
write.csv(results2$Enrichment.Results$A_C, file = "~/data/DBH/CD11C_aggr/suppl.Table_ms/enrichment_A_C.csv")

# Export enrichment results for B_C
write.csv(results2$Enrichment.Results$B_C, file = "~/data/DBH/CD11C_aggr/suppl.Table_ms/enrichment_B_C.csv")
####


heatmap_differential(results2, angle_col=0,fontsize_col=14,cluster_cols = F, cluster_rows=F)

###

cp <- stats::na.omit(results2[["DA.results"]][["cp"]])

mat.heatmap <- scale_rows(cp)

breaksList <- seq(min(mat.heatmap), max(mat.heatmap), by = .01)
plot.pheatmap <- pheatmap::pheatmap(mat.heatmap, scale = "none",
                                    color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n=7,"RdBu")))(length(breaksList)),
                                    show_rownames = FALSE)

write.csv(mat.heatmap,"~/data/DBH/CD11C_aggr/suppl.Table_ms/mat.heatmap_A_B_C.csv")

########
pca_plot(results2,overlaid.info=contrasts)

##############
PC1 <- PC2 <- NULL

cp <- results2[["DA.results"]][["cp"]]

# eliminate NaN values before-hand if there is any.

pca <- stats::prcomp(t(stats::na.omit(cp)), center = TRUE)

d  <- round(pca$sdev^2/sum(pca$sdev^2)*100, digits=1)

xl <- sprintf("PC 1: %.1f %%", d[1])

yl <- sprintf("PC 2: %.1f %%", d[2])


plot.df <- data.frame(names = rownames(pca$x),PC1 = as.numeric(pca$x[,1]),
                      PC2 = as.numeric(pca$x[,2])) 


write.csv(plot.df, file = "~/data/DBH/CD11C_aggr/suppl.Table_ms/pca_scores.csv")


# Export eigenvalues and variance explained
write.csv(data.frame(PC=c("PC1","PC2"),
                     VarianceExplained = d[1:2]),
          file = "~/data/DBH/CD11C_aggr/suppl.Table_ms/pca_eigenvalues_variance.csv")



###
heatmap_var_peaks(results2, heatmap.peak.count = 78, cluster_cols = F,angle_col=0, cluster_rows=F)


#################done#####




`


