library(som)
library(beanplot)
library(RColorBrewer)
library(pheatmap)
library(dendextend)

setwd("/NextGenSeqData/project-data/chirag/projects/zebrafishCage/alternativeTranscript/githubDataCode")

# Heatmap of reference and alternative promtoers
A = read.table ("altPromoterTable.txt", header=TRUE)
A1= A[c(1:33)]
M1 = as.matrix(log2(A1[,-1]+1))
rownames(M1)=paste(A$GeneID,A$P1TranscriptID,A$P2TranscriptID, sep="|")
#rownames(M1)=A1$GeneID

# Scale data in the range of 0-1
S1=apply(M1, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))
S2=t(S1)
S2=na.omit(S2)

# Add column annotations
promoter <- factor( c ( rep("Reference",16), rep("Alternative", 16) ) )
col_annotation  <- data.frame(row.names=colnames(S2), promoter)

# Add row annotations
# Do hierarchical clustering and get genes in order
my_hclust_gene <- hclust(dist(S2), method = "complete")
geneOrder =  data.frame(my_hclust_gene$labels, my_hclust_gene$order) 

# cut the tree to get clusters
cluster <- cutree(tree = as.dendrogram(my_hclust_gene), k = 5)
#row_annotation= data.frame(paste("cluster", cluster, sep=""))
#row_annotation= data.frame(cluster)
row_annotation= data.frame(as.factor(cluster))

new_order= cbind(row_annotation, as.factor(A$CorType), as.factor(A$CDSType) )
row_annotation= data.frame(new_order)
colnames(row_annotation)=c("cluster","cor","cds")

# color
my_colour = list( promoter = c(Alternative = "#5977ff", Reference = "#f74747"),
    cor =c ("0" = "grey", "1" = "cyan", "2"= "brown"), 
    cds =c ("0" = "grey", "1" = "cyan", "2"= "brown"),    
    cluster = c("1" = "red", "2" = "blue", "3" = "grey", "4" = "orange", "5"= "purple") )


pdf (file="FigAltPromoterHeatmap.pdf", useDingbats=FALSE )
pheatmap(S2,cluster_cols=0, fontsize=8, annotation_col = col_annotation, annotation_row = row_annotation, cutree_rows = 5, show_rownames = FALSE,  show_colnames = TRUE, annotation_colors = my_colour )
dev.off()



# Supplementary figure
# Supplementary figure
# Compare expression levels of reference and alternative promoters
A = read.table ("altPromoterTable.txt", header=TRUE)
A1= A[c(1:33)]
A1=(log2(A1[,-1]+1))
library(beanplot)

pdf (file="FigExpLevelComparision.pdf", useDingbats=FALSE)
ylim=c(0,14)                           
ylab="log2(TPM)"
col=c("red","cyan")
pos=c(1,2, 4,5, 7,8, 10,11, 13,14, 16,17, 19,20, 22,23, 25,26, 28,29, 31,32, 34,35, 37,38, 40,41, 43,44, 46,47 )

boxplot(A1$P1_fert_egg,A1$P2_fert_egg, A1$P1_s1_cell, A1$P2_s1_cell, A1$P1_s16_cell, A1$P2_s16_cell, A1$P1_s64_cell, A1$P2_s64_cell, A1$P1_s128_cell, A1$P2_s128_cell, A1$P1_s512_cell, A1$P2_s512_cell, A1$P1_high, A1$P2_high, A1$P1_oblong, A1$P2_oblong, A1$P1_sphere, A1$P2_sphere, A1$P1_s30pc_epi, A1$P2_s30pc_epi, A1$P1_shield, A1$P2_shield, A1$P1_s4_somites, A1$P2_s4_somites, A1$P1_s14_19_somites, A1$P2_s14_19_somites,  A1$P1_prim5, A1$P2_prim5, A1$P1_prim25, A1$P2_prim25, A1$P1_longPec, A1$P2_longPec,  outline=T, col=c("red","cyan"), las=2,  ylim=ylim, cex=0.3, notch=T, at=pos, ylab=ylab)

legend("topleft", legend=c("Canonical promoter","Alternative promoter"), fill=col, col=col)
s=t.test(A1$P1_fert_egg, A1$P2_fert_egg)
mtext(s$p.value, side=3, line=-4, adj=0.05)
s=t.test(A1$P1_longPec, A1$P2_longPec)
mtext(s$p.value, side=3, line=-2, adj=0.95)

dev.off()



