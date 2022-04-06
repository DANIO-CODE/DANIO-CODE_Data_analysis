

setwd("/NextGenSeqData/project-data/chirag/projects/zebrafishCage/alternativeTranscript/githubDataCode")

A = read.table("KEGGPathways.txt", sep="\t", header=TRUE)
A1=A[order(A$P.value),]
A1=A1[A1$P.value <=0.05,]

pdf (file="FigKEGGPathway.pdf", useDingbats=FALSE)
par (mar=c(5,17,2,2))
barplot(rev(-log2(A1$P.value)), horiz=T, las=1, names.arg=rev(A1$Term),xlim=c(0,10), xlab="-log2(P-value)")
abline(v=-log2(0.05), lty=2)
box()
dev.off()







