This directory contains 3 R scripts that should be executed in the following order:

1) **constructPromoterome.Rmd** promoterome construction and Supplementary Figure 4b
2) **exprClustering.Rmd**: promoterome expression clustering by self-organising maps (part of Figure 4f panel) 
3) **promoterClassification.Rmd** promoter chromatin architecture clustering with k-means, and the creation of panels for Figure 4 (a, b, d, e, f, g) (figure 4c is generated with figure 3 code) and all panels for Supplementary Figure 13

MyLibs directory contains R helper functions in R.
../../Code/Scripts contains shell and Perl scripts for processing ATAC and ChIP-seq data from bam file to input files reqiured in stage 3.
