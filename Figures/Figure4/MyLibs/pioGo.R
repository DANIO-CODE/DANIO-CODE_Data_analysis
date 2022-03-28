
# Find enriched gene ontologies:
# arguments
#   - thisGenes - a vector of gene ids, eg ENSDARG000000001
#   - allGenes  - background of all genes. Set to all genes you included in your analysis:
#                 differential accessibility, etc
#   - txDb      - a txDb object for a given genome

# For annotating distal elements (not promoters) you need to run annotation first with pioAnnotateRegions
# set a cutoff on distance to TSS or so and provide the list of transcripts and alltranscript set
# for instance:
# foo = pioGO(thisGenes=peaksUp.anno$GENEID)
pioGO = function(thisGenes, allGenes = NULL, txDb = NULL, genome='danRer10')
{
  require(goseq)
  #require(org.Dr.eg.db)
  require(GenomicFeatures)   # not sure if necessery
  
  #TxDb.Drerio.ucsc = loadDb("/mnt/biggles/csc_home/piotr/Work/GenomeAnnotations/TxDb/TxDb.Drerio.ucsc.sqlite")
  #TxDb.Drerio.ens = loadDb("/mnt/biggles/csc_home/piotr/Work/GenomeAnnotations/TxDb/TxDb.Drerio.ensembl.sqlite")
  if(is.null(allGenes))
    allGenes = genes(TxDb)$gene_id
  
  thisGenes = thisGenes[!is.na(thisGenes)]
  pwf = data.frame(row.names=allGenes[!is.na(allGenes)])
  pwf["DEgenes"] = 0L
  pwf["bias.data"] = 1L
  pwf["pwf"] = 1.0
  
  #thisTx = unique(sort(anno[anno$K27me3_128C==0,"geneId", drop=T]))
  pwf[thisGenes, "DEgenes"] = 1L
  aaa = goseq(pwf, genome, 'ensGene', method = "Wallenius", use_genes_without_cat = F)
  aaa$bh_adj_over_represented_pvalue = p.adjust(aaa$over_represented_pvalue, method="BH")
  aaa$bh_adj_under_represented_pvalue = p.adjust(aaa$under_represented_pvalue, method="BH")
  return(aaa)
}

pioAssignGo = function(thisGenes, genome = "danRer10", id="ensGene")
{
  require(goseq)
  getgo(thisGenes, genome, id)
}

pioFindGenesFromGO = function(thisGenes, cat)
{
  require(goseq)
  gos = pioAssignGo(thisGenes)
  ret = lapply(gos, FUN=
  function(g)
  {
    cat %in% g
  })
  ret = unlist(ret)
  ret = names(ret)[ret]
  return(ret)
}
