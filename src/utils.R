library(dplyr)
library(tidyverse)
library(DESeq2)
library(tximport)
library(biomaRt)
library(ggplot2)
library(SummarizedExperiment)
library(stringr)
library(ggfortify)
library(data.table)
library(VennDiagram)


# metadata extraction:
treatment_getter = function(x){
  if(str_detect(x[[1]], 'cnt')){
    return('cnt')
  } else{
    return('ascr')
  }
}

time_getter = function(x){
  if(str_detect(x[[1]], 'Bef') | str_detect(x[[1]], '50')){
    return('50')
  } else{
    return('58')
  }
}

genotype_getter = function(x){
  if(str_detect(x[[1]], 'tph')){
    return('tph1')
  } else{
    return('wt')
  }
}

# DESeq wrappers:
loadDDS = function(rds, which, design, column='Genotype'){
  ensembl=useMart("ensembl")
  worm = useDataset("celegans_gene_ensembl",mart=ensembl)
  wormProteinCoding = getBM(attributes=c("ensembl_gene_id", "external_gene_name",
                                         "description"), filters='biotype',
                            values=c('protein_coding'), mart=worm,
                            useCache = FALSE)
  common = names(rds)[names(rds) %in% wormProteinCoding$ensembl_gene_id]
  
  if(which != ''){
    rds = rds[common, rds[[column]] == which]
  } else{
    rds = rds[common,]
  }
  assay = round(as.matrix(assay(rds)))
  se = SummarizedExperiment(list(raw=assay),
                            colData=colData(rds),
                            rowData=rowData(rds))

  dds = DESeqDataSet(se, design=design)
  
  if(which == 'both'){
    dds$Genotype <- relevel(dds$Genotype, ref = "wt")
  }
  dds$Treatment <- relevel(dds$Treatment, ref = "cnt")
  dds$Age <- relevel(dds$Age, ref = "50")
  
  keep <- rowSums(counts(dds)) >= 20
  dds <- dds[keep,]
  dds = DESeq(dds)
  
  return(dds)
}

get_results = function(dds, factor, numerator, denominator,
                       name='', contrast=TRUE, filter=TRUE){
  alpha = 0.01
  if(contrast){
    res = results(dds, alpha=alpha, contrast = c(factor,
                                                 paste(numerator, sep=''),
                                                 paste(denominator, sep='')))
  } else{
    res = results(dds, alpha=alpha, name=name)
  }
  resOrdered <- res[order(res$pvalue),]
  resOrdered['gene_name'] = rowData(dds)[rownames(resOrdered),]$gene_name
#  if(sum(is.na(resOrdered$padj)) > 0){
#    resOrdered[is.na(resOrdered$padj),]$padj = 1
#  }
  if(filter == TRUE){
    resOrdered = resOrdered[resOrdered$padj < alpha,]
  }
 resOrdered['logq'] = -log10(resOrdered$padj)
  return(resOrdered)
}

# plotting functions:
plotSubsetPCA = function(vst, dds, cond, intgroup=c('Genotype', 'Treatment'),
                         title){
  plotPCA(vst[,rownames(colData(dds)[cond,])],
          intgroup=intgroup) + ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
}

cdf = function(res, title, save=FALSE, path=''){
  p = ggplot(as.data.frame(res), aes(x=-logq)) +
      geom_histogram(aes(y = cumsum(..count..)), bins=100) +
      scale_y_continuous(trans='log10') +
      geom_vline(xintercept=log10(0.01), color='red', size=2) +
      xlim(-45, 0) + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

scatter = function(which, a, b, alpha, X, Y, rank=FALSE){
  a = a[order(rownames(a)),]
  b = b[order(rownames(b)),]
  
  apvals = rownames(a[a$padj < 0.01,])
  bpvals = rownames(b[b$padj < 0.01,])
  rpvals = apvals[(apvals %in% bpvals)]
  
  tmp = cbind(a[rpvals,]$log2FoldChange,
              b[rpvals,]$log2FoldChange) %>%
    as.data.frame %>%  dplyr::rename(X = 'V1', Y = 'V2')
  rownames(tmp) = rpvals
  if(rank == TRUE){
    p = ggplot(tmp, aes(rank(X), rank(Y)),  ggplot(tmp, aes(X, Y))) +
      geom_point(alpha=0.3, size=2) +
      geom_abline(intercept=0, slope=1) + xlab(paste('rank', X, sep=' ')) +
      ylab(paste('rank', Y, sep=' '))
  } else{
    p = ggplot(tmp, aes(X, Y)) + geom_point(alpha=0.3, size=2) +
      geom_abline(intercept=0, slope=1) + xlab(X) + ylab(Y)
  }
  return(p)
}

venn = function(r50, r58, genotype, alpha=0.01){
  r50 = r50[order(rownames(r50)),]
  r58 = r58[order(rownames(r58)),]
  
  r50name = rownames(r50[r50$padj < alpha,])
  r58name = rownames(r58[r58$padj < alpha,])
  
  venn.diagram(
    x = list(r50name, r58name),
    category.names = c("50h" , "58h"), filename=paste('../output/venn/',
                                                      paste(genotype, 'png',
                                                            sep='.'),
                                                      sep='/'),
    output=TRUE
  ) 
}

