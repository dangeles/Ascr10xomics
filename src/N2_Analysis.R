library(dplyr)
library(tidyverse)
library(DESeq2)
library(tximport)
#library(biomaRt)
library(ggplot2)
library(SummarizedExperiment)
library(stringr)
library(ggfortify)
library(data.table)
source('utils.R')

# get tph1 + n2 data
rds = readRDS('../nxt/salmon/salmon.merged.gene_counts.rds')
rds$Names = rds$names %>% gsub('\\.[0-9].*', '', .)

# annotate metadata:
colData(rds)['Treatment'] = sapply(rds$names, treatment_getter) %>% as.factor
colData(rds)['Genotype'] = sapply(rds$names, genotype_getter) %>% as.factor
colData(rds)['Age'] = sapply(rds$names, time_getter) %>% as.factor
colData(rds)['MultiLevel'] = rds$Treatment + '-' + colData(rds)$Age

# put into DESeq
dds50 = loadDDS(rds, which='pqm1', design=~Treatment)
res50 = get_results(dds50, 'Treatment', 'ascr', 'cnt', filter=FALSE)
s50 = lfcShrink(dds50, 'Treatment_ascr_vs_cnt')


res50 %>% as.data.frame %>% write.csv('../data/diff_exp/DE_pqm1_50.csv', quote = FALSE)
s50 %>% as.data.frame %>% write.csv('../data/diff_exp/DE_pqm1_50_shrunken.csv', quote = FALSE)