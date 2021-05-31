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
source('utils.R')

# get tph1 + n2 data
rds = readRDS('../nxt/salmon/salmon_merged_gene_counts.rds')
rds$Names = rds$names %>% gsub('\\.[0-9].*', '', .)

# annotate metadata:
colData(rds)['Treatment'] = sapply(rds$Names, treatment_getter) %>% as.factor
colData(rds)['Genotype'] = sapply(rds$Names, genotype_getter) %>% as.factor
colData(rds)['Age'] = sapply(rds$Names, time_getter) %>% as.factor

# put into DESeq
dds50 = loadDDS(rds[,rds$Age == '50'], which='wt', design=~Treatment)
res50 = get_results(dds50, 'Treatment', 'ascr', 'cnt', filter=FALSE)

# at time 58
dds58 = loadDDS(rds[,rds$Age == '58'], which='wt', design=~Treatment)
res58 = get_results(dds58, 'Treatment', 'ascr', 'cnt', filter=FALSE)

# compare gene expression thru time
ddsTime = loadDDS(rds[,rds$Genotype == 'wt'], which='wt', design=~Age)
resTime = get_results(ddsTime, 'Age', '58', '50', filter=FALSE)


res50 %>% as.data.frame %>% write.csv('../data/diff_exp/DE_N250.csv', quote = FALSE)
res58 %>% as.data.frame %>% write.csv('../data/diff_exp/DE_N258.csv', quote = FALSE)
resTime %>% as.data.frame %>% write.csv('../data/diff_exp/DE_Time.csv', quote = FALSE)