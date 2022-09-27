"""
This script preps a master table containing all data in an easily accesible
format.

First, the data is loaded. The genes are subset only to those genes that are
commonly identified in ALL datasets (pqm1, 50 and 58hrs, with and without
ascr#10).

Next, I add a few columns with useful information.

Next, I add the TPM information for 50 and 58hrs.

run this script from `python/prep_master_table.py`

Author: David Angeles-Albores
email: davidaalbores@gmail.com
"""

import pandas as pd
import numpy as np
import scipy

################################################################################
################################################################################
################################################################################
# ensembl data
ens = pd.read_csv('../data/ensembl.tsv', sep='\t')
# counts matrix and metadata:
mat = pd.read_csv('../data/matrix.csv')
# DESeq2 results:
n50 = pd.read_csv('../data/diff_exp/DE_N250.csv', index_col=0)
n58 = pd.read_csv('../data/diff_exp/DE_N258.csv', index_col=0)
# pqm1 = pd.read_csv('../data/diff_exp/DE_pqm1_50.csv', index_col=0)
# tph1_50 = pd.read_csv('../data/diff_exp/DE_tph1_50.csv', index_col=0)
# tph1_58 = pd.read_csv('../data/diff_exp/DE_tph1_58.csv', index_col=0)
# tapper datasets
df = pd.read_excel('../data/mmc1.xls')  # supplementary data from Tapper et al

treatment = {c: 'control' if 'cnt' in c else 'ascr' for c in mat.columns}
time = {c: '50' if '50' in c else '58' for c in mat.columns}

n50.columns = [c + '-50' for c in n50.columns]
n58.columns = [c + '-58' for c in n58.columns]
# pqm1.columns = [c + '-pqm1' for c in pqm1.columns]
# tph1_50.columns = [c + '-tph1-50' for c in tph1_50.columns]
# tph1_58.columns = [c + '-tph1-58' for c in tph1_58.columns]

# fix tapper columns:
df.columns = ['Rank', 'Transcript', 'Gene', 'Locus', 'ensembl_gene_id',
              'AvailableExperiments', 'ExpUp', 'ExpDown', 'VoteScore', 'tvalue',
              'pvalue', 'fdr', 'DBE', 'DAE', 'daf16', 'pqm1']

# restrict analyses to detected genes
common = np.intersect1d(n50.index, n58.index)
# common = np.intersect1d(pqm1.index, common)
# common = np.intersect1d(tph1_50.index, common)
# common = np.intersect1d(tph1_58.index, common)

# keep common genes
n50 = n50.reindex(common)
n58 = n58.reindex(common)
# pqm1 = pqm1.reindex(common)
# tph1_50 = tph1_50.reindex(common)
# tph1_58 = tph1_58.reindex(common)
ens = ens[ens.ensembl_gene_id.isin(common)]

# join datasets and add ensembl info:
res = n50.join(n58)
#res = res.join(pqm1)
# res = res.join(tph1_50)
# res = res.join(tph1_58)
res = res.join(ens.set_index('ensembl_gene_id'
              ).drop(columns='ensembl_transcript_id').drop_duplicates())
res = res.join(df.set_index('ensembl_gene_id'), rsuffix='_tapper')

# ignore mito genes:
res = res[res.chromosome_name != 'MtDNA']
res = res.dropna(subset=['padj-58', 'padj-50'])#, 'padj-pqm1'])
################################################################################
################################################################################
################################################################################

def sign_wt(x):
    if x['log2FoldChange-50'] * x['log2FoldChange-58'] < 0:
        return 'Different'
    else:
        return 'Same'

# def sign_pqm1(x):
#     if x['log2FoldChange-50'] * x['log2FoldChange-pqm1'] < 0:
#         return 'Different'
#     else:
#         return 'Same'
#
#
# def sign_tph1(x, time='50'):
#     if x['log2FoldChange-' + time] * x['log2FoldChange-tph1-' + time] < 0:
#         return 'Different'
#     else:
#         return 'Same'

def color(x):
    if x > 0:
        return 'Positive'
    else:
        return 'Negative'


def sig(x):
    if (x['padj-58'] < 0.05) & (x['padj-50'] < 0.05):
        return 'DE in both'
    elif (x['padj-58'] < 0.05):
        return 'DE at 58hrs'
    elif (x['padj-50'] < 0.05):
        return 'DE at 50hrs'
    else:
        return 'Not DE'


# def sig_pqm1(x):
#     if (x['padj-58'] < 0.05) & (x['padj-50'] < 0.05) & (x['padj-pqm1'] < 0.05):
#         return 'DE in all'
#     elif (x['padj-pqm1'] < 0.05):
#         return 'DE in pqm-1'
#     else:
#         return 'Not DE in pqm-1'

################################################################################
################################################################################
################################################################################

# useful columns
res.external_gene_name = res.external_gene_name.astype(str)
res['Size'] = (res.start_position - res.end_position).abs()

res['Sign-50'] = res['log2FoldChange-50'].map(color)
res['Sign-58'] = res['log2FoldChange-58'].map(color)
# res['Sign-pqm1'] = res['log2FoldChange-pqm1'].map(color)
# res['Sign-tph1-50'] = res['log2FoldChange-tph1-50'].map(color)
# res['Sign-tph1-58'] = res['log2FoldChange-tph1-58'].map(color)

cat_type = pd.CategoricalDtype(categories=['Same', 'Different'], ordered=True)
res['Sign-WT'] = res.apply(sign_wt, axis=1).astype(cat_type)
# res['Sign-pqm1'] = res.apply(sign_pqm1, axis=1).astype(cat_type)
# res['Sign-tph1-50'] = res.apply(sign_tph1, axis=1).astype(cat_type)
# res['Sign-tph1-58'] = res.apply(sign_tph1, args=('58',), axis=1).astype(cat_type)

cat_type = pd.CategoricalDtype(categories=['Not DE', 'DE at 50hrs',
                                           'DE at 58hrs', 'DE in both'],
                               ordered=True)
res['Significance-WT'] = res.apply(sig, axis=1).astype(cat_type)

# cat_type = pd.CategoricalDtype(categories=['Not DE in pqm-1', 'DE in pqm-1',
                                           # 'DE in all'],
                               # ordered=True)
# res['Significance-pqm1'] = res.apply(sig_pqm1, axis=1).astype(cat_type)

res['Ratio'] = (res['log2FoldChange-58']) / res['log2FoldChange-50']

# order your chromosomes as categories, kids!
cat_type = pd.CategoricalDtype(categories=['I', 'II', 'III', 'IV', 'V', 'X'],
                               ordered=True)
res['chromosome_name'] = res.chromosome_name.astype(cat_type)

################################################################################
################################################################################
################################################################################

# count matrix formatting:
c = [c for c in mat.columns if 'tph1' not in c]
mat = mat[c]
mat = mat[sorted(mat.columns)[-6:-3] +
          sorted(mat.columns)[0:3] +
          sorted(mat.columns)[-3:] +
          sorted(mat.columns)[3:6]]
mat = mat / mat.sum(axis=0) * 10 ** 6  # transform to CPM

# melt the matrix into tidy format
melted = mat.reset_index().melt(id_vars='index', var_name='Sample',
                                value_name='counts')
melted.rename({'index': 'ensembl_gene_id'}, axis=1, inplace=True)
melted['treatment'] = melted.Sample.map(treatment)
melted['time'] = melted.Sample.map(time).astype('str')

# add a column for minimum counts in any one of the samples
agged = melted.groupby(['ensembl_gene_id']).agg({'counts': np.min}).reset_index()
agged = agged[agged.ensembl_gene_id.isin(res.index)]

res = res.join(agged.rename(columns={'counts': 'MinCountsDetected'}
              ).set_index('ensembl_gene_id'))

res.columns = [c.replace('_', '').replace('chromosomename', 'chromosome')
               for c in res.columns]

logBM = lambda x: np.log10(x + 10 ** -6)
res['logBM-50'] = res['baseMean-50'].apply(logBM)
res['logBM-58'] = res['baseMean-58'].apply(logBM)
# res['logBM-pqm1'] = res['baseMean-pqm1'].apply(logBM)

res.to_csv('../data/master_table.tsv', sep='\t')
