import json
import pandas as pd
import numpy as np
import scipy
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import fdrcorrection

from matplotlib import rc
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{cmbright}')
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

rc = {'lines.linewidth': 2,
      'axes.labelsize': 20,
      'axes.titlesize': 20,
      'axes.facecolor': 'DFDFE5'}
sns.set_context('notebook', rc=rc)
sns.set_style("ticks")

mpl.rcParams['xtick.labelsize'] = 23
mpl.rcParams['ytick.labelsize'] = 23
mpl.rcParams['legend.fontsize'] = 20



def load_tissues(bckgrnd, n_annotations, n_tissues,
                 path='../data/tissue_annotations.tsv' ):
    """
    A function to load C. elegans tissues and annotations. The function loads
    all the annotations, then removes all genes not found in the provided
    background list. Tissues are renamed so that '-' are replaced by '_'.

    Tissues are further subset by removing any tissues with < some threshold of
    annotated genes; and genes are removed if they appear in more than some
    threshold of tissues.

    Params:
    -------
    bckgrnd:   list-like. Background list of genes
    n_annotations:  int. Minimum number of genes a tissue must have to remain
                         in the list
    n_tissues:      int. Maximum number of tissues a gene can be annotated to
                    before it gets removed.

    Output:
    tissues:    pd.DataFrame. Trimmed tissue-annotation dataframe.
    """
    tissues = pd.read_csv(path, sep='\t', header=None,
                          names=['wbid', 'external_gene_name', 'target_id',
                                 'species', 'tissue'])
    tissues = tissues[tissues.wbid.isin(bckgrnd)]
    tissues.tissue = tissues.tissue.str.replace('_', '-')

    print('tissues originally:', tissues.tissue.nunique())
    # remove tissues that have too few annotations:
    counts = tissues.groupby('tissue').wbid.count()
    tissues = tissues[tissues.tissue.isin(counts.index[counts > n_annotations])]
    # remove genes that appear everywhere:
    promiscuous = tissues.groupby('wbid').tissue.nunique()
    keep = promiscuous[promiscuous < n_tissues].index
    tissues = tissues[tissues.wbid.isin(keep)]
    print('tissues afterwards:', tissues.tissue.nunique())

    return tissues


def fdr_correct(data, alpha=0.05, method='indep', col='pval'):
    """A wrapper to perform FDR correction and append results to a dataframe"""
    fdr = fdrcorrection(data[col], alpha=alpha, method=method)
    data['fdr'] = fdr[1]
    data['neglogq'] = -data.fdr.apply(np.log10)
    data['sig'] = fdr[0]
    data.neglogq = data.neglogq.replace(np.inf, 350)

    return data

def test_tissue_direction(sig, tissues, col='log2FoldChange58'):
    """
    A function to test the probability that so many genes are changing in the
    same direction.

    Params:
    -------
    sig:   pd.DataFrame. Data, containing only significantly altered genes.
    tissues:    pd.DataFrame. Dataframe containing tissue-gene annotations.
    col:    str. Name of the column where the log2FC are kept in `sig`

    Output:
    -------
    data:    pd.DataFrame. Results of statistical tests.
    """
    # parameters:
    size = sig.shape[0]
    p_pos = (sig[col] > 0).sum() / size
    p_neg = 1 - p_pos

    # count # of genes annotated to tissues:
    tissues = tissues[tissues.wbid.isin(sig.index)]
    data = []

    # go through each tissue
    for n, g in tissues.groupby('tissue'):
        # find genes annotated to current tissue:
        found = sig.reindex(g.wbid.unique()).dropna()
        # count the positive genes:
        pos = (found[col] > 0).sum()

        # only test if you found at least 5 genes
        if len(found) < 5:
            data += [[n, np.nan,  pos / g.wbid.nunique()]]
            continue

        # test two-sided p-value that positive skew is this large or larger,
        # given average skew
        pval = scipy.stats.binom_test(pos, len(found), p_pos)
        data += [[n, pval, pos / len(found)]]

    # store in dataframe:
    data = pd.DataFrame(data, columns=['tissue', 'pval', 'FracPos']
                       ).sort_values('pval', ascending=True)
    data = data.dropna()
    data['FracPosExpected'] = p_pos
    return data


def pretty_GSEA_plots(binom_data, fc_data, alpha=0.05, x='neglogq', y='tissue',
                      hue='data', size=10,
                      palette={'50': 'tab:blue', '58': 'tab:green'}):
    """
    GSEA Plots
    """

    fig, ax = plt.subplots(figsize=(25, 12), ncols=3, sharey = True)
    # generate stripplots of q-values and fraction of positive log2FC
    sns.stripplot(x=x, y=y, data=binom_data, jitter=False,
                  size=size, hue=hue,
                  palette=palette, ax=ax[0])

    sns.stripplot(x='FracPos', y=y, data=binom_data, jitter=False,
                  size=size, hue=hue,
                  palette=palette, ax=ax[1])

    # plot log2FC in the same order as the stripplots above:
    tissue_order = pd.CategoricalDtype(categories=binom_data.tissue.unique(),
                                       ordered=True)
    tmp = fc_data[fc_data.tissue.isin(binom_data.tissue)].copy()
    tmp.tissue = tmp.tissue.astype(tissue_order)
    tmp.sort_values('tissue', inplace=True)

    # plot points:
    sns.stripplot(x='log2FoldChange', y='tissue', size=6, hue=hue, data=tmp,
                  palette=palette, ax=ax[2], dodge=True)

    # significance threshold:
    ax[0].axvline(-np.log10(alpha), ls='--', color='red',
                  label='alpha threshold')
    # set a horizontal line at 0 for the log2FC:
    ax[2].axvline(0, ls='-', color='black')

    # add vertical lines at Mean Positive Fraction (null expectation)
    if hue is not None:
        for n, g in binom_data.groupby(hue):
            ax[1].axvline(g.FracPosExpected.unique()[0],
                  ls='--', label='Expected Positive Fraction for ' + n,
                  color=palette[n])

    # remove legends, labels, gray:
    for ai in ax:
        ai.set_ylabel('')
        ai.legend([])
        ai.yaxis.grid(color='gray', linewidth=.5)

    # add one legend back in
    legend = ax[1].legend(loc=(0.05, 0.05), fontsize=20)
    legend.get_title().set_fontsize('14') #legend 'Title' fontsize

    return fig, ax, legend
