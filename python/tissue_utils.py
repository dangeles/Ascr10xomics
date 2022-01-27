import pandas as pd
import numpy as np
import scipy
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import fdrcorrection

from matplotlib import rc
rc('text', usetex=False)
# rc('text.latex', preamble=r'\usepackage{cmbright}')
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
    data.neglogq = data.neglogq.replace(np.inf, np.nanmax(data.neglogq) + 1)
    return data

def test_tissue_direction(sig, tissues, col='log2FoldChange-58'):
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
        found = sig.reindex(g.wbid.unique()).dropna(subset=[col])
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


def similarity_trimming(considered_tissues, tissues, remove=[]):
    """
    Checks terms pairwise for similarity of annotations, and in the event that
    the similarity is too great, adds the tissue with fewer annotations to a
    list.

    Params:
    -------
    considered_tissues:    list. Tissues to compare among.
    tissues:    pd.DataFrame. Dataframe of terms and annotations
    remove:     list. List of tissues to remove.

    Output:
    remove:     list. List of tissues that are highly similar and should be
                      removed.
    """
    for i, t in enumerate(considered_tissues):
        for t2 in considered_tissues[i + 1:]:
            if t == t2:
                continue
            current = tissues.tissue.isin([t, t2])
            counted = tissues[current].groupby('wbid').tissue.count()
            both = (counted == 2).sum()
            min_size = tissues[current].groupby('tissue').wbid.count().min()
            min_tissue = tissues[current].groupby('tissue').wbid.count().idxmin()

            if (both / min_size > 0.75) & (min_tissue.lower() != 'pharynx'):
                remove += [min_tissue]
    return remove


def pretty_GSEA_plots(binom_data, fc_data, alpha=0.05, x='neglogq',
                      y='tissue', hue='data', size=10,
                      palette={'50hrs': 'tab:blue', '58hrs': 'tab:green'}):
    """
    GSEA Plots
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    # generate stripplots of q-values and fraction of positive log2FC
    sns.stripplot(x=x, y=y, data=binom_data, jitter=False,
                  size=size, hue=hue,
                  palette=palette, ax=ax)

    # significance threshold:
    ax.axvline(-np.log10(alpha), ls='--', color='red',
               label='alpha threshold')

    return fig, ax


def corr_plots(data, x, y, hue, size, tissue, qval=0.05, sig_col='padj-50',
               kind='pearson', **kwargs):
    """
    A wrapper to plot scatter plots or rank scatter plots.

    Params:
    -------
    data     pd.DataFrame. Data
    x        str. Column containing x-coordinates
    y        str. Column containing y-coordinates
    hue      str. Column that specifies point color
    size     str. Column that specifies point size.
    tissue:  str. Tissue genes should be expressed in
    qval     float. significance threshold
    kind:    str. If `rank`, then plots a rank-plot
    """
    cond = (data.tissue.str.contains(tissue)) & (data[sig_col] < qval)
    tmp = data[cond].copy()

    if kind == 'rank':
        tmp['rank-' + x] = tmp[x].rank()
        tmp['rank-' + y] = tmp[y].rank()
        x = 'rank-' + x
        y = 'rank-' + y

    sns.relplot(data=tmp, x=x, y=y, size=size, hue=hue, kind="scatter",
                **kwargs)
    plt.plot((tmp[x].min(), tmp[x].max()), (tmp[x].min(), tmp[x].max()),
             color='black', lw=3)

    if kind != 'rank':
        plt.axvline(0, ls='--', lw=1, color='black')
        plt.axhline(0, ls='--', lw=1, color='black')

    plt.title(tissue)
