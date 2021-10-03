import pandas as pd
import numpy as np
import scipy
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rc
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{cmbright}')
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})


rc = {'lines.linewidth': 2,
      'axes.labelsize': 18,
      'axes.titlesize': 18,
      'axes.facecolor': 'DFDFE5'}
sns.set_context('notebook', rc=rc)
sns.set_style("dark")

mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['legend.fontsize'] = 14


def calculate_pearson(dataframe, x, pval, shuffle=False, window_size=5):
    """
    Calculates a running pearson correlation arranged by genomic position.

    Params:
    -------
    dataframe:    pd.DataFrame. Dataset
    x:            str. One of `50`, `58` or `pqm1`
    shuffle:      Bool. If true, shuffles the dataset to try to
                        generate a null distribution
    window_size:  int. Number of genes to use for correlation calculation.

    Output:
    -------
    corr:         pd.DataFrame. A dataframe with the pearson correlation values
                                at each genomic position.
    """
    df = dataframe.copy()
    df.sort_values(['chromosome', 'startposition'], inplace=True)

    # select a subset of genes
    if pval is None:
        subset = df['padj-' + x] < 1
    else:
        subset = df['padj-' + x] < pval

    corrs = []
    avg_poss = []
    chroms = []
    for n, g in df[subset].groupby('chromosome'):
        if shuffle is True:
            g = g.sample(frac=1)

        # calculate running correlation:
        correlation = g.startposition.rolling(window_size).corr(g['log2FoldChange-' + x])
        avg_pos = g.startposition.rolling(window_size).mean()

        # store values:
        corrs += correlation.values.tolist()
        avg_poss += avg_pos.values.tolist()
        chroms += [n] * len(correlation)

    # place into dataframe:
    corr = pd.DataFrame([chroms, avg_poss, corrs]).T
    corr.columns = ['chromosome', 'pos', 'pearson']
    corr.pearson = corr.pearson.astype(np.float64)
    corr.pos = corr.pos.astype(np.float64)
    corr = corr.dropna()
    corr.set_index('pos', inplace=True)
    return corr


def get_corr_and_shuffle(res, x, pval, window=6):
    """
    Calculates running window correlation + null distribution.

    Params:
    -------
    x:      str. One of `50`, `58` or `pqm1`
    pval:   float. Significance threshold.
    window: int. Defaults to 6. Number of genes to pool for correlation.

    Output:
    -------
    corr:   pd.DataFrame. Dataframe with running pearson correlations.
    shuffled: pd.DataFrame. Dataframe with shuffled results (null distribution)
    """
    corr = calculate_pearson(res, x, pval, shuffle=False, window_size=window)

    # shuffled
    reps = 10 ** 3
    randoms = []
    for i in range(reps):
        tmp = calculate_pearson(res, x, pval, shuffle=True, window_size=window)
        tmp['rep'] = [i] * len(tmp)
        tmp.index = corr.index
        randoms +=  [tmp]
    shuffled = pd.concat(randoms)
    return corr, shuffled


def plot_corrs(corr, shuffled):
    """
    Given a dataframe `corr` containing sliding windows of pearson correlations,
    and a `shuffled` dataframe of correlations with N replicates, plot the
    correlations.
    """
    fig, ax = plt.subplots(figsize=(42, 16), ncols=3, nrows=2, sharey=True,
                           sharex=True)

    i = 0
    for n, g in corr.groupby('chromosome'):
        # axis coordinates:
        y = i % 3
        x = int((i - y) / 3)

        # 5 and 95 percentiles of random distribution:
        rand5 = shuffled[shuffled.chromosome == n].dropna().groupby('pos').pearson.quantile(0.05).mean()
        rand95 = shuffled[shuffled.chromosome == n].dropna().groupby('pos').pearson.quantile(0.95).mean()

        # select correlation points above the 95 percentile
        sel = (g.pearson > rand95)
        # count number of peaks selected:
        n_corr_peaks = sel.sum()

        #plot lines
        ax[x, y].plot(g.pearson)
        # plot points above 95%
        ax[x, y].scatter(g[sel].index, g[sel].pearson, color='red', s=90, zorder=np.inf)
        # color in the null interval:
        ax[x, y].fill_between(g.index, rand5, rand95, color='gray', alpha=0.5)
        ax[x, y].set_title(n, fontsize=35)
        ax[x, y].annotate('No. Corr Peaks: {0}'.format(n_corr_peaks),
                          (1 * 10 ** 7, .8), fontsize=30)
        ax[x, y].tick_params(axis="x", labelsize=30)
        ax[x, y].tick_params(axis="y", labelsize=30)
        ax[x, y].yaxis.get_offset_text().set_fontsize(26)

        i += 1

    ax[0, 0].set_ylabel('Correlation', fontsize=30)
    ax[1, 0].set_ylabel('Correlation', fontsize=30)
    ax[1, 1].set_xlabel('Genomic Position', fontsize=30)

    return fig, ax


def calculate_enrichment(selection=None, alpha=0.05, suppress=False, test='pqm1'):
    if selection is None:
        selection = (res['padj-50'] < 0.05) & (res['padj-58'] < 0.05)

    selection = selection & (res['Sign-WT'] == 'Same')

    group = test
    classII = [("Y")]
    null = res.groupby(group).Locus.count()
    obs = res[selection].groupby(group).Locus.count()

    fold = (obs / null * null.sum() / obs.sum())
    results = pd.DataFrame([null, obs, fold],
                           index=['Null', 'Observed', 'Fold Change'])
    results = results.T
    results.Null = results.Null.astype(int)
    results.Observed = results.Observed.fillna(0).astype(int)

    if results.Observed.sum() < 3:
        return None, None, [np.inf] * 6

    # hypergeometric test:
    Total = null.sum()
    ClassII = null.reindex(classII).values[0]
    draws = obs.sum()
    ObsClassII = obs.reindex(classII).values[0]
    pval = scipy.stats.hypergeom.sf(ObsClassII, Total, ClassII, draws)
    antipval = scipy.stats.hypergeom.cdf(ObsClassII, Total, ClassII, draws)

    # print message:
    message = 'P-value associated with observing {0} class `{1}` genes out of {2} possible class 2 genes in this dataset, given {3} draws: {4:.2g}'
    antimessage = 'P-value associated with observing {0} class `{1}` genes out of {2} possible class 2 genes in this dataset, given {3} draws: {6:.2g}'
    return message, antimessage, (ObsClassII, classII[0], ClassII, draws, pval, Total, antipval)
