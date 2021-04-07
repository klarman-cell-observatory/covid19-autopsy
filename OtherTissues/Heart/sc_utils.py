# Utility functions for single cell RNA-seq analysis in scanpy.

import pandas as pd
import numpy as np

import anndata
import scanpy as sc

from scipy.io import mmwrite

import matplotlib.pyplot as plt

from typing import Tuple, Union, List, Dict, Callable, Iterable, Optional
from functools import partial
import os
import gzip
import shutil
import inspect
import string
import gc


def get_marker_gene_summary_table(adata: anndata.AnnData,
                                  marker_gene_key: str = 'rank_genes_groups1.0',
                                  expression_cutoffs: List[int] = [0]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Create an excel-ready dataframe of cluster marker genes and their relevant
    statistics including AUC and % cells in cluster expressing above a given level
    and % cells outside the cluster expressing above a given level. This grabs
    genes from a scanpy test... will calculate AUC if the test was Wilcoxon.

    Args:
        adata: The scanpy adata object
        marker_gene_key: The marker gene info must be stored in adata.uns[marker_gene_key]
        expression_cutoffs: List of expression cutoffs to use to determine "expressing"
            versus "non-expressing" cells in the percentage calculation.

    Returns:
        df: Dataframe with all the information
        df_stacked: Long format dataframe with all genes in column 0 and cluster labels in column 1

    """

    assert marker_gene_key in adata.uns.keys(), f'{marker_gene_key} is not in adata.uns'
    assert 'names' in adata.uns[marker_gene_key].keys(), \
        f'adata.uns["{marker_gene_key}"] must have the key "names"'
    cluster_key = adata.uns[marker_gene_key]['params']['groupby']
    marker_method = adata.uns[marker_gene_key]['params']['method']
    assert cluster_key in adata.obs.keys(), \
        f'adata.uns["{marker_gene_key}"] was grouped by {cluster_key}, but {cluster_key} is not in adata.obs'

    frames = []

    # get a dataframe of gene names
    gene_df = pd.DataFrame(np.array(adata.uns[marker_gene_key]['names'].tolist(), dtype=str),
                           columns=adata.uns[marker_gene_key]['names'].dtype.names)

    frames.append(gene_df)

    # get a dataframe of auc
    if marker_method == 'wilcoxon':

        assert 'scores' in adata.uns[marker_gene_key].keys(), \
            f'Ran a Wilcoxon test, but "scores" key in adata.uns["{marker_gene_key}"] is missing.'

        auc_df = pd.DataFrame(np.array(adata.uns[marker_gene_key]['scores'].tolist(), dtype=float),
                              columns=adata.uns[marker_gene_key]['scores'].dtype.names)

        for k in auc_df.columns:
            # k is cluster label in louvain

            # conversion courtesy of Mark, modify auc_df in place
            n1 = adata.obs[cluster_key].isin([str(k)]).sum()
            n2 = adata.X.shape[0] - n1
            auc_df[k] = (auc_df[k] * np.sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)) + ((n1 * n2) / 2)  # Wilcoxon U
            auc_df[k] = auc_df[k] / (n1 * n2)  # AUC

        auc_df = auc_df.rename(columns=dict(zip(auc_df.columns,
                                                [i + ': AUC' for i in auc_df.columns])))

        # remove nan genes that were filtered
        auc_df[(gene_df[:] == 'nan').values] = None

        frames.append(auc_df)

    else:

        print(f'AUC data not included, since the key "scores" is not in adata.uns["{marker_gene_key}"]. '
              f'  Did you run a Wilcoxon test?')

    # remove nan genes that were filtered
    gene_df[gene_df[:] == 'nan'] = None

    for expression_cutoff in expression_cutoffs:

        # set up dataframes for percentage calculations
        in_df = gene_df.copy()
        in_df = in_df.rename(columns=dict(zip(in_df.columns,
                                              [i + f': pct. expressing > {expression_cutoff}'
                                               for i in gene_df.columns])))
        out_df = gene_df.copy()
        out_df = out_df.rename(columns=dict(zip(out_df.columns,
                                                ['Not ' + i + f': pct. expressing > {expression_cutoff}'
                                                 for i in gene_df.columns])))

        print(f'Expression cutoff = {expression_cutoff}')
        print('Working on cluster ', end='')

        # NOTE: it is important that the cluster order in the gene test dataframe
        # is in order by population of cluster (decreasing)
        for k, gene_df_col, in_df_col, out_df_col in zip(adata.obs[cluster_key].value_counts().index,
                                                         gene_df.columns,
                                                         in_df.columns,
                                                         out_df.columns):
            # taking some code from scanpy.tl.filter_rank_genes_groups to make code faster

            var_names = pd.DataFrame(adata.uns[marker_gene_key]['names'])[k].values  # gene names

            adata.obs['__is_in_cluster__'] = pd.Categorical(adata.obs[cluster_key] == k)

            # obs_tidy has rows='__is_in_cluster__', columns=var_names
            matrix = adata[:, var_names].X.toarray()
            obs_tidy = pd.DataFrame(matrix, columns=var_names)
            obs_tidy.set_index(adata.obs['__is_in_cluster__'], '__is_in_cluster__', inplace=True)

            # compute fraction of cells having value > expression_cutoff
            # transform obs_tidy into boolean matrix
            obs_bool = (obs_tidy > expression_cutoff).astype(bool)

            # compute the fraction expressing per group
            fraction_obs = obs_bool.groupby(level=0).sum() / obs_bool.groupby(level=0).count()

            # store in dataframes
            in_df[in_df_col] = fraction_obs.loc[True].values * 100.  # percent
            out_df[out_df_col] = fraction_obs.loc[False].values * 100.

            print(str(k) + ', ', end='')

        print('done.')

        frames.append(in_df)
        frames.append(out_df)

    # concatenate columns
    df = pd.concat(frames, axis=1)

    # interleave columns appropriately
    df = df[df.columns[[i + j * gene_df.shape[1]
                        for i in range(gene_df.shape[1])
                        for j in range(len(frames))]]]

    # create a stacked dataframe with a column for cluster
    try:
        tmp = []
        for cluster_colname in gene_df.columns:

            # subset to this cluster's columns
            sub = df[df.columns[[(cluster_colname == d
                                  or d.startswith(cluster_colname + ': ')
                                  or d.startswith('Not ' + cluster_colname + ': '))
                                 for d in df.columns]].values]

            # add a "cluster" label in column 1
            sub.insert(1, 'cluster', cluster_colname)

            # rename columns and append
            new_col_names = ['Gene', 'Cluster', 'AUC']
            for e in expression_cutoffs:
                new_col_names.append('pct.cells.expr>' + str(e) + '.target')
                new_col_names.append('pct.cells.expr>' + str(e) + '.other')

            sub.columns = new_col_names
            tmp.append(sub)

        # concatenate it together
        df_stacked = pd.concat(tmp)

    except Exception:
        print('Failed to make stacked dataframe.')
        df_stacked = None

    return df, df_stacked


def cell_count_dotplot(adata: anndata.AnnData,
                       key: str,
                       cluster_key: str = 'louvain',
                       normalize: str = 'column',
                       row_normalization_factors: Optional[List[float]] = None,
                       show: bool = True,
                       xrotation: float = 0,
                       yrotation: float = 0,
                       size_factor: float = 1.,
                       color: Union[str, None] = None,
                       specified_key_order: Optional[List[str]] = None,
                       axis_labels: bool = True,
                       axis_pad: float = 0.8,
                       dot_size_legend_sizes: Optional[List[int]] = [0.05, 0.4, 1.0],
                       dot_size_legend_title: Optional[str] = None,
                       dot_size_legend_labelspacing: Optional[float] = None,
                       figsize: Optional[List[float]] = None):
    """Makes a dotplot that reflects the number of cells in a given condition.
    The condition could be tissue and cluster, for example.  One part of the
    two-part condition is meant to be cluster.  The dot sizes can be normalized
    per row or per column or both.
    Args:
        adata: AnnData object
        key: Group by this adata.obs[key] to tabulate cells
        cluster_key: Clusters are in adata.obs[cluster_key]
        normalize: How to normalize the dot sizes {'row', 'column', 'row_then_column'}
        row_normalization_factors: If normalization involves rows, you can input
            normalization factors here.  Must be in plotted order, bottom to top. This
            is useful if doing "row_then_column" normalization on a subset of data.
        show: True to show plot
        xrotation: Rotation of x-axis (cluster) labels
        yrotation: Rotation of y-axis (groupby) labels
        size_factor: Scale the dot sizes by a multiplicative factor
        color: If None, uses colors in adata.uns[cluster_key + '_colors'], otherwise
            uses specified color for all dots
        specified_key_order: Control the order of the groupby y-axis features
        axis_labels: True to label axes
        axis_pad: Padding between axes and dots
        dot_size_legend_sizes: If not None, create a dot size legend using markers
            with these sizes
        dot_size_legend_title: Title for dot size legend, if creating one
        dot_size_legend_labelspacing: Labelspacing param for pyplot.legend()
        figsize: Override default figure size
    """

    assert key in adata.obs.keys(), f'Input key {key} is not in adata.obs'
    assert cluster_key in adata.obs.keys(), f'Input cluster_key {cluster_key} is not in adata.obs'
    assert normalize in {'row', 'column', 'row_then_column'}, \
        f'normalize must be in ["row", "column", "row_then_column"] but was {normalize}'

    counts_tissue_cluster_df = pd.crosstab(adata.obs[cluster_key], adata.obs[key])

    if specified_key_order is not None:
        for item in specified_key_order:
            assert item in counts_tissue_cluster_df.columns, \
                f'Tried to re-order key values, but "{item}" is not in adata.obs["{key}"]'
        counts_tissue_cluster_df = counts_tissue_cluster_df[specified_key_order]

    if normalize == 'row':
        if row_normalization_factors is None:
            counts_tissue_cluster_df = counts_tissue_cluster_df.div(
                counts_tissue_cluster_df.sum(axis=0), axis=1)  # row normalize
        else:
            counts_tissue_cluster_df = counts_tissue_cluster_df.div(
                row_normalization_factors, axis=1)
    elif normalize == 'column':
        counts_tissue_cluster_df = counts_tissue_cluster_df.div(
            counts_tissue_cluster_df.sum(axis=1), axis=0)  # column normalize
    elif normalize == 'row_then_column':
        if row_normalization_factors is None:
            counts_tissue_cluster_df = counts_tissue_cluster_df.div(
                counts_tissue_cluster_df.sum(axis=0), axis=1)  # row normalize
        else:
            counts_tissue_cluster_df = counts_tissue_cluster_df.div(
                row_normalization_factors, axis=1)
        counts_tissue_cluster_df = counts_tissue_cluster_df.div(
            counts_tissue_cluster_df.sum(axis=1), axis=0)  # column normalize
    else:
        raise ValueError(f'normalize must be in ["row", "column", "row_then_column"] but was {normalize}')

    scatter_df = counts_tissue_cluster_df.stack().rename_axis(['y', 'x']).reset_index(name='val')

    counts = pd.crosstab(adata.obs[cluster_key],
                         adata.obs[key]).sum(axis=1).values.tolist()

    y = scatter_df['x'].values
    x = scatter_df['y'].values
    s = scatter_df['val'].values * 500 * size_factor

    xvals = scatter_df['y'].cat.codes
    yvals = scatter_df['x'].cat.codes
    xlim = [xvals.min() - axis_pad, xvals.max() + axis_pad]
    ylim = [yvals.min() - axis_pad, yvals.max() + axis_pad]

    if color is None:
        color = np.tile(np.expand_dims(adata.uns[cluster_key + '_colors'], 1),
                        adata.obs[key].nunique()).flatten()

    if figsize is None:
        figsize = (adata.obs[cluster_key].unique().size / 2,
                   adata.obs[key].unique().size / 2)
    plt.figure(figsize=figsize)
    plt.scatter(x, y, s=s, c=color)
    if axis_labels:
        plt.ylabel(key)
        plt.xlabel('Cluster label')
    plt.xticks(rotation=xrotation)
    plt.yticks(rotation=yrotation)
    ax1 = plt.gca()
    ax1.grid(False)
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)

    ax2 = ax1.twiny()

    ax2.scatter(x, y, s=s, c=color)  # I can't figure out a way around re-plotting

    ax2.set_xticklabels([str(c) for c in counts])
    plt.xticks(rotation=90)
    plt.xlabel('Cells per cluster')
    ax2.grid(False)
    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim)

    # Optional dot size legend
    if dot_size_legend_sizes is not None:
        _dot_size_legend(sizes=dot_size_legend_sizes,
                         display_scale_fcn=lambda siz: siz * 500 * size_factor,
                         marker='o',
                         labelspacing=dot_size_legend_labelspacing,
                         title=dot_size_legend_title)

    if show:
        plt.show()


def gene_lookup_table(adata: anndata.AnnData,
                      gene_names: List[str],
                      groupby_key: str,
                      expression_cutoffs: List[int] = [0],
                      use_raw: bool = False) -> pd.DataFrame:
    """Get a long-format table of expression inside and outside groups, for given genes."""

    df_dict = get_expressing_fraction(adata, gene_names, groupby_key,
                                      expression_cutoffs, use_raw)
    return exp_fraction_to_long_form_table(df_dict)


def get_expressing_fraction(adata: anndata.AnnData,
                            gene_names: List[str],
                            groupby_key: str,
                            expression_cutoffs: List[int] = [0],
                            use_raw: bool = False) -> Dict[int, Dict[str, pd.DataFrame]]:
    """Return fraction of cells expressing a gene(s) above a cutoff(s), for a given grouping of cells.

    out[Cutoff 1]['target']:

            group1  |  group2  | ...
           ---------o----------o----
    gene1 | target  |  target  | ...
    ----------------o----------o----
    gene2 | target  |  target  | ...
    ----------------o----------o----
    ...   |         |          |

    out[Cutoff 1]['other']:

            group1  |  group2  | ...
           ---------o----------o----
    gene1 |  other  |  other   | ...
    ----------------o----------o----
    gene2 |  other  |  other   | ...
    ----------------o----------o----
    ...   |         |          |

    out[Cutoff 2]['target']:

            group1  |  group2  | ...
           ---------o----------o----
    gene1 | target  |  target  | ...
    ----------------o----------o----
    gene2 | target  |  target  | ...
    ----------------o----------o----
    ...   |         |          |

    out[Cutoff 2]['other']:

            group1  |  group2  | ...
           ---------o----------o----
    gene1 |  other  |  other   | ...
    ----------------o----------o----
    gene2 |  other  |  other   | ...
    ----------------o----------o----
    ...   |         |          |

    etc...

    """

    assert len(gene_names) < 1000, f'You are trying to calculate for {len(gene_names)} genes.  Keep it to < 1000.'

    # taking some code from scanpy.tl.filter_rank_genes_groups to make code faster

    groupings = adata.obs[groupby_key].unique()

    # set up output dicts and dataframes
    out = {}
    for exp_cut in expression_cutoffs:
        out[exp_cut] = {}
        out[exp_cut]['target'] = pd.DataFrame(index=gene_names, columns=groupings)
        out[exp_cut]['other'] = pd.DataFrame(index=gene_names, columns=groupings)

    # dense matrix of count data
    if not use_raw:
        if len(gene_names) == 1:
            matrix = adata[:, gene_names].X
        else:
            matrix = adata[:, gene_names].X.toarray()
    else:
        if len(gene_names) == 1:
            matrix = adata.raw.X[:, [g in gene_names for g in adata.var.index]].todense()
        else:
            matrix = adata.raw.X[:, [g in gene_names for g in adata.var.index]].toarray()

    for k in groupings:
        adata.obs['__is_in_group__'] = pd.Categorical(adata.obs[groupby_key] == k)

        # obs_tidy has rows='__is_in_cluster__', columns=var_names
        obs_tidy = pd.DataFrame(matrix, columns=gene_names)
        obs_tidy.set_index(adata.obs['__is_in_group__'], '__is_in_group__', inplace=True)

        for exp_cut in expression_cutoffs:
            # compute fraction of cells having value > expression_cutoff
            # transform obs_tidy into boolean matrix
            obs_bool = (obs_tidy > exp_cut).astype(bool)

            # compute the fraction expressing per group
            fraction_obs = obs_bool.groupby(level=0).sum() / obs_bool.groupby(level=0).count()

            # get fractions
            in_frac = fraction_obs.loc[True].values
            out_frac = fraction_obs.loc[False].values

            out[exp_cut]['target'][k] = in_frac
            out[exp_cut]['other'][k] = out_frac

    return out


def exp_fraction_to_long_form_table(exp_df_dict: Dict[int, Dict[str, pd.DataFrame]]):
    """Go from output of get_expressing_fraction() to a single, long-format table."""

    # make a list of dataframes
    df_list = []
    for exp_cut in exp_df_dict.keys():
        # rename columns to contain the relevant information
        target_df = exp_df_dict[exp_cut]['target']
        target_df = target_df.rename(columns=dict(zip(target_df.columns,
                                                      [c + f'::pct.exp>{exp_cut}.target'
                                                       for c in target_df.columns])))
        df_list.append(target_df)
        other_df = exp_df_dict[exp_cut]['other']
        other_df = other_df.rename(columns=dict(zip(other_df.columns,
                                                    [c + f'::pct.exp>{exp_cut}.other'
                                                     for c in other_df.columns])))
        df_list.append(other_df)

    # make a big concatenated dataframe with all the data
    df = pd.concat(df_list, axis=1).copy()

    df = long_form_table(df, col_name_delim='::', new_col_name='cluster')

    return df.rename(columns=dict(zip(['index'], ['gene'])))


def long_form_table(df: pd.DataFrame,
                    col_name_delim: str,
                    new_col_name: str):
    """Go from a wide table, with information in column names, to a long table.

    Args:
        df: Input wide dataframe, with two pieces of information in the column
            names, delimited by col_name_delim
        col_name_delim: Delimiter for the two pieces of information in the
            column names
        new_col_name: Name of new column


    Note: The first piece of column name information becomes a new column in
    the dataframe.  The second piece stays with the data.

    See second answer here:
    https://stackoverflow.com/questions/38862832/pandas-melt-several-groups-of-columns-into-multiple-target-columns-by-name

    """

    # use a multi-index, splitting column names on col_name_delim
    df.columns = pd.MultiIndex.from_tuples(tuple(df.columns.str.split(col_name_delim)))

    df = df.stack(level=0).reset_index(level=1, drop=False).reset_index()

    return df.rename(columns=dict(zip(['level_1'], [new_col_name])))


def adata_to_single_cell_portal_files(
        adata: anndata.AnnData,
        cluster_keys: List[str],
        umap_key: str,
        save_directory: str = './',
        filename_prefix: str = '',
        expression_normalization: Callable[[np.float32], np.float32] = lambda x: np.log2(x + 1),
        use_gzip: bool = True) -> List[str]:
    """Prepare files that can be uploaded to the Single Cell Portal.

    Args:
        adata: The scanpy anndata object, annotated.
        cluster_keys: adata.obs keys for the clusterings to be uploaded as cluster files.
        umap_key: adata.obsm key specifying the UMAP (or tSNE) coordinates (e.g. 'X_scvi_umap')
        expression_normalization: Lambda function used to normalize count data.
        filename_prefix: Prefix to put on every output file.
        save_directory: Where to write output files.
        use_gzip: True to use gzip compression for matrix, barcodes, and genes.

    Returns:
        filenames: List of saved files.

    """

    assert os.access(save_directory, os.W_OK), 'Cannot write to directory ' + save_directory
    for k in cluster_keys:
        assert k in adata.obs.keys(), f'{k} is not a valid clustering in adata.obs.keys()'
    assert umap_key in adata.obsm.keys(), umap_key + ' is not a valid umap key in adata.obsm.keys()'

    print('adata.X data should be non-normalized count data (can be CellBender corrected).')
    pieces = inspect.getsourcelines(expression_normalization)[0][0].split(',')
    lam_str = ''
    for piece in pieces:
        if 'expression_normalization' in piece:
            lam_str = piece
    print('Using normalization:', lam_str.split('=')[1].rstrip())

    # normalize data
    sparse_mat = adata.X.copy()
    sparse_mat.data = expression_normalization(sparse_mat.data)

    # write matrix to sparse matrix market format
    matrix_file = os.path.join(save_directory, filename_prefix + '_matrix.mtx')
    mmwrite(matrix_file, sparse_mat.transpose())

    # sort matrix market format file by gene

    # this is taken from the Single Cell Portal at
    # https://raw.githubusercontent.com/broadinstitute/single_cell_portal/master/scripts/SortSparseMatrix.py
    headers = []
    with open(matrix_file) as matrix:
        line = next(matrix)
        while line.startswith("%"):
            headers = headers + [line]
            line = next(matrix)
        headers = headers + [line]
        df = pd.read_table(matrix, sep='\s+', names=['genes', 'barcodes', 'expr'])

    # sort sparse matrix
    print('Sorting sparse matrix')
    df = df.sort_values(by=['genes', 'barcodes'])

    # save sparse matrix
    with open(matrix_file, 'w+') as output:
        output.write(''.join(headers))
    df.to_csv(matrix_file, sep=' ', index=False, header=0, mode='a')

    # use gzip compression
    if use_gzip:
        with open(matrix_file, 'rb') as f_in:
            with gzip.open(matrix_file + '.gz', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(matrix_file)  # clean up un-compressed file
        matrix_file = matrix_file + '.gz'
    print('Saved sparse matrix to:', matrix_file)

    # save barcodes file
    if use_gzip:
        barcodes_file = os.path.join(save_directory, filename_prefix + '_barcodes.tsv.gz')
    else:
        barcodes_file = os.path.join(save_directory, filename_prefix + '_barcodes.tsv')
    pd.Series(adata.obs.index).to_csv(barcodes_file, sep='\t', index=False)
    print('Saving barcodes to:', barcodes_file)

    # save genes file
    if use_gzip:
        genes_file = os.path.join(save_directory, filename_prefix + '_genes.tsv.gz')
    else:
        genes_file = os.path.join(save_directory, filename_prefix + '_genes.tsv')
    pd.Series(adata.var.index).to_csv(genes_file, sep='\t', index=False)
    print('Saving genes to:', genes_file)

    # metadata file
    metadata_file = os.path.join(save_directory, filename_prefix + '_metadata.txt')
    header_line_1 = '\t'.join(['NAME'] + adata.obs.columns.tolist()) + '\n'

    def redefine_types(t):
        if str(t) == 'category':
            return 'group'
        elif str(t).startswith('int') or str(t).startswith('float'):
            return 'numeric'
        else:
            raise ValueError(f'unknown dtype {t}')

    header_line_2 = '\t'.join(['TYPE'] + adata.obs.dtypes.apply(redefine_types).values.tolist()) + '\n'

    with open(metadata_file, 'w+') as file:
        file.write(header_line_1)
        file.write(header_line_2)
    adata.obs.to_csv(metadata_file, sep='\t', header=0, mode='a')
    print('Saving cell metadata to:', metadata_file)

    # cluster file(s): cluster names and visualization coordinates
    cluster_files = []
    for k in cluster_keys:
        cluster_df = pd.DataFrame(data={'X': adata.obsm[umap_key][:, 0],
                                        'Y': adata.obsm[umap_key][:, 1],
                                        'Category': adata.obs[k].values},
                                  index=adata.obs[k].index)
        cluster_files.append(os.path.join(save_directory, filename_prefix + '_cluster_' + k + '.txt'))
        header_line_1 = '\t'.join(['NAME', 'X', 'Y', 'Category']) + '\n'
        header_line_2 = '\t'.join(['TYPE', 'numeric', 'numeric', 'group']) + '\n'

        # save file
        with open(cluster_files[-1], 'w+') as file:
            file.write(header_line_1)
            file.write(header_line_2)
        cluster_df.to_csv(cluster_files[-1], sep='\t', header=0, mode='a')
        print(f'Saving cluster {k} to:', cluster_files[-1])

    print('Done.')
    return [matrix_file, barcodes_file, genes_file, metadata_file] + cluster_files


def grouping_pca(adata: anndata.AnnData,
                 summed_counts: np.ndarray,
                 layer: Optional[str] = None,
                 row_labels: Optional[np.ndarray] = None,
                 color_labels: Optional[np.ndarray] = None,
                 color_label_order: Optional[np.ndarray] = None,
                 label_samples: bool = False,
                 n_hvgs: int = 2000,
                 max_pcs: int = 4,
                 which_pcs: List[int] = [0, 1],
                 title: Optional[str] = 'PCA of summed counts',
                 ms: int = 8,
                 marker_rotation: List[str] = ['o', 'p', 's', 'D', '<', '>', 'P', 'X', '*', 'H', 'd'],
                 alpha: float = 0.8,
                 figsize: Tuple[float] = (5, 5),
                 show: bool = True,
                 **kwargs):
    """Perform PCA at the level of counts summed by a grouping.  Create plots.

    Args:
        adata: AnnData object containing all relevant metadata.
        summed_counts: Matrix of summed counts.  Groupings rows by n_genes columns.
        layer: Layer of adata that contains raw counts (used for highly variable
            gene selection only).  If None, uses adata.X
        row_labels: Labels which name each row of summed_counts, and will name
            each point in the PCA plot.
        color_labels: Labels for summed_counts rows, where the same label gets
            the same color in the PCA plot.
        color_label_order: Optionally specify the order of the legend labels
            used as colors in the plot.
        label_samples: Display a label for each sample on the PCA plot.
        n_hvgs: Number of highly variable genes to use during PCA.
        max_pcs: How many PCs to compute.
        which_pcs: Principal components to plot.
        title: Plot title of choice.
        ms: Marker size for plot.
        marker_rotation: List of markers to be rotated through in the plot.
        alpha: Transparency alpha for pyplot in [0, 1]
        figsize: Figure size.
        show: True to show the plot.
        **kwargs: Passed to pyplot.plt()

    """

    if color_label_order is not None:
        assert len(set(np.unique(color_labels)) - set(color_label_order)) == 0, \
            'Input color_label_order must contain all the unique labels in color_labels'

    from sklearn.decomposition import PCA

    assert max(which_pcs) <= max_pcs, f'Select which_pcs with values less than {max_pcs}'

    # re-compute highly-variable genes
    try:
        if layer is not None:
            sc.pp.highly_variable_genes(adata, layer=layer, n_top_genes=n_hvgs, flavor='seurat_v3')
        else:
            sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs, flavor='seurat_v3')
        adata_tmp = adata
    except Exception:
        print('PCA plot warning: unable to use "seurat_v3" highly variable gene '
              'selection method.  Using old way.  It is recommended to upgrade '
              'to scanpy 1.6.0+')
        adata_tmp = adata.copy()
        sc.pp.normalize_total(adata_tmp, target_sum=1e4)
        sc.pp.log1p(adata_tmp)
        sc.pp.highly_variable_genes(adata_tmp, n_top_genes=n_hvgs)

    # restrict to highly-variable genes and reasonably-expressed genes
    gene_logic = (adata_tmp.var['highly_variable'].values
                  & (np.array(adata_tmp.X.sum(axis=0)).squeeze() > 10))
    summed_counts_hvg = summed_counts[:, gene_logic]
    del adata_tmp
    gc.collect()

    # log scale
    summed_counts_hvg = np.log1p(summed_counts_hvg) * 10000

    # normalize summed counts per grouping
    norm_summed_counts_hvg = summed_counts_hvg / summed_counts_hvg.sum(axis=1, keepdims=True)

    # z score genes
    means = norm_summed_counts_hvg.mean(axis=0)
    norm_summed_counts_hvg = norm_summed_counts_hvg - means
    z_gene = (np.power(norm_summed_counts_hvg, 2).mean(axis=0, keepdims=True)
              - np.power(norm_summed_counts_hvg.mean(axis=0, keepdims=True), 2) + 1e-5)
    z_norm_summed_sample_hvg = norm_summed_counts_hvg / z_gene

    # # plot a disgnostic for z-scoring of genes
    # plt.plot(np.argsort(means),
    #          np.sqrt(np.power(z_norm_summed_sample_hvg, 2).mean(axis=0)
    #                  - np.power(z_norm_summed_sample_hvg.mean(axis=0), 2)), '.', ms=2, alpha=0.1)
    # plt.xlabel('rank(mean)')
    # plt.ylabel('stdev')
    # plt.title('genes after z-scoring')
    # plt.show()

    # run PCA
    pca_obj = PCA(n_components=max_pcs)
    pca = pca_obj.fit_transform(z_norm_summed_sample_hvg)

    # create PCA plot
    plt.figure(figsize=figsize)

    if color_label_order is None:
        color_label_order = np.unique(color_labels)

    for i, c in enumerate(color_label_order):
        if np.sum(color_labels == c) == 0:
            # provide a way to skip missing legend labels
            # useful to share markers even when using different data subsets
            plt.plot(pca[color_labels == c, which_pcs[0]],
                     pca[color_labels == c, which_pcs[1]],
                     marker_rotation[i % len(marker_rotation)],
                     ms=ms, alpha=alpha, **kwargs)
        else:
            plt.plot(pca[color_labels == c, which_pcs[0]],
                     pca[color_labels == c, which_pcs[1]],
                     marker_rotation[i % len(marker_rotation)],
                     ms=ms, label=c, alpha=alpha, **kwargs)

    plt.title(title)
    x_variance = pca_obj.explained_variance_ratio_[which_pcs[0]]
    y_variance = pca_obj.explained_variance_ratio_[which_pcs[1]]
    plt.xlabel(f'PC {which_pcs[0]}: ({x_variance * 100:.0f}% variance)')
    plt.ylabel(f'PC {which_pcs[1]}: ({y_variance * 100:.0f}% variance)')
    plt.xticks([])
    plt.yticks([])
    plt.gca().legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

    # create PCA plot with labels
    if label_samples:
        from adjustText import adjust_text
        texts = [plt.text(pca[i, which_pcs[0]], pca[i, which_pcs[1]],
                          row_labels[i], color='black', fontsize=12)
                 for i in range(pca.shape[0])]
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray'))

    if show:
        plt.show()


def sum_adata_per_group(adata: anndata.AnnData,
                        groupby: str = 'sample',
                        layer: Optional[str] = None,
                        categorical_keys: List[str] = ['tissue'],
                        numeric_keys: List[str] = [],
                        min_cells: int = 10,
                        verbose: bool = True) -> Dict[str, np.ndarray]:
    """Return a dense matrix whose columns are genes and rows are groupings.

    Args:
        adata: AnnData object with un-scaled counts in adata.X
        groupby: adata.obs key used to group data over which sums will be performed
        layer: Layer of adata with count data to sum.
        categorical_keys: List of adata.obs keys which are constant within groupby
            and which will be kept track of when grouping
        numeric_keys: List of adata.obs keys which are numeric values which will
            be summed when grouping
        min_cells: Minimum required number of cells per group.  If the number of
            cells is less than min_cells, the grouping for that key value is
            omitted.
        verbose: Print status updates.

    """

    import warnings
    warnings.simplefilter(action='ignore', category=FutureWarning)

    # input checking
    check_keys = [groupby] + categorical_keys + numeric_keys
    for key in check_keys:
        assert key in adata.obs.keys(), \
            f"sum_adata_per_group got {key} as an adata.obs key, but the valid " \
            f"keys are: {adata.obs.keys()}"

    # number of groupings
    unique_groups = adata.obs[groupby].unique()
    nrows = unique_groups.size
    ngenes = adata.X.shape[1]

    if verbose:
        print(f'Summing adata.X per "{groupby}": {nrows} unique groups')

    # create output data structure
    out = dict()
    out['summed_counts'] = np.zeros((nrows, ngenes), dtype=np.int)
    out['groupby'] = np.empty(nrows, dtype='<U100')
    for key in categorical_keys:
        out[key] = np.empty(nrows, dtype='<U100')
    for key in numeric_keys:
        out[key] = np.zeros(nrows, dtype=np.float32)

    # go through each unique group
    for ind, group_name in enumerate(unique_groups):

        # subset data to this group
        adata_subset = adata[adata.obs[groupby] == group_name]

        # sum the counts in this group
        if len(adata_subset) >= min_cells:
            if layer is None:
                summed_counts = np.array(adata_subset.X.sum(axis=0),
                                         dtype=np.float).squeeze()
            else:
                summed_counts = np.array(adata_subset.layers[layer].sum(axis=0),
                                         dtype=np.float).squeeze()
            out['summed_counts'][ind, :] = summed_counts
        if verbose:
            print('.', end='')

        # keep track of the sample and cluster names
        out['groupby'][ind] = group_name
        for key in categorical_keys:
            # there should be only one unique value here, but if there's more
            # than one, concatenate them with commas (np.unique = alphabetical)
            out[key][ind] = ', '.join(np.unique(adata_subset.obs[key].unique().tolist()).tolist())
        for key in numeric_keys:
            out[key][ind] = adata_subset.obs[key].mean()  # mean so that it's intensive

        del adata_subset
        gc.collect()

    # take care of the empty groupings
    empty_groups = out['groupby'][out['summed_counts'].sum(axis=1) == 0]

    if verbose:
        print(f'\nThe following groupings yield fewer than {min_cells} cells:')
        print(empty_groups)

    if empty_groups.size > 0:
        not_empty_logic = [r not in empty_groups for r in out['groupby']]

        # remove empty groupings from counts and from metadata
        out['summed_counts'] = out['summed_counts'][not_empty_logic, :]
        out['groupby'] = out['groupby'][not_empty_logic]
        for key in categorical_keys + numeric_keys:
            out[key] = out[key][not_empty_logic]

    if verbose:
        print('Done summing counts.\n')

    return out


def pseudobulk_pca_plot(adata: anndata.AnnData,
                        groupby: str,
                        label: str,
                        label_order: Optional[List[str]] = None,
                        layer: Optional[str] = None,
                        n_hvgs: int = 2000,
                        min_cells: int = 25,
                        which_pcs: List[int] = [0, 1],
                        max_pcs: int = 4,
                        title: str = 'PCA: pseudobulk per sample',
                        label_each_point: bool = False,
                        ms: int = 5,
                        marker_rotation: List[str] = ['o', 'p', 's', 'D', '<', '>', 'P', 'X', '*', 'H', 'd'],
                        figsize: Tuple[float] = (5, 5),
                        show: bool = True,
                        **kwargs):
    """Create a PCA plot after summing counts over a given grouping of cells.

    Args:
        adata: AnnData object
        groupby: Key of adata.obs used to group cells for pseudobulk summation
        label: Key of adata.obs used to apply labels to pseudobulk PCA points.
            Can be different from "groupby", but should there should probably
            only be one unique "groupby" condition per "label".
        label_order: Order of labels in legend.
        layer: Layer of adata that contains raw counts.
        n_hvgs: Number of highly-variable genes to use when doing PCA.
        min_cells: Minimum number of cells that must exist to constitute a valid
            group from "groupby".  If there are fewer cells in a group, there
            will be no point plotted on the PCA plot for that group.
        which_pcs: Which PCs to plot.
        max_pcs: Number of PCs computed (can plot PCs up to this number).
        title: Title of figure
        label_each_point: True to add a text annotation label to every point
            (this usually looks like a mess but it can be done).
        ms: Marker size for points on the plot
        marker_rotation: Order in which marker labels are rotated through in
            the plot.
        figsize: Figure size
        show: True to show plot, False to not show it yet (helpful for by-hand
            additions / changes to labels, etc.)
        **kwargs: Passed to pyplot.plt()

    Example:

        >>> pseudobulk_pca_plot(
        >>>     adata,
        >>>     groupby='sample',
        >>>     layer='cellbender_0.01',
        >>>     label='tissue',
        >>>     label_order=adata.obs['tissue'].cat.categories,
        >>>     ms=8,
        >>>     markeredgewidth=0,  # helpful for PDF rendering
        >>> )

    """

    assert groupby in adata.obs.keys(), f'Specified groupby "{groupby}" but the ' \
        f'only columns in adata.obs are {adata.obs.keys()}'
    assert label in adata.obs.keys(), f'Specified label "{label}" but the ' \
        f'only columns in adata.obs are {adata.obs.keys()}'
    for pc in which_pcs:
        assert (pc >= 0) and (pc < max_pcs)

    # normalize counts per cell
    if layer is not None:
        assert layer in adata.layers.keys(), f'Specified layer "{layer}" but the ' \
            f'only layers in adata are {adata.layers.keys()}'
        adata.X = adata.layers[layer].copy()
    else:
        print('Assuming raw count data is stored in adata.X (if it is not, '
              'specify the location of raw count data using the "layer" input).')
    sc.pp.normalize_total(adata)

    # sum counts over groupings to create pseudobulk expression vectors
    summation = sum_adata_per_group(adata,
                                    groupby=groupby,
                                    layer=layer,
                                    categorical_keys=[label],
                                    min_cells=min_cells)

    # create PCA plot
    grouping_pca(
        adata,
        layer=layer,
        n_hvgs=n_hvgs,
        summed_counts=summation['summed_counts'],
        row_labels=summation['groupby'],
        color_labels=summation[label],
        color_label_order=label_order,
        label_samples=label_each_point,
        which_pcs=which_pcs,
        max_pcs=max_pcs,
        title=title,
        ms=ms,
        marker_rotation=marker_rotation,
        figsize=figsize,
        show=show,
        **kwargs,
    )


def pseudobulk_pca_plot_deseq2(
        adata: anndata.AnnData,
        groupby: str,
        label: str,
        label_order: Optional[List[str]] = None,
        layer: Optional[str] = None,
        n_hvgs: int = 500,
        min_cells: int = 25,
        which_pcs: List[int] = [0, 1],
        title: str = 'PCA: pseudobulk per sample',
        label_each_point: bool = False,
        ms: int = 5,
        marker_rotation: List[str] = ['o', 'p', 's', 'D', '<', '>', 'P', 'X', '*', 'H', 'd'],
        alpha: float = 0.8,
        figsize: Tuple[float] = (5, 5),
        temp_dir: str = './',
        show: bool = True,
        **kwargs):
    """Mark Chaffin's implementation of the DESeq2 pipeline for creating a
    pseudobulk PCA plot, summing counts over a given grouping of cells.

    NOTE: This is the "tried-and-true" alternative to `pseudobulk_pca_plot()`.

    Args:
        adata: AnnData object
        groupby: Key of adata.obs used to group cells for pseudobulk summation
        label: Key of adata.obs used to apply labels to pseudobulk PCA points.
            Can be different from "groupby", but should there should probably
            only be one unique "groupby" condition per "label".
        layer: Layer of adata that contains raw counts.
        n_hvgs: Number of highly-variable genes to use when doing PCA.
        min_cells: Minimum number of cells that must exist to constitute a valid
            group from "groupby".  If there are fewer cells in a group, there
            will be no point plotted on the PCA plot for that group.
        which_pcs: Which PCs to plot.
        title: Title of figure
        label_each_point: True to add a text annotation label to every point
            (this usually looks like a mess but it can be done).
        ms: Marker size for points on the plot
        marker_rotation: Order in which marker labels are rotated through in
            the plot.
        alpha: Transparency alpha for pyplot in [0, 1]
        figsize: Figure size
        temp_dir: Temporary directory (tsv will be written here) to pass data
            from R to python.
        show: True to show plot, False to not show it yet (helpful for by-hand
            additions / changes to labels, etc.)
        **kwargs: Passed to pyplot.plt()

    """

    assert os.access(temp_dir, os.W_OK), 'Cannot write to temp_dir, which was ' \
        f'specified as: {temp_dir}'
    assert groupby in adata.obs.keys(), f'Specified groupby "{groupby}" but the ' \
        f'only columns in adata.obs are {adata.obs.keys()}'
    assert label in adata.obs.keys(), f'Specified label "{label}" but the ' \
        f'only columns in adata.obs are {adata.obs.keys()}'
    for pc in which_pcs:
        assert (pc >= 0) and (pc < 8)

    # normalize counts per cell
    if layer is not None:
        assert layer in adata.layers.keys(), f'Specified layer "{layer}" but the ' \
            f'only layers in adata are {adata.layers.keys()}'
        adata.X = adata.layers[layer].copy()
    else:
        print('Assuming raw count data is stored in adata.X (if it is not, '
              'specify the location of raw count data using the "layer" input).')
    sc.pp.normalize_total(adata)

    # sum counts over groupings to create pseudobulk expression vectors
    summation = sum_adata_per_group(adata,
                                    groupby=groupby,
                                    categorical_keys=[label],
                                    min_cells=min_cells)

    # use rpy2 for communication between python and R
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()  # allow passing of numpy arrays to R

    # import R packages
    from rpy2.robjects.packages import importr
    DESeq2 = importr('DESeq2')
    datatable = importr('data.table')

    # it's too hard to translate all the R to rpy2
    # so wrap up most of the R in a function
    from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

    Rfunctionstring = """
        runPCA <- function(all_counts,
                           coldata,
                           nhvg,
                           tmpdir) {

            # create DESeq dataset
            dds <- DESeqDataSetFromMatrix(countData = all_counts,
                                          colData = coldata,
                                          design= ~ 1)

            # eliminate genes with counts below 10
            keep <- rowSums(counts(dds)) >= 10
            dds <- dds[keep,]

            dds <- DESeq(dds)

            # run a variance stabilizing transformation
            vsd <- vst(dds, blind=FALSE)

            # pick top n.hvg highly-variable genes
            rv <- rowVars(assay(vsd))
            var_gene <- as.data.frame(cbind(rownames(assay(vsd)), rv))
            select <- order(rv, decreasing = TRUE)[seq_len(min(nhvg, length(rv)))]

            # run PCA
            pca <- prcomp(t(assay(vsd)[select, ]))

            # calculate percentage of variation coming from each PC
            percentVar <- pca$sdev^2 / sum(pca$sdev^2)

            # make this data into a nice data.frame
            intgroup.df <- as.data.frame(colData(vsd))
            d <- data.frame(PC1 = pca$x[, 1], 
                            PC2 = pca$x[, 2], 
                            PC3 = pca$x[, 3], 
                            PC4 = pca$x[, 4], 
                            PC5 = pca$x[, 5], 
                            PC6 = pca$x[, 6],
                            PC7 = pca$x[, 7], 
                            PC8 = pca$x[, 8], 
                            intgroup.df, 
                            name = colData(vsd)[, 1])

            # save table
            write.table(d, paste0(tmpdir, 'pca_tmp.tsv'), sep='\t', row.names=F, quote=F)

            return(percentVar)
        }
        """

    # create the R function
    Rfunc = SignatureTranslatedAnonymousPackage(Rfunctionstring, 'Rfunc')

    # run the R function
    percent_variance = Rfunc.runPCA(all_counts=np.transpose(summation['summed_counts']),
                                    coldata=summation['groupby'],
                                    nhvg=n_hvgs,
                                    tmpdir=temp_dir)
    percent_variance = np.array(list(percent_variance), dtype=float)

    # read the R output into python
    df = pd.read_csv(os.path.join(temp_dir, 'pca_tmp.tsv'), sep='\t')
    pca = df[[f'PC{i}' for i in range(1, 9)]].to_numpy()

    # create PCA plot
    plt.figure(figsize=figsize)

    color_labels = summation[label]

    if label_order is None:
        label_order = np.unique(summation[label])

    for i, c in enumerate(label_order):
        if np.sum(color_labels == c) == 0:
            # provide a way to skip missing legend labels
            # useful to share markers even when using different data subsets
            plt.plot(pca[color_labels == c, which_pcs[0]],
                     pca[color_labels == c, which_pcs[1]],
                     marker_rotation[i % len(marker_rotation)],
                     ms=ms, alpha=alpha, **kwargs)
        else:
            plt.plot(pca[color_labels == c, which_pcs[0]],
                     pca[color_labels == c, which_pcs[1]],
                     marker_rotation[i % len(marker_rotation)],
                     ms=ms, label=c, alpha=alpha, **kwargs)

    plt.title(title)
    x_variance = percent_variance[which_pcs[0]]
    y_variance = percent_variance[which_pcs[1]]
    plt.xlabel(f'PC {which_pcs[0]}: ({x_variance * 100:.0f}% variance)')
    plt.ylabel(f'PC {which_pcs[1]}: ({y_variance * 100:.0f}% variance)')
    plt.xticks([])
    plt.yticks([])
    plt.gca().legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

    # create PCA plot with labels
    if label_each_point:
        from adjustText import adjust_text
        texts = [plt.text(pca[i, which_pcs[0]], pca[i, which_pcs[1]],
                          row_labels[i], color='black', fontsize=12)
                 for i in range(pca.shape[0])]
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray'))

    if show:
        plt.show()


def ppv_table(adata: anndata.AnnData,
              separate_groupings_key: Union[str, None],
              separate_groupings: Optional[Iterable[Union[int, str]]] = None,
              expression_cutoff: int = 0,
              layer: Optional[str] = None,
              verbose: bool = True) -> pd.DataFrame:
    """Create a table of positive predictive values of each gene for each grouping.

    Args:
        adata: AnnData object.
        separate_groupings_key: Cell groupings defined in adata.obs[separate_groupings_key]
        separate_groupings: If groupings are not adata.obs[separate_groupings_key].unique(),
            they can be defined here as desired.  E.g.:
            separate_groupings=([0, 1, 2], 3, 4, 5, [6, 8], 7).  None means the
            whole dataset is to be considered at once.
        expression_cutoff: Counts above which a cell is considered "expressing"
        layer: Layer of the AnnData object to use for count data.  If None, uses
            adata.X as count data.
        verbose: True to print a dot each time a grouping completes.

    Returns:
        df: Pandas dataframe indexed by gene, with groupings as column names.
            Each entry is PPV for the given gene in the given grouping.

    """

    # TODO: multiprocessing

    if layer is not None:
        assert layer in adata.layers.keys(), \
            f'Specified layer "{layer}" is not a valid adata.layers key: {adata.layers.keys()}'

    df = pd.DataFrame(index=adata.var.index)  # genes are the index

    if separate_groupings_key is not None:
        if separate_groupings is None:
            separate_groupings = adata.obs[separate_groupings_key].unique()
    else:
        separate_groupings = (['everything'],)

    for group in separate_groupings:

        if (type(group) is int) or (type(group) is str):
            group = [group]

        column_name = ', '.join([g for g in group])

        # determine cells in this grouping
        if separate_groupings_key is None:
            adata_group = adata
        else:
            adata_group = adata[adata.obs[separate_groupings_key].isin(group)]

        # determine cells outside this grouping
        logic_out_group = [ind not in adata_group.obs.index for ind in adata.obs.index]
        adata_out_group = adata[logic_out_group]

        # fraction of cells in this group
        prevalence = adata_group.shape[0] / adata.shape[0]

        # locate the count data we're using
        if layer is None:
            X_in = adata_group.X
            X_out = adata_out_group.X
        else:
            X_in = adata_group.layers[layer]
            X_out = adata_out_group.layers[layer]

        # true positives (in-group cells expressing) and false positives
        frac_exp_target = (np.array((X_in > expression_cutoff).sum(axis=0)).squeeze()
                           / adata_group.shape[0])
        frac_exp_other = (np.array((X_out > expression_cutoff).sum(axis=0)).squeeze()
                          / adata_out_group.shape[0])

        # definitions for clarity
        sensitivity = frac_exp_target  # true positive rate
        false_positive_rate = frac_exp_other
        specificity = 1. - false_positive_rate  # true negative rate

        # calculate PPV in two steps to avoid NaNs and their warnings for division by zero
        positive_predictive_value = sensitivity * prevalence
        nnz_logic = (positive_predictive_value != 0.)
        positive_predictive_value[nnz_logic] = (positive_predictive_value[nnz_logic]
                                                / (sensitivity[nnz_logic] * prevalence
                                                   + (1. - specificity[nnz_logic]) * (1. - prevalence)))

        # add the data as a column to the dataframe
        df[column_name] = pd.Series(positive_predictive_value, index=adata.var.index)

        # and try our best to free up memory (possible bug with anndata or something...?)
        gc.collect()  # garbage collect in python

        if verbose:
            print('.', end='')

    if verbose:
        print('\n', end='')

    return df


def calculate_bkg_prob(adata: anndata.AnnData,
                       testing_key: str,
                       df: pd.DataFrame,
                       bcs_included: Optional[List[str]] = None,
                       separate_groupings_column: Optional[str] = None,
                       separate_groupings_key: Optional[str] = None,
                       test_group_column: str = 'test.group',
                       comparison_group_column: str = 'comparison',
                       comparison_separator: Optional[str] = None,
                       sample_key: str = 'sample',
                       layer: Optional[str] = None,
                       gene_column: str = 'gene',
                       verbose: bool = True) -> pd.DataFrame:
    """Calculate a heuristic estimate of the probability that each gene in a
    DE test is showing up as significant due to background RNA contamination.
    Args:
        adata: AnnData object
        testing_key: adata.obs key used to define testing groups
        df: A dataframe obtained from a differential expression test, e.g.
            limma_voom_DE()
        bcs_included: An alternative to specifying `separate_groupings_column`
             and `separate_groupings_key`, this is a way to specify which barcodes
             were included during DE testing, if it was not all of adata.  For
             example, if the DE test was done just on the cardiomyocytes, this would
             be a list of all the barcodes of cardiomyocytes.
        separate_groupings_column: If there are several DE tests in one table,
            this column is used to group the dataframe into entirely separate tests
        separate_groupings_key: The key to adata.obs that corresponds to the
            groupings in separate_groupings_column
        test_group_column: Name of the column that specifies the test group
        comparison_group_column: Name of the column that specifies the
            comparison group
        comparison_separator: Separator (such as '_') that separates labels in
            the df[comparison_group_column] values
        sample_key: Key of adata.obs that specifies sample (single 10x run, where
            background RNA can be shared)
        layer: Specify the layer of the AnnData object that contains the counts
            to be used.  None defaults to adata.X
        gene_column: 'gene' unless otherwise specified
        verbose: True to print out intermediate updates
    Returns:
        df: The full DE table, augmented with additional columns
    Example:
    Imagine a DE test was conducted using a subset of the full adata: adata_subset,
    and the testing_key in the DE model was 'cluster'.
    The resulting DE dataframe is 'df'.
    df = calculate_bkg_prob(adata,
                            testing_key='cluster',
                            df=df,
                            bcs_included=adata_subset.obs.index)
    """

    # TODO: multiprocessing

    # check inputs
    assert test_group_column in df.columns, \
        f'test_group_column "{test_group_column}" is not in df.columns:\n{df.columns}'
    assert comparison_group_column in df.columns, \
        f'comparison_group_column "{comparison_group_column}" is not in df.columns:\n{df.columns}'
    assert testing_key in adata.obs.keys(), \
        f'testing_key "{testing_key}" is not a key of adata.obs:\n{adata.obs.keys()}'
    if separate_groupings_key is not None:
        assert separate_groupings_key in adata.obs.keys(), \
            f'separate_groupings_key "{separate_groupings_key}" is not a key of adata.obs:\n{adata.obs.keys()}'
    if separate_groupings_column is not None:
        assert separate_groupings_key is not None, \
            'Included separate_groupings_column, so you must include a separate_groupings_key to ' \
            'designate where in adata.obs the corresponding group information is stored.'
    assert sample_key in adata.obs.keys(), \
        f'sample_key should specify a categorical column of adata.obs that contains ' \
        f'information about which 10x run each cell comes from.\nThe columns of ' \
        f'adata.obs are {adata.obs.keys()}'

    # try to approximate the background RNA profile
    approx_background_profile = np.array(adata.X.sum(axis=0)).squeeze() / adata.X.sum()  # normalized
    gene_order = np.argsort(approx_background_profile)
    gene_to_cdf_x = np.argsort(gene_order)
    cdf = np.cumsum(approx_background_profile[gene_order])
    adata.var['gene_bkg_prob'] = cdf[[gene_to_cdf_x[i] for i in range(len(adata.var))]]

    if verbose:
        plt.figure(figsize=(12, 3))
        plt.plot(cdf)
        plt.title('Empirical CDF of "background" gene probability')
        plt.xlabel('gene index, sorted')
        plt.ylabel('background probability')
        plt.show()
        print(f'Gene most probable to be background RNA is \n{adata.var.iloc[gene_order[-1]]}\n\n')

    df_out = df.copy()
    new_columns = ['bkg.prob', 'PPV.expr>0', 'PPV.expr>1',
                   'frac.cells.in.group.expr>0', 'frac.cells.out.group.expr>0',
                   'frac.cells.in.group.expr>1', 'frac.cells.out.group.expr>1']
    for col in new_columns:
        try:
            del df_out[col]
        except Exception:
            pass

    original_column_order = df_out.columns.tolist()

    # create new columns, empty
    for col in new_columns:
        df_out[col] = np.nan

    # figure out the unique groupings of cells that have different bkg.prob
    if comparison_separator is None:
        if np.any(np.array([', ' in g for g in df[comparison_group_column].unique()])):
            joiner = ', '
        else:
            joiner = '_'
    else:
        joiner = comparison_separator

    if separate_groupings_column is not None:
        test_and_comparison = df.apply(
            lambda x: ':::'.join(x[separate_groupings_column].split(joiner)
                                 + [joiner.join(sorted(x[test_group_column].split(joiner)
                                                       + x[comparison_group_column].split(joiner)))]), axis=1)
    else:
        test_and_comparison = df.apply(
            lambda x: ':::'.join([joiner.join(sorted(x[test_group_column].split(joiner)
                                                     + x[comparison_group_column].split(joiner)))]), axis=1)

    df_out['tmp'] = test_and_comparison
    if verbose:
        print(f'Unique conditions under which bkg.prob will be computed:\n{df_out["tmp"].unique()}')

    # subset to relevant genes
    relevant_genes = set(df[gene_column].unique())
    adata = adata[:, [g in relevant_genes for g in adata.var_names]]

    # go through each unique (grouping + test.group + comparison)
    for group in df_out['tmp'].unique():

        # subset adata to grouping
        if bcs_included is not None:
            adata_subset = adata[bcs_included]
        else:
            adata_subset = adata

        group_labels = group.split(':::')[0].split(joiner)
        if verbose:
            print(f'Working on grouping {group_labels}')
        if separate_groupings_key is not None:
            adata_subset = adata_subset[adata_subset.obs[separate_groupings_key].isin(group_labels)]

        # subset adata to cells involved in test or comparison
        test_and_comparison_labels = group.split(':::')[-1].split(joiner)
        if verbose:
            print(f'Working on test conditions that include {test_and_comparison_labels}')
        adata_subset = adata_subset[adata_subset.obs[testing_key]
            .isin(test_and_comparison_labels)]
        if verbose:
            print(f'{adata_subset.shape[0]} cells')

        # cells inside and outside the group (group includes test and comparison)
        # the out-group only contains cells from the same samples as the in-group
        sample_bcs = adata[adata.obs[sample_key].isin(adata_subset.obs[sample_key].unique())].obs.index
        bc_out_group = list(set(sample_bcs) - set(adata_subset.obs.index))

        # locate the count data we're using
        if layer is None:
            X_in = adata_subset.X
        else:
            X_in = adata_subset.layers[layer]

        # calculate fraction of cells in group expressing > 0 counts
        tp_0 = np.array((X_in > 0).sum(axis=0)).squeeze()
        frac_exp_gr_0_target = tp_0 / adata_subset.shape[0]

        # calculate fraction of cells in group expressing > 1 counts
        tp_1 = np.array((X_in > 1).sum(axis=0)).squeeze()
        frac_exp_gr_1_target = tp_1 / adata_subset.shape[0]

        # things ill-defined without an "out-group"
        frac_exp_gr_0_other = np.nan
        frac_exp_gr_1_other = np.nan
        ppv_0 = np.nan
        ppv_1 = np.nan
        bkg_prob = 0.

        if len(bc_out_group) == 0:
            # this is the case where all cells are in the "in" group
            # but it might not be the case for other comparisons (contrasts)
            print('Skipping this grouping + test condition, since all cells are included '
                  '(and we cannot estimate background without out-of-group cells).')

        else:
            adata_out_group = adata[bc_out_group]

            # locate the count data we're using
            if layer is None:
                X_out = adata_out_group.X
            else:
                X_out = adata_out_group.layers[layer]

            # calculate fraction of cells out group expressing > 0 counts
            fp_0 = np.array((X_out > 0).sum(axis=0)).squeeze()
            frac_exp_gr_0_other = fp_0 / X_out.shape[0]

            # calculate fraction of cells out group expressing > 1 counts
            fp_1 = np.array((X_out > 1).sum(axis=0)).squeeze()
            frac_exp_gr_1_other = fp_1 / X_out.shape[0]

            # TODO: amend PPV calculation somehow...
            # TODO: maybe just do tp / (tp + fp) ...

            ppv_0 = frac_exp_gr_0_target / (frac_exp_gr_0_target + frac_exp_gr_0_other + 1e-5)
            ppv_1 = frac_exp_gr_1_target / (frac_exp_gr_1_target + frac_exp_gr_1_other + 1e-5)

            #         ppv_0 = tp_0 / (tp_0 + fp_0 + 1.)
            #         ppv_1 = tp_1 / (tp_1 + fp_1 + 1.)
            mean_ppv = (ppv_0 + ppv_1) / 2.

            # create an estimated "background probability"
            # that this DE result is influenced by background RNA
            prob_from_another_cluster = 1. - mean_ppv  # false positive rate
            bkg_prob = prob_from_another_cluster * adata.var['gene_bkg_prob']

        # add gene count information
        gene_grouping_df = pd.DataFrame(data={gene_column: adata.var.index.values,
                                              'tmp': group,
                                              'bkg.prob': bkg_prob,
                                              'PPV.expr>0': ppv_0,
                                              'PPV.expr>1': ppv_1,
                                              'frac.cells.in.group.expr>0': frac_exp_gr_0_target,
                                              'frac.cells.out.group.expr>0': frac_exp_gr_0_other,
                                              'frac.cells.in.group.expr>1': frac_exp_gr_1_target,
                                              'frac.cells.out.group.expr>1': frac_exp_gr_1_other})

        # add data to table: merge with groupby
        # https://stackoverflow.com/questions/52397657/merge-pandas-dataframe-with-overwrite-of-columns
        df_out = (pd.merge(left=df_out, right=gene_grouping_df, how='left',
                           on=[gene_column, 'tmp'], suffixes=('___x', '___y'))
                  .groupby(lambda x: x.split('___')[0], axis=1).last())

        gc.collect()

    # clean up and re-order columns to match input order
    del df_out['tmp']
    df_out = df_out[original_column_order + new_columns]

    return df_out


def limma_voom_DE(adata: anndata.AnnData,
                  summation_key: str,
                  separate_groupings_key: Optional[str],
                  testing_key: str,
                  model: str,
                  model_keys: List[str],
                  one_versus_all_contrasts: bool,
                  min_cells_per_test_group: int = 10,
                  duplicate_correlation_key: Optional[str] = None,
                  separate_groupings: Optional[Iterable[Union[int, str]]] = None,
                  separate_groupings_for_background_calc: Optional[Iterable[Union[int, str]]] = None,
                  fdr: float = 1.,
                  working_directory: str = '.',
                  additional_contrasts: List[str] = [],
                  low_expression_mean_threshold: float = 1.,
                  voom_function: str = 'voomWithQualityWeights',
                  voom_lowess_span: float = 0.2,  # this is < 0.5, the actual default in voom
                  calculate_background_estimate: bool = True,
                  full_adata_for_background_estimate: Optional[anndata.AnnData] = None,
                  gene_id_key: str = 'gene_id',
                  create_PCA_plot: bool = True,
                  label_PCA_plot: bool = False,
                  verbose: bool = True) -> pd.DataFrame:
    """Perform formal differential expression tests using the limma-voom framework.

    NOTE: This follows the count-summation recommendation in Lun and Marioni 2017,
    see https://academic.oup.com/biostatistics/article/18/3/451/2970368

    Args:
        adata: AnnData object containing counts matrix in adata.X (un-scaled
            integer counts) along with per-cell metadata in adata.obs
        summation_key: Limma-voom testing according to the approach of Lun and
            Marioni recommends summation of counts over samples considered
            replicates, in order to maintain FDR control.  Cells with the same
            adata.obs[summation_key] will be summed as a replicate group.
            For example, this could be "individual" depending on the design.
        separate_groupings_key: adata.obs key defining groups of data to be
            tested separately.  In lots of single-cell work, it is best to
            perform DE tests for one cell type at a time.  For example, this
            could be a louvain clustering key in adata.obs.  If on the other
            hand, you are testing for marker genes between clusters, this can
            be left None.
        separate_groupings: Groupings of adata.obs[separate_groupings_key] if
            you do not want to use the default, which is that everything in
            adata.obs[separate_groupings_key].unique() is a group.
        separate_groupings_for_background_calc: Groupings of adata to use for
            calculating background in each DE test (the genes highly expressed
            in some other group may be background RNA).  By default, this is
            separate_groupings.
        testing_key: adata.obs key that defines groupings between which testing
            contrasts are created.
        min_cells_per_test_group: Minimum number of cells allowed to form a
            group used for testing. Groupings with fewer cells will be excluded
            from tests.
        model: Model string, e.g., '~ treatment + individual'.  The words used
            in this string must exactly match model_keys, and it must start
            with the tilde character '~'.
        model_keys: adata.obs keys referred to by the model.
        one_versus_all_contrasts: Auto-generate contrasts by doing one group
            versus all others lumped together, for each group in testing_key.
        duplicate_correlation_key: adata.obs key that denotes samples which are
            correlated, but whose correlation is not captured by the model.
            Uses the duplicateCorrelation function in limma.
        fdr: False discovery rate.  Only genes with Benjamini-Hochberg adjusted
            p-values below fdr will be included in the final table.
        working_directory: Directory for intermediate output files.
        additional_contrasts: List of user-specified contrasts that do not follow
            the one-versus-all model.  Each contrast is a string, e.g.,
            'LV.vs.RV = chamberLV - chamberRV', where the names on the
            right-hand side of the equals sign are exact matches to values in
            adata.obs[testing_key] with the string testing_key prepended (here
            "chamber") without a space.  NOTE: There must be whitespace between
            level names ("chamberLV") and mathematical symbols.
        low_expression_mean_threshold: If mean expression of a gene over the
            summed-groups is less than this value, it will not be included
            in DE testing.
        voom_function: Must be 'voom' or 'voomWithQualityWeights'.  The quailty
            weights have to do with an attempt to be robust to outlier samples.
        voom_lowess_span: Voom parameter: width of lowess smoothing window as
            a proportion.  Default in limma's voom function is 0.5.  Testing in
            PCL has shown that a lower value (near 0.2) produces results that
            depend less on the details of how genes are filtered.
        calculate_background_estimate: True to compute extra columns including
            PPV per gene per cell type and background contaminant probability.
        full_adata_for_background_estimate: The full adata object, only needed
            if calculate_background_estimate is True, and it could be the same
            as adata if you do not subset adata at all before running this
            function.
        gene_id_key: Gene IDs should be in adata.var[gene_id_key], if available.
        create_PCA_plot: True to make a PCA plot of summed counts per grouping,
            to look for sample outliers or similar samples.
        label_PCA_plot: Add labels to every dot on the PCA plot.  Can be busy!
        verbose: True to print intermediate messages.

    Returns:
        df: A pandas dataframe with all the results included.

    """

    if verbose:
        from IPython.display import display, HTML

    # input checks
    for name, key in {'summation_key': summation_key,
                      'separate_groupings_key': separate_groupings_key,
                      'duplicate_correlation_key': duplicate_correlation_key,
                      'testing_key': testing_key,
                      'model_keys': model_keys}.items():
        if name == 'model_keys':
            for subkey in model_keys:
                assert subkey in adata.obs.keys(), \
                    f'Input "{name}" contained "{subkey}", which is not one of the ' \
                    f'keys of adata.obs: {adata.obs.keys()}'
        else:
            if key is None:
                continue
            assert key in adata.obs.keys(), \
                f'Input "{name}" was "{key}", which is not one of the keys ' \
                f'of adata.obs: {adata.obs.keys()}'
    if separate_groupings is not None:
        for grouping in separate_groupings:
            # if the grouping itself is a list, look at each element
            if not hasattr(grouping, '__iter__'):
                grouping = [grouping]
            for element in grouping:
                assert element in adata.obs[separate_groupings_key].unique(), \
                    f'The identifier "{element}" specified in separate_groupings ' \
                    f'was not found in adata.obs["{separate_groupings_key}"].unique(): ' \
                    f'{adata.obs[separate_groupings_key].unique()}'
    working_directory = working_directory.rstrip('/')
    assert os.access(working_directory, os.W_OK), \
        f'Cannot write to working_directory:\n{working_directory}'
    assert (voom_lowess_span > 0) and (voom_lowess_span < 1), \
        'voom_lowess_span must be between zero and one, exclusive.'
    assert voom_function in ['voom', 'voomWithQualityWeights'], \
        f'voom_function must be either "voom" or "voomWithQualityWeights", but ' \
        f'the specified input was "{voom_function}"'
    if calculate_background_estimate:
        assert full_adata_for_background_estimate is not None, \
            'calculate_background_estimate is True, and so ' \
            'full_adata_for_background_estimate must not be None (even if it ' \
            'is the same as adata, then make this input adata as well)'

    # ensure that the first term in the model is the testing key
    # currently this is necessary due to R's factor-naming convention of prepending
    whitelist = string.ascii_letters + '_' + ' '
    model_vars = np.array(''.join(char for char in model.split('~')[1]  # RHS of equation
                                  if char in whitelist).split())  # don't keep the math
    if model_vars[0] != testing_key:
        model_rhs = ' + '.join([testing_key] + model_vars.tolist().remove(testing_key))
        if '0' in model:
            model = '~ 0 + ' + model_rhs
        else:
            model = '~ ' + model_rhs
        print(f'Warning: The testing_key was specified as {testing_key}, but the '
              f'first term in the model is {model_vars[0]}. In the current '
              f'implementation, the first term in the model must be the '
              f'testing_key. Will run using\nmodel="{model}"')

    # use rpy2 for communication between python and R
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()  # allow passing of numpy arrays to R

    # import R packages
    data_table = importr('data.table')
    edgeR = importr('edgeR')
    DESeq2 = importr('DESeq2')
    limma = importr('limma')
    base = importr('base')

    # limit the analysis to one grouping at a time
    if separate_groupings_key is not None:
        if separate_groupings is None:
            separate_groupings = adata.obs[separate_groupings_key].unique()
    else:
        separate_groupings = (['everything'],)

    if separate_groupings_for_background_calc is None:
        separate_groupings_for_background_calc = separate_groupings

    # create a PPV per grouping per gene table to do later contamination calculations
    if calculate_background_estimate:
        if verbose:
            print('Calculating PPV based on counts > 0')
        ppv_0_df = ppv_table(adata,
                             separate_groupings_key=separate_groupings_key,
                             separate_groupings=separate_groupings_for_background_calc,
                             expression_cutoff=0,
                             verbose=verbose)
        if verbose:
            display(HTML(ppv_0_df.head().to_html()))
            print('Calculating PPV based on counts > 1')
        ppv_1_df = ppv_table(adata,
                             separate_groupings_key=separate_groupings_key,
                             separate_groupings=separate_groupings_for_background_calc,
                             expression_cutoff=1,
                             verbose=verbose)

        if verbose:
            display(HTML(ppv_1_df.head().to_html()))

    # the list that will contain all the output dataframes
    dfs = []

    for group in separate_groupings:

        # the list that will contain the output dataframes from this group
        group_dfs = []
        duplicate_correlation_is_not_degenerate = True  # assume the best at first

        if verbose:
            if separate_groupings_key is not None:
                print(f'Working on {separate_groupings_key} group "{group}" ', end='')
            else:
                print('Working on the entire dataset at once ', end='')
            print('========================================\n')

        # copy a subset of adata corresponding to this group
        if separate_groupings_key is None:
            adata_group = adata
        else:
            if (type(group) is int) or (type(group) is str):
                group = [group]  # list
            adata_group = adata[adata.obs[separate_groupings_key].isin(group)].copy()

        # some stats on the experimental design if using duplicateCorrelation
        if duplicate_correlation_key is not None:
            dc_table = pd.crosstab(adata_group.obs[duplicate_correlation_key],
                                   adata_group.obs[testing_key])
            testing_keys_per_dc_key = np.sum((dc_table.to_numpy() > min_cells_per_test_group),
                                             axis=1, keepdims=False)

        if verbose:
            print(f'total cells in group = {adata_group.shape[0]}\n')
            test_key_dtype = adata_group.obs[testing_key].dtype.name
            test_key_is_numeric = test_key_dtype.startswith('float') or test_key_dtype.startswith('int')
            if test_key_is_numeric:
                print(f'numeric values in "{testing_key}"\n')
            else:
                print(f'cells per "{testing_key}":')
                print(adata_group.obs[testing_key].value_counts())
                print(' ')
            print(f'cells per "{summation_key}":')
            print(adata_group.obs[summation_key].value_counts())
            print(' ')
            if duplicate_correlation_key is not None:
                print(f'number of unique "{testing_key}"s per "{duplicate_correlation_key}":')
                print(f'{duplicate_correlation_key}\tnum unique {testing_key}')
                print('\n'.join(['\t\t\t'.join([str(a), str(b)])
                                 for a, b in zip(dc_table.index.tolist(), testing_keys_per_dc_key)]))
                print(' ')

        # check for issues running duplicateCorrelation: degeneracy
        # meaning: there is no duplicate_correlation group that has a sample in
        # more than one test group
        if duplicate_correlation_key is not None:
            if np.all(testing_keys_per_dc_key <= 1):
                # duplicateCorrelation will be degenerate, so don't use it
                print(f'WARNING: duplicateCorrelation is totally degenerate, since there is '
                      f'no "{duplicate_correlation_key}" in more than one "{testing_key}".\n'
                      f'Analysis will proceed without duplicateCorrelation!\n')
                duplicate_correlation_is_not_degenerate = False

        # some keys have numeric data and others are categoricals
        numeric_keys = [key for key in model_keys
                        if (adata_group.obs[key].dtype.name.startswith('float')
                            or adata_group.obs[key].dtype.name.startswith('int'))]
        nonnumeric_keys = [key for key in model_keys
                           if not (adata_group.obs[key].dtype.name.startswith('float')
                                   or adata_group.obs[key].dtype.name.startswith('int'))]

        # sum adata per group
        cat_keys = nonnumeric_keys
        if (duplicate_correlation_key is not None) and duplicate_correlation_is_not_degenerate:
            cat_keys = nonnumeric_keys + [duplicate_correlation_key]
        out = sum_adata_per_group(adata=adata_group,
                                  groupby=summation_key,
                                  categorical_keys=cat_keys,
                                  numeric_keys=numeric_keys,
                                  min_cells=min_cells_per_test_group,
                                  verbose=verbose)

        # check and see if we are left with more than one summed sample
        if out['summed_counts'].shape[0] < 2:
            # gracefully give up, and skip testing this grouping
            rowname = out['groupby']
            print(f'WARNING: No tests to be run for group {group}, which has only {rowname}')
            print('Skipping this grouping! *****************************\n')
            continue

        # optionally create PCA plot
        if create_PCA_plot:
            try:
                grouping_pca(adata=adata_group,
                             summed_counts=out['summed_counts'],
                             row_labels=out['groupby'],
                             color_labels=out[testing_key],
                             label_samples=label_PCA_plot,
                             title='PCA of summed counts\n[' + ', '.join(group) + '] by ' + testing_key)
            except Exception:
                pass  # don't fail because of a PCA plot

            # use R for another version: MDS plot
            # limma.plotMDS(out['summed_counts'], col=out[testing_key])  # for some reason this dies
            # perhaps here:
            # https://stackoverflow.com/questions/43110228/how-to-plot-inline-with-rpy2-in-jupyter-notebook

        # set up variables for DE testing in R
        gene_ids = np.array(adata.var.index)
        tmp_keys = nonnumeric_keys
        if (duplicate_correlation_key is not None) and duplicate_correlation_is_not_degenerate:
            tmp_keys = tmp_keys + [duplicate_correlation_key]
        for key in tmp_keys:
            # squish names by eliminating spaces if they exist
            out[key] = np.array([s.replace(" ", "") for s in out[key]])

        # if the testing key is numeric values, no contrasts need to be defined
        if testing_key in numeric_keys:
            contrast_list = [testing_key]
            use_contrasts = False  # is this correct in general?
            if '0' in model:
                model = f'~ {model.split("~ ")[-1].split("0 + ")[-1]}'
                print(f'Warning: Numerical testing factor. Coercing model to "{model}"')

        else:

            # if the model contains an intercept, no contrasts need to be defined
            if '0' not in model:
                use_contrasts = False
                whitelist = string.ascii_letters + '_' + ' '
                model_vars = np.array(''.join(char for char in model.split('~')[1]  # RHS of equation
                                              if char in whitelist).split())  # don't keep the math
                # the first column after the intercept
                first_term = model_vars[0]
                assert testing_key == first_term, \
                    f'You are using a model with an intercept.  The specified model ' \
                    f'must have the testing_key (in this case, "{testing_key}") as its ' \
                    f'first term, but the first term in this model is "{first_term}"'

                # testing key is categorical
                testing_levels = adata.obs[testing_key].unique()
                contrast_list = [testing_levels[-1] + '.vs.not = '
                                 + first_term + testing_levels[-1]]

            else:
                use_contrasts = True

                # set up list of "contrasts" for tests
                contrast_list = []
                levels = set([t for t in np.unique(out[testing_key])])

                # define contrasts programmatically
                if one_versus_all_contrasts:
                    if len(levels) <= 1:
                        pass  # no contrasts if there's just one level
                    else:
                        for t in levels:
                            levels_without_prefix = [s for s in levels - {t}]
                            levels_with_prefix = [testing_key + s for s in levels - {t}]
                            contrast = (t + '.vs.' + '_'.join(levels_without_prefix) + ' = '
                                        + testing_key + t + ' - (' + ' + '.join(levels_with_prefix) + ') / '
                                        + str(len(levels) - 1))
                            contrast_list.append(contrast)

                            # just do one contrast if there are two levels, since the other one is symmetric
                            if len(levels) == 2:
                                break  # just do one

            # add user-defined contrasts to list
            for contrast in additional_contrasts:

                # check to make sure this contrast is applicable for this grouping
                # (some groupings could lack levels needed for some contrasts, due
                # to a lack of cells in that condition)
                levels_with_prefix = [testing_key + s for s in levels]
                whitelist = string.ascii_letters + ' ' + '_' + string.digits  # needed to allow numerals in names
                contrast_vars = np.array(''.join(char for char in contrast.split('=')[1]  # RHS of equation
                                                 if char in whitelist).split())
                contrast_vars = contrast_vars[[not s.isnumeric() for s in contrast_vars]]  # remove pure numbers
                missing_levels = contrast_vars[[v not in levels_with_prefix
                                                for v in contrast_vars]]
                if missing_levels.size > 0:
                    print(f'Warning: Levels {missing_levels} were missing from the '
                          f'levels in "{separate_groupings_key}" {group} summed over '
                          f'"{summation_key}".\nTherefore the user-defined contrast "{contrast}" '
                          f'is being excluded from this test.')
                else:
                    # ensure this is not redundant with existing one-versus-all contrasts
                    existing_contrast_names = [eqn.split('=')[0].strip() for eqn in contrast_list]
                    this_contrast_name = contrast.split('=')[0].strip()
                    contrast_is_redundant = False
                    for name in existing_contrast_names:
                        if ({this_contrast_name.split('.')[0], this_contrast_name.split('.')[-1]}
                                == {name.split('.')[0], name.split('.')[-1]}):
                            contrast_is_redundant = True
                    if not contrast_is_redundant:
                        contrast_list.append(contrast)
                    else:
                        print(f'Warning: Dropping redundant contrast "{contrast}"')

        # if there are no contrasts left, then we are done here: there's no test
        if len(contrast_list) < 1:
            print(f'No contrasts remaining!  Skipping DE testing for "{separate_groupings_key}" {group}\n')
            continue

        # turn contrasts into one large string command for R
        make_contrast_command = 'makeContrasts(' + ', '.join(contrast_list) + ', levels=colnames(design))'

        if verbose:
            print('Contrasts: ', end='')
            if not use_contrasts:
                print('(not fit using contrasts.fit(); fit using an intercept)', end='')
            print('', end='\n')
            print('\n'.join(contrast_list) + '\n')

        # handle duplicateCorrelation
        if (duplicate_correlation_key is not None) and duplicate_correlation_is_not_degenerate:
            use_duplicatecorr = True
            duplicatecorr = out[duplicate_correlation_key]
            duplicatecorr_name = duplicate_correlation_key
            corfit_command = 'corfit <- duplicateCorrelation(v.all, design, block=duplicatecorr)'
            vfit_command = 'vfit <- lmFit(v.all, design, block=duplicatecorr, correlation=corfit$consensus)'
        else:
            use_duplicatecorr = False
            duplicatecorr = 'none'
            duplicatecorr_name = 'none'
            corfit_command = '\n'
            vfit_command = 'vfit <- lmFit(v.all, design)'

        # it's too hard to translate all the R to rpy2
        # so wrap up most of the R in a function
        from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

        Rfunctionstring = """
        runLimma <- function(summed_sample, 
                             gene_ids, 
                             groups, 
                             model,
                             verbose, 
                             use_contrasts,
                             use_duplicatecorr,
                             duplicatecorr,
                             working_directory,
                             """ + ', '.join(model_keys) + """) {

            # python boolean to R boolean
            verbose <- as.logical(verbose)
            use_duplicatecorr <- as.logical(use_duplicatecorr)
            use_contrasts <- as.logical(use_contrasts)

            # transpose
            counts <- t(summed_sample)

            # name the rows with gene names
            rownames(counts) <- gene_ids
            colnames(counts) <- groups
            counts <- as.matrix(counts)

            if (isTRUE(verbose)) {
                message(paste0("After trimming genes, summed counts matrix: [", 
                               dim(counts)[1], ", ", dim(counts)[2], "]"))
            }

            # keep track of sample and tissue 
            sample <- as.factor(groups)
        """ + ''.join(['\t' + key + ' <- ' + key + '\n' for key in numeric_keys]) + """
            if (isTRUE(verbose)) {
        """ + ''.join(['\t\tmessage("\n' + key + ' is numeric")\n'
                       for key in numeric_keys]) + """
            }
        """ + ''.join(['\t' + key + ' <- as.factor(' + key + ')\n' for key in nonnumeric_keys]) + """
            if (isTRUE(verbose)) {
        """ + ''.join(['\t\tmessage(paste0("\nLevel of ' + key + ': ", as.character(levels(' + key + '))))\n'
                       for key in nonnumeric_keys]) + """
            }

            # check for errors coming downstream due to too few levels
            smallest_model_levels <- min(c(""" + ', '.join(['length(levels(' + key + '))'
                                                            for key in nonnumeric_keys]) + """))
            if (smallest_model_levels < 2) {
                # this will cause an error in model.matrix()
                message(paste0("You are getting the following error from R limma ", 
                               "due to the fact that there are fewer than two levels ", 
                               "in one of the factors in your model.  I will mercifully ", 
                               "overlook the error and skip this test.\nError:\n", 
                               "Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) :",
                               " contrasts can be applied only to factors with 2 or more levels\n\n", 
                               "WARNING: skipping this test!"))
                return(c('aborted', 'aborted'))
            }

            # the model
            design <- model.matrix(""" + model + """)
            if (isTRUE(verbose)) {
                message("design <- model.matrix(""" + model + """)")
                message(paste0("Dimension of design matrix: [", 
                               dim(design)[1], ", ", dim(design)[2], "]"))
            }

            # check for degeneracy
            if (dim(design)[1] < dim(design)[2] + 1) {
                message(paste0("WARNING: experimental design is degenerate, since there are ",
                               "fewer rows than columns+1 in the design matrix.\n", 
                               "https://support.bioconductor.org/p/59168/ \n",
                               "Aborting this test!\n"))
                return(c('aborted', 'aborted'))
            }

            # deseq2 normalization
            sf <- estimateSizeFactorsForMatrix(counts)
            eff.lib <- sf * mean(colSums(counts))

            # voom
            if (isTRUE(verbose)) {
                message("Running """ + voom_function + """")
            }
            y <- DGEList(counts, lib.size=eff.lib)
            v.all <- """ + voom_function + """(y, design, span=""" + str(voom_lowess_span) + """, 
                                                plot=TRUE, save.plot=TRUE)

            # save components of the plot to pass to python
            data_voom <- data.frame(v.all$voom.xy, v.all$voom.line)
            voom_filename <- file.path(working_directory, "R_voom_plot.tsv")
            if (isTRUE(verbose)) {
                message(paste0("Writing voom plot information to ", voom_filename))
            }
            write.table(data_voom, file=voom_filename, append=FALSE,
                        row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

            # duplicateCorrelation, if called for
            if (isTRUE(use_duplicatecorr)) {
                if (isTRUE(verbose)) {
                    message("Running duplicateCorrelation, blocking on '""" + duplicatecorr_name + """'")
                }
            """ + corfit_command + """

                # as per the limma User Guide p. 126, we re-run voom and then duplicateCorrelation again
                if (isTRUE(verbose)) {
                    message("Re-running """ + voom_function + """, blocking on '""" + duplicatecorr_name + """'")
                }
                v.all <- """ + voom_function + """(y, design, span=""" + str(voom_lowess_span) + """, 
                                                    plot=TRUE, save.plot=TRUE, 
                                                    block=duplicatecorr, correlation=corfit$consensus)

                # save components of the plot to pass to python
                data_voom <- data.frame(v.all$voom.xy, v.all$voom.line)
                voom_filename <- file.path(working_directory, "R_voom_plot.tsv")
                if (isTRUE(verbose)) {
                    message(paste0("Writing voom plot information to ", voom_filename))
                }
                write.table(data_voom, file=voom_filename, append=FALSE,
                            row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

                # second round of duplicateCorrelation
                if (isTRUE(verbose)) {
                    message("Re-running duplicateCorrelation after voom, blocking on '""" + duplicatecorr_name + """'")
                }
            """ + corfit_command + """
            }

            # run limma DE
            if (isTRUE(verbose)) {
                message("Running lmFit")
            }
        """ + vfit_command + """

            # fit contrasts
            contr_matrix <- eval(parse(text='""" + make_contrast_command + """'))
            contrast_list <- colnames(contr_matrix)

            if (isTRUE(use_contrasts)) {
                if (isTRUE(verbose)) {
                    message("Fitting contrasts")
                }
                vfit <- contrasts.fit(vfit, contrasts=contr_matrix)

                contr_matrix_filename <- file.path(working_directory, "R_limma_contrasts.tsv")

                # write out the contrast matrix
                if (isTRUE(verbose)) {
                    message(paste0("Writing contrast matrix to ", contr_matrix_filename))
                }
                write.table(contr_matrix, file=contr_matrix_filename, append=FALSE,
                            row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
            } else {
                contr_matrix_filename <- 'not_applicable'
            }

            vfit <- eBayes(vfit, robust=TRUE)

            vres_filenames <- c()
            for (contrast in contrast_list) {

                # get the final output from limma with topTable
                if (isTRUE(use_contrasts)) {
                    vres <- topTable(vfit, coef=as.character(contrast), n=Inf, sort.by="P")
                } else {
                    vres <- topTable(vfit, coef=2, n=Inf, sort.by="P")  # is this right?
                }

                vres_filename <- file.path(working_directory, paste0(as.character(contrast), ".tsv"))
                vres_filenames <- c(vres_filenames, vres_filename)

                # write to disk: it's the only reliable way to transfer to python
                if (isTRUE(verbose)) {
                    message(paste0("Writing topTable to ", vres_filename))
                }
                write.table(vres, file=vres_filename, append=FALSE,
                            row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
            }

            return(c(voom_filename, contr_matrix_filename, vres_filenames))

        }
        """

        # create the R function
        Rfunc = SignatureTranslatedAnonymousPackage(Rfunctionstring, 'Rfunc')

        # pass a variable number of arguments as keyword arguments
        kwargs = dict()
        for key in model_keys:
            kwargs[key] = out[key]

        # TODO: filtering strategy...
        # only pre-filter things that are really zero
        # run limma to get "raw" p values
        # try several filtering thresholds, and for each:
        #     do BH correction based on the number of genes
        # plot the number of significant genes as a function of threshold
        # pick the min threshold that maximizes the number of significant genes

        # choose which genes to test (eliminate very low expression stuff)
        # use deseq2 normalization strategy
        # https://rdrr.io/bioc/DESeq2/src/R/core.R#sym-estimateSizeFactorsForMatrix
        counts = out['summed_counts'].copy()
        gene_exp_always_greater_than_zero = np.min(counts, axis=0, keepdims=False) > 0
        counts = counts[:, gene_exp_always_greater_than_zero]  # keep genes with no zeros
        log_counts = np.log2(counts)
        log_geometric_mean_per_gene = np.mean(log_counts, axis=0, keepdims=True)
        median_z_per_sample = np.median(log_counts - log_geometric_mean_per_gene, axis=1)
        scale_factors_per_sample = np.exp(median_z_per_sample)
        normalized_counts = out['summed_counts'] / np.expand_dims(scale_factors_per_sample, axis=1)

        # gene filtering logic
        include_gene_logic = (normalized_counts.mean(axis=0) >= low_expression_mean_threshold)

        # apply gene filter
        out['summed_counts'] = out['summed_counts'][:, include_gene_logic]
        gene_ids = gene_ids[include_gene_logic]

        # run differential expression testing in R limma
        if verbose:
            print('Running limma-voom in R...')

        outs = Rfunc.runLimma(summed_sample=out['summed_counts'],
                              gene_ids=gene_ids,
                              groups=out['groupby'],
                              model=model,
                              use_contrasts='TRUE' if use_contrasts else 'FALSE',
                              use_duplicatecorr='TRUE' if use_duplicatecorr else 'FALSE',
                              duplicatecorr=duplicatecorr if use_duplicatecorr else 'N/A',
                              working_directory=working_directory,
                              verbose='TRUE' if verbose else 'FALSE',
                              **kwargs)
        voom_file = outs[0]
        contrast_file = outs[1]
        vres_files = outs[2:]

        # chance to realize that R aborted because of degeneracy
        if outs[0] == 'aborted':
            continue

        # print the contrast matrix in jupyter
        if verbose and use_contrasts:
            print('Contrast matrix used by Limma:')
            contrast_df = pd.read_csv(contrast_file, sep='\t', index_col=0)
            display(HTML(contrast_df.to_html()))

        # display the voom plot
        if verbose:
            print('Voom plot:')
            voom_df = pd.read_csv(voom_file, sep='\t', index_col=0)
            plt.figure(figsize=(6, 6))
            plt.plot(voom_df['x'].values, voom_df['y'].values, '.', ms=1, color='black')
            plt.plot(voom_df['x.1'].values, voom_df['y.1'].values, '-', color='red')
            plt.xlabel('log2 ( count size + 0.5 )')
            plt.ylabel('Sqrt ( standard deviation )')
            plt.title('voom: Mean-variance trend')
            plt.show()

        # for each contrast, grab the output DE table
        for file in vres_files:

            # table from limma
            table = (pd.read_csv(file, sep='\t', index_col=0)
                     .reset_index()  # make gene its own column
                     .rename(columns={'index': 'gene'}))  # and name it

            # limit to specified FDR
            table = table[table['adj.P.Val'] < fdr]

            # add in cell type and cluster columns
            table.insert(loc=1,
                         column=(separate_groupings_key if (separate_groupings_key
                                                            is not None) else 'clusters.included'),
                         value=', '.join(group))

            # get gene count information
            if testing_key not in numeric_keys:

                # this is not well-defined if the testing key is not a categorical grouping

                # add in DE contrast information
                contrast = os.path.basename(os.path.splitext(file)[0])
                test_group = contrast.split('.')[0]
                comparison_group = contrast.split('.')[-1]
                table.insert(loc=2,
                             column='test.group',
                             value=test_group)
                table.insert(loc=3,
                             column='comparison',
                             value=comparison_group)

                testing_keys = adata_group.obs[testing_key].apply(lambda s: s.replace(" ", ""))  # match R
                if test_group not in testing_keys.unique():

                    # check if it's a compound contrast on the left-hand side of the equation
                    test_groups = test_group.split('_')
                    if not np.all([g in testing_keys.unique() for g in test_groups]):
                        print(f'WARNING: Cannot find a group called {test_group}, or groups '
                              f'called {test_groups} in adata, where the '
                              f'groups are {testing_keys.unique()}')
                else:
                    test_groups = [test_group]  # wrap in a list

                if comparison_group not in testing_keys.unique():

                    # check if it's a compound contrast on the left-hand side of the equation
                    comparison_groups = comparison_group.split('_')
                    if not np.all([g in testing_keys.unique() for g in comparison_groups]):
                        print(f'WARNING: Cannot find a group called {comparison_group}, or groups '
                              f'called {comparison_groups} in adata, where the '
                              f'groups are {testing_keys.unique()}')
                else:
                    comparison_groups = [comparison_group]  # wrap in a list

                # grab cells in the test group and in the comparison group
                test_group_logic = [k in test_groups for k in testing_keys]
                comparison_group_logic = [k in comparison_groups for k in testing_keys]
                adata_test = adata_group[test_group_logic]
                adata_comparison = adata_group[comparison_group_logic]

                # get mean count information, in test and comparison
                mean_test = np.array(adata_test.X.mean(axis=0)).squeeze()
                mean_comparison = np.array(adata_comparison.X.mean(axis=0)).squeeze()
                test_group_frac_expr0 = np.array((adata_test.X > 0).mean(axis=0)).squeeze()
                comparison_frac_expr0 = np.array((adata_comparison.X > 0).mean(axis=0)).squeeze()
                test_group_frac_expr1 = np.array((adata_test.X > 1).mean(axis=0)).squeeze()
                comparison_frac_expr1 = np.array((adata_comparison.X > 1).mean(axis=0)).squeeze()

                # add gene count information
                gene_count_df = pd.DataFrame(data={'gene': adata_group.var.index,
                                                   'test.group.cell.mean.counts': mean_test,
                                                   'comparison.cell.mean.counts': mean_comparison,
                                                   'test.group.frac.expr>0': test_group_frac_expr0,
                                                   'comparison.frac.expr>0': comparison_frac_expr0,
                                                   'test.group.frac.expr>1': test_group_frac_expr1,
                                                   'comparison.frac.expr>1': comparison_frac_expr1})
                table = pd.merge(left=table, right=gene_count_df, how='left', on=['gene'])

                if verbose:
                    print('Contrast ' + test_group + ' vs. ' + comparison_group)
                    print(f'number of cells in {test_group} = {adata_test.shape[0]}')
                    print(f'number of cells in {comparison_group} = {adata_comparison.shape[0]}')

                # clean up
                del adata_test
                del adata_comparison

            else:

                # numeric values in testing_key

                # add in DE contrast information
                contrast = os.path.basename(os.path.splitext(file)[0])
                test_group = contrast.split('.')[0]
                comparison_group = '(numeric)'  # just a placeholder so plotting functions work
                table.insert(loc=2,
                             column='test.group',
                             value=test_group)
                table.insert(loc=3,
                             column='comparison',
                             value=comparison_group)

                if verbose:
                    print('Effect: ' + test_group)

            # add table to list of dataframes
            dfs.append(table)
            group_dfs.append(table)

            if verbose:
                display(HTML(table.head().to_html()))

            # clean up
            gc.collect()

        # due to memory limitations, save these groupings TSVs
        tmp_filename = os.path.join(working_directory, 'grouping_' + "_".join(group) + '_DE.tsv')
        pd.concat(group_dfs).to_csv(tmp_filename, sep='\t', na_rep='-')
        print(f'Saved grouping {group} data to {tmp_filename}\n')

        # and try our best to free up memory
        del adata_group
        gc.collect()  # garbage collect in python
        robjects.r('rm(list=ls())')  # clear R memory
        robjects.r('gc()')  # garbage collect in R
        gc.collect()  # garbage collect in python again.  oh yes.

    if len(dfs) == 0:
        print(f'No tests gave any results below the specified FDR of {fdr}!')
        return None

    # concatenate into one massive table
    df = pd.concat(dfs)

    try:
        gene_id_df = pd.DataFrame(data={'gene': adata.var.index,
                                        'gene_id': adata.var[gene_id_key]})
        df = pd.merge(left=df, right=gene_id_df, how='left', on=['gene'])
    except KeyError:
        print('Warning: unable to add gene IDs, since they were not present in adata.var["gene_id"]')

    if calculate_background_estimate:
        try:
            df = calculate_bkg_prob(full_adata_for_background_estimate,
                                    testing_key=testing_key,
                                    df=df,
                                    separate_groupings_column=separate_groupings_key,
                                    separate_groupings_key=separate_groupings_key,
                                    bcs_included=adata.obs.index)
        except Exception:
            print('Warning: unable to run calculate_bkg_prob() successfully!  '
                  'You may want to try running it separately.  Returning the '
                  'dataframe without running calculate_bkg_prob()')

    else:
        if verbose:
            print('Background RNA probabilities have not been calculated since '
                  'calculate_background_estimate=False. If you want to try to '
                  'calculate bkg.prob on your own, it is recommended to try:\n\n'
                  'df = calculate_bkg_prob(full_adata_with_all_cells_included,\n'
                  '                        testing_key=testing_key,\n'
                  '                        df=df,  # this is from limma_voom_DE()\n'
                  '                        separate_groupings_column=separate_groupings_key,  '
                  '# remove this line if no separate_groupings_key used\n'
                  '                        separate_groupings_key=separate_groupings_key,  '
                  '# remove this line if no separate_groupings_key used\n'
                  '                        bcs_included=adata.obs.index)\n')

    return df


def volcano_plot(vres: pd.DataFrame,
                 cluster_column_label: str,
                 cluster_value: str,
                 test_group_value: str,
                 comparison_value: str,
                 test_group_column_label: str = 'test.group',
                 comparison_group_label: str = 'comparison',
                 FDR: float = 0.01,
                 min_AveExpr: Optional[float] = None,
                 num: int = 5,
                 label_genes_list: List[str] = [],
                 title: str = '',
                 xlim: Optional[List[float]] = None,
                 bkg_prob_vmin: float = 0.4,
                 bkg_prob_vmax: float = 0.6,
                 bkg_prob_label_threshold: float = 0.5,
                 fontsize: float = 12,
                 figsize: Tuple[float] = (8, 5),
                 rasterized: bool = True,
                 show: bool = True):
    """Create one volcano plot from data in a DE table."""

    assert len(vres) > 0, 'Input dataframe for volcano_plot() is empty'

    from adjustText import adjust_text

    # limit to cluster of interest
    vres = vres[vres[cluster_column_label] == cluster_value]

    # limit to test group of interest
    vres = vres[vres[test_group_column_label] == test_group_value]

    # limit to comparison group of interest
    vres = vres[vres[comparison_group_label] == comparison_value]

    # limit to certain minimum AveExpr
    if min_AveExpr is not None:
        vres = vres[vres['AveExpr'] >= min_AveExpr]

    # set gene to the index
    if 'gene' in vres.columns:
        vres = vres.set_index(keys='gene')

    assert len(vres) > 0, \
        'There are no entries left in the dataframe after limiting to cluster, ' \
        'test group, comparison group' + ('' if min_AveExpr is None else ', min_AveExpr')

    # find p-value threshold
    vres = vres.sort_values(by='P.Value')
    threshold_ind = np.where(vres['adj.P.Val'] > FDR)[0][0]
    threshold = -1 * np.log10(vres.iloc[threshold_ind]['P.Value'])

    plt.figure(figsize=figsize)

    plt.plot(vres['logFC'].values,
             vres['P.Value'].apply(lambda x: -1 * np.log10(x)).values,
             '.', ms=3, alpha=0.2, color='black', rasterized=rasterized)
    plt.plot([vres['logFC'].min() - 1, vres['logFC'].max() + 1],
             [threshold] * 2,
             ':', color='red')

    logic = (vres['P.Value'].apply(lambda x: -1 * np.log10(x)).values > threshold)

    if 'bkg.prob' in vres.columns:
        plt.scatter(vres['logFC'].values[logic],
                    vres['P.Value'].apply(lambda x: -1 * np.log10(x)).values[logic],
                    c=vres['bkg.prob'].values[logic],
                    cmap='cool', s=12, alpha=1.,
                    rasterized=rasterized,
                    vmin=bkg_prob_vmin, vmax=bkg_prob_vmax)
    else:
        plt.plot(vres['logFC'].values[logic],
                 vres['P.Value'].apply(lambda x: -1 * np.log10(x)).values[logic],
                 '.', ms=8, alpha=0.5, color='cyan',
                 rasterized=rasterized)
    plt.title(title)
    plt.xlabel('Effect size\n[log2 fold change]')
    plt.ylabel('Significance\n[-log10 p-value]')

    if 'bkg.prob' in vres.columns:
        # only label non-background
        gene_list = vres[logic & (vres['bkg.prob'] <= bkg_prob_label_threshold)]
    else:
        gene_list = vres[logic]
    gene_list_pos = gene_list[gene_list['logFC'] > 0]
    gene_list_neg = gene_list[gene_list['logFC'] < 0]

    # top num by logFC and top num by P.Value, on each side of logFC = 0
    gene_list = []
    gene_list_pos = gene_list_pos.sort_values(by='logFC', ascending=False)
    gene_list.extend(gene_list_pos.index[:num].tolist())

    gene_list_neg = gene_list_neg.sort_values(by='logFC', ascending=True)
    gene_list.extend(gene_list_neg.index[:num].tolist())

    gene_list_pos = gene_list_pos.sort_values(by='P.Value', ascending=True)
    gene_list.extend(gene_list_pos.index[:num].tolist())

    gene_list_neg = gene_list_neg.sort_values(by='P.Value', ascending=True)
    gene_list.extend(gene_list_neg.index[:num].tolist())

    gene_list = np.unique(gene_list).tolist()

    # check and make sure we don't have a gene showing up multiple times
    # because that will break the code
    for g in gene_list:
        if not isinstance(vres.index.get_loc(g), int):
            print('WARNING: gene ' + g + ' has multiple entries in the table.'
                  + '\nLabeling of this gene cannot proceed.  Please check '
                  + 'to ensure the table is alright.')

    texts = [plt.text(vres['logFC'][i],
                      vres['P.Value'].apply(lambda x: -1 * np.log10(x))[i],
                      vres.index[i], color='green', fontsize=fontsize)
             for i in [vres.index.get_loc(g) for g in gene_list]]

    for g in label_genes_list:
        try:
            i = vres.index.get_loc(g)
            texts.extend([plt.text(vres['logFC'][i],
                                   vres['P.Value'].apply(lambda x: -1 * np.log10(x))[i],
                                   vres.index[i], color='green', fontsize=fontsize)])
        except Exception:
            print(f'gene {g} not in the filtered data table')
            pass

    if xlim is not None:
        plt.xlim(xlim)

    adjust_text(texts,
                x=vres['logFC'].values,
                y=vres['P.Value'].apply(lambda x: -1 * np.log10(x)).values,
                expand_points=(1.1, 1.25),
                arrowprops=dict(arrowstyle='-', color='lightgray'))
    plt.grid(False)

    if show:
        plt.show()


def all_volcano_plots(df: pd.DataFrame,
                      cluster_column_name: str,
                      FDR: float,
                      test_group_column_name: str = 'test.group',
                      comparison_group_column_name: str = 'comparison',
                      label_genes_list: List[str] = [],
                      min_AveExpr: Optional[float] = None,
                      title_prefix: str = '',
                      **kwargs):
    """Produce all possible volcano plots, given a dataframe that contains DE output."""

    assert cluster_column_name in df.columns, \
        f'cluster_column_name "{cluster_column_name}" is not a column of the dataframe.'

    for k in np.unique(df[cluster_column_name].values):
        print(f'cluster {k}')

        for group in np.unique(df[df[cluster_column_name] == k][test_group_column_name].values):

            for comparison in np.unique(df[(df[cluster_column_name] == k)
                                           & (df[test_group_column_name]
                                              == group)][comparison_group_column_name].values):
                print(f'{group} vs. {comparison}')

                volcano_plot(df,
                             cluster_column_label=cluster_column_name,
                             cluster_value=k,
                             test_group_column_label=test_group_column_name,
                             test_group_value=group,
                             comparison_group_label=comparison_group_column_name,
                             comparison_value=comparison,
                             label_genes_list=label_genes_list,
                             FDR=FDR,
                             min_AveExpr=min_AveExpr,
                             title=f'{title_prefix} {group} vs. {comparison}\ncluster {k}: FDR {FDR}',
                             **kwargs)


def gsea(de_table: pd.DataFrame,
         gsea_gmt_file: str,
         gene_key: str = 'gene',
         ranking_key: str = 't',
         genome_mapping_tsv: Optional[str] = None,
         genome_mapping_mapfrom_column: Optional[str] = None,
         genome_mapping_human_genenames_column: Optional[str] = None,
         n_permutations: int = 100000,
         min_size: int = 0,
         max_size: int = 500,
         working_directory: str = '.',
         verbose: bool = True) -> pd.DataFrame:
    """Perform Gene Set Enrichment Analysis using the fgsea package in R.
    Can be used to perform GO enrichment analysis by choosing the appropriate
    GMT file, such as c5.all.v7.0.symbols.gmt.  Be aware that this version of
    GO enrichment does not account for relatedness of the GO terms, as GOstats
    does.

    NOTE: GSEA is run using human gene names.  If using a non-human species, a
    gene mapping file is required, along with column names that specify where
    the information is contained.

    For fgsea examples, see: https://stephenturner.github.io/deseq-to-fgsea/

    Args:
        de_table: Table that includes gene names and/or gene IDs as well as some
            numerical value that can be used to rank genes.  T-statistic
            probably works best.  This table can be, for example, the output of
            limma_voom_DE(), limited to a certain DE test of interest.
        gsea_gmt_file: GMT file for GSEA analysis, for example, a file like
            c2.cp.biocarta.v7.0.symbols.gmt
        gene_key: Column name for either gene names (human or non-human) or
            gene IDs (non-human only).
        ranking_key: Column name used to rank genes for GSEA algorithm.
        genome_mapping_tsv: TSV file that contains a mapping between species, if
            working with a non-human species.
        genome_mapping_mapfrom_column: If a genome_mapping_tsv is supplied, this
            contains the column name of either gene names or IDs for the current
            species.
        genome_mapping_human_genenames_column: If a genome_mapping_tsv is
            supplied, this contains the column name of human gene names.
        n_permutations: Number of permutations used for GSEA.  (I think) the
            smallest achievable p-value is ~ 1 / n_permutations.  This is the
            nperm parameter of fgsea.
        min_size: minSize parameter for fgsea.
        max_size: maxSize parameter for fgsea.
        working_directory: Directory where temporary files are saved.  They will
            be removed after (assuming there are no errors).
        verbose: True to print out intermediate information.

    """

    if verbose:
        from IPython.display import display, HTML

    # input checks
    assert ranking_key in de_table.columns, \
        f'ranking_key "{ranking_key}" is not a column name in the de_table.'
    if 'test.group' in de_table.columns:
        if de_table['test.group'].nunique() > 1:
            print('WARNING: multiple values in the "test.group" column of '
                  'de_table.  Did you intend to filter the DE table to one '
                  'comparison before running GSEA?')
    if 'comparison' in de_table.columns:
        if de_table['comparison'].nunique() > 1:
            print('WARNING: multiple values in the "comparison" column of '
                  'de_table.  Did you intend to filter the DE table to one '
                  'comparison before running GSEA?')
    working_directory = working_directory.rstrip('/')
    assert os.access(working_directory, os.W_OK), \
        f'Cannot write to working_directory:\n{working_directory}'
    assert os.access(gsea_gmt_file, os.R_OK), \
        f'Cannot find the pathway file:\n{gsea_gmt_file}'

    # map genes to human (if necessary)
    if genome_mapping_tsv is not None:
        translator = pd.read_csv(genome_mapping_tsv, sep='\t')

        lookup = dict(zip(translator[genome_mapping_mapfrom_column].values,
                          translator[genome_mapping_human_genenames_column].values))

    def translate(g):
        try:
            out = lookup[g]
        except KeyError:
            out = 'none'
            if verbose:
                print('Warning: ' + g + ' could not be translated to a human gene')
        return out

    # create ranking RNK file
    rnk = pd.DataFrame({'gene': de_table[gene_key].values,
                        't': de_table[ranking_key].values})
    if genome_mapping_tsv is not None:
        rnk['gene'] = rnk['gene'].apply(translate)
        rnk = rnk[rnk['gene'] != 'none']  # eliminate genes that couldn't be translated

    # save RNK file
    rnk.to_csv(os.path.join(working_directory, 'tmp.rnk'),
               sep='\t', index=False, header=['# gene', 't'])  # the pound sign is needed

    # import R packages
    from rpy2.robjects.packages import importr
    fgsea = importr('fgsea')
    tidyverse = importr('tidyverse')
    datatable = importr('data.table')

    # it's too hard to translate all the R to rpy2
    # so wrap up most of the R in a function
    from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

    Rfunctionstring = """
    runGSEA <- function(rnk.file,
                        gmt.file,
                        out.file,
                        minSize,
                        maxSize,
                        nperm) {

        ranks <- read.table(rnk.file, header=FALSE, colClasses = c("character", "numeric"))
        ranks <- deframe(ranks)  # needs tidyverse
        pathways <- gmtPathways(gmt.file)

        # run fgsea
        fgseaRes <- fgsea(pathways, ranks, minSize=minSize, maxSize=maxSize, nperm=nperm)

        # arrange the table
        fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))

        # save the output
        fwrite(fgseaResTidy, file = out.file)  # needs data.table

        return(file)

    }
    """

    # create the R function
    Rfunc = SignatureTranslatedAnonymousPackage(Rfunctionstring, 'Rfunc')

    # run fgsea
    if verbose:
        print('Running GSEA...')

    out_file = os.path.join(working_directory, 'tmp.gsea')
    Rfunc.runGSEA(os.path.join(working_directory, 'tmp.rnk'),
                  gsea_gmt_file,
                  out_file,
                  minSize=min_size,
                  maxSize=max_size,
                  nperm=n_permutations)

    # read in the output table created by fgsea
    df = pd.read_csv(out_file)

    try:
        # clean up temporary files
        os.system('rm ' + os.path.join(working_directory, 'tmp.rnk'))
        os.system('rm ' + out_file)
    except Exception:
        print('Unable to clean up temporary files')

    return df


def plot_discrete_scatter(adata: anndata.AnnData,
                          gene: str,
                          basis: str = 'umap',
                          layer: Optional[str] = None,
                          color_map: str = 'Oranges',
                          sort: bool = True,
                          show: bool = True,
                          **kwargs):
    """Display a discrete scatter plot.

    Args:
        adata: The AnnData object
        gene: Plot expression for this gene
        basis: Information for the 2D representation of gene expression is taken
            from adata.obsm['X_' + basis]
        layer: If you use layers in adata, specify the appropriate layer for counts
        color_map: A matplotlib colormap
        sort: True to plot the points with most counts on top
        show: True to show figure
        kwargs: Any arguments to matplotlib.pyplot.scatter(), the most useful one
            being vmax, which allows you to specify a max color value. Can also
            specify alpha and s (size of dots).

    Example:
        plot_discrete_scatter(adata_myocytes, basis='umap', gene='Hcn4',
                              layer='cellbender_0.01',
                              s=2, vmax=5, alpha=0.5)

    """

    assert gene in adata.var_names, f'Gene "{gene}" is not in adata.var_names'

    import matplotlib

    if layer is None:
        print('Data in adata.X should be integer counts for this discrete plot')
        counts = np.array(adata.X[:, adata.var_names == gene].todense()).squeeze()
    else:
        counts = np.array(adata.layers[layer][:, adata.var_names == gene].todense()).squeeze()

    if sort:
        order = np.argsort(counts)
    else:
        order = np.arange(counts.size)

    fig = plt.figure(figsize=(5.25, 5))
    if 'vmax' in kwargs.keys():
        vmax = kwargs['vmax']
    else:
        vmax = np.max(counts)
    cmap = plt.get_cmap(color_map, vmax + 1)

    norm = matplotlib.colors.BoundaryNorm(np.arange(-0.5, vmax + 1.5), cmap.N)
    scatter = plt.scatter(adata.obsm['X_' + basis][:, 0][order],
                          adata.obsm['X_' + basis][:, 1][order],
                          c=counts[order], cmap=cmap, norm=norm, **kwargs)
    plt.ylabel(basis.replace('_', ' ').upper() + ' 2')
    plt.xlabel(basis.replace('_', ' ').upper() + ' 1')
    plt.title(gene)
    plt.xticks([])
    plt.yticks([])
    plt.grid(False)

    cbar = fig.colorbar(scatter, ticks=np.arange(0, vmax + 1, np.ceil(vmax / 20)),
                        fraction=0.05, aspect=35)
    cbar.set_label('counts')
    cbar.set_alpha(1)
    cbar.draw_all()

    if show:
        plt.show()

    return fig


def gsea_from_DE(df: pd.DataFrame,
                 gmt_files: List[str],
                 bkg_prob_cutoff: Optional[float] = None,
                 bkg_prob_column: str = 'bkg.prob',
                 separate_groupings_column: str = 'clusters.included',
                 test_group_column: str = 'test.group',
                 comparison_column: str = 'comparison',
                 gene_key: str = 'gene_id',
                 ranking_key: str = 't',
                 genome_mapping_tsv: Optional[str] = None,
                 genome_mapping_mapfrom_column: str = 'rat.gene_id',
                 genome_mapping_human_genenames_column: str = 'human.gene_name',
                 n_permutations: int = 1000000,
                 working_directory: str = '.',
                 verbose: bool = True) -> pd.DataFrame:
    """Perform GSEA using a table created by limma_voom_DE(), or other DE table.

    Args:
        df: Results dataframe obtained from a differential expression analysis.
            Expects same outputs as produced by limma.
        gmt_files: List of paths to GMT tiles for running GSEA.  GSEA will be
            run for each GMT file.
        bkg_prob_cutoff: Subset DE results to non-background-contamination
            results. This is very important. This threshold will be used to
            subset the full data frame. If None, it will skip this subsetting.
        bkg_prob_column: Column that specifies background contamination probability.
        separate_groupings_column: Column that indicates different subgroupings
            of cells used to perform the indicated DE test.
        test_group_column: Column that specifies the test group for the DE test.
        comparison_column: Column that specifies the comparison group for the
            DE test.
        gene_key: Name of column to be used to specify the gene. For human data,
            use gene names. For other species, use Ensembl IDs and use a
            gene_mapping_tsv that links IDs to human gene names.
        ranking_key: Column to use for ranking DE results for GSEA.
        genome_mapping_tsv: If using non-human species, we need to link gene IDs
            to human gene names. This TSV files specifies how that is done.
        genome_mapping_mapfrom_column: Column in the genome_mapping_tsv file that
            specifies the non-human species' gene IDs.
        genome_mapping_human_genenames_column: Column in the genome_mapping_tsv
            file that specifies the human gene names.
        n_permutations: Used by fgsea, the larger this number, the longer it
            takes, and the smaller the achievable p-value.
        working_directory: Directory in which to write temporary files.
        verbose: True to print intermediate outputs.

    """

    assert separate_groupings_column in df.columns, \
        f'The input separate_groupings_column is not a column of the input dataframe'

    # subset the dataframe
    if bkg_prob_cutoff is not None:
        df = df[df[bkg_prob_column] < bkg_prob_cutoff]
    else:
        print('Warning: DE results that may be affected by background contamination '
              'should be removed. The input parameter bkg_prob_cutoff has been set '
              'to None, so this function is not subsetting the DE results. Ensure '
              'that the input dataframe is not contaminated with background DE results.')

    biglist = []

    for group in df[separate_groupings_column].unique():

        df_subset = df[df[separate_groupings_column] == group]

        for testgroup in df_subset[test_group_column].unique():

            for compgroup in df_subset[df_subset[test_group_column] == testgroup][comparison_column].unique():

                if verbose:
                    print(f'Within "{group}", working on {testgroup} vs. {compgroup}', end='')

                vres = df_subset[(df_subset[test_group_column] == testgroup)
                                 & (df_subset[comparison_column] == compgroup)]

                dfs = []

                for gmt in gmt_files:
                    dfs.append(gsea(de_table=vres,
                                    gsea_gmt_file=gmt,
                                    gene_key=gene_key,
                                    ranking_key=ranking_key,
                                    genome_mapping_tsv=genome_mapping_tsv,
                                    genome_mapping_mapfrom_column=genome_mapping_mapfrom_column,
                                    genome_mapping_human_genenames_column=genome_mapping_human_genenames_column,
                                    n_permutations=n_permutations,
                                    working_directory=working_directory,
                                    verbose=False)
                               )
                    if verbose:
                        print('.', end='')
                if verbose:
                    print('\n', end='')

                # concat
                gdf = pd.concat(dfs)
                gdf.insert(loc=0, column=comparison_column, value=compgroup)
                gdf.insert(loc=0, column=test_group_column, value=testgroup)
                gdf.insert(loc=0, column=separate_groupings_column, value=group)

                # keep track
                biglist.append(gdf.copy())

    gsea_df = pd.concat(biglist)

    return gsea_df


def plot_pathway_heatmap(gsea_df: pd.DataFrame,
                         fdr: float = 0.05,
                         max_per_group: Optional[int] = None,
                         sort_by: Optional[Union[str, List[str]]] = None,
                         test_group_order: Optional[List[str]] = None,
                         test_group_mapping: Optional[Dict[str, str]] = None,
                         test_group_column: str = 'test.group',
                         comparison_column: str = 'comparison',
                         p_value_column: str = 'padj',
                         effect_size_column: str = 'NES',
                         pathway_column: str = 'pathway',
                         color_map: str = 'PiYG',
                         test_name_fcn: Callable[[str], str] = lambda x: x,
                         p_value_style: str = 'both',
                         min_marker_size: float = 20,
                         max_marker_size: float = 100,
                         size_legend_values: List[float] = [1, 2, 4, 6],
                         colorbar_aspect: int = 10,
                         fontsize: int = 10,
                         marker: str = 's',
                         vmin: float = -2.5,
                         vmax: float = 2.5,
                         show: bool = True):
    """Plot a graphical visualization of the output of multiple GSEA analyses.

    The plot shows a grid of colored squares, with pathways on the y-axis and
    DE tests on the x-axis.  It only really makes sense to plot this for
    one-versus-all DE tests, where each condition on the x-axis is tested
    against all others.  The color is proportional to normalized effect size,
    and significant pathway-tests are given a black border.  The only pathways
    that are included on the y-axis are those that show at least one significant
    p-value in any test condition.

    Args:
        gsea_df: Dataframe that contains the output of GSEA runs, as from
            gsea_from_DE()
        fdr: False discovery rate for which GSEA results are considered significant
            for the purposes of this plot
        max_per_group: Maximum number of significant results to display per group.
            If None, all significant results with the specified fdr are shown.
        sort_by: Sort the pathways using this gsea_df column
        test_group_order: In case the x-axis values need to be in a specific order,
            that can be specified as a list here. Values must match values of
            gsea_df[test_group_column].unique()
        test_group_mapping: In case the x-axis values need to have different names
            than the values in gsea_df[test_group_column].unique(), the mapping
            from gsea_df[test_group_column].unique() to something else can be
            specified here as a dictionary
        test_group_column: Column that contains test group designations
        comparison_column: Column that contains comparison designations
        p_value_column: Column that contains adjusted p-values
        effect_size_column: Column that contains normalized effect sizes
        pathway_column: Column that contains pathway names
        color_map: Valid matplotlib colormap name
        test_name_fcn: Function to apply to pathway names prior to using as
            x-axis labels.  For example `lambda s: s.lower()` to make everything
            lowercase.
        p_value_style: Controls how the p-value is displayed: ['border', 'size', 'both'].
            'border' puts a black border around FDR significant results.
            'size' makes the marker sizes proportional to -log10(p-val) ** 2
        min_marker_size: The smallest marker is this size
        max_marker_size: The largest marker is this size
        size_legend_values: If not using 'border' style, then this series of values
            becomes the dot sizes shown in the p-value legend
        colorbar_aspect: Aspect ratio of the colorbar
        fontsize: Font size for axis labels
        marker: Marker for the plot: default is 's', a square box
        vmin: Lower limit for the effect size colorbar
        vmax: Upper limit for the effect size colorbar
        show: True to show plot, otherwise return figure without showing

    Returns:
        matplotlib figure object

    """

    # all the pathways that have significant p-values in any association

    tests = gsea_df[test_group_column].unique()

    # input checks
    assert min_marker_size > 0, 'min_marker_size must be > 0'
    assert p_value_style in ['border', 'size', 'both'], \
        f'p_value_style must be one of ["border", "size", "both"] but was "{p_value_style}"'
    if test_group_order is not None:
        for label in test_group_order:
            assert label in tests, f'test_group_order label "{label}" is not one of ' \
                f'the tests in the table: {tests}'
    for test in tests:
        unique_comparisons = gsea_df[gsea_df[test_group_column] == test][comparison_column].unique()
        assert len(unique_comparisons) == 1, \
            f'For test group {test}, there are multiple comparisons {unique_comparisons}.' \
                f'  This function is for one-versus-all comparisons, and does not make sense ' \
                f'for other types of tests.'
    if sort_by is not None:
        sort_by = list(sort_by)  # ensure it's a list
        for col in sort_by:
            assert col in gsea_df.columns, f'gsea_df.columns does not contain specified ' \
                f'sort_by "{col}": {gsea_df.columns}'
    else:
        sort_by = effect_size_column

    # order the xlabels
    if test_group_order is not None:
        test_labels = test_group_order
    else:
        # lexicographic order if no order is specified
        test_labels = np.sort(list(tests))

    # sort the dataframe
    if sort_by is not None:
        gsea_df = gsea_df.sort_values(by=sort_by)

    # subset DE table to relevant pathways
    gsea_df = gsea_df[[x in test_labels for x in gsea_df[test_group_column]]]
    dfs = []
    for test in test_labels:
        # this sorts pathways in the order the x-label tests will appear, by effect size
        tmp_df = gsea_df[(gsea_df[test_group_column] == test)
                         & (gsea_df[p_value_column] <= fdr)].copy().sort_values(by=sort_by,
                                                                                ascending=False)
        tmp_df = tmp_df[tmp_df[effect_size_column] > 0]  # without this, changing x order changes pathways
        dfs.append(tmp_df if max_per_group is None else tmp_df.head(max_per_group))
    pathway_list = pd.concat(dfs)[pathway_column].unique()

    if len(pathway_list) == 0:
        print('Warning: no pathways found to be significant!')
        return None

    if len(test_labels) == 1:
        # special case of one test: we need to infer the result of the opposite test
        gsea_df_sym = gsea_df.copy()
        gsea_df_sym[test_group_column] = gsea_df[comparison_column].copy()
        gsea_df_sym[comparison_column] = gsea_df[test_group_column].copy()
        gsea_df_sym[effect_size_column] = gsea_df[effect_size_column].apply(lambda x: -1 * x).copy()
        gsea_df = pd.concat([gsea_df, gsea_df_sym])
        test_labels = np.sort(list(gsea_df[test_group_column].unique()))

    # grab NES and pvalue info for those pathways in the celltypes tested
    mat = np.zeros((len(pathway_list), len(test_labels)))
    matsig = np.zeros((len(pathway_list), len(test_labels)))
    matp = np.zeros((len(pathway_list), len(test_labels)))

    for i, pathway in enumerate(pathway_list):
        for j, test in enumerate(test_labels):
            #             print(pathway, test)
            #             from IPython.display import display, HTML
            #             display(HTML(gsea_df[(gsea_df[pathway_column] == pathway)].to_html()))

            mat[i, j] = gsea_df[(gsea_df[pathway_column] == pathway)
                                & (gsea_df[test_group_column] == test)][effect_size_column].values.item()
            matsig[i, j] = 1 if (gsea_df[(gsea_df[pathway_column] == pathway)
                                         & (gsea_df[test_group_column] == test)][p_value_column].values.item()
                                 <= fdr) else 0
            matp[i, j] = -np.log10(gsea_df[(gsea_df[pathway_column] == pathway)
                                           & (gsea_df[test_group_column] == test)]
                                   [p_value_column].values.item())

    # make sure no mat values are actually zero
    mat[mat == 0] = 1e-50  # all points need to be caught by np.nonzero

    # figsize
    fig = plt.figure(figsize=(max(1., len(test_labels) / 2.), len(pathway_list) * 7 / 20.))

    # dot sizes vary depending on plotting style
    border_padding = np.power(max_marker_size / 7.2, 1.75)  # a heuristic
    if p_value_style == 'border':
        s = max_marker_size
    else:
        # max maps to max_marker_size
        conversion = lambda x: (min_marker_size
                                + ((x - np.square(matp).min()) / np.square(matp).max()
                                   * (max_marker_size - min_marker_size)))
        s = conversion(np.square(matp))

    # surrounding boxes for p values
    if p_value_style != 'size':
        plt.scatter(np.nonzero(matsig)[1], np.nonzero(matsig)[0],
                    c='white', s=max_marker_size + border_padding,
                    marker='s', edgecolors='black')

    # plot heatmap using NES values
    plt.scatter(np.nonzero(mat)[1], np.nonzero(mat)[0],
                c=mat[mat != 0], s=s, marker=marker, linewidths=0,  # edgecolors=None,
                cmap=color_map, alpha=0.6)

    # colorbar
    cb = plt.colorbar(aspect=colorbar_aspect, alpha=0.6, pad=0.2)
    cb.ax.tick_params(labelsize=fontsize)
    cb.set_ticks([-2, -1, 0, 1, 2])
    cb.ax.set_title('NES', fontsize=fontsize)
    plt.clim([vmin, vmax])

    ax = plt.gca()
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.spines['left'].set_color('white')

    plt.grid(False)
    plt.yticks(range(mat.shape[0]), pathway_list, fontsize=fontsize)

    # dot size legend
    if p_value_style != 'border':
        _dot_size_legend(sizes=size_legend_values,
                         display_scale_fcn=lambda x: conversion(np.square(x)),
                         marker=marker,
                         title='-log10(pval)')

    # optionally rename x-axis labels
    if test_group_mapping is not None:
        test_labels = [test_group_mapping[t] for t in test_labels]

    plt.xticks(ticks=range(max(2, len(test_labels))),
               labels=[test_name_fcn(t) for t in test_labels],
               fontsize=fontsize, rotation=90)
    plt.xlim([-0.5, max(2, len(test_labels)) - 0.5])
    plt.ylim([-0.5, max(2, len(pathway_list)) - 0.5])
    plt.title(f'FDR {fdr}')

    if show:
        plt.show()

    return fig


def _dot_size_legend(sizes: List[float],
                     marker: str,
                     title: str,
                     display_scale_fcn: Callable[[float], float] = lambda x: x,
                     color: str = 'k',
                     vertical: bool = True,
                     loc: str = 'center left',
                     bbox_to_anchor: Tuple[float] = (1, 0.5),
                     borderpad: float = 6,
                     labelspacing: Optional[float] = None,
                     fontsize: int = 12):
    """Simple way to create a legend for scatter plot dot sizes"""
    # https://blogs.oii.ox.ac.uk/bright/2014/08/12/point-size-legends-in-matplotlib-and-basemap-plots/

    ghosts = []
    for s in sizes:
        ghosts.append(plt.scatter([], [],
                                  s=display_scale_fcn(s),
                                  marker=marker,
                                  color=color,
                                  edgecolors='none'))
    labels = [str(s) for s in sizes]

    leg = plt.legend(ghosts, labels,
                     ncol=1 if vertical else len(sizes),
                     frameon=False, fontsize=fontsize,
                     handlelength=2, borderpad=borderpad,
                     loc=loc, bbox_to_anchor=bbox_to_anchor,
                     labelspacing=labelspacing,
                     handletextpad=1, title=title, scatterpoints=1)
    plt.setp(leg.get_title(), fontsize=fontsize)

    return leg


def plot_pathway_bars(gsea_df: pd.DataFrame,
                      fdr: float = 0.05,
                      test_groups: Optional[List[str]] = None,
                      colors: Optional[List[str]] = None,
                      n_pathways: int = 3,
                      sort_first_by: str = 'effect',
                      test_group_column: str = 'test.group',
                      p_value_column: str = 'padj',
                      effect_size_column: str = 'NES',
                      pathway_column: str = 'pathway',
                      enriched_only: bool = True,
                      labels_on_top: bool = True,
                      figsize: Tuple[float] = (30, 2),
                      fontsize: int = 16,
                      asterisk_size: float = 8,
                      pathway_name_fcn: Callable[[str], str] = lambda x: x,
                      show: bool = True):
    """Plot a bar graph of several GSEA pathway effect sizes for each cell cluster.
    Significant pathways are marked with an asterisk.
    Args:
        gsea_df: Dataframe that contains the output of GSEA runs, as from
            gsea_from_DE()
        fdr: False discovery rate used to mark significant pathways
        test_groups: Names of test group DE test GSEA results to plot, in order.
            If None, it defaults to gsea_df[test_group_column].unique()
        colors: Colors, as hex strings, in order.  Number should be number of
            test groups.  As long as colors[i] is a legitimate color, it doesn't
            have to be a hex string.
        n_pathways: Number of pathways to show per test_group
        sort_first_by: One of ['effect', 'pval']. Bars will be sorted first by
            this attribute, then by the other.
        test_group_column: Name of column with test group designations
        p_value_column: Column that contains adjusted p-values
        effect_size_column: Column that contains normalized effect sizes
        pathway_column: Column that contains pathway names
        enriched_only: True to plot only gene sets enriched in the test group.
            False to plot n_pathways enriched gene sets followed by
            n_pathways depleted gene sets (i.e., NES < 0).
        labels_on_top: True to show pathway labels on top axis
        figsize: Figure size
        fontsize: Size of labels
        asterisk_size: Size of asterisks that denote significant pathways
        pathway_name_fcn: Function to apply to pathway names prior to using as
            x-axis labels.  For example `lambda s: s.lower()` to make everything
            lowercase.
        show: True to show figure
    """

    assert sort_first_by in ['effect', 'pval'], f'sort_first_by must be in ["effect", "pval"]' \
        f' but was "{sort_first_by}"'

    # alter the location of labels
    if labels_on_top:
        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
        plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True

    fig = plt.figure(figsize=figsize)

    pathways = []
    significant = []

    if test_groups is None:
        test_groups = gsea_df[test_group_column].unique()
        print(f'x-axis groupings are, in order: {test_groups}')
    if colors is None:
        colors = ['black'] * len(test_groups)

    for i, test in enumerate(test_groups):

        if test not in gsea_df[test_group_column].unique():
            print(f'Warning: {test} not in gsea_df["{test_group_column}"]. Skipping...')
            continue

        if enriched_only:
            sub_df = gsea_df[(gsea_df[test_group_column] == test)
                             & (gsea_df[effect_size_column].apply(lambda x: x > 0))]

            # sort
            sub_df['logp'] = sub_df[p_value_column].apply(lambda x: -1 * np.log10(x))
            if sort_first_by == 'effect':
                sub_df = sub_df.sort_values(by=[effect_size_column, 'logp'], ascending=False)
            elif sort_first_by == 'pval':
                sub_df = sub_df.sort_values(by=['logp', effect_size_column], ascending=False)
            else:
                raise ValueError(f'sort_first_by is {sort_first_by}')

            plt.gca().bar(np.arange(n_pathways) + (n_pathways + 1) * i,
                          sub_df[effect_size_column][:n_pathways],
                          color=colors[i])
            y_vals = np.array(sub_df[effect_size_column][:n_pathways]) + 0.2
            y_vals[np.array(sub_df[p_value_column][:n_pathways]) > fdr] = None
            significant.extend(y_vals.tolist() + [None])
            pathways.extend(list(sub_df[pathway_column][:n_pathways].values) + [' '])

            plt.xticks(ticks=range((n_pathways + 1) * (len(test_groups) + 1) - 1),
                       labels=[pathway_name_fcn(p) for p in pathways],
                       rotation=90, fontsize=fontsize)
            plt.xlim([-1, (n_pathways + 1) * len(test_groups) - 1])

        else:

            # NES > 0
            sub_df = gsea_df[(gsea_df[test_group_column] == test)
                             & (gsea_df[effect_size_column].apply(lambda x: x > 0))]

            # sort
            sub_df['logp'] = sub_df[p_value_column].apply(lambda x: -1 * np.log10(x))
            if sort_first_by == 'effect':
                sub_df = sub_df.sort_values(by=[effect_size_column, 'logp'], ascending=False)
            elif sort_first_by == 'pval':
                sub_df = sub_df.sort_values(by=['logp', effect_size_column], ascending=False)
            else:
                raise ValueError(f'sort_first_by is {sort_first_by}')

            plt.gca().bar(np.arange(n_pathways) + (2 * n_pathways + 1) * i,
                          sub_df[effect_size_column][:n_pathways],
                          color=colors[i])
            y_vals = np.array(sub_df[effect_size_column][:n_pathways]) + 0.2
            y_vals[np.array(sub_df[p_value_column][:n_pathways]) > fdr] = None
            significant.extend(y_vals.tolist())
            pathways.extend(list(sub_df[pathway_column][:n_pathways].values))

            # NES < 0
            sub_df = gsea_df[(gsea_df[test_group_column] == test)
                             & (gsea_df[effect_size_column].apply(lambda x: x < 0))]

            # sort
            sub_df['logp'] = sub_df[p_value_column].apply(lambda x: -1 * np.log10(x))
            if sort_first_by == 'effect':
                sub_df = sub_df.sort_values(by=[effect_size_column, 'logp'], ascending=True)
            elif sort_first_by == 'pval':
                sub_df = sub_df.sort_values(by=['logp', effect_size_column], ascending=False)
            else:
                raise ValueError(f'sort_first_by is {sort_first_by}')

            plt.gca().bar(np.arange(n_pathways) + (2 * n_pathways + 1) * i + n_pathways,
                          sub_df[effect_size_column][:n_pathways],
                          color=colors[i])
            y_vals = np.array(sub_df[effect_size_column][:n_pathways]) - 0.2
            y_vals[np.array(sub_df[p_value_column][:n_pathways]) > fdr] = None
            significant.extend(y_vals.tolist() + [None])
            pathways.extend(list(sub_df[pathway_column][:n_pathways].values) + [' '])

            plt.plot([-1, (2 * n_pathways + 1) * (len(test_groups) + 1) - 1], [0, 0],
                     color='gray', alpha=0.1)

            plt.xticks(ticks=range((2 * n_pathways + 1) * (len(test_groups) + 1) - 1),
                       labels=[pathway_name_fcn(p) for p in pathways],
                       rotation=90, fontsize=fontsize)
            plt.xlim([-1, (2 * n_pathways + 1) * len(test_groups) - 1])

    plt.plot(significant, 'k*', ms=asterisk_size)
    plt.ylabel('NES', fontsize=fontsize)
    plt.yticks(rotation=90, fontsize=fontsize)
    plt.grid(False)

    # return the location of labels to the usual places
    if labels_on_top:
        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
        plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False

    if show:
        plt.show()

    return fig


def plot_DE_matrix(df: pd.DataFrame,
                   label_order: Optional[List[str]] = None,
                   fun: Callable[[pd.DataFrame], float] =
                   lambda d: d['logFC'].apply(lambda x: np.abs(x) > 2).sum(),
                   colorbar_label: str = 'Genes with |logFC| > 2',
                   bkg_max: float = 0.4,
                   leave_off: bool = True,
                   grid: bool = True,
                   vmin: Optional[float] = None,
                   vmax: Optional[float] = None,
                   label_fcn: Callable[[str], str] = lambda s: s,
                   cmap: str = 'Greys',
                   show: bool = True):
    """Display the results of many one vs. one DE tests as a heat matrix.

    Matrix shows number of DE genes detected above various thresholds for each
    comparisons.

    "One vs. one" comparison means that we are not looking at one vs. all.  For
    example, for chambers {LA, LV, RA, RV}, this graphic would use the results
    of tests for RA.vs.RV, LA.vs.RA, etc., but not RA.vs.LA_RV_LV, etc.

    Args:
        df: DE result dataframe (as from limma_voom_DE), containing one vs. one
            comparisons
        label_order: Order of labels (df['test.group'].unique() values)
        fun: Function to apply to dataframe to calculate a value for plotting.
            Input a lambda function, for example `lambda d: d['logFC'].std()`
        colorbar_label: Explain `fun` with this label for the color values
        bkg_max: Max bkg.prob of interest
        leave_off: True to leave off axis labels that are not used
        grid: True to show a grid around comparisons the plot is displaying
        vmin: Min color axis limit
        vmax: Max color axis limit
        label_fcn: Function to apply to test.group values before displaying
        cmap: Matplotlib colormap name
        show: True to show plot

    Returns:
        matplotlib figure object

    """

    if label_order is None:
        tests = np.unique(df['test.group'].unique().tolist()
                          + df['comparison'].unique().tolist())
        label_order = tests
    else:
        tests = label_order

    for label in label_order:
        assert label in (df['test.group'].unique().tolist()
                         + df['comparison'].unique().tolist()), \
            f'Label "{label}" from label_order is not in ' \
            f'df["test.group"].unique(): {df["test.group"].unique()} or ' \
            f'df["comparison"].unique(): {df["comparison"].unique()}'

    mat = np.zeros((len(tests), len(tests)))
    mat[:] = np.nan

    for i, test in enumerate(tests):

        for j in range(i):

            # Subset to the relevant comparison by creating two dataframes
            df_subset = df[(((df['test.group'] == test) & (df['comparison'] == tests[j]))
                            | ((df['test.group'] == tests[j]) & (df['comparison'] == test)))
                           & (df['bkg.prob'] <= bkg_max)]

            if len(df_subset) == 0:
                print(f'WARNING: missing comparison {test} and {tests[j]}... '
                      f'maybe not enough cells in condition?')
                continue

            # Calculate the desired value.
            mat[i, j] = fun(df_subset)

    fig = plt.figure(figsize=(4, 4))
    plt.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax)
    if leave_off:
        plt.xticks(np.arange(len(label_order) - 1),
                   labels=[label_fcn(x) for x in label_order[:-1]], rotation=90)
        plt.yticks(np.arange(len(label_order) - 1) + 1,
                   labels=[label_fcn(x) for x in label_order[1:]])
    else:
        plt.xticks(np.arange(len(label_order)),
                   labels=[label_fcn(x) for x in label_order], rotation=90)
        plt.yticks(np.arange(len(label_order)),
                   labels=[label_fcn(x) for x in label_order])

    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.ax.get_yaxis().labelpad = 25
    cbar.ax.set_ylabel(colorbar_label, rotation=270)
    plt.grid(False)
    plt.ylim([len(label_order) - 0.5, -0.5])
    plt.xlim([-0.5, len(label_order) - 0.5])

    if grid:
        for x in np.arange(0.5, len(label_order), 1):
            plt.plot([x] * 2, [x, len(label_order)], color='black', alpha=0.5, lw=0.5)
            plt.plot([-0.5, x], [x] * 2, color='black', alpha=0.5, lw=0.5)

    if show:
        plt.show()

    return fig


def markers_from_DE(df: pd.DataFrame,
                    n_max: int,
                    logFC_min: float = 1.5,
                    pval_max: float = 0.01,
                    bkg_prob_max: float = 0.5,
                    test_group_mean_count_min: float = 0.2,
                    test_group_frac_expr_min: float = 0.1,
                    sort_fcn: Union[str, Callable[[pd.DataFrame], pd.Series]] \
                            = lambda x: x['logFC'] * x['adj.P.Val'].apply(lambda x: -np.log10(x)),
                    test_group_key_order: Optional[Iterable[str]] = None,
                    test_group_name_fcn: Callable[[str], str] = lambda s: s):
    """Get a dictionary of marker genes by prioritizing genes from a DE table.

    NOTE: Subset the input DE table to one-versus-all comparisons.

    Args:
        df: Dataframe from a differential expression test.
        n_max: Max number of genes per group.
        logFC_min: Exclude genes with logFC below this value.
        pval_max: Exclude genes with adj.P.Val above this value.
        bkg_prob_max: Exclude genes with bkg.prob >= this value.
        test_group_mean_count_min: Exclude genes with test.group.cell.mean.counts
            below this value.
        test_group_frac_expr_min: Exclude genes with test.group.frac.expr>0
            below this value.
        sort_fcn: Name of a pre-defined sort function, or an arbitrary lambda
            function that operates on the dataframe and returns a series of
            floats, with genes being prioritized by highest values.
        test_group_key_order: Order of test.group keys in the output dict.
        test_group_name_fcn: Rename test.group values (keys in the output dict)
            according to this function.

    """

    required_columns = ['logFC', 'adj.P.Val', 'test.group.cell.mean.counts',
                        'test.group.frac.expr>0', 'test.group.frac.expr>1',
                        'comparison.frac.expr>0', 'comparison.frac.expr>1']
    assert all([c in df.columns for c in required_columns]),\
        f'The dataframe needs columns {required_columns}'

    unique_markers = {}

    if 'bkg.prob' not in df.columns:
        print('WARNING: no "bkg.prob" column in the input dataframe... '
              'will skip the bkg.prob cutoff condition.')
        df['bkg.prob'] = 0

    if test_group_key_order is None:
        test_group_key_order = df['test.group'].unique()

    # pre-defined sorting functions
    if type(sort_fcn) == str:

        # set up some quantities
        frac_in_key = 'test.group.frac.expr>0'
        frac_out_key = 'comparison.frac.expr>0'
        if '_1' in sort_fcn:
            frac_in_key = 'test.group.frac.expr>1'
            frac_out_key = 'comparison.frac.expr>1'
        tp = lambda x: x[frac_in_key]
        fp = lambda x: x[frac_out_key]
        precision = lambda x: tp(x) / (tp(x) + fp(x))
        recall = lambda x: tp(x)
        f = lambda x, beta: ((1 + beta ** 2) * precision(x) * recall(x)
                             / (beta ** 2 * precision(x) + recall(x)))

        if (sort_fcn == 'ppv') or (sort_fcn == 'ppv_1'):
            sort_fcn = precision
        elif (sort_fcn == 'f1') or (sort_fcn == 'f1_1'):
            sort_fcn = partial(f, beta=1.)
        elif (sort_fcn == 'f0.5') or (sort_fcn == 'f0.5_1'):
            sort_fcn = partial(f, beta=0.5)
        elif (sort_fcn == 'f2') or (sort_fcn == 'f2_1'):
            sort_fcn = partial(f, beta=2.)
        elif sort_fcn == 'logFC':
            sort_fcn = lambda x: x['logFC']
        else:
            raise ValueError('sort_fcn is implemented for '
                             '["logFC", "ppv", "ppv_1", "f1", "f1_1", '
                             '"f0.5", "f0.5_1", "f2", "f2_1"].  '
                             'Custom functions can be input as lambda '
                             'functions that operate on the dataframe.')

    for k in test_group_key_order:

        # subset dataframe to genes of interest
        genes = df[(df['test.group'] == k)
                   & (df['adj.P.Val'] <= pval_max)
                   & (df['logFC'] >= logFC_min)
                   & (df['bkg.prob'] < bkg_prob_max)
                   & (df['test.group.cell.mean.counts'] >= test_group_mean_count_min)
                   & (df['test.group.frac.expr>0'] >= test_group_frac_expr_min)]

        # warn if input dataframe has multiple comparisons
        if genes['comparison'].nunique() > 1:
            print('WARNING: Dataframe needs to be subset to one-versus-all comparisons '
                  'for "marker gene" determination to make sense. There are multiple '
                  f'"comparison"s {genes["comparison"].unique()} for "test.group" {k}')

        # ensure genes are more highly-expressed in a different cluster (yes this is needed)
        genes_all_tests = df[df['gene'].apply(lambda g: g in genes['gene'].values)]
        max_mean_groups = (genes_all_tests[['gene', 'logFC', 'test.group.cell.mean.counts']]
                           .groupby('gene').max().reset_index())
        max_mean_groups.rename(columns={'logFC': 'max.logFC',
                                        'test.group.cell.mean.counts': 'max.counts'},
                               inplace=True)
        genes = pd.merge(left=genes,
                         right=max_mean_groups[['gene', 'max.logFC', 'max.counts']],
                         how='left', on='gene')

        # eliminate genes with higher logFC in other tests, or higher expression elsewhere
        genes = genes[np.logical_not(genes['max.logFC'] > genes['logFC'])
                      & np.logical_not(genes['max.counts'] > genes['test.group.cell.mean.counts'])]

        if len(genes) == 0:
            continue

        # prioritize genes
        genes['sort'] = sort_fcn(genes)
        genes = genes.sort_values(by='sort', ascending=False)['gene'].values

        # add genes to list, checking for uniqueness (probably redundant)
        i = 0
        genelist = []
        for gene in genes:
            if i < n_max:
                if gene not in [a for b in list(unique_markers.values()) for a in b]:
                    genelist.append(gene)
                    i = i + 1

        if len(genelist) > 0:
            unique_markers.update({test_group_name_fcn(k): genelist})

    return unique_markers
