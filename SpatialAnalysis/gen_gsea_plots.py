import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

figid2name = {
    # For WTA S vs C
    'Figure_4c_gsea': 'WTA PanCK+ alveolar',
    'ED_Figure_8j_gsea': 'WTA PanCK- alveolar',
    # For CTA S vs C
    'ED_Figure_8h_gsea': 'CTA PanCK+ alveolar',
    'ED_Figure_8i_gsea': 'CTA PanCK- alveolar',
}

def translate_pathway(old_name):
    s_list = list(map(lambda s: s.capitalize(), old_name.split('_')[1:]))
    return ' '.join(s_list)

def plot_gsea(df, dirname, fig_id, gene_set, suffix):
    outname = '_'.join(fig_id.split('_')[:-1] + ['right'])
    df_plot = df.sort_values(['Log Q', 'NES Abs'], ascending=False)
    fig = plt.figure()
    ax = sns.barplot(x='Log Q', y='pathway', data=df_plot, color='green')
    ax.set_xlabel('-log10(q-value)')
    plt.title(f"{figid2name[fig_id]} {suffix}")
    plt.tight_layout()
    fig.savefig(f"{dirname}/{outname}.{gene_set}.{suffix}.png", dpi=500)

def gen_gsea(filename):
    dir_name = os.path.dirname(filename)
    data_name = os.path.basename(filename).split('.')[0]
    gene_set = os.path.basename(filename).split('.')[1]
    df = pd.read_csv(filename)
    df = df.loc[df['padj']<0.05].copy()
    df['pathway'] = df['pathway'].apply(translate_pathway)
    df['Log Q'] = -np.log10(df['padj'])
    df['NES Abs'] = np.abs(df['NES'])

    df_up = df.loc[df['NES']>0]
    plot_gsea(df_up, dir_name, data_name, gene_set, 'up')

    df_down = df.loc[df['NES']<0]
    plot_gsea(df_down, dir_name, data_name, gene_set, 'down')


if __name__ == '__main__':
    gen_gsea("Figure_4/Figure_4c_gsea.H.csv")
    gen_gsea("ED_Figure_8/ED_Figure_8j_gsea.H.csv")
    gen_gsea("ED_Figure_8/ED_Figure_8h_gsea.H.csv")
    gen_gsea("ED_Figure_8/ED_Figure_8i_gsea.H.csv")
