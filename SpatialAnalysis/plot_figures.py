import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

if not os.path.exists("ED_Figure_8"):
    os.mkdir('ED_Figure_8')

df_segment = pd.read_csv("data/wta_cta.csv")
df_segment.dropna(inplace=True)
df_segment.reset_index(drop=True, inplace=True)
df_segment['CTA_ID'] = df_segment['CTA_ID'].apply(lambda s: s.replace('-', '.'))
df_segment['WTA_ID'] = df_segment['WTA_ID'].apply(lambda s: s.replace('-', '.'))

count_wta = pd.read_csv("data/Broad-COVID_WTA_Q3Norm_TargetCountMatrix.txt", sep='\t')
count_cta = pd.read_csv("data/Broad-COVID_CTA_Q3Norm_TargetCountMatrix.txt", sep='\t')

wta_genes = count_wta['Gene'].values
cta_genes = count_cta['Gene'].values
common_genes = np.intersect1d(wta_genes, cta_genes)
common_genes = np.delete(common_genes, np.where(common_genes=='SARS-CoV-2 Neg')[0])
print(f"Consider {common_genes.size} genes in common.")

filt_count_wta = count_wta.loc[np.isin(count_wta['Gene'], common_genes),].copy()
filt_count_wta.set_index('Gene', inplace=True)
filt_count_cta = count_cta.loc[np.isin(count_cta['Gene'], common_genes),].copy()
filt_count_cta.set_index('Gene', inplace=True)

cta_seg_idx = []
wta_seg_idx = []
for idx, row in df_segment.iterrows():
    if row['WTA_ID'] in filt_count_wta.columns:
        wta_seg_idx.append(idx)
    if row['CTA_ID'] in filt_count_cta.columns:
        cta_seg_idx.append(idx)

# AOIs with data in both CTA and WTA.
filt_segments = df_segment.iloc[np.intersect1d(np.array(cta_seg_idx), np.array(wta_seg_idx)),].copy()
filt_segments['patient'] = filt_segments['Slide'].apply(lambda s: s.split('-')[0])
print(f"Consider {filt_segments.shape[0]} AOIs with data in both CTA and WTA.")

bwhID2broadID = {'C01': 'D22', 'C02': 'D23', 'C03': 'D24', 'S01': 'D18', 'S02': 'D19', 'S03': 'D20', 'S09': 'D21', 'S10': 'D8', 'S11': 'D9', 'S16': 'D10', 'S18': 'D11', 'S28': 'D12'}
filt_segments['broad_ID'] = filt_segments['patient'].apply(lambda s: bwhID2broadID[s])

distance_dict = {
    'S01': '40mm', 'S02': '40mm', 'C03': '40mm',
    'C01': '20mm', 'C02': '20mm', 'S03': '20mm',
    'S09': '10mm', 'S10': '10mm', 'S11': '10mm', 'S16': '10mm', 'S18': '10mm', 'S28': '10mm',
}
filt_segments['distance'] = pd.Categorical(filt_segments['patient'].apply(lambda s: distance_dict[s]), ordered=True, categories=['10mm', '20mm', '40mm'])
filt_segments['broad_ID'] = pd.Categorical(filt_segments['broad_ID'], ordered=True, categories=['D21', 'D8', 'D9', 'D10', 'D11', 'D12', 'D22', 'D23', 'D20', 'D18', 'D19', 'D24'])

# Gene correlation
df_corr = pd.DataFrame({'gene': common_genes})
df_corr['corr'] = 0.0
for _, row in df_corr.iterrows():
    gene = row['gene']
    vec_cta = filt_count_cta.loc[gene, filt_segments['CTA_ID'].values]
    vec_wta = filt_count_wta.loc[gene, filt_segments['WTA_ID'].values]
    df_corr.loc[df_corr['gene']==gene, 'corr'] = pearsonr(vec_cta, vec_wta)[0]
df_corr.sort_values('corr', ascending=False, inplace=True)
df_corr.reset_index(drop=True, inplace=True)

fig = plt.figure(dpi=500)
ax = fig.subplots()
ax = sns.violinplot(y='corr', data=df_corr, cut=0)
plt.title("Gene expression in WTA and CTA data")
ax.set_xlabel('Common genes')
ax.set_ylabel('Pearson correlation')
fig.savefig('gene_correlation.pdf')

# AOI correlation
assert np.sum(filt_count_cta.index!=filt_count_wta.index) == 0
aoi_corr = []
for _, row in filt_segments.iterrows():
    vec_cta = filt_count_cta[[row['CTA_ID']]].values.flatten()
    vec_wta = filt_count_wta[[row['WTA_ID']]].values.flatten()
    cor, _ = pearsonr(vec_cta, vec_wta)
    if np.isnan(cor):
        print(f"AOI with WTA_ID {row['WTA_ID']} has at least one contant gene expression.")
    aoi_corr.append(cor)

filt_segments['corr'] = aoi_corr
filt_segments.dropna(inplace=True)
print(f"Consider {filt_segments.shape[0]} AOIs after removing those with NaN correlation.")

fig = plt.figure(dpi=500)
ax = fig.subplots()
ax = sns.violinplot(y='corr', data=filt_segments, cut=0)
ax.set_xlabel('AOIs')
ax.set_ylabel('Pearson Correlation')
fig.savefig('ED_Figure_8/ED_Figure_8f.pdf')

# Get whisker data
q1, q3 = np.percentile(filt_segments['corr'], [25, 75])
fig = plt.figure(dpi=500)
ax = plt.boxplot(x='corr', data=filt_segments)
whisker_list = [item.get_ydata() for item in ax['whiskers']]
lower_whisker = whisker_list[0][1]
upper_whisker = whisker_list[1][1]
fig.savefig('aoi_corr_boxplot.pdf')
print("AOI correlation numeric summary:")
print(f"  Minimum: {filt_segments['corr'].min()}")
print(f"  Lower whisker: {lower_whisker}")
print(f"  Q1: {q1}")
print(f"  Median: {filt_segments['corr'].median()}")
print(f"  Q3: {q3}")
print(f"  Upper whisker: {upper_whisker}")
print(f"  Maximum: {filt_segments['corr'].max()}")


# Strip plot.
fig = plt.figure(dpi = 500)
ax = fig.subplots()
ax = sns.stripplot(x='broad_ID', y='corr', hue='distance', data=filt_segments, size=3)
ax.set_xlabel('Broad ID')
ax.set_ylabel('Pearson Correlation of AOI')
fig.savefig('ED_Figure_8/ED_Figure_8g.pdf')
