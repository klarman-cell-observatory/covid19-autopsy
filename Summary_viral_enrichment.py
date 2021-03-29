# Contents of this file:
# * Functions
# 1) Enrichment scores - general UMAP
# 2) Enrichment per patient
# 3) Enrichment in sub-clusters
#
# Author: Rachelly Normand, March 2021

import anndata as an
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import scanpy as sc
import os
import numpy as np
from numpy import Inf
import random
import re
import seaborn as sns
from statsmodels.stats.multitest import fdrcorrection


################################################################
#                       Functions
################################################################

def viral_enrichment_pipeline(data, cluster_col, title, save_str, patient_specific, all_or_subclustering):
    viral_enrichment_score(data, cluster_col, "Viral+")
    viral_enrichment_umap(data, enrichment_score="Enrichment_"+enrichment_name, save_str=save_str, title= title, all_or_subclustering=all_or_subclustering)
    FDRs = viral_enrichment_significance(data, cluster_col, enrichment_name, save_str=save_str, perm_no=1000)
    bar_infected_per_cluster(data = data, cluster_col=cluster_col, FDR_per_cluster=FDRs,title=title, save_str=save_str)
    if not patient_specific:
        bar_infected_per_cluster_stacked(data, cluster_col, save_str=save_str)



# This function computes viral enrichment score (inspired by Bost et. al.)
# data - should include the columns clusters_name and "Viral_counts" in its obs object, that contains the total
#        number of viral UMIs per cell
# clusters_name is the column name that contains the relevant clusters name by which the score should be grouped
# infected_threshold - the number of viral UMIs a cell should have to be counted as infected
#
# The enrichment score the function computes is = log(observed/expected), where:
#    observed = are the number of infected cells per cluster (cells with vUMI>9)
#    expected = the number of total infected cells * cell frequency per cluster
#
# Returns: the data object with the score column per cell in its obs object, named by "type"
def viral_enrichment_score(data, cluster_col, viral_count_col, suffix=""):
    viral_count_df = data.obs[[cluster_col,viral_count_col]].copy()

    grouped_clust = viral_count_df.groupby(data.obs[cluster_col])
    n_perClust = grouped_clust.size()
    n_total = sum(n_perClust)

    obs_perClust = grouped_clust.sum()[viral_count_col]
    total_infCells = sum(obs_perClust)
    exp_perClust = (n_perClust / n_total) * total_infCells

    enrichment_perClust = pd.DataFrame(np.log(0.0001+obs_perClust/exp_perClust))
    enrichment_perClust[enrichment_perClust==-Inf]=0
    new_col_name = "Enrichment_"+suffix
    enrichment_perClust = enrichment_perClust.rename(columns={0:new_col_name})

    data.obs = data.obs.merge(enrichment_perClust, left_on=cluster_col, right_index= True, how="left")


# Computing the significance of given enrichment scores by creating 100 permutations of the given counts per cell and computing a p-value by comparing the real enrichment score to those observed in a distribution of enrichment scores per the computer permutations
# data - should include the following columns in obs:
#           - <clusters_col> - that contains the grouping onto clusters
#           - 'Viral_counts' - that contains the total number of viral UMIs per cell
# clusters_col - the column in data.obs in which the cluster labels are saved
# tit_prefix - the title prefix for the file that will save the ann_data with the permutation
# perm_no - number of permutations
def viral_enrichment_significance(data, clusters_col, save_str, perm_no):
    perm_file = save_str+"_viral_perm.h5ad"
    fdr_file = save_str + "_FDRs.csv"

    if not os.path.exists(perm_file):
        print("The file "+perm_file+" doesn't exist.")
        print("(A) Shuffling viral data...")
        # A. Shuffle the viral UMI counts across cells 100 times and compute the enrichment per cluster for each permutation
        cell_no = data.shape[0]
        data_perm = data.copy()
        data_perm.obs.reset_index(inplace=True)
        for i in range(0,perm_no):
            rand_indx = random.sample(range(0,cell_no), cell_no)
            data_perm.obs["ViralPerm_" + str(i)] = data_perm.obs.loc[rand_indx,"Viral+"].reset_index(drop=True)
            viral_enrichment_score(data_perm, clusters_col, viral_count_col="ViralPerm_" + str(i), suffix="_"+str(i))

        data_perm.write_h5ad(perm_file)
    else:
        if not os.path.exists(fdr_file):
            print("(A) Reading existing viral-shuffeling file: "+perm_file)
            data_perm = an.read_h5ad(perm_file)

    # B. Distribution -> p-val
    if not os.path.exists(fdr_file):
        print("(B) Computing FDR from distribution of enrichment scores...")
        # Collecting enrichment score per permutation
        curr_cols = [clusters_col]
        col_name = "ViralEnrichment"
        curr_cols.extend([col for col in data_perm.obs.columns if re.search(col_name+"_", col)])
        df_for_summary = data_perm.obs[curr_cols]
        enrich_df = pd.DataFrame(index=data_perm.obs[clusters_col].drop_duplicates().sort_values())
        for i in range(0,perm_no):
            x = df_for_summary[["Enrichment_"+str(i), clusters_col]].drop_duplicates()
            enrich_df = enrich_df.merge(x, left_index=True, right_on=clusters_col).set_index(clusters_col)
        pd.DataFrame.transpose(enrich_df).hist(bins=15)

        # comparing real enrichments per cluster to permutations
        real_enrichments = data.obs[[clusters_col, "Enrichment_"]].drop_duplicates().set_index(clusters_col)
        FDR= {}
        for clust in set(data.obs[clusters_col]):
            FDR[clust] = sum(float(real_enrichments.loc[clust]) < enrich_df.loc[clust,:]) / perm_no
        FDR_df = pd.DataFrame.from_dict(FDR, orient="index").reset_index().rename(columns={"index":"clusters",0:"FDR"})
        FDR_df.to_csv(save_str+"_FDRs.csv")
    else:
        print("(B) Reading existing FDR file: " + fdr_file)
        FDR_df = pd.read_csv(fdr_file, header=0, index_col=[0])

    adj_FDR = fdrcorrection(FDR_df["FDR"])
    FDR_df["adj_FDR"] = adj_FDR[1]
    FDR_df.to_csv(save_str + "_FDRs.csv")

    return FDR_df


# Plots a umap of the data with the given enrichment score colored
# data - an anndata object
# enrichment_score - this column must be present in data.obs.columns
def viral_enrichment_umap(adata, enrichment_score, save_str="", title="", all_or_subclustering="all"):
    # setting colors
    adata = adata.copy()
    data = adata.obs[enrichment_score]
    num_levels = 20
    vmin, midpoint, vmax = data.min(), 0, data.max()
    levels = np.linspace(vmin, vmax, num_levels)
    midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
    vals = np.interp(midp, [vmin, midpoint, vmax], [0, 0.5, 1])
    colors = plt.cm.coolwarm(vals)
    cmap, norm = matplotlib.colors.from_levels_and_colors(levels, colors)
    adata.obs["Viral+"] = pd.Categorical(adata.obs["Viral+"])
    plt.clf()
    sc.set_figure_params(figsize=(5, 5))
    fig, ax = plt.subplots(1, 1)
    sc.pl.umap(adata, color=enrichment_score, color_map=cmap, title=title +" "+ enrichment_score, show=False, ax=ax)
    size=23
    if all_or_subclustering!="all":
        size=50
    sc.pl.umap(adata[adata.obs["Viral+"]==True,:], color="Viral+", palette=["black"], size=size, show=False, ax=ax, legend_loc=None, title=title)
    plt.tight_layout()
    if save_str == "":
        plt.show()
    else:
        plt.savefig("figures/"+save_str+"_viral_umap.pdf")



# Bar plots of infected cells per cluster
# data must contain "Viral_counts" in obs and cluster_col which contains the cluster groups
# cluster_col - the column that annotated the clusters, by which to group the summary
# FDR_per_cluster - a dataframe stating the FDR value per cluster, will be used for coloring the bar plots
def bar_infected_per_cluster(data, cluster_col, FDR_per_cluster, title, save_str=""):
    grouped_clust = data.obs.groupby(data.obs[cluster_col])
    infClust = pd.DataFrame(grouped_clust["Viral+"].sum()).reset_index()
    infClust.columns=["clusters","Infected_cell_count"]
    infClust=infClust.merge(FDR_per_cluster, left_on="clusters", right_on="clusters")

    plt.clf()
    ax = sns.barplot(data=infClust, x="clusters", y="Infected_cell_count", palette = cm.viridis(infClust['adj_FDR']), dodge=False)
    plt.colorbar(cm.ScalarMappable(cmap=cm.viridis))
    ax.grid(False)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    # Adding asterik in significant bars
    for p, sig in zip(ax.patches, infClust["adj_FDR"]):
        if sig <=0.05:
            x = p.get_x()
            w = p.get_width()
            y = p.get_height()
            if pd.isna(y):
                y=0
            ax.text(x + w / 2., y, '*', ha='center')
    plt.title(title, size=16)
    plt.xticks(fontsize=8)  # myeloid
    plt.yticks(fontsize=8)  # myeloid
    plt.tight_layout()

    if save_str=="":
        plt.show()
    else:
        plt.savefig("figures/"+save_str+"_bars.pdf")

# bar plot of viral+ cells per cluster - stacked by patient.
def bar_infected_per_cluster_stacked(data, cluster_col, save_str=""):
    grouped_clust = data.obs.groupby([data.obs[cluster_col], data.obs["patient"]])
    infClust = pd.DataFrame(grouped_clust["Viral+"].sum()).reset_index()
    infClust.columns=["clusters","patient","Infected_cell_count"]
    infClust.loc[infClust["Infected_cell_count"].isna(), "Infected_cell_count"]=0

    infClust_dense = infClust.pivot_table(index = "clusters", columns="patient", values="Infected_cell_count")
    plt.figure(figsize=(7, 15)) #Myeloid, for All->(15,15)
    infClust_dense.plot.bar(stacked=True)
    plt.grid(None)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5),fontsize=8)
    plt.xticks(fontsize=8) #myeloid
    plt.yticks(fontsize=8)  # myeloid
    plt.tight_layout()
    if save_str != "":
        plt.savefig("figures/" + save_str + "_stacked_bar.pdf")
















########################################
# 1) Enrichment scores - general umap
########################################
all_data = an.read_h5ad("LungData.h5ad") # Use the full file of data available in the single Cell Portal: https://singlecell.broadinstitute.org/single_cell/study/SCP1052/covid-19-lung-autopsy-samples



# Focusing on the 7 patients with viral+ cells
#all_data.obs[["Viral+","patient"]].value_counts()
focus_patients = ["P054921","P079042","P348762","P230638","P852049","P166169","P334354"]
data = all_data[all_data.obs["patient"].isin(focus_patients),:].copy()

viral_enrichment_pipeline(data=data, cluster_col="Cluster", title="All", save_str="All_7patients", patient_specific=False, all_or_subclustering="all")

########################################
# 2) Enrichment per patient
########################################
# Re-read the data to avoid duplicate enrichment scores
for p in focus_patients:
    print("***"+p+"+***")
    curr_data = all_data[all_data.obs["patient"]==p,:].copy()
    viral_enrichment_pipeline(data=curr_data, cluster_col="Cluster", title=p, save_str=p, patient_specific=True, all_or_subclustering="all")




####################################
# 3) Enrichment in sub-clusters
####################################
# Myeloid
#============
Mac_cells = all_data[all_data.obs["Cluster"]=="Myeloid"]
Mac_cells_focused = Mac_cells[Mac_cells.obs["patient"].isin(focus_patients),:].copy()
viral_enrichment_pipeline(data=Mac_cells_focused, cluster_col="ClusterName", title="Myeloid 7 patients", save_str="Myeloid_7patients", patient_specific=False,all_or_subclustering="subcluster")

# Per patient
# Re-read the data to avoid duplicate enrichment scores
for p in focus_patients:
    curr_data = Mac_cells[Mac_cells.obs["patient"]==p,:].copy()
    viral_enrichment_pipeline(data=curr_data, cluster_col="ClusterName", title=p+" Myeloid", save_str=p+"myeloid", patient_specific=True, all_or_subclustering="subcluster")





# Endothelial
#============
Endo_cells = all_data[all_data.obs["Cluster"]=="Endothelail"]
Endo_cells_focused = Endo_cells[Endo_cells.obs["patient"].isin(focus_patients),:].copy()
viral_enrichment_pipeline(data=Endo_cells_focused, cluster_col="ClusterName", title="Endothelial 7 patients", save_str="Endo_7patients", patient_specific=False, all_or_subclustering="subcluster")

for p in focus_patients:
    curr_data = Endo_cells[Endo_cells.obs["patient"]==p,:].copy()
    viral_enrichment_pipeline(data=curr_data, cluster_col="ClusterName", title=p+" Endothelial", save_str=p+"endothelial", patient_specific=True, all_or_subclustering="subcluster")


