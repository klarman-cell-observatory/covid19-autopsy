# Summary of manual clustering and sub-clustering.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (1) Clustering main lineages (ED figure 2i, ED figure 4a)
# (2) Sub-clustering - Myeloid, endothelial, T+NK, B+Plasma (ED figure 2j-m)
# (3) Sub-clustering stats (ED figure 4c)
# (4) Creating heatmaps of top genes per sub-clustering (for each lineage) (ED figure 4b)
# (5) Comparing UMAPs manual and predicted (ED figure 4d,e)

import os
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt
import pegasus as pg
import scanpy as sc
import anndata as an
import pandas as pd
import re
import seaborn as sb
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import subprocess
import scipy.io
import pickle




####################################################
# (1) Clustering main lineages (ED figure 2i, ED fig4a)
####################################################
def signature_score_per_cell(data, gene_set) :
    # Get rid of genes that aren't in data
    orig_len = len(gene_set)
    gene_set = [gene for gene in gene_set if gene in data.var_names]
    print(str(len(gene_set)) + "/" + str(orig_len)  + " of the gene set genes are measured" )

    # Limit the data to just those genes
    dat = data[:,gene_set].X
    dat = dat.toarray()
    mean = dat.mean(axis=0)
    var = dat.var(axis=0)
    std = np.sqrt(var)

    with np.errstate(divide="ignore", invalid="ignore"):
        dat = (dat - mean) / std
    dat[dat < -5] = -5
    dat[dat > 5] = 5

    scores = dat.mean(axis = 1)
    return(scores)

def plot_zscore_signature(markers, data, pdf, pop_name):
    sc.set_figure_params(figsize=(10, 8))
    gene_set = [gene for gene in markers if gene in data.var_names]
    data.obs["zscore_" + pop_name] = signature_score_per_cell(data, gene_set)
    gene_list = "".join([x + "\n" for x in gene_set])
    fig = sc.pl.umap(data, color="zscore_"+pop_name, cmap="Purples", size=40, show=False, return_fig=True)
    plt.suptitle(pop_name+" Z-score")
    plt.annotate(gene_list, xy=(1.2, 0), xycoords=('axes fraction', 'axes fraction'))
    pdf.savefig(fig, bbox_inches="tight")


data = an.read_h5ad("LungData.h5ad") # Use the full file of data available in the single Cell Portal: https://singlecell.broadinstitute.org/single_cell/study/SCP1052/covid-19-lung-autopsy-samples
General_immune = ["PTPRC"]
Myeloid = ["C1QA","C1QB","C1QC", "MRC1","CD163","MARCO","HLA-DPB1","HLA-DRB1", "HLA-DRB5", "HLA-DQB1","HLA-DQA1","HLA-DMA","HLA-DMB","HLA-DRA","LYZ","FCER1G","FCER2A","FCER2A"]
MAST=["CPA3"]
B_cells = ["CD19","MS4A1","MZB1","CD79B","JCHAIN"]
T_NK_cells= ["TRAC","TRBC1","TRBC2","TRDC","TRGC2", "CD8A","CD8B","CD3G","CD3E","CD3D","GZMA","GZMB","GZMH","GZMK","GZMM", "GNLY","PFR1","CTLA4","KLRB1","KLRC1","NKG7"]

lin_list = pd.Series(data = [General_immune, Myeloid, T_NK_cells, B_cells, MAST], index = ["General_immune", "Myeloid", "T_NK_cells", "B_cells", "MAST"])

with PdfPages("manual_annot_main_lineages.pdf") as pdf:
    for lin_name in lin_list.index:
        plot_zscore_signature(lin_list[lin_name], data, pdf, lin_name)







####################################################
# (2) Sub-clustering - Myeloid, endothelial, T+NK, B+Plasma (ED figure 2j-m)
####################################################
# Filters markers made by Pegasus for those with a given AUC or larger
# Considers only "up" markers
# feature should be set to True in case the gene names are in "feature" column rather than in the index
def filter_markers(marker_dict, roc_per, feature=False):
    markers_per_clust = {}
    for clust in marker_dict.keys():
        curr = marker_dict[str(clust)]["up"]
        if feature:
            markers_per_clust[str(clust)] = curr[curr.auroc > roc_per]["feature"]
        else:
            markers_per_clust[str(clust)] = curr[curr.auroc>roc_per].index.values
        print("Cluster " + str(clust) + " : " + str(len(markers_per_clust[str(clust)])))
    return(markers_per_clust)




def rand_index(adata, title):
    if not os.path.exists('RandInd_dictionaries'):
        os.makedirs('RandInd_dictionaries')

    resamp_perc = 0.9
    adata = adata.copy()
    indx_array = adata.obs.index.values
    n_cells = range(adata.shape[0])
    resamp_size = round(adata.shape[0] * resamp_perc)

    for res in [0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.1]:#, 1.3,1.5, 1.7, 1.9]:
        print(res)
        rand_indx_dict = {}
        pg.neighbors(adata,rep="pca_harmony")
        pg.leiden(adata, rep="pca_harmony", resolution = res)
        rand_list = []
        for iter in range(20) :
            samp_indx = random.sample(n_cells, resamp_size)
            samp_indx = indx_array[samp_indx]
            samp_data = adata[samp_indx]
            true_class = samp_data.obs["leiden_labels"]

            pg.neighbors(samp_data, rep="pca_harmony")
            pg.leiden(samp_data, rep = "pca_harmony", resolution = res)
            new_class = samp_data.obs["leiden_labels"]

            rand_list.append(adjusted_rand_score(true_class, new_class))

        rand_indx_dict[str(res)] = rand_list
        file_name = "RandInd_dictionaries/Dict_"+ title +"_"+str(res)+".pckl"
        filehandler = open(file_name,"wb")
        pickle.dump(rand_indx_dict, filehandler)
        filehandler.close()



def rand_index_summary(outdir ="RandInd_dictionaries", title=""):
    final_dict = {}
    for res in [0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.1]:#, 1.3,1.5, 1.7, 1.9]:
        final_dict[res] = pickle.load(file=open(outdir + "/Dict_" + title +"_"+ str(res) + ".pckl", 'rb'))[str(res)]

    plot_df = pd.DataFrame(final_dict).T
    plot_df = plot_df.reset_index()
    plot_df = plot_df.melt(id_vars="index")
    ax = sb.boxplot(x="index", y="value", data=plot_df)
    ax.set(xlabel='Leiden resolution', ylabel='Adjusted rand_index')
    plt.title(title, fontsize= 20)
    plt.show()


# Runs the 1st step of sub-clustering:
# Pre-processing of the partial data: qc_metrics, HVG, PCA, harmony, neighbors, umap
# adata - the partial data being sub-clustered. This data should be copied in the
#         slicing process, as this function will update fields in the adata object
# hvg_no - number of highly variable genes to consider for PCS
# PCs_no - number of PCs to compute for the rest of the analysis
# tit - title for the plot
# ** adata should have a column "sample" in obs, by which the data will be harmonized
def subcluster_1_preprocess(data, hvg_no, PCs_no, tit, in_place=False):
    #adata.obsm = None # Removing former PCA and umap values
    data.uns.clear()  # To eliminate problems with size of uns (specifically "fmat_highly_variable_features")

    if not in_place:
        adata=data.copy()
    else:
        adata = data

    # Preprocessing
    pg.qc_metrics(adata, min_genes=200, min_umis=400)
    pg.highly_variable_features(adata, consider_batch=False, n_top=hvg_no)
    pg.pca(adata, n_components=PCs_no)
    adata.obs['Channel'] = adata.obs['sample']
    pg.run_harmony(adata)
    pg.neighbors(adata, rep="pca_harmony")
    pg.umap(adata, rep="pca_harmony")
    pg.leiden(adata, rep="pca_harmony", resolution=0.5)
    sc.pl.umap(adata, size=15, title=tit + " (#HVG: " + str(hvg_no) + ", #PCs: " + str(PCs_no) + ")",
               color="leiden_labels")

# Runs the 2nd step of sub-clustering:
# Run rand-index on resolutions 0.1-1.3 and plots a summary plot that allows picking
# the optimal clustering resolution
# Before rand_index the pre-processing is run once more with the in_place mode
def subcluster_2_rand_index(adata, dict_tit, hvg_no, PCs_no):
    subcluster_1_preprocess(adata, hvg_no, PCs_no, dict_tit, in_place=True)

    rand_index(adata, dict_tit)
    rand_index_summary(outdir = "RandInd_dictionaries", title= dict_tit)

# Runs the 3rd step of sub-clustering:
# 1) Clusters the data with the given resolution
# 2) Plots QC plots
# 3) Computes DE genes (saves the full list in an excel file and the filtered list as a dictionary)
# ** adata should have a column "patient" in obs
def subcluster_3_cluster_QC_DE(adata, res, tit):
    # 1)clustering
    leiden_tit = "leiden_res_"+str(res)
    pg.leiden(adata, rep="pca_harmony", resolution= res, class_label=leiden_tit)

    # # 2) QC plots
    sc.set_figure_params(figsize=(10,10), fontsize=30)
    sc.pl.umap(adata, color=[leiden_tit, "method"], size=35, show=False)
    plt.suptitle(tit, fontsize=50)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()

    sc.set_figure_params(figsize=(16, 12), fontsize=30)
    sc.pl.umap(adata, color=["patient", "sample"], size=95, show=False)
    plt.suptitle(tit, fontsize=50)
    plt.show()

    sc.pl.umap(adata, size=95, color=['n_genes', 'n_counts', 'percent_mito'], show=False)
    plt.suptitle(tit, fontsize=50)
    plt.show()
    composition_barplot(adata, leiden_tit, "sample", tit+" - Sample count per cluster")
    composition_barplot(adata, leiden_tit, "patient", tit+" - Patient count per cluster")

    # 3) DE genes
    pg.de_analysis(adata, leiden_tit)
    markers = pg.markers(adata)
    pg.write_results_to_excel(markers, "DE_genes_"+ tit +".xlsx")
    filt_markers = filter_markers(markers, 0.75, False)
    filehandler = open("FiltMarkers"+tit+".pckl", "wb")
    pickle.dump(filt_markers, filehandler)
    filehandler.close()
    return filt_markers



# The an_data object should include the following columns in the obs part: "condition","Patient"
# an_data - the AnnData object for which a psuedo-bulk is needed
# out_dir - the directory in which the ourput files for the R execution will be saved
# prefix - prefix for the file names
# with_patient_cor - a boolean value indicating wither the patient should be considered as a covariate
# an_data.obs should contian the columns "leiden" and "patient" that contains the groups of cell we are comparing and the patients IDs which will be considered as covariates if with_patient_cor is True
def psuedo_bulk_wrapper(an_data, out_dir, prefix, with_patient_cor=True):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    os.chdir(out_dir)

    # 1. Create files for R code
    #===========================
    # count matrix
    scipy.io.mmwrite(os.path.join(prefix+"_data.mtx"), an_data.raw.X)
    # meta data
    meta = an_data.obs[["leiden","patient"]]
    meta.to_csv(os.path.join(prefix+"_meta.csv"))
    # gene list:
    pd.DataFrame(an_data.var_names).to_csv(prefix+"_genes.csv", index=False)

    # 2. Run pseudo-bulk in R
    #===========================
    # Build subprocess command
    # [<command>, <path to script> + <arguments>
    args = [prefix+"_data.mtx", prefix+"_genes.csv", prefix+"_meta.csv", prefix, str(with_patient_cor)]
    cmd = ['Rscript', 'Summary_pseudobulk.R'] + args

    # check_output will run the command and store to result
    x = subprocess.check_output(cmd, universal_newlines=True)
    print(x)
    os.chdir("..")



# Defining main lineages by marker expression
Mac_cells = data[data.obs["leiden_res_1.3"].isin(["3","7","18","15"]),:].copy()
T_NK_cells = data[data.obs["leiden_res_1.3"].isin(["5"]),:].copy()
B_Plasma_cells = data[data.obs["leiden_res_1.3"].isin(["12"]),:].copy()
Endo_cells = data[data.obs["leiden_res_1.3"].isin(["2","8","15"]),:].copy()



# Finding optimal #PCs and #HVG per sub-population
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for hvg_no in [500,1000,1500, 2000]:
    for PCs_no in [10,15,25,50]:
        print("*** hvg_no:" + str(hvg_no) + ", PCs_no: " + str(PCs_no) + " ***")
        subcluster_1_preprocess(Mac_cells, hvg_no, PCs_no, "Myeloid")

for hvg_no in [500,1000,1500, 2000]:
    for PCs_no in [10,15,25,50]:
        print("*** hvg_no:" + str(hvg_no) + ", PCs_no: " + str(PCs_no) + " ***")
        subcluster_1_preprocess(T_NK_cells, hvg_no, PCs_no, "T-NK")

for hvg_no in [500,1000,1500, 2000]:
    for PCs_no in [4,6,8]:
        print("*** hvg_no:" + str(hvg_no) + ", PCs_no: " + str(PCs_no) + " ***")
        subcluster_1_preprocess(B_Plasma_cells, hvg_no, PCs_no, "B-Plasma")

for hvg_no in [500,1000,1500, 2000]:
    for PCs_no in [10,15,25,50]:
        print("*** hvg_no:" + str(hvg_no) + ", PCs_no: " + str(PCs_no) + " ***")
        subcluster_1_preprocess(Endo_cells, hvg_no, PCs_no, "Endo_cells")


# Making rand index plots for all subsets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subcluster_2_rand_index(Mac_cells, "myeloid", 2000, 10)
subcluster_2_rand_index(T_NK_cells, "t_nk", 2000, 10)
subcluster_2_rand_index(B_Plasma_cells, "b_plasma", 1000, 4)
subcluster_2_rand_index(Endo_cells, "endo", 500, 15)


# Computing DE genes by AUC and pseudo-bulk
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subcluster_3_cluster_QC_DE(Mac_cells, res=0.9, tit="Myeloid") # I tried a range of leiden resolutions
Mac_cells.write_h5ad("Myeloid.h5ad")
Mac_cells.obs["leiden"] =Mac_cells.obs["leiden_res_0.9"]
psuedo_bulk_wrapper(an_data = Mac_cells, out_dir= "pseudo_myeloid_0.9", prefix="Myeloid_0.9",with_patient_cor=True)

subcluster_3_cluster_QC_DE(T_NK_cells, res=0.7, tit="T-NK") # I tried a range of leiden resolutions
T_NK_cells.write_h5ad("T_NK.h5ad")
T_NK_cells.obs["leiden"] =T_NK_cells.obs["leiden_res_0.7"]
psuedo_bulk_wrapper(an_data = T_NK_cells, out_dir= "pseudo_T_NK_0.7", prefix="T_NK_0.7",with_patient_cor=True)

subcluster_3_cluster_QC_DE(B_Plasma_cells, res=0.3, tit="B_Plasma") # I tried a range of leiden resolutions
B_Plasma_cells.write_h5ad("B_Plasma.h5ad")
B_Plasma_cells.obs["leiden"] =B_Plasma_cells.obs["leiden_res_0.7"]
psuedo_bulk_wrapper(an_data = B_Plasma_cells, out_dir= "pseudo_B_Plasma_0.3", prefix="B_Plasma_0.3",with_patient_cor=True)

subcluster_3_cluster_QC_DE(Endo_cells, res=0.5, tit="B_Plasma") # I tried a range of leiden resolutions
Endo_cells.write_h5ad("Endothelial.h5ad")
Endo_cells.obs["leiden"] =Endo_cells.obs["leiden_res_0.5"]
psuedo_bulk_wrapper(an_data = Endo_cells, out_dir= "pseudo_Endo_0.5", prefix="Endo_0.5",with_patient_cor=True)



#######################################################################
# (3) Sub-clustering stats (ED figure 4c)
#######################################################################

# Plots a composition bar plot of the frequency of a given attribute per cluster
def composition_barplot(adata, xattr, yattr, title, save_str="", cmap=None, fig_size=None):
    if fig_size==None:
        figsize = (10, 7)
    else:
        figsize = fig_size
    df = pd.crosstab(adata.obs.loc[:, xattr], adata.obs.loc[:, yattr])
    df = df.div(df.sum(axis=1), axis=0) * 100.0
    if (cmap==None):
        ax = df.plot(kind = "bar", stacked = True, figsize=figsize,legend=True, grid=False)
    else:
        ax = df.plot(kind="bar", stacked=True, figsize=figsize, legend=True, grid=False, cmap=cmap)
    df.to_excel("DataForFigures/ED_Fig4c"+save_str+".xlsx")
    ax.figure.subplots_adjust(right=0.9)
    ax.set_title(title,fontsize= 25)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    #plt.rc('xtick',labelsize=18)
    plt.xticks(fontsize=15)
    plt.tight_layout()
    if save_str=="":
        plt.show()
    else:
        plt.savefig(save_str+".pdf")

Mac_cells = an.read_h5ad("Myeloid.h5ad")
T_NK_cells = an.read_h5ad("T_NK.h5ad")
B_Plasma_cells = an.read_h5ad("B_Plasma.h5ad")
Endo_cells = an.read_h5ad("Endothelial.h5ad")

composition_barplot(Mac_cells, "ClusterName", "Sample_ID", "Myeloid - sample per cluster",save_str="Myeloid_stats_sample")
composition_barplot(Mac_cells, "ClusterName", "method", "Myeloid - method per cluster",save_str="Myeloid_stats_method")
composition_barplot(T_NK_cells, "ClusterName", "Sample_ID", "T+NK - sample per cluster",save_str="T_NK_stats_sample")
composition_barplot(T_NK_cells, "ClusterName", "method", "T+NK - method per cluster",save_str="T_NK_stats_method")
composition_barplot(B_Plasma_cells, "ClusterName", "Sample_ID", "B+Plasma - sample per cluster",save_str="B_Plasma_stats_sample")
composition_barplot(B_Plasma_cells, "ClusterName", "method", "B+Plasma - method per cluster", save_str="B_Plasma_stats_method")
composition_barplot(Endo_cells, "ClusterName", "Sample_ID", "Endothelial - sample per cluster",save_str="Endo_stats_sample")
composition_barplot(Endo_cells, "ClusterName", "method", "Endothelial - method per cluster", save_str="Endo_stats_method")



#######################################################################
# (4) Creating heatmaps of top genes per sub-clustering (for each lineage) (ED figure 4b)
#######################################################################
# Reading excel file into a dictionary since pg.markers doesn't work with Pegasus1.2 on older AnnData
def excel_markers_to_dict(adata, cluster_label, DE_file_path):
    all_DEGs = {}
    for clust in list(set(adata.obs[cluster_label])):
        print(clust)
        all_DEGs[clust] = {}
        all_DEGs[clust]["up"] = pd.read_excel(DE_file_path, sheet_name="up " + clust, index_col=0)
        all_DEGs[clust]["down"] = pd.read_excel(DE_file_path, sheet_name="down " + clust, index_col=0)
    return all_DEGs

# Merging 3 DEG dictionaries: two are based on AUC DEGs - one defined by AUC threshold, one defined by top X genes fromt he top. the thirsd dictionary is based on pseudo-bulk analysis of OVA - One (cluster) Vs. All (other clusters).
def merging_3_DEG_dict(top_genes_OVA, DEGs_AUC_thresh, DEGs_AUC_top):
    for k, v in list(top_genes_OVA.items()): # Making the keys AUC and OVA dictionaries similar
        clust = k.split(" ")[0] # removing "vs all" string in tab titles
        top_genes_OVA[clust] = top_genes_OVA.pop(k)
    DEGs_AUC_for_merge={}
    for k in DEGs_AUC_thresh.keys():
        DEGs_AUC_for_merge[k] = list(DEGs_AUC_thresh[k])
    full_dict = merge_dicts(merge_dicts(top_genes_OVA, DEGs_AUC_for_merge), DEGs_AUC_top) # Merging dictionaries
    return(full_dict)


# Input: AnnData and a dictionary with up-regulated genes per cluster
# Output: A data frame of the mean gene expression per cluster of each gene. The data frame is sorted by the clusters, so it is ready for ploting using seaborn clustermap
def get_exp_df_for_heatmap(adata, gene_dict):
    #top_genes_final = pd.DataFrame.from_dict(gene_dict).melt()
    # The next line is done with orient=index and then transpose and remove NAs to allow using dictionaries with different number of genes per cluster
    top_genes_final = pd.DataFrame.from_dict(gene_dict, orient='index').transpose().melt().dropna()
    top_genes_final.rename(columns={"variable": "Cluster", "value": "Gene"}, inplace=True)
    mean_ex_mat = get_mean_exp_per_cluster(adata, top_genes_final["Gene"])
    mean_ex_mat = mean_ex_mat.merge(top_genes_final, left_index=True, right_on="Gene").set_index("Gene")
    mean_ex_mat = mean_ex_mat.drop("Cluster", axis=1)
    mean_ex_mat = mean_ex_mat.drop_duplicates()  # genes may be upregulated for more than one cluster
    mean_ex_mat.columns = [x.split(":")[1] for x in mean_ex_mat.columns] # leave only cluster name as column name in the heatmap
    return(mean_ex_mat)

# This function gets an AnnData and a list of genes and returns a dataframe with the mean expression per cluster of each gene. The AnnData is required to have a varm["de_res"] attribute created by Pegasus
def get_mean_exp_per_cluster(data, gene_list):
    gene_names = data.var_names.copy()
    de_df = pd.DataFrame(data.varm["de_res"])
    de_df["Gene"]= gene_names
    de_df = de_df.set_index("Gene")
    cols = list()
    for col in de_df.columns:
        if re.search("mean_logExpr:", col):
            cols.append(col)
    gene_list_rel = [x for x in gene_list if x in gene_names]
    return(de_df.loc[gene_list_rel,cols])


# Plots heatmap to show the genes that characterize each sub-cluster
# - Reads relevant DEG files (top X DEGs + by AUC threshold + pseudo-bulk DEGs), filters them, combines them and plots a heatmap.
# - Writes the data the plot is made of
def heatmap_for_subclustering(adata, cluster_label, DE_file_path, OVA_file_path, logFC_OVA_thresh, title, pdf_out_path,data_out_path, AUC_thresh):
    print("AUC DEGs")
    all_DEGs = excel_markers_to_dict(adata=adata, cluster_label=cluster_label, DE_file_path=DE_file_path)
    DEGs_AUC_thresh = filter_markers(all_DEGs, AUC_thresh)  # DEGs by AUC threshold
    DEGs_AUC_top = top_AUC_DEGs(all_DEGs, 10)  # DEGs by AUC - 10 from the top

    print("AUC pseudo-bulk - OVA")
    DEGs_OVA = pd.read_excel(OVA_file_path, sheet_name=None)  # DEGs from pseudo-bulk OVA
    top_genes_OVA = filter_markers_PB(DEGs_OVA, logFC_thresh=logFC_OVA_thresh, adj_pval_thresh=0.05)

    print("Mergint the 3 tables")
    full_dict = merging_3_DEG_dict(top_genes_OVA, DEGs_AUC_thresh, DEGs_AUC_top)
    mean_ex_mat_all = get_exp_df_for_heatmap(adata, full_dict)

    print("Plotting")
    cells_fig = sb.clustermap(mean_ex_mat_all, col_cluster=False, row_cluster=False, z_score=0, cmap="bwr",center=0, yticklabels=True)
    cells_fig.fig.set_figwidth(6)
    cells_fig.fig.set_figheight(25)
    plt.title(title, fontsize=20, horizontalalignment="left")
    cells_fig.savefig(pdf_out_path, bbox_inches="tight")

    # Lets print the amount of DEGs from each table
    def gene_no(dict):
        genes = [dict[x] for x in dict.keys()]
        genes_list = [item for sublist in genes for item in sublist]
        return (len(set(genes_list)))
    print("top_genes_OVA: " + str(gene_no(top_genes_OVA)))
    print("DEGs_AUC_thresh: " + str(gene_no(DEGs_AUC_thresh)))
    print("DEGs_AUC_top: " + str(gene_no(DEGs_AUC_top)))

    print("Saving data")
    mean_ex_mat_all.to_csv(data_out_path)


# ~~~~~~~~~~~~~~~~~ Endothelial cells  ~~~~~~~~~~~~~~~~~ #
Endo_cells = an.read_h5ad("Endothelial.h5ad")
heatmap_for_subclustering(adata= Endo_cells,
                          cluster_label = 'leiden_res_0.5',
                          DE_file_path="Input/Endothelial_DE_results.xlsx",
                          OVA_file_path = "Input/Endo_final_0.5_de_OVA.xlsx",
                          logFC_OVA_thresh=4.2,
                          title = "Endothelial cells",
                          pdf_out_path = "Heatmap_endothelial2.pdf",
                          data_out_path = "Endothelial_heatmap_data.csv",
                          AUC_thresh=0.8)



# ~~~~~~~~~~~~~~~~~ T/NK cells  ~~~~~~~~~~~~~~~~~ #
T_NK_filt = an.read_h5ad("T_NK.h5ad")
heatmap_for_subclustering(adata= T_NK_filt,
                          cluster_label = 'leiden_res_0.7',
                          DE_file_path="Input/T_NK_DE_results.xlsx",
                          OVA_file_path = "Input/TNK_0.7_de_OVA.xlsx",
                          logFC_OVA_thresh=3,
                          title = "T+NK cells",
                          pdf_out_path = "Heatmap_T_NK.pdf",
                          data_out_path = "Data_heatmap_T+NK.csv",
                          AUC_thresh=0.8)


# ~~~~~~~~~~~~~~~~~ Myeloid cells  ~~~~~~~~~~~~~~~~~ #
Myeloid = an.read_h5ad("Myeloid.h5ad")
heatmap_for_subclustering(adata= Myeloid,
                          cluster_label = 'leiden_final',
                          DE_file_path="Input/Myeloid_DE_results.xlsx",
                          OVA_file_path = "Input/Myeloid_final_de_OVA.xlsx",
                          logFC_OVA_thresh=4,
                          title = "Myeloid cells",
                          pdf_out_path = "Heatmap_Myeloid.pdf",
                          data_out_path = "Data_heatmap_Myeloid.csv",
                          AUC_thresh=0.8)


# ~~~~~~~~~~~~~~~~~ Myeloid cells  ~~~~~~~~~~~~~~~~~ #
B_plasma = an.read_h5ad("B_Plasma.h5ad")
heatmap_for_subclustering(adata= B_plasma,
                          cluster_label = 'leiden_res_0.3',
                          DE_file_path="Input/B_Plasma_DE_results.xlsx",
                          OVA_file_path = "Input/B_Plasma_final_de_OVA.xlsx",
                          logFC_OVA_thresh=2,
                          title = "Myeloid cells",
                          pdf_out_path = "Heatmap_B_Plasma.pdf",
                          data_out_path = "Data_heatmap_B_Plasma.csv",
                          AUC_thresh=0.75)


#######################################################################
# (5) Comparing UMAPs manual and predicted (ED figure 4d,e)
#######################################################################
Myeloid = an.read_h5ad("Myeloid.h5ad")
T_NK = an.read_h5ad("T_NK.h5ad")

def umap_compare_manual_predicted(data, preds):
    for p in preds:
        data.obs[p] = False
        data.obs.loc[data.obs['predictions']==p,p] = True
        data.obs[p] = pd.Categorical(data.obs[p])
        sc.pl.umap(data, color=p, size=30, palette=["blue", "orange"])
        plt.show()

umap_compare_manual_predicted(data=Myeloid, preds=["cDC","macrophage","monocyte","neutrophil"])
umap_compare_manual_predicted(data=T_NK, preds=["CD4+ T cell","CD8+ T cell","NK cell","nkt","Treg"])

