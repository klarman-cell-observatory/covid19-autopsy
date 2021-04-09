#!/usr/bin/env Rscript
# To save a sparse *count matrix* from Python in a Markert Matrix sparse format run:
#    import scipy.io
#    scipy.io.mmwrite("P3_data.mtx", P3_data.raw.X)
# For meta data, in python:
#    meta = P3_data.obs[["leiden_labels","Patient"]]
#    meta = meta.rename(columns={"leiden_labels":"leiden","Patient":"donor"})
#    meta.to_csv("P3_meta.csv")
# For gene list:
#    pd.DataFrame(P3_data.var_names).to_csv("genes.csv", index=False)

# The arguments given to this script should be are:
# (1) data file name, (2) gene file name, (3) meta-data file name, (4) prefix for output files

source("Summary_pseudobulk_functions.R")
args = commandArgs(trailingOnly=TRUE)
data_file = args[1]
genes_file = args[2]
meta_file = args[3]
prefix = args[4]
with_patient_cor = args[5]

print("** Args: **")
print(paste0("data_file=",data_file))
print(paste0("genes_file=",genes_file))
print(paste0("meta_file=",meta_file))
print(paste0("prefix=",prefix))

data = t(readMM(data_file))
genes = read.table(genes_file, sep=",", header = T)
meta = read.table(meta_file, sep=",", header = T, row.names = 1)
CalcPseudoBulk(data, genes, meta, prefix, with_patient_cor)





