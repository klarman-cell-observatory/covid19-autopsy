# Heart data analyses

These python notebooks comprise the portion of the work on the heart data that 
made it directly into the paper.

## `08_heart_annotation`
Manual annotation of clusters in the heart data, as well as how these manual 
annotations align with auto-generated cell type labels.  Includes the breakdown 
of the heart data by cell type per donor.

- Extended Data Figure 10a
- Extended Data Figure 10b
- Extended Data Figure 10c
- Extended Data Figure 10j
- Extended Data Figure 11c

## `17_heart_DE_metastudy`
Limma-voom differential expression analysis of COVID-19 versus healthy heart 
data from other studies.

- Extended Data Figure 11d
- Extended Data Figure 11g
- Extended Data Figure 11h
- Extended Data Figure 11i

## `18_heart_DE_metastudy_umaps_and_extraQC`
Extra cell QC on the heart data (COVID-19 plus healthy control data), which was 
carried out in order to create a UMAP that is easier to interpret visually.  
Without the extra QC, the UMAP is a bit chaotic, due to the inclusion of 
relatively low-quality cells from the various studies.

- Extended Data Figure 11f
- Extended Data Figure 11j
- Extended Data Figure 11k

## `20_heart_DE_metastudy_extraQC`
A companion notebook to `17_heart_DE_metastudy`: this is the same differential 
expression analysis, but run only on the aggresively-QCed cells left after the 
cell QC procedure in `18_heart_DE_metastudy_umaps_and_extraQC`.  The results 
are very much in agreement with the results in `17_heart_DE_metastudy`, and so 
this information was not included separately in the paper.

## `22_heart_DE_meta_cell_proportions`
An examination of the proportions of each cell type present in samples from the 
metastudy of COVID-19 and healthy heart.  Contains the scCODA analysis.

- Extended Data Figure 11e

## `23_heart_single_cell_portal_files`
Creation of the heart `.h5ad` file (and other associated files) uploaded to the 
Single Cell Portal.

## `sc_utils.py`
A collection of useful functions written and maintained by Stephen Fleming and 
the Precision Cardiology Lab at the Broad Institute.  This may be released as 
its own package in the future if there is interest from the community.  One of 
its main features is the simplification of using python to run R pipelines for 
summation + voom + limma differential expression testing in a uniform way using 
an AnnData object.
