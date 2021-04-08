# A Spatial Atlas of COVID-19 Lung

## Overview

1. [Packages in Use](#1-packages-in-use)

    1.1. [R packages](#11-r-packages)

    1.2. [Python packages](#12-python-packages)
2. [Data](#2-data)
3. [WTA and CTA Expression](#3-wta-and-cta-expression)
4. [DE and GSEA Analysis](#4-de-and-gsea-analysis)

    4.1. [Between subjects and controls for WTA and CTA data](#41-between-subjects-and-controls-for-wta-and-cta-data)

    4.2. [Between SARS-CoV-2 high and low AOIs](#42-between-sars-cov-2-high-and-low-aois)

## 1. Packages in Use

### 1.1. R packages

**R 4.0.3** and **Bioconductor 3.11** are used in analysis.

|Package|Version|Release Date|
|---|---|---|
|[limma](https://bioconductor.org/packages/3.11/bioc/html/limma.html)|3.44.3|2020/04/28|
|[edgeR](https://bioconductor.org/packages/3.11/bioc/html/edgeR.html)|3.30.3|2020/04/28|
|[EnhancedVolcano](https://bioconductor.org/packages/3.11/bioc/html/EnhancedVolcano.html)|1.6.0|2020/04/28|
|[fgsea](https://bioconductor.org/packages/3.11/bioc/html/fgsea.html)|1.14.0|2020/04/28|
|[stringr](https://cran.r-project.org/web/packages/stringr/index.html)|1.4.0|2019/02/09|
|[rlist](https://cran.r-project.org/web/packages/rlist/index.html)|0.4.6.1|2016/04/04|
|[plyr](https://cran.r-project.org/web/packages/plyr/index.html)|1.8.6|2020/03/03|
|[qusage](https://bioconductor.org/packages/3.11/bioc/html/qusage.html)|2.22.0|2020/04/28|


### 1.2. Python packages

**Python 3.7** is used in analysis.

|Package|Version|Release Date|
|---|---|---|
|[Pegasus](https://pypi.org/project/pegasuspy/1.3.0/)|1.3.0|2021/02/02|
|[PegasusIO](https://pypi.org/project/pegasusio/0.2.10/)|0.2.10|2021/02/02|
|[pandas](https://pypi.org/project/pandas/1.2.2/)|1.2.2|2021/02/09|
|[numpy](https://pypi.org/project/numpy/1.20.1/)|1.20.1|2021/02/07|
|[scipy](https://pypi.org/project/scipy/1.6.0/)|1.6.0|2020/12/30|
|[seaborn](https://pypi.org/project/seaborn/0.11.1/)|0.11.1|2020/12/20|
|[matplotlib](https://pypi.org/project/matplotlib/3.3.3/)|3.3.3|2020/11/11|

## 2. Data

Nanostring GeMx Q3-normalized WTA and CTA count matrices are used in analysis, which are available on the [GEO](https://www.ncbi.nlm.nih.gov/geo/) database under accession no. GSE162911. After downloading the data, they are put as the following:

|File|Description|
|---|---|
|`data/Broad-COVID_WTA_Q3Norm_TargetCountMatrix.txt`|Q3-normalized WTA count matrix in tab-separated txt format, with gene names as rows and AOIs as columns.|
|`data/Broad-COVID_WTA_SegmentProperties.txt`|Segment attributes of WTA AOIs in tab-separated txt format.|
|`data/Broad-COVID_CTA_Q3Norm_TargetCountMatrix.txt`|Q3-normalized CTA count matrix in tab-separated txt format, with gene names as rows and AOIs as columns.|
|`data/Broad-COVID_CTA_SegmentProperties.txt`|Segment attributes of CTA AOIs in tab-separated txt format.|
|`data/Broad-COVID_WTA_BioProbeCountMatrix.txt`|Raw WTA count matrix in tab-separated txt format, only used for extracting WTA targets.|
|`data/Broad-COVID_CTA_BioProbeCountMatrix.txt`|Raw CTA count matrix in tab-separated txt format, only used for extracting CTA targets.|

Besides, there are other files put into repo in advance for user's convenience:

|File|Description|
|---|---|
|`data/annotation_file_wta.txt`|Primary Morph type of each WTA AOI, which is available on GEO.|
|`data/annotation_file_cta.txt`|Primary Morph type of each CTA AOI, which is available on GEO.|
|`data/wta_cta.csv`|Correspondence between WTA and CTA AOIs, which is available on GEO.|
|`utils/lung_spatial_markers.gmt`|Gene markers for cell type deconvolution of WTA data in gmt format, which is in Extended Data Table 8 of the paper.|
|`utils/lung_spatial_markers_CTA.gmt`|Gene markers for cell type deconvolution of CTA data in gmt format, which is in Extended Data Table 8 of the paper.|
|`utils/h.all.v7.1.symbols.gmt`|MSigDB Hallmark v7.1 gene set in gmt format ([download](https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.1/h.all.v7.1.symbols.gmt))|
|`utils/c2.cp.v7.1.symbols.gmt`|MSigDB C2-CP v7.1 gene set in gmt format ([download](https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.1/c2.cp.v7.1.symbols.gmt))|
|`utils/c5.bp.v7.1.symbols.gmt`|MSigDB C5-BP v7.1 gene set in gmt format ([download](https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.1/c5.bp.v7.1.symbols.gmt))|

Among the three MSigDB gene sets, only Hallmark is considered in the analysis.

## 3. WTA and CTA Expression

Run
```
python classify_targets.py wta
```
to generate lists of WTA targets (`wta_target.csv`) and SARS-CoV-2 related targets (`wta_sars_cov_2.csv`). Run
```
python classify_targets.py cta
```
to generate lists of CTA targets (`cta_targets.csv`) and SARS-CoV-2 related targets (`cta_sars_cov_2.csv`). They are used for generating **Figure 8D** and **Extended Data Table 6(a)**.

Run
```
python plot_figures.py
```
to generate **Extended Data Figure 8E** (`ED_Figure_8/ED_Figure_8e.pdf`) and **Extended Data Figure 8F** (`ED_Figure_8/ED_Figure_8f.pdf`).

Run `heatmap.ipynb` Jupyter notebook to generate **Figure 4B** and **Extended Data Figure 8B, 8C, 8G, 8K** in folders `Figure_4` and `ED_Figure_8`, respectively.

## 4. DE and GSEA Analysis

### 4.1. Between subjects and controls for WTA and CTA data

The following command executes the DE and GSEA analysis between subjects and controls for WTA data:

```
Rscript WTA_Patient_vs_Control.R
```
When finished, it generates the following output:

|File|Description|
|---|---|
|`Figure_4/Figure_4c_left.pdf`|Left part of **Figure 4C** in paper. DE volcano plot on WTA PanCK+ AOIs.|
|`Figure_4/Figure_4c_left_DE.up.csv`|Up-regulated genes for subjects on WTA PanCK+ AOIs. Part of **Extended Data Table 6(b)** in paper.|
|`Figure_4/Figure_4c_left_DE.down.csv`|Up-regulated genes for controls on WTA PanCK+ AOIs. Part of **Extended Data Table 6(b)** in paper.|
|`ED_Figure_8/ED_Figure_8h_left.pdf`|Left part of **Extended Data Figure 8H** in paper. DE volcano plot on WTA PanCK- AOIs.|
|`ED_Figure_8/ED_Figure_8h_left_DE.up.csv`|Up-regulated genes for subjects on WTA PanCK- AOIs. Part of **Extended Data Table 6(b)** in paper.|
|`ED_Figure_8/ED_Figure_8h_left_DE.down.csv`|Up-regulated genes for controls on WTA PanCK+ AOIs. Part of **Extended Data Table 6(b)** in paper.|

The following command executes for CTA data:

```
Rscript CTA_Patient_vs_Control.R
```
When finished, it generates the following output:

|File|Description|
|---|---|
|`ED_Figure_8/ED_Figure_8i_left.pdf`|Left part of **Extended Data Figure 8I** in paper. DE volcano plot on CTA PanCK+ AOIs.|
|`ED_Figure_8/ED_Figure_8i_left_DE.up.csv`|Up-regulated genes for subjects on CTA PanCK+ AOIs. Part of **Extended Data Table 6(c)** in paper.|
|`ED_Figure_8/ED_Figure_8i_left_DE.down.csv`|Up-regulated genes for controls on CTA PanCK+ AOIs. Part of **Extended Data Table 6(c)** in paper.|
|`ED_Figure_8/ED_Figure_8j_left.pdf`|Left part of **Extended Data Figure 8J** in paper. DE volcano plot on CTA PanCK- AOIs.|
|`ED_Figure_8/ED_Figure_8j_left_DE.up.csv`|Up-regulated genes for subjects on CTA PanCK- AOIs. Part of **Extended Data Table 6(c)** in paper.|
|`ED_Figure_8/ED_Figure_8j_left_DE.down.csv`|Up-regulated genes for controls on CTA PanCK- AOIs. Part of **Extended Data Table 6(c)** in paper.|

When both analyses are finished, run the following command to generate GSEA bar plots:

```
python gen_gsea_plots.py
```
which gives the following output:

|File|Description|
|---|---|
|`Figure_4/Figure_4c_right.H.up.png`|Significant pathways for subjects on WTA PanCK+ AOIs. Right part of **Figure 4C** in paper, but in bar plot.|
|`Figure_4/Figure_4c_right.H.down.png`|Significant pathways for controls on WTA PanCK+ AOIs. Right part of **Figure 4C** in paper, but in bar plot.|
|`ED_Figure_8/ED_Figure_8h_right.H.up.png`|Significant pathways for subjects on WTA PanCK- AOIs. Right part of **Extended Data Figure 8H** in paper, but in bar plot.|
|`ED_Figure_8/ED_Figure_8h_right.H.down.png`|Significant pathways for controls on WTA PanCK- AOIs. Right part of **Extended Data Figure 8H** in paper, but in bar plot.|
|`ED_Figure_8/ED_Figure_8i_right.H.up.png`|Significant pathways for subjects on CTA PanCK+ AOIs. Right part of **Extended Data Figure 8I** in paper, but in bar plot.|
|`ED_Figure_8/ED_Figure_8i_right.H.down.png`|Significant pathways for controls on CTA PanCK+ AOIs. Right part of **Extended Data Figure 8I** in paper, but in bar plot.|
|`ED_Figure_8/ED_Figure_8j_right.H.up.png`|Significant pathways for subjects on CTA PanCK- AOIs. Right part of **Extended Data Figure 8J** in paper, but in bar plot.|
|`ED_Figure_8/ED_Figure_8j_right.H.down.png`|Significant pathways for controls on CTA PanCK- AOIs. Right part of **Extended Data Figure 8J** in paper, but in bar plot.|

### 4.2. Between SARS-CoV-2 high and low AOIs

First, run `Signature_Score_SARS-CoV-2_High_vs_Low.ipynb` Jupyter notebook to generate **Figure 4D** and **Extended Data Figure 8L**, which are scatter plots on COVID-19 signature scores of WTA AOIs. Then they are categorized into SARS-CoV-2 *High*, *Medium*, and *Low* groups based on these scores. This categorization information is generated as `WTA_SARS-CoV-2_PanCK.sig.csv` (for WTA PanCK+ AOIs) and `WTA_SARS-CoV-2_Syto13.sig.csv` (for WTA PanCK- AOIs), which will be used in the next step.

Next, run the following command to perform DE analysis:

```
Rscript SARS-CoV-2_High_vs_Low.R
```
which gives the following output:

|File|Description|
|---|---|
|`Figure_4/Figure_4e.pdf`|DE volcano plot on WTA PanCK+ AOIs. **Figure 4E** in paper.|
|`Figure_4/Figure_4e_DE.up.csv`|Up-regulated genes for SARS-CoV-2 High group on WTA PanCK+ AOIs. Part of **Extended Data Table 9(b)** in paper.|
|`Figure_4/Figure_4e_DE.down.csv`|Up-regulated genes for SARS-CoV-2 Low group on WTA PanCK+ AOIs. Part of **Extended Data Table 9(b)** in paper.|
|`ED_Figure_8/ED_Figure_8m.pdf`|DE volcano plot on WTA PanCK- AOIs. **Extended Data Figure 8M** in paper.|
|`ED_Figure_8/ED_Figure_8m_DE.up.csv`|Up-regulated genes for SARS-CoV-2 High group on WTA PanCK- AOIs. Part of **Extended Data Table 9(b)** in paper.|
|`ED_Figure_8/ED_Figure_8m_DE.down.csv`|Up-regulated genes for SARS-CoV-2 Low group on WTA PanCK- AOIs. Part of **Extended Data Table 9(b)** in paper.|
