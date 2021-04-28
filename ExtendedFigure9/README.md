# Extended Figure 9

Scripts to generate Extended Figure 9: b,c,d

Extended Data Figure 9. GeoMx WTA DSP analysis of lung biopsies reveals region- and inflammation-specific expression programs
a. Region selection. Serial sections of lung biopsies (five donors, D13-17; image depicts serial sections of D14) processed with GeoMx WTA-DSP with 4-color staining (DNA, CD45, CD68, PanCK), RNAscope with probes against (SARS-CoV-2 S-gene (utilized to derive semi-quantitative viral load scores), ACE2, TMPRSS2), H&E staining, and immunohistochemistry with anti-SARS-CoV-2 S-protein. Scale bar: 100 µm. b-d. Regions and inflammation specific expression programs. b. The first two principal components (PCs, x and y axes) from lung ROI gene expression profiles from donors D13-17, spanning normal-appearing alveoli (green; D14=6 AOIs, D15=2 AOIs, D16=5 AOIs, D17=4 AOIs); inflamed alveoli (magenta; D13=14 AOIs, D14=18 AOIs, D15=7 AOIs, D16=3 AOIs ,D17=8 AOIs); bronchial epithelium (blue;  D14 =2 AOIs, D15 =1 AOI, D16 =2 AOIs, D17 =3 AOIs), and arterial blood vessels (black; D13=2 AOIs, D15=3 AOIs). c. GSEA score (circle size, legend) of the enrichment of the interferon-γ pathway in each normal-appearing (green; 6 AOIs) and inflamed (magenta; 18 AOIs) alveolar AOIs (dot) from the section of donor D14 (in a), placed in their respective physical coordinates on the tissue section (as in a). d. Expression (color bar, log2(counts per million)) of IFNγ pathway genes (rows) from normal-appearing (green, n=6) and inflamed alveoli (magenta, n=18) AOIs (columns) from D14 lung biopsy.

## Packages in Use


**R version 3.6.0** was used to generate the plots.

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Package </th>
   <th style="text-align:left;"> Loaded version </th>
   <th style="text-align:left;"> Date </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> annotate </td>
   <td style="text-align:left;"> 1.64.0 </td>
   <td style="text-align:left;"> 2019-10-29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AnnotationDbi </td>
   <td style="text-align:left;"> 1.48.0 </td>
   <td style="text-align:left;"> 2019-10-29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Biobase </td>
   <td style="text-align:left;"> 2.46.0 </td>
   <td style="text-align:left;"> 2019-10-29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BiocGenerics </td>
   <td style="text-align:left;"> 0.32.0 </td>
   <td style="text-align:left;"> 2019-10-29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cowplot </td>
   <td style="text-align:left;"> 1.1.1 </td>
   <td style="text-align:left;"> 2020-12-30 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> devtools </td>
   <td style="text-align:left;"> 2.3.2 </td>
   <td style="text-align:left;"> 2020-09-18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dplyr </td>
   <td style="text-align:left;"> 1.0.3 </td>
   <td style="text-align:left;"> 2021-01-15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> edgeR </td>
   <td style="text-align:left;"> 3.28.1 </td>
   <td style="text-align:left;"> 2020-02-26 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> forcats </td>
   <td style="text-align:left;"> 0.5.1 </td>
   <td style="text-align:left;"> 2021-01-27 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> genefilter </td>
   <td style="text-align:left;"> 1.68.0 </td>
   <td style="text-align:left;"> 2019-10-29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggplot2 </td>
   <td style="text-align:left;"> 3.3.3 </td>
   <td style="text-align:left;"> 2020-12-30 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> graph </td>
   <td style="text-align:left;"> 1.64.0 </td>
   <td style="text-align:left;"> 2019-10-29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GSEABase </td>
   <td style="text-align:left;"> 1.48.0 </td>
   <td style="text-align:left;"> 2019-10-29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GSVA </td>
   <td style="text-align:left;"> 1.34.0 </td>
   <td style="text-align:left;"> 2019-10-29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IRanges </td>
   <td style="text-align:left;"> 2.20.2 </td>
   <td style="text-align:left;"> 2020-01-13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> kableExtra </td>
   <td style="text-align:left;"> 1.3.1 </td>
   <td style="text-align:left;"> 2020-10-22 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> knitr </td>
   <td style="text-align:left;"> 1.31 </td>
   <td style="text-align:left;"> 2021-01-27 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> limma </td>
   <td style="text-align:left;"> 3.42.2 </td>
   <td style="text-align:left;"> 2020-02-03 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pheatmap </td>
   <td style="text-align:left;"> 1.0.12 </td>
   <td style="text-align:left;"> 2019-01-04 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> purrr </td>
   <td style="text-align:left;"> 0.3.4 </td>
   <td style="text-align:left;"> 2020-04-17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> readr </td>
   <td style="text-align:left;"> 1.4.0 </td>
   <td style="text-align:left;"> 2020-10-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S4Vectors </td>
   <td style="text-align:left;"> 0.24.4 </td>
   <td style="text-align:left;"> 2020-04-09 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stringr </td>
   <td style="text-align:left;"> 1.4.0 </td>
   <td style="text-align:left;"> 2019-02-10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tibble </td>
   <td style="text-align:left;"> 3.0.5 </td>
   <td style="text-align:left;"> 2021-01-15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyr </td>
   <td style="text-align:left;"> 1.1.2 </td>
   <td style="text-align:left;"> 2020-08-27 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyverse </td>
   <td style="text-align:left;"> 1.3.0 </td>
   <td style="text-align:left;"> 2019-11-21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> usethis </td>
   <td style="text-align:left;"> 2.0.0 </td>
   <td style="text-align:left;"> 2020-12-10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> XML </td>
   <td style="text-align:left;"> 3.99-0.3 </td>
   <td style="text-align:left;"> 2020-01-20 </td>
  </tr>
</tbody>
</table>

## Data

- `LungDSP_Q3Norm_counts.tsv`: Q3-normalized WTA count matrix in tab-separated txt format, with gene names as rows and AOIs as columns.
- `LungDSP_sample_annotation.tsv`: Segment attributes of WTA AOIs in tab-separated txt format.
- `c2.cp.v7.1.symbols.gmt`: Canonical Pathways Annotation from MSigDB



