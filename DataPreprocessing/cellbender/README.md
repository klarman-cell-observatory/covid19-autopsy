# Analyses related to CellBender QC

These python notebooks comprise some of the work done to QC and perform checks 
on the outputs of the 
[`cellbender remove-background`](https://cellbender.readthedocs.io/en/latest/) 
pipeline.

## `10_cellbender_before_and_after`
Exploration of a dataset in terms of raw data from `CellRanger` and the output of 
`CellBender`.  Examines changes in gene expression, changes in cells kept, and 
resulting changes in the look of the final UMAP, among other things.

- Extended Data Figure 1e
- Extended Data Figure 1f
- Extended Data Figure 1g

## `15_decontx_script`
Details the way we ran `DecontX`, as a baseline to compare against `CellBender 
remove-background`, at the request of a reviewer.

## `16_decontx_cellbender_before_and_after`
A parallel notebook to `10_cellbender_before_and_after`, but this time also 
including the `DecontX` output data.

- Extended Data Figure 1h
