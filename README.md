# A framework for cross-species cluster annotation transfer
This framework aims to transfer existing single-cell cluster annotations to a new dataset by detecting cluster similarities.

## Input
This framework takes as input two CSV files, containing the differentially expressed genes (DEGs), or marker genes, of each cluster. Inference of these DEGs should be done prior, and is not part of the framework. We recommend to filter the input data to only retain high-quality DEGs. This is currently also not part of the framework and should be done prior.

To infer DEGs, we suggest using [Seurat](https://satijalab.org/seurat/)'s [`FindAllMarkers()`](https://www.rdocumentation.org/packages/Seurat/versions/5.0.3/topics/FindAllMarkers) function, although other methods might be employed. To filter the input DEGs, we recommend to either set a cutoff (on log fold change and/or adjusted p-value), or by taking the top N best DEGs (sorting on q-value, and using fold change to break ties).

Each input CSV file should have the following column names:
- "gene ID": The gene ID of the inferred marker gene
- "cluster": The cluster name of which the gene is a marker
- "avg_log2FC": The average log fold change in expression of the cells in the given cluster, compared to all the other clusters 
- "p_val_adj": The adjusted p-value calculated when computationally inferring the marker gene
- "Orth_group" (only needed when making a cross-species comparison): the ID of the orthologous group the marker gene is part of

NOTE: The orthologous group IDs were calculated using [OrthoFinder](https://github.com/davidemms/OrthoFinder). Precomputed orthology files, including several commonly used plant species, are available through [PLAZA](https://bioinformatics.psb.ugent.be/plaza/).

## Parameters

