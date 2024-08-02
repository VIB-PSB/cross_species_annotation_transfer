# A framework for cross-species cluster annotation transfer
This framework aims to transfer existing single-cell cluster annotations to a new dataset by detecting cluster similarities.

## Input
This framework takes as input two CSV files, containing the differentially expressed genes (DEGs), or marker genes, of each cluster. Inference of these DEGs should be done prior, and is not part of the framework. We recommend to filter the input data to only retain high-quality DEGs. This is currently also not part of the framework and should be done prior.

To infer DEGs, we suggest using [Seurat](https://satijalab.org/seurat/)'s [`FindAllMarkers()`](https://www.rdocumentation.org/packages/Seurat/versions/5.0.3/topics/FindAllMarkers) function, although other methods might be employed. To filter the input DEGs, we recommend to either set a cutoff (on log fold change and/or adjusted p-value), or by taking the top N best DEGs (sorting on q-value, and using fold change to break ties).

Each input CSV file should have the following column names:
- "gene ID": the gene ID of the inferred marker gene
- "cluster": the cluster name of which the gene is a marker
- "avg_log2FC": the average log fold change in expression of the cells in the given cluster, compared to all the other clusters 
- "p_val_adj": the adjusted p-value calculated when computationally inferring the marker gene
- "Orth_group" (only needed when making a cross-species comparison): the ID of the orthologous group the marker gene is part of

NOTE: The orthologous group IDs were calculated using [OrthoFinder](https://github.com/davidemms/OrthoFinder). Precomputed orthology files, including several commonly used plant species, are available through [PLAZA](https://bioinformatics.psb.ugent.be/plaza/).

## Parameters
- *separator* (default = ","): the column separator of the input CSV files
- *nb_of_background_sets* (default = 1000): the number of background sets (determines the minimal p-value that can be calculated)
  <br/>NOTE: this parameter strongly affects the runtime
- *minimal_real_matches* (default = 2): the minimal number of background sets that has to show an overlap with the query cluster DEGs, larger than or equal to the real overlap between the query and reference cluster. E.g., if this parameter is set to 2, and there is only 1 background set of which the overlap with the query cluster is larger than or equal to the real overlap, the p-value will automatically be set to 1.
- *compare_orthogroups* (default = `False`): use the orthogroups to compare the clusters of the different dataset. This should be set to `True` if a cross-species is being performed, it should be set to `False` if a comparison within the same species is being performed.
