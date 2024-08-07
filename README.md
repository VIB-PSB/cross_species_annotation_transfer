# A framework for cross-species cluster annotation transfer
This framework aims to transfer existing single-cell cluster annotations to a new dataset by detecting cluster similarities.

## Input
This framework takes as input minimally two CSV files, containing the differentially expressed genes (DEGs), or marker genes, of each cluster. Inference of these DEGs should be done prior, and is not part of the framework. We recommend filtering the input data to only retain high-quality DEGs. This is currently also not part of the framework and should be done prior.

To infer DEGs, we suggest using [Seurat](https://satijalab.org/seurat/)'s [`FindAllMarkers()`](https://www.rdocumentation.org/packages/Seurat/versions/5.0.3/topics/FindAllMarkers) function, although other methods might be employed. To filter the input DEGs, we recommend to either set a cutoff (on log fold change and/or adjusted p-value), or by taking the top N best DEGs (sorting on q-value, and using fold change to break ties).

Each input CSV file should have the following column names:
- "gene ID": the gene ID of the inferred marker gene
- "cluster": the cluster name of which the gene is a marker
- "avg_log2FC": the average log fold change in expression of the cells in the given cluster, compared to all the other clusters 
- "p_val_adj": the adjusted p-value calculated when computationally inferring the marker gene
- "Orth_group" (only needed when making a cross-species comparison): the ID of the orthologous group the marker gene is part of

NOTE: The orthologous group IDs were calculated using [OrthoFinder](https://github.com/davidemms/OrthoFinder). Precomputed orthology files, including several commonly used plant species, are available through [PLAZA](https://bioinformatics.psb.ugent.be/plaza/).

When more than two input files are given, all pairwise combinations of the given datasets will be evaluated.

## Parameters
- *separator* (default = ","): the column separator of the input CSV files
- *nb_of_background_sets* (default = 1000): the number of background sets (determines the minimal p-value that can be calculated)
- *minimal_real_matches* (default = 2): the minimal number of background sets that has to show an overlap with the query cluster DEGs, larger than or equal to the real overlap between the query and reference cluster. E.g., if this parameter is set to 2, and there is only 1 background set of which the overlap with the query cluster is larger than or equal to the real overlap, the p-value will automatically be set to 1.
- *compare_orthogroups* (default = `False`): use the orthogroups (orthologous groups) to compare the clusters of the different dataset. This should be set to `True` if a cross-species is being performed, it should be set to `False` if a comparison within the same species is being performed.
- *significance_threshold* (default = 0.05): the adjusted p-value threshold at which the DEG overlap between two clusters is considered significant

NOTE: the *nb_of_background_sets* parameter strongly affects the runtime

## Output
The framework produces two output files, and a log file in which the progress can be followed during execution.

- `cluster_deg_overlap_statistics.xlsx`: this output file is the most verbose and contains information on the all-versus-all cluster comparisons. For each comparison, two tabs are present where the query and reference dataset is reversed. I.e., in the first comparisons, backgrounds sets are being generated of the first dataset, while in the second, background sets are generated for the other dataset. This output file contains the following fields:
  - \<name query dataset\>: the cluster in the query dataset
  - \<name reference dataset\>: the cluster in the reference dataset
  - Number of matches: the number of (real) DEG matches between the query and reference cluster
  - Enrichment fold: the enrichment fold between the number of real matches and the average number of matches in the background
  - P-value: the p-value of the enrichment
  - Q-value: the adjusted p-value (i.e. q-value), calculated by adjusting the p-value for multiple testing using the [Benjamini-Hochberg procedure](https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html)
  - Significance: does the q-value pass the significance threshold?
  - Match IDs: the IDs of the matching DEGs (either the orthologous group IDs in case of a cross-species comparison, or the gene IDs in case of a comparison within-species)
  - Matching query genes: the IDs of the matching query DEGs (always gene IDs, thus identical to "Match IDs" in case of a comparison within-species)
  - Matching reference genes: the IDs of the matching reference DEGs (always gene IDs, thus identical to "Match IDs" in case of a comparison within-species)
- `cluster_deg_overlap_statistics_combined_best_hits.xlsx`: this output file is a subset of the `cluster_deg_overlap_statistics.xlsx` file, but has merged the two-way comparison together by selecting the one with the lowest q-value. Then, for each cluster, only the best (significant) hit in the other dataset is retained. This is done for both datasets.
  - This output file contains the same output fields as `cluster_deg_overlap_statistics.xlsx`.
- `<dataset1>-<dataset2>_overlap_significance.png`: this is an additional output file, visualizing the results of  `cluster_deg_overlap_statistics.xlsx`. The figure contains 3 x 2 figures. The two columns left and right correspond to the two-way comparison between the two datasets. The three rows show for every cluster-cluster comparison the log<sub>10</sub>(q-value), the significance, and the enrichment fold respectively.
- `log.txt`: the log file, tracking the progress of the run during execution

## Running
Execute the Python script in this folder.
```
python annotation_transfer_analysis.py
```

## Dependencies
This framework was developed using Python 3.8.0, with dependencies listed in [`requirements.txt`](https://github.com/VIB-PSB/cross_species_annotation_transfer/blob/main/requirements.txt).
