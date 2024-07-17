#!/usr/bin/env python3

# Importing modules
import sys
sys.path.append('/group/transreg/jasst/01_utilities/bin/python3')

from itertools import combinations
from significance_deg_overlap import * # custom module

# Define input file locations
file_locs = {"tae_sc":"input/tae_log2FC0.5_FDR0.05.csv", "tae_sn":"input/tae_sn_top_100_markers_original.csv"}

# Collect all results
overlap_statistics_all = dict()

# Run for all species combinations
for species1, species2 in combinations(file_locs.keys(), 2):
    
    # Read in the input
    species2df_dict = read_files(file_locs, sep=";", include_orthogroups=False)
    
    # Arrange the input in a good format
    species2cluster2genes_and_groups = convert_df_to_dict(species2df_dict, include_orthogroups=False)
    
    # Calculate all statistics for a species-species comparison
    print("[STARTING] {} - {} comparison".format(species1, species2))
    overlap_statistics = calculate_two_way_statistics(species1, species2, \
                                                      species2cluster2genes_and_groups, \
                                                      gene2group_dict=None, \
                                                      nb_of_background_sets=1000, \
                                                      minimal_real_matches=2, \
                                                      compare_orthogroups=False, \
                                                      all_clusters_as_background=True)
    
    # Combine all results
    overlap_statistics_all.update(overlap_statistics)
    
    # Generate plot
    fig = generate_visualization(overlap_statistics, species1, species2)

# Write everything to Excel
write_statistics_to_excel(overlap_statistics_all)
