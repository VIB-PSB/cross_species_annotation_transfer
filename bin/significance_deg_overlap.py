#!/usr/bin/env python3

##########################################
### SIGNIFICANCE OF DEG OVERLAP MODULE ###
##########################################

import math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from natsort import natsorted
from statistics import median, mean
from collections import Counter
from statsmodels.stats.multitest import multipletests

def read_files(file_locs, sep=",", include_orthogroups=True):
    
    """
    Read DEG (Differentially Expressed Genes) files for multiple species and create a dictionary.

    Parameters:
        file_locs (dict): A dictionary with keys as species and values as the location of DEG files for each species.
        sep (str): The delimiter used in the DEG files. Default is ",".

    Returns:
        dict: A dictionary with keys as species and values as dataframes containing selected columns
              (cluster, gene ID, Orthogroup, logFC, q-value) from the corresponding DEG files.
    """
    
    species2df_dict = dict()
    
    # For each given species, read in the DEG file
    for species, file_loc in file_locs.items():
        deg_df = pd.read_csv(file_loc, sep=sep, dtype=str)
        
        # Remove NAs (due to no orthogroup, or no gene conversion)
        if include_orthogroups:
            deg_df = deg_df.dropna(subset=["Orth_group", "gene ID"])
        else:
            deg_df = deg_df.dropna(subset=["gene ID"])
        
        # Select the cluster, gene and orthogroup columns and add to dict
        if include_orthogroups:
            species2df_dict[species] = deg_df[["cluster", "gene ID", "Orth_group", "avg_log2FC", "p_val_adj"]]
        else:
            species2df_dict[species] = deg_df[["cluster", "gene ID", "avg_log2FC", "p_val_adj"]]
    
    return species2df_dict


def convert_df_to_dict(species2df_dict, include_orthogroups=True, remove_duplicate_genes=False):
    
    """
    Convert DEG (Differentially Expressed Genes) dataframes to a structured dictionary format.

    Parameters:
        species2df_dict (dict): A dictionary with keys as species and values as dataframes containing
                                DEG information with columns "cluster", "gene ID", and "Orth_group".
        remove_duplicate_genes (bool): A flag indicating whether to remove duplicate gene IDs within the same cluster.
                                       Default is False.

    Returns:
        tuple: A tuple containing two dictionaries:
            - species2gene2group: A dictionary with keys as species and values as dictionaries mapping gene IDs to Orthogroups.
            - species2cluster2genes_and_groups: A dictionary with keys as species and values as dictionaries
                                               mapping cluster names to sets of genes and Orthogroups.
    """
    
    species2cluster2genes_and_groups = dict()
    if include_orthogroups:
        species2gene2group = dict()
    
    # For each given species, convert the DEG dataframe
    for species, df in species2df_dict.items():
        
        # Put all the gene-to-group relations in a dict
        if include_orthogroups:
            species2gene2group[species] = dict(zip(df["gene ID"], df["Orth_group"]))
        
        # Structure the data per cluster in a dict
        cluster2genes_and_groups = dict()

        # Add the genes and groups of the individual clusters to the dict
        for cluster in set(df.cluster):

            # Select only the rows of a given cluster
            if include_orthogroups:
                cluster_df = df.query("cluster == @cluster")[["cluster", "gene ID", "Orth_group"]]
            else:
                cluster_df = df.query("cluster == @cluster")[["cluster", "gene ID"]]
            
            # Remove duplicate genes in the same cluster (might be introduced when changing gene IDs to another version)
            cluster_df = cluster_df.drop_duplicates()
            
            # If explicitly stated, allow removal of duplicates
            if remove_duplicate_genes:
                cluster_df = cluster_df.drop_duplicates("gene ID")

            # Extract the genes and groups as sets
            genes = set(cluster_df["gene ID"])
            assert len(genes) == len(cluster_df["gene ID"]) # check for duplicated DEGs in one cluster
            if include_orthogroups:
                groups = set(cluster_df["Orth_group"])

            # Add them to the dict
            if include_orthogroups:
                cluster2genes_and_groups[cluster] = {"genes": genes, "groups": groups}
            else:
                cluster2genes_and_groups[cluster] = {"genes": genes}

        # Add the genes and groups of all clusters together to the dict
        if include_orthogroups:
            cluster2genes_and_groups["all_clusters"] = {"genes": set(df["gene ID"]), "groups": set(df["Orth_group"])}
        else:
            cluster2genes_and_groups["all_clusters"] = {"genes": set(df["gene ID"])}
        
        # Add to species dict
        species2cluster2genes_and_groups[species] = cluster2genes_and_groups
    
    if include_orthogroups:
        return species2gene2group, species2cluster2genes_and_groups
    else:
        return species2cluster2genes_and_groups


def get_background(background_degs, cluster_size, nb_of_background_sets=10000, seed=0):
    
    """
    Generate a background model by randomly sampling genes from the background DEGs.

    Parameters:
        background_degs (list): A list of genes from which random samples will be drawn to create background sets.
        cluster_size (int): The desired size of each background set (number of genes in each set).
        nb_of_background_sets (int): The number of background sets to generate. Default is 10000.
        gene2group_dict (dict): A dictionary mapping gene IDs to orthogroups. If provided, the function
                                converts randomly sampled genes to their corresponding orthogroups. Default is None.

    Returns:
        list: A list containing "nb_of_background_sets" sets, each representing a background set of genes.
              If gene2group_dict is provided, the genes are converted to orthogroups in each background set.
    """

    # Calculate the probability of finding a given DEG in the background (as it could be included in multiple clusters)
    deg2probability = {deg: (nb_of_deg_occurences/len(background_degs)) for deg, nb_of_deg_occurences in Counter(background_degs).items()}
    background_degs_unique, probabilities = zip(*deg2probability.items())

    # Initialize random generator with seed
    rng = np.random.default_rng(seed)

    # Generate a background model
    background_list = [None] * nb_of_background_sets
    for i in range(nb_of_background_sets):

        # Draw a set of random genes from the background DEGs, given the calculated probabilities
        background = rng.choice(background_degs_unique, cluster_size, p=probabilities, replace=False)

        # Add to background list
        background_list[i] = list(background)

    return background_list


def calculate_overlap_statistics(species1, species2, cluster1, cluster2, species2cluster2genes_and_groups, gene2group_dict=None, nb_of_background_sets=10000, compare_orthogroups=True, minimal_real_matches=1, all_clusters_as_background=False):
    
    """
    Calculate statistical measures of overlap between two clusters (gene sets or orthogroup sets).

    Parameters:
        species1 (str): The name of the reference species.
        species2 (str): The name of the species being compared to the reference species (query species).
        cluster1 (str): The name of the cluster in the reference species for comparison.
        cluster2 (str): The name of the cluster in the query species.
        species2cluster2genes_and_groups (dict): A dictionary with keys as species and values as dictionaries
                                                 mapping cluster names to sets of genes and Orthogroups.
        gene2group_dict (dict): A dictionary mapping gene IDs to orthogroups. Required if compare_families is False.
        nb_of_background_sets (int): The number of background sets to generate for statistical comparison. Default is 10000.
        minimal_real_matches (int): The minimum number of real matches required for meaningful statistical analysis. Default is 1.
        compare_orthogroups (bool): Flag to indicate whether to compare orthogroups (across species), or genes IDs (within species). 
                                    Default is True.
        all_clusters_as_background (bool): Flag to indicate whether to sample from all clusters as background, or all clusters,
                                           except for the cluster of interest. Default is False.

    Returns:
        tuple: A tuple containing three values:
            - real_matches (int): The number of matches between the real cluster and the query cluster.
            - enrichment_fold (float): The enrichment fold, indicating the degree of enrichment of real matches compared to background.
            - p_val (float): The p-value indicating the significance of the enrichment.

    """
    
    # Check that there is orthology info given when comparing orthogroups
    if compare_orthogroups:
        assert gene2group_dict is not None

    elements_to_compare = "groups" if compare_orthogroups else "genes"
    
    # Define all the DEGs (can be either DE genes, or DE groups) to draw from for the background, allowing duplicates
    all_cluster_degs = list()
    for cluster in natsorted(species2cluster2genes_and_groups[species1].keys() - {"all_clusters"}):
        
        # Don't include the reference cluster if all_clusters_as_background == False
        if all_clusters_as_background or (cluster != cluster1): 
            all_cluster_degs.extend(species2cluster2genes_and_groups[species1][cluster][elements_to_compare])
    
    # Size of cluster of interest (in reference species)
    cluster_size = len(species2cluster2genes_and_groups[species1][cluster1][elements_to_compare])
    
    # Calculate the matches of the real cluster, with the query cluster
    real_degs = species2cluster2genes_and_groups[species1][cluster1][elements_to_compare]
    query_degs = species2cluster2genes_and_groups[species2][cluster2][elements_to_compare]
    real_matches = sum(1 for group in real_degs if group in query_degs)

    assert nb_of_background_sets >= 100

    # Initialize 10 background sets of DEGs
    current_nb_of_background_sets = 10
    background_degs_list = get_background(all_cluster_degs, cluster_size, nb_of_background_sets=current_nb_of_background_sets)

    p_val_estimate = 1.0
    while current_nb_of_background_sets < nb_of_background_sets:
        
        # Add 9 times more background sets (scale up x10) unless that exceeds the given limit
        current_nb_of_background_sets = 10 * current_nb_of_background_sets
        if current_nb_of_background_sets > nb_of_background_sets: # above the limit
            additional_background_sets = int(nb_of_background_sets - current_nb_of_background_sets/10)
            current_nb_of_background_sets = nb_of_background_sets
        else:
            additional_background_sets = int(0.9 * current_nb_of_background_sets)
        background_degs_list.extend(get_background(all_cluster_degs, cluster_size, nb_of_background_sets=additional_background_sets, seed=len(background_degs_list)))

        # Calculate the matches of all background groups with the query cluster
        background_matches = list()
        for background_degs in background_degs_list:
            background_matches.append(sum(1 for deg in background_degs if deg in query_degs))
        
        # Calculate p-value
        times_background_above_real = sum(1 for matches in background_matches if matches >= real_matches)
        if times_background_above_real == 0:
            times_background_above_real = 0.9
        p_val_estimate = times_background_above_real / len(background_degs_list)

        # Overwrite p-value to 1 if there is not enough information
        if real_matches < minimal_real_matches:
            p_val_estimate = 1.0

        # Stop scaling up if p-value is already at sufficient precision (100 times higher than the lowest possible value)
        if p_val_estimate >= 100/len(background_degs_list):
            break
    
    p_val = p_val_estimate

    # Calculate enrichment fold
    mean_background_matches = mean(background_matches)
    if mean_background_matches == 0:
        mean_background_matches = 1 / nb_of_background_sets # make sure background never becomes zero
    enrichment_fold = round(real_matches / mean_background_matches, 2)

    # Gather the matching DEGs and corresponding genes
    matching_degs = [deg for deg in real_degs if deg in query_degs]
    if compare_orthogroups:
        matching_genes_species1 = [gene for gene in species2cluster2genes_and_groups[species1][cluster1]["genes"] if gene2group_dict[gene] in matching_degs]
        matching_genes_species2 = [gene for gene in species2cluster2genes_and_groups[species2][cluster2]["genes"] if gene2group_dict[gene] in matching_degs]
    else:
        # If the matching DEGs are not orthogroups but gene IDs, simply use those
        matching_genes_species1 = matching_degs
        matching_genes_species2 = matching_degs

    # Convert to lists
    matching_degs = ",".join(sorted(matching_degs))
    matching_genes_species1 = ",".join(sorted(matching_genes_species1))
    matching_genes_species2 = ",".join(sorted(matching_genes_species2))
    
    return real_matches, enrichment_fold, p_val, matching_degs, matching_genes_species1, matching_genes_species2


def correct_p_values(p_value_df):
    
    """
    Correct p-values for multiple testing using the Benjamini-Hochberg procedure.

    Parameters:
        p_value_df (pd.DataFrame): A DataFrame containing p-values for cluster comparisons.
                                   Rows represent clusters from species 1, columns represent clusters from species 2.

    Returns:
        pd.DataFrame: A DataFrame with adjusted p-values for multiple testing.
                      Rows and columns represent clusters from species 1 and 2, respectively.
    """
    
    # Convert wide format (clusters species 1 X clusters species 2) 
    # into long format (col1 = cluster1, col2 = cluster2, col3 = p-value)
    df_reset = p_value_df.reset_index()
    df_long_format = pd.melt(df_reset, id_vars=['index'], var_name='Column', value_name='Value')
    df_long_format.columns = ['cluster1', 'cluster2', 'p_val']
    
    # Calculate the adjusted p-values and add to long dataframe
    df_long_format["adj_p_val"] = multipletests(df_long_format['p_val'], method = 'fdr_bh')[1]
    
    # Put the adjusted p-values in similar dataframe as the original p-values
    adjusted_p_value_df = p_value_df.copy()
    for i in df_long_format.index:
        c1, c2, adj_p_val = df_long_format.loc[i,"cluster1"], df_long_format.loc[i,"cluster2"], df_long_format.loc[i,"adj_p_val"]
        adjusted_p_value_df.loc[c1,c2] = adj_p_val
    
    return adjusted_p_value_df


def get_all_statistics(species1, species2, species2cluster2genes_and_groups, gene2group_dict=None, nb_of_background_sets=10000, minimal_real_matches=1, compare_orthogroups=True, all_clusters_as_background=False, verbose=True):
    
    """
    Calculate overlap statistics for all combinations of clusters between two species.

    Parameters:
        species1 (str): The name of the reference species.
        species2 (str): The name of the species being compared to the reference species (query species).
        species2cluster2genes_and_groups (dict): A dictionary with keys as species and values as dictionaries
                                                 mapping cluster names to sets of genes and Orthogroups.
        gene2group_dict (dict): A dictionary mapping gene IDs to orthogroups. 
        nb_of_background_sets (int): The number of background sets to generate for statistical comparison. Default is 10000.
        minimal_real_matches (int): The minimum number of real matches required for meaningful statistical analysis. Default is 1.
        compare_families (bool): Flag to indicate whether to collapse all genes into gene families (orthogroups) instead of 
                                 matching all individual genes. Default is False.
        all_clusters_as_background (bool): Flag to indicate whether to sample from all clusters as background, or all clusters,
                                           except for the cluster of interest. Default is False.
        verbose (bool): Flag to indicate whether to print progress information during the calculation. Default is True.

    Returns:
        tuple: A tuple containing four DataFrames:
            - all_matches_df (pd.DataFrame): Number of matches between clusters from species 1 and 2.
            - all_enrichment_folds_df (pd.DataFrame): Enrichment fold for each cluster pair.
            - all_p_vals_df (pd.DataFrame): P-values for each cluster pair.
            - all_adj_p_vals_df (pd.DataFrame): Adjusted p-values (multiple testing correction) for each cluster pair.
    """
    
    # Extract the cluster names of both species
    all_reference_clusters = natsorted(species2cluster2genes_and_groups[species1].keys() - {"all_clusters"})
    all_query_clusters = natsorted(species2cluster2genes_and_groups[species2].keys() - {"all_clusters"})
    
    # Define the empty dataframes with cluster of species 1 as rows and clusters of species 2 as columns
    empty_df = pd.DataFrame(0, index=all_reference_clusters, columns=all_query_clusters)
    all_matches_df, all_enrichment_folds_df, all_p_vals_df = empty_df.copy(), empty_df.copy(), empty_df.copy()
    matching_degs_df, matching_genes_species1_df, matching_genes_species2_df = empty_df.copy(), empty_df.copy(), empty_df.copy()

    # Loop over each combination of clusters and calculate statistics for that combination
    for cluster1 in all_reference_clusters:
        if verbose:
            with open("log.txt", "a") as out_file:
                out_file.write("[RUNNING] cluster {}".format(cluster1)+"\n")
            print("[RUNNING] cluster {}".format(cluster1))
        for cluster2 in all_query_clusters:
            nb_of_matches, enrichment_fold, p_val, matching_degs, \
            matching_genes_species1, matching_genes_species2 = \
                calculate_overlap_statistics(species1, species2, cluster1, cluster2, \
                                             species2cluster2genes_and_groups, \
                                             gene2group_dict=gene2group_dict, \
                                             nb_of_background_sets=nb_of_background_sets, \
                                             minimal_real_matches=minimal_real_matches, \
                                             compare_orthogroups=compare_orthogroups, \
                                             all_clusters_as_background=all_clusters_as_background)
            # Fill in statistics in the results dataframe
            all_matches_df.loc[cluster1, cluster2] = nb_of_matches
            all_enrichment_folds_df.loc[cluster1, cluster2] = enrichment_fold
            all_p_vals_df.loc[cluster1, cluster2] = p_val
            matching_degs_df.loc[cluster1, cluster2] = matching_degs
            matching_genes_species1_df.loc[cluster1, cluster2] = matching_genes_species1
            matching_genes_species2_df.loc[cluster1, cluster2] = matching_genes_species2
            
    # Calculate adjusted p-values (multiple testing correction)
    all_adj_p_vals_df = correct_p_values(all_p_vals_df)
    
    return all_matches_df, all_enrichment_folds_df, all_p_vals_df, all_adj_p_vals_df, matching_degs_df, matching_genes_species1_df, matching_genes_species2_df


def calculate_two_way_statistics(species1, species2, species2cluster2genes_and_groups, gene2group_dict=None, nb_of_background_sets=1000, minimal_real_matches=2, compare_orthogroups=True, all_clusters_as_background=False, verbose=True):
    
    """
    Calculate two-way overlap statistics between two species, considering both directions of comparison.

    Parameters:
        species1 (str): The name of the first species for comparison.
        species2 (str): The name of the second species for comparison.
        species2cluster2genes_and_groups (dict): A dictionary with keys as species and values as dictionaries
                                                 mapping cluster names to sets of genes and Orthogroups.
        gene2group_dict (dict): A dictionary mapping gene IDs to orthogroups. Required if compare_families is False.
        nb_of_background_sets (int): The number of background sets to generate for statistical comparison. Default is 1000.
        minimal_real_matches (int): The minimum number of real matches required for meaningful statistical analysis. Default is 2.
        compare_families (bool): Flag to indicate whether to compare gene families instead of orthogroups. Default is True.
        all_clusters_as_background (bool): Flag to indicate whether to sample from all clusters as background, or all clusters,
                                           except for the cluster of interest. Default is False.
        verbose (bool): Flag to indicate whether to print progress information during the calculation. Default is True.

    Returns:
        dict: A dictionary containing four DataFrames for each direction of comparison (e.g. species 1 to species 2):
            key (str) = "species1-to-species2:
            value (tuple) = A tuple containing four DataFrames:
                - all_matches_df_mean (pd.DataFrame): Mean number of matches between clusters from both species.
                - all_enrichment_folds_df_mean (pd.DataFrame): Mean enrichment fold for each cluster pair.
                - all_p_vals_df_mean (pd.DataFrame): Mean p-values for each cluster pair.
                - all_adj_p_vals_df_mean (pd.DataFrame): Mean adjusted p-values (multiple testing correction) for each cluster pair.
    """
    
    # Calculate all statistics from species 1 to species 2
    if verbose:
        with open("log.txt", "a") as out_file:
            out_file.write("Calculating background sets for {}".format(species1)+"\n")
        print("Calculating background sets for {}".format(species1))
    all_matches_df1, all_enrichment_folds_df1, all_p_vals_df1, all_adj_p_vals_df1, \
    matching_degs_df1, matching_genes_species1_df1, matching_genes_species2_df1 = \
    get_all_statistics(species1, species2, \
                       species2cluster2genes_and_groups, \
                       gene2group_dict=gene2group_dict, \
                       nb_of_background_sets=nb_of_background_sets, \
                       minimal_real_matches=minimal_real_matches, \
                       compare_orthogroups=compare_orthogroups, \
                       all_clusters_as_background=all_clusters_as_background, \
                       verbose=verbose)
    
    # Calculate all statistics from species 2 to species 1
    if verbose:
        with open("log.txt", "a") as out_file:
            out_file.write("Calculating background sets for {}".format(species2)+"\n")
        print("Calculating background sets for {}".format(species2))
    all_matches_df2, all_enrichment_folds_df2, all_p_vals_df2, all_adj_p_vals_df2, \
    matching_degs_df2, matching_genes_species1_df2, matching_genes_species2_df2 = \
    get_all_statistics(species2, species1, \
                       species2cluster2genes_and_groups, \
                       gene2group_dict=gene2group_dict, \
                       nb_of_background_sets=nb_of_background_sets, \
                       minimal_real_matches=minimal_real_matches, \
                       compare_orthogroups=compare_orthogroups, \
                       all_clusters_as_background=all_clusters_as_background, \
                       verbose=verbose)
    
    # Put all statistics in a dict
    all_statistics_two_way = dict()
    all_statistics_two_way["{}-to-{}".format(species1, species2)] = (all_matches_df1, all_enrichment_folds_df1, \
                                                                     all_p_vals_df1, all_adj_p_vals_df1, \
                                                                     matching_degs_df1, matching_genes_species1_df1, \
                                                                     matching_genes_species2_df1)
    all_statistics_two_way["{}-to-{}".format(species2, species1)] = (all_matches_df2.T, all_enrichment_folds_df2.T, \
                                                                     all_p_vals_df2.T, all_adj_p_vals_df2.T, \
                                                                     matching_degs_df2.T, matching_genes_species1_df2.T, \
                                                                     matching_genes_species2_df2.T)
    
    return all_statistics_two_way


def generate_visualization(overlap_statistics, species1, species2, cmap="viridis"):
    
    """
    Generate a visualization of overlap statistics between clusters from two species.

    Parameters:
        overlap_statistics (dict): A dictionary containing four DataFrames for each direction of comparison (e.g. species 1 to species 2):
            key (str) = "species1-to-species2:
            value (tuple) = A tuple containing four DataFrames:
                - all_matches_df_mean (pd.DataFrame): Mean number of matches between clusters from both species.
                - all_enrichment_folds_df_mean (pd.DataFrame): Mean enrichment fold for each cluster pair.
                - all_p_vals_df_mean (pd.DataFrame): Mean p-values for each cluster pair.
                - all_adj_p_vals_df_mean (pd.DataFrame): Mean adjusted p-values (multiple testing correction) for each cluster pair.
        species1 (str): The name of the first species for comparison.
        species2 (str): The name of the second species for comparison.

    Returns:
        matplotlib.figure.Figure: The generated matplotlib figure.
    """
    
    # Define subplots
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(10, 12))
    
    # Plot first -log(q-value) plot
    sns.heatmap(-np.log(overlap_statistics["{}-to-{}".format(species1, species2)][3]), ax=ax1, cmap=cmap, cbar_kws={"label": "- log(q-value)"})
    ax1.set(xlabel="Cluster in {}".format(species2), ylabel="Cluster in {}".format(species1))
    ax1.set_title("Significance of {} clusters enrichment\nfor {} clusters".format(species1, species2))
    
    # Plot second -log(q-value) plot
    sns.heatmap(-np.log(overlap_statistics["{}-to-{}".format(species2, species1)][3]), ax=ax2, cmap=cmap, cbar_kws={"label": "- log(q-value)"})
    ax2.set(xlabel="Cluster in {}".format(species2), ylabel="Cluster in {}".format(species1))
    ax2.set_title("Significance of {} clusters enrichment\nfor {} clusters".format(species2, species1))
    
    # Plot first binary significance plot
    sns.heatmap(overlap_statistics["{}-to-{}".format(species1, species2)][3] < 0.05, ax=ax3, cmap=cmap, cbar=True)
    ax3.set(xlabel="Cluster in {}".format(species2), ylabel="Cluster in {}".format(species1))
    ax3.set_title("Is overlap significant?")
    
    # Plot first binary significance plot
    sns.heatmap(overlap_statistics["{}-to-{}".format(species2, species1)][3] < 0.05, ax=ax4, cmap=cmap, cbar=True)
    ax4.set(xlabel="Cluster in {}".format(species2), ylabel="Cluster in {}".format(species1))
    ax4.set_title("Is overlap significant?")
    
    # Plot first enrichment plot
    sns.heatmap(overlap_statistics["{}-to-{}".format(species1, species2)][1], ax=ax5, cmap=cmap, vmin=0, vmax=10, cbar_kws={"label": "Enrichment fold"})
    ax5.set(xlabel="Cluster in {}".format(species2), ylabel="Cluster in {}".format(species1))
    ax5.set_title("Enrichment fold of {} clusters\nfor {} clusters".format(species1, species2))
    
    # Plot second enrichment plot
    sns.heatmap(overlap_statistics["{}-to-{}".format(species2, species1)][1], ax=ax6, cmap=cmap, vmin=0, vmax=10, cbar_kws={"label": "Enrichment fold"})
    ax6.set(xlabel="Cluster in {}".format(species2), ylabel="Cluster in {}".format(species1))
    ax6.set_title("Enrichment fold of {} clusters\nfor {} clusters".format(species2, species1))
    
    # Clean up and save plot
    fig.tight_layout()
    fig.savefig("{}-{}_overlap_significance.png".format(species1, species2), bbox_inches="tight")
    
    return fig


def remove_two_way_redundancy(all_comparisons):
    
    """
    Remove the two-way redundancy in the comparisons (remove the reverse comparison).

    Parameters:
        all_comparisons (iterable): An iterable containing all comparisons.

    Returns:
        unique_comparison_couples (list): A list of unique comparison couples with redundancy removed.
    """
    
    remaining_comparisons = set(all_comparisons)
    unique_comparison_couples = list()
    
    for comparison in all_comparisons:

        # If in remaining comparisons, add to list and remove from remaining comparisons
        if comparison in remaining_comparisons:
            unique_comparison_couples.append(comparison)
            remaining_comparisons.remove(comparison)

        # Remove the reverse comparison from remaining comparisons
        species1, species2 = comparison.split("-to-")
        comparison_reverse = species2 + "-to-" + species1
        if comparison_reverse in remaining_comparisons:
            remaining_comparisons.remove(comparison_reverse)
    
    return unique_comparison_couples


def combine_and_write_statistics_to_excel(overlap_statistics_long_dict, output_file):
    
    """
    Combine (take the best of each two-way comparison) and write overlap statistics to an Excel file with multiple sheets.

    Parameters:
        overlap_statistics_long_dict (dict): A dictionary with the long format DataFrames for each comparison.
        output_file (str): Output file name

    Returns:
        None: The function writes the Excel file directly.
    """
    
    # Initialize combined long table dict
    overlap_statistics_long_combined_dict = dict()

    # Remove the two-way redundency in the comparisons (remove the reverse comparison)
    unique_comparison_couples = remove_two_way_redundancy(overlap_statistics_long_dict.keys())

    for comparison in unique_comparison_couples:

        # Get both dataframes 
        statistics_df1 = overlap_statistics_long_dict[comparison]
        species1, species2 = comparison.split("-to-")
        comparison_reverse = species2+"-to-"+species1
        statistics_df2 = overlap_statistics_long_dict[comparison_reverse]

        # Define result dataframe
        statistics_combined = statistics_df1.copy()

        for i in statistics_combined.index:

            # Get the q-values and enrichment folds
            q_val1, q_val2 = statistics_df1.loc[i,"q-value"], statistics_df2.loc[i,"q-value"]
            enrich1, enrich2 = statistics_df1.loc[i,"enrichment fold"], statistics_df2.loc[i,"enrichment fold"]

            if (q_val1 < q_val2) or ((q_val1 == q_val2) and (enrich1 > enrich2)):
                pass # Do not overwrite the values originating from statistics_df1
            else:
                # Overwrite the statistics_df1 values with the statistics_df2 values (only for the statistics)
                statistics_combined.loc[i,"q-value"] = q_val2
                statistics_combined.loc[i,"enrichment fold"] = enrich2
                statistics_combined.loc[i,"p-value"] = statistics_df2.loc[i,"p-value"]
                statistics_combined.loc[i,"significant"] = statistics_df2.loc[i,"significant"]

        # Add combined table to the dict
        overlap_statistics_long_combined_dict[comparison] = statistics_combined
    
    # Open Excel writer to write to multiple Excle sheets
    with pd.ExcelWriter(".".join(output_file.split(".")[:-1])+"_best_hits.xlsx") as writer:
        
        # Write all combined comparisons to an Excel
        for comparison, overlap_statistics_long_combined in overlap_statistics_long_combined_dict.items():
            
            # Select significant hits and sort them
            statistics_df = overlap_statistics_long_combined \
                                .query("significant") \
                                .sort_values(["q-value", "enrichment fold"], ascending=(True, False))
            
            # Get the best match for each cluster, both ways, then merge them again
            statistics_df1 = statistics_df.drop_duplicates(statistics_df.columns[0], keep="first")
            statistics_df2 = statistics_df.drop_duplicates(statistics_df.columns[1], keep="first")
            statistics_df = pd.concat([statistics_df1, statistics_df2])
            statistics_df = statistics_df.drop_duplicates()
            
            # Write the final combined dataframe to a sheet in the Excel file
            statistics_df.to_excel(writer, index=False, sheet_name=comparison)

            
def write_statistics_to_excel(overlap_statistics, significance_cutoff=0.05):
    
    """
    Write overlap statistics to an Excel file with multiple sheets.

    Parameters:
        overlap_statistics (dict): A dictionary containing four DataFrames for each direction of comparison (e.g. species 1 to species 2):
            key (str) = "species1-to-species2:
            value (tuple) = A tuple containing four DataFrames:
                - all_matches_df_mean (pd.DataFrame): Mean number of matches between clusters from both species.
                - all_enrichment_folds_df_mean (pd.DataFrame): Mean enrichment fold for each cluster pair.
                - all_p_vals_df_mean (pd.DataFrame): Mean p-values for each cluster pair.
                - all_adj_p_vals_df_mean (pd.DataFrame): Mean adjusted p-values (multiple testing correction) for each cluster pair.
        significance_cutoff: q-value cutoff at which an overlap is considered significant. Default is q-value < 0.05.

    Returns:
        None: The function writes the Excel file directly.

    Excel File Structure:
        The Excel file will have multiple sheets, each corresponding to a species-to-species comparison.
        Each sheet contains a long-format DataFrame with columns:
            - "{} cluster" (str): Cluster name from the first species.
            - "{} cluster" (str): Cluster name from the second species.
            - "number of matches" (float): Mean number of matches between clusters.
            - "enrichment fold" (float): Mean enrichment fold for each cluster pair.
            - "p-value" (float): Mean p-values for each cluster pair.
            - "q-value" (float): Mean adjusted p-values (multiple testing correction) for each cluster pair.
            - "significant" (boolean): True if q-value is significant, false otherwise
    """
    
    # Open Excel writer to write to multiple Excle sheets
    with pd.ExcelWriter("cluster_deg_overlap_statistics.xlsx") as writer:
        
        # Initialize long table dict
        overlap_statistics_long_dict = dict()
        
        # For each species-to-species comparison in the input dict, make a dataframe to write to Excel
        for comparison in overlap_statistics.keys():
            
            # Get species from species-to-species comparison (str)
            species1, species2 = comparison.split("-to-")
            
            # Convert all individual dataframes from wide to long format
            overlap_statistics_long = None
            for i, statistic_name in enumerate(["number of matches", "enrichment fold", "p-value", "q-value", "match IDs", "matching genes 1", "matching genes 2"]):
                df = overlap_statistics[comparison][i]
                df_reset = df.reset_index()
                df_long_format = pd.melt(df_reset, id_vars=["index"], var_name="Column", value_name="Value")
                df_long_format.columns = ["{} cluster".format(species1), "{} cluster".format(species2), statistic_name]
                
                # Combine different statistics together in the final dataframe as extra columns
                if i == 0:
                    overlap_statistics_long = df_long_format
                else:
                    overlap_statistics_long[statistic_name] = df_long_format[statistic_name]
                
                # Add a column with a boolean value, showing if the overlap is significant
                if statistic_name == "q-value":
                    overlap_statistics_long["significant"] = overlap_statistics_long[statistic_name] < significance_cutoff
            
            # Add long format table to dict
            overlap_statistics_long_dict[comparison] = overlap_statistics_long
            
            # Write the final dataframe to a sheet in the Excel file
            overlap_statistics_long.to_excel(writer, index=False, sheet_name=comparison)
    
    # Combine the statistics of both two-way comparisons and write results to Excel
    combine_and_write_statistics_to_excel(overlap_statistics_long_dict, "cluster_deg_overlap_statistics_combined.xlsx")

            