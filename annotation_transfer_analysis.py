#!/usr/bin/env python3

# Define input file locations
file_locs = {"wheat":"input/wheat_log2fc0.5_fdr0.05.csv", "rice":"input/rice_log2fc0.5_fdr0.05.csv"}

# Define parameters
compare_orthogroups=True # True when comparing between species, False when comparing within the same species
separator = ","
nb_of_background_sets = 1000
minimal_real_matches = 2
significance_cutoff = 0.05

# Import functions from bin folder
import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bin'))
from significance_deg_overlap import *

# Run the analysis
run_overlap_analysis(
    file_locs, 
    separator, 
    compare_orthogroups, 
    nb_of_background_sets, 
    minimal_real_matches, 
    significance_cutoff
)
