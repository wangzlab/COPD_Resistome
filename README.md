# COPD_Resistome

This repository contains pipeline for the resistome analysis on sputum metagenomic data for COPD patients, as well as related R script for data visualization.

This pipeline integrates the blast6out results from ARGsOAP v2.0 stage two analysis and taxonomic annotation (by Kraken2) of the extracted ARG-related reads from ARGsOAP v2.0 stage one analysis. 

The pipeline generates read counts of each ARG, and calculates the read counts contributed by from each microbial taxa.
