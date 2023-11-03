# IPMAV Protocol IBt, UNAM
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![DOI](https://zenodo.org/badge/712122995.svg)](https://zenodo.org/doi/10.5281/zenodo.10070393)

## Description: 
This code includes the implementation of our methodology for analyzing SARS-CoV-2 allelic datasets for the identification of Intra Patient Minor Allelic Variants (in English IPMAVs, in Spanish VAMIPs) to detect coinfections in high-throughput sequencing samples
The manuscript is currently under review.

# IMPORTANT NOTE: Most of the time, the term VAMIPs is used in the scripts

# Table of contents:
Preprocessing (bash):
1) Process iVar tables and build allelic map
 * 01_cluster_commands.sh
Actual analyses (R):
2) Explore the IPMAVs in the allelic map
 * 02_explore_vamips.R
3) Single mutation analyses
 * 03_Analyze_vamip_contingency_Bar_n_boxplots.R
4) Haplotype construction
 * 04_2023-05-03_VAMIPS_per_genome_comparison.R
5) Detect coinfections
 * 05_Extract_putative_coinfections.R
6) Calculate co-occurring variants in CoViGen-Mex report
 * 06_Variant_Analysis_extra_cooccurrence_info_for_VAMIPs.R
7) Filter putative coinfections and plot
 * 07_Analyze_variant_completeness_per_genome.R


## How to run the scripts

To run the code:
- Download the full content of the directory containing this README file.
- Standard bash (we developed the code under GNU bash, version 5.1.16(1) 2020 on Arch Linux, kernel 5.15.133)
- R v4.0 or higher is required (we developed the code under version:  4.2.3 (2023-03-15)
- Most of the code was written on Rbase for maximum compatibility, except for R libraries "vioplot" and "stringr"

## Credits
The R scripts in this directory have been written by Rodrigo García-López for Carlos Arias's Virology Group at IBt, UNAM, Cuernavaca, Mexico as part of the Mexican Consortium for Genomic surveillance (CoViGen-Mex).
We acknowledge the Mexican Consortium for Genomic surveillance Mexicano de Vigilancia Genómica (CoViGen-Mex), the Mexican Institute of Social Security (IMSS) and the authors from the originating laboratories that collected the samples, the sequencing laboratories that produced the genomic sequences, available at GISAIDs Epicov database (EPI_SET ID: EPI_SET_231031gz for full list of authors). We appreciate the computer assistance provided by Jerome Verleyen, and Juan Manuel Hurtado, as well as the project LANCAD-UNAM-DGTIC-396 of the Dirección General de Cómputo y Tecnologías de la Información (DGTIC-UNAM) which provided supercomputing resources in MIZTLI.

### FULL CITATION
Coming soon
