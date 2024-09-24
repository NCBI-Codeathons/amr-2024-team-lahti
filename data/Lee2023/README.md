18.4.2024 Chouaib & Leo


## Using the TreeSummarizedExperiment R data object

Just download Lee2023.rds

This contains the data in a readily processed format.

## Constructing the TreeSummarizedExperiment R data obhect

In order to create the Lee2023.rds data yourself, you can do the
following.

Download these files:
- data.zip
- make_TreeSE.R

unzip the data.zip

make_TreeSE.R creates TSE object
- NAs have been replaced with zeroes
- zero variance rows and cols have been removed


## Misc

analysis.R contains additional, unfinished scripts to replicate analyses in the original
publication (Lee et al., 2023, see below). This is to be finalized.

## Data

This folder includes data downloaded from the following paper

Lee et al. (2023):
https://www.nature.com/articles/s41467-023-36633-7

The data files are stored in gzipped format in the data/ subfolder

And these URL sources:
 
1) Sample metadata
https://github.com/kihyunee/gut_resistotype/blob/data/sfile-microbiome-sample-detail.tsv

2) read-based normalized ARG abundance table
https://zenodo.org/records/7383076/files/cpg_table.20220516_scg_valid_only.tsv

3) assembly-based normalized ARG abundance table
https://zenodo.org/records/7383076/files/DS3.SCG_normalized_ARG_abund.columns_CARD_ref.n_6104.tsv

**The following files were received from Kihyun Lee on April 19, 2024 (personal communication).**

4) taxonomic abundance table
samples_with_1k_core_COG_ORFs.adult_stool.SGB_level.tsv

5) rowData for taxonomic abundances
SGB_info.taxonomy.tsv

6) rowData for ARG info 
ARG_family_info.tsv

7) Country-level data. Country-level antibiotic usage rate statistics from CDDEP and WHO. - citation for CDDEP = ResistanceMap 2019: Antibiotic resistance. The Center for Disease, Dynamics Economics & Policy. Accessed: 30-07-2019 (2019). Citation for WHO = WHO report on surveillance of antibiotic consumption: 2016-2018 early implementation. Geneva: World Health Organization (2018).
country_AMU_data.tsv
