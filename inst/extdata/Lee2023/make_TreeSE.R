library(mia)
library(dplyr)

# Load country-level data
countrydata <- read.table(
  "data/country_AMU_data.tsv",
  header = TRUE,
  sep = "\t",
  # row.names = 1,
  check.names = FALSE)
countrydata$Country <-  factor(countrydata$Country)

## loading metadata
col_data <- read.table(
  "data/sfile-microbiome-sample-detail.tsv",
  header = TRUE, sep = "\t", check.names = FALSE)
# Making a new ID variable to match the assembly data with the rest of the data
col_data <- col_data %>% mutate(tmp_ID = paste0(Study, "__", Sample_ID))
rownames(col_data) <- col_data$tmp_ID
col_data$tmp_ID <- NULL

## Loading the taxonomic abundance table
assay <- read.table(
  "data/samples_with_1k_core_COG_ORFs.adult_stool.SGB_level.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
# helper list for subsetting and renaming rownames
new_ids <- lapply(rownames(assay), function(old_id) {
  rownames(col_data[which(old_id==col_data$Sample_ID),]) }) %>%
  setNames(rownames(assay)) %>% unlist()
# subsetting
assay <- assay[names(new_ids),]
# renaming
rownames(assay) <- new_ids %>% unname()

# Convert to numeric matrix
assay <- as.matrix(assay)
# Replace NAs with zeroes
assay[is.na(assay)] <- 0
# Remove rows and cols with zero variance
assay <- assay[rowSds(assay)>0, colSds(assay)>0]



## Loading rowData
library(data.table)
row_data <- fread("data/SGB_info.taxonomy.tsv", sep="\t")
rownames(row_data) <- as.character(row_data$sgbid)
# renaming taxonomy ranks 
colnames(row_data)[3:9] <- c("Kingdom", "Phylum", "Class", "Order",
                             "Family", "Genus", "Species")

row_data <- as.data.frame(row_data)

## common taxa
common_taxa <- intersect(rownames(row_data), colnames(assay))

## Loading Read-based normalized ARG abundance table
assay_read <- read.table(
  "data/cpg_table.20220516_scg_valid_only.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
# helper list for subsetting and renaming rownames
new_ids <- lapply(rownames(assay_read), function(old_id) {
  rownames(col_data[which(old_id==col_data$Sample_ID),]) }) %>%
  setNames(rownames(assay_read)) %>% unlist()
# subsetting
assay_read <- assay_read[names(new_ids),]
# renaming
rownames(assay_read) <- new_ids %>% unname()


# Convert to numeric matrix
assay_read <- as.matrix(assay_read)
# Replace NAs with zeroes
assay_read[is.na(assay_read)] <- 0
# Remove rows and cols with zero variance
assay_read <- assay_read[rowSds(assay_read)>0, colSds(assay_read)>0]



## Loading  assembly-based normalized ARG abundance table
assay_assembly <- read.table(
  "data/DS3.SCG_normalized_ARG_abund.columns_CARD_ref.n_6104.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)


# Convert to numeric matrix
assay_assembly <- as.matrix(assay_assembly)
# Replace NAs with zeroes
assay_assembly[is.na(assay_assembly)] <- 0
# Remove rows and cols with zero variance
assay_assembly <- assay_assembly[rowSds(assay_assembly)>0, colSds(assay_assembly)>0]


## Loading rowData for arg read data
row_data_arg <- read.table(
  "data/ARG_family_info.tsv",
  header = TRUE, sep = "\t",
  row.names = 1, # ARG_family
  check.names = FALSE)
## common arg family
common_arg <- Reduce(intersect, list(rownames(row_data_arg),
                                     colnames(assay_read),
                                     colnames(assay_assembly)))

## Common ids
common_ids <- Reduce(intersect, list(rownames(assay_read),
                                     rownames(assay),
                                     rownames(assay_assembly),
                                     rownames(col_data)))





## Constructing the TreeSE object with the read assay and metadata, with an
# alternative experiment object for assembly data.

tse <- TreeSummarizedExperiment(
  assays=SimpleList(abundances=t(assay[common_ids, common_taxa])),
  rowData=row_data[common_taxa, ],
  colData=col_data[common_ids,],
  metadata=list(Country=countrydata)
)
altExp(tse, "read") <- TreeSummarizedExperiment(
  assays = SimpleList(abundances=t(assay_read[common_ids, common_arg])),
  rowData = row_data_arg[common_arg, ])
altExp(tse, "assembly") <- TreeSummarizedExperiment(
  assays = SimpleList(abundances=t(assay_assembly[common_ids, common_arg])),
  rowData = row_data_arg[common_arg, ])



# Some manipulations
tse$antibiotic_current_use_binary <- factor(tse$antibiotic_current_use_binary, levels=c("no", "yes"))

# Convert taxa abundances to relative abundances and remove the original version
tse <- transformAssay(tse, assay.type="abundances", method="relabundance")
assay(tse, "abundances") <- NULL


## Saving the object
saveRDS(tse, "Lee2023.rds")
#rm(common_ids, col_data, assay_assembly, assay_read, new_ids, assay, row_data,
#   row_data_arg, common_arg, common_taxa, countrydata)


