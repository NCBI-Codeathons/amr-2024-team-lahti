#' \code{teamlahti} Package.
#'
#' \code{teamlahti} aims to develop [R/Bioconductor](https://www.bioconductor.org/) 
#' tools to analyze AMR data by integrating publicly available human gut 
#' metagenomes with ARG profiles (e.g. 
#' [metagenomics](https://www.nature.com/subjects/metagenomics)) and 
#' epidemiological, transmission data to study AMR gene co-migration.
#'
#' @name teamlahti-package
"_PACKAGE"
NULL


#' teamlahti datasets
#'
#' teamlahti provides Lee2023 dataset derived from independent experimental 
#' studies.The dataset represents instances of the TreeSummarizedExperiment to
#' serve as tools to practice the teamlahti functionality.
#' 
#' Currently, the following datasets are available:
#' \itemize{
#'   \item{\code{\link{Lee2023}}: A TreeSummarizedExperiment with 9,251 samples}
#' }
#' 
#' @name teamlahti-datasets
#' @docType data
#' @keywords datasets
#' 
#' @examples
#' # Load dataset from teamlahti
#' library(teamlahti)
#' data("Lee2023", package = "teamlahti")
#' 
#' # In this case, the dataset is a TreeSE, so it is renamed as tse
#' tse <- Lee2023
#' 
#' # Print summary
#' tse
NULL

#'  Metagenome assembly data
#' 
#' Lee2023 includes a total of 9,251 human microbiome samples from multiple body 
#' sites (7718 stool, 783 oral cavity, 410 skin, 150 airway, 93 nasal cavity, 
#' 88 vagina, and 9 milk) which had been sequenced in various studies.
#'
#' @format A TreeSummarizedExperiment with 4182 features and 5195 samples. The
#' rowData contains taxonomic information at Kingdom, Phylum, Class, Order,
#' Family and Genus level. See [dataset summary](https://ncbi-codeathons.github.io/amr-2024-team-lahti/articles/data_summary.html).
#' 
#' @name Lee2023
#' @docType data
#' @keywords datasets
#' @usage data(Lee2023)
#' @seealso \code{\link{teamlahti-datasets}}
#' @author Lee, K., Raguideau, S., Sirén, K. et al.
#' @references
#' Lee, K., Raguideau, S., Sirén, K. et al.
#' Population-level impacts of antibiotic usage on the human gut microbiome. 
#' Nat Commun 14, 1191 (2023). \url{https://doi.org/10.1038/s41467-023-36633-7} 
#' 
NULL
