
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Antimicrobial resistome data analysis with R/Bioconductor

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/NCBI-Codeathons/amr-2024-team-lahti/branch/main/graph/badge.svg)](https://codecov.io/gh/NCBI-Codeathons/amr-2024-team-lahti?branch=main)
[![Slack](https://img.shields.io/badge/Slack-Join-blue.svg?logo=slack)](https://ncbiamrcodeathon.slack.com/x-p7700037081187-7765806765012-7763213820995/archives/C07L5HA3UFR)
[![Watch on
GitHub](https://img.shields.io/github/watchers/NCBI-Codeathons/amr-2024-team-lahti.svg?style=social)](https://github.com/NCBI-Codeathons/amr-2024-team-lahti/watchers)
[![Star on
GitHub](https://img.shields.io/github/stars/NCBI-Codeathons/amr-2024-team-lahti.svg?style=social)](https://github.com/NCBI-Codeathons/amr-2024-team-lahti/stargazers)
<!-- badges: end -->

## Team roles

- **Team Lead** / [Leo Lahti](https://github.com/antagomir) Conveys
  scientific objectives to the team, and coordinates work.

- **Tech Lead** / [Akewak Jeba](https://github.com/ake123) & [Muluh
  Geraldson](https://github.com/Daenarys8) Coordinates software
  installation and data acquisition, manages version control and the
  team’s GitHub repository, troubleshoots technical issues with tech
  support

- **Writer** / [Dattaray Mongad](https://github.com/microDM) & [Nitin
  Bayal](https://github.com/nixonbyl) Ensures that all work is
  documented, manages GitHub README and Team Report

- **Flex** / [Mahkameh Salehi](https://github.com/mahkamehsalehi) &
  [Shivang ](https://github.com/ShivangPB) Fills various roles and
  assumes responsibilities for tasks as needed

## Project goals

*This project aims to develop
[R/Bioconductor](https://www.bioconductor.org/) tools to analyze AMR
data by integrating publicly available human gut metagenomes with ARG
profiles (e.g.,
[metagenomics](https://www.nature.com/subjects/metagenomics)) and
epidemiological, transmission data to study AMR gene co-migration.*

## Approach

Reproducible data science workflow and dissemination material to
demonstrate the utilization of the latest Bioconductor multi-assay data
integration framework in the context of resistome profiling in large
population studies. This covers dedicated data structures and the
complementary data analysis and visualization methods. The outcomes are
reported via an automated pkgdown website.

## Overall analysis plan

1.  *Prepare data*
    - Prepare the [Lee et
      al. (2023)](https://github.com/microbiome/data/tree/main/Lee2023)
      data set for analysis in R/Bioc and provide sufficient
      documentation in the Github repository. This can be done by
      adapting and refining the previously created data set.
2.  *Prepare documentation*
    - Create minimal website for this project using pkgdown (perhaps we
      need to convert the repo into a minimal R package to do this?
      Should be easy, a pkg template is readily available via RStudio)

    - Prepare summary of the data set to be presented on a slide set;
      for examples, you can see here for examples (“Example data”) and
      their source code -\> render the slide set as html, we can later
      see how to publish it.

    - Create one or more tutorials on resistome data analysis using the
      latest R/Bioconductor methods and release them via the Github
      repository (vignettes/ subfolder); these can be browsed via the
      github homepage (created in point (2))

    - Motivate the use of Bioconductor methods.
3.  *Analyse the resistome data set*
    - Provide exploratory summaries of the data set using R/Bioc
      methods.
    - Replicate of some key results from the manuscript using
      R/Bioconductor methods to verify the data.
    - Add new analyses (TBD)
4.  *Develop methodology*
    - Consider automation of specific analysis or visualization methods
      as functions that can be share via new R pkg, or existing
      packages.
5.  *Finalize workflow*
    - Polish the data analysis workflow/s & documentation so that they
      can be used as training material on resistome analyses using
      R/Bioconductor.

## Results

The main outcome of this work is the reproducible project site including
data science workflows for resistome analysis. This forms the basis for
an extendable collection of resistome analysis workflows for reanalysis
of data for new discoveries, methods development and educational
purposes.

## Future work

Future work will standardize and automate some of the key tools, and
release them through open source software libraries such as
*[mia](https://microbiome.github.io/mia/)* in R/Bioconductor.

## NCBI codeathon disclaimer

This software was created as part of an NCBI codeathon(See
[details](https://www.nlm.nih.gov/ncbi/workshops/2024-09_AMR-Codeathon/codeathon-details.html))
, a hackathon-style event focused on rapid innovation. While we
encourage you to explore and adapt this code, please be aware that NCBI
does not provide ongoing support for it.

For general questions about NCBI software and tools, please visit: [NCBI
Contact Page](https://www.ncbi.nlm.nih.gov/home/about/contact/)
