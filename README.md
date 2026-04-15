 <!-- badges: start -->
  ![Static Badge](https://img.shields.io/badge/MBNeuroLab-darkblue)
  ![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
  ![Python version](https://img.shields.io/badge/python-3.11%20%7C%203.12-blue)
  ![R package check](https://github.com/eriklarsen4/Endo/actions/workflows/r-package-check.yml/badge.svg)
  ![Python package check](https://github.com/eriklarsen4/Endo/actions/workflows/python-tests.yml/badge.svg)
  
 <!-- badges: end -->

**Endo**

Repository for the [TMEM184B Endoylysosomal Acidification paper](<https://journals.biologists.com/jcs/article/doi/10.1242/jcs.263908/368440/Transmembrane-protein-184B-TMEM184B-modulates>)

Adapted from source code used to analyze data in **Dr. Martha Bhattacharya's lab**
at the **University of Arizona**

  + The lab studies transmembrane protein 184b (TMEM184B), a protein involved in
  axon degeneration and cancer, across multiple neuroscience model systems

The [Journal of Cell Science publication](<https://journals.biologists.com/jcs/article/138/15/jcs263908/368852/TMEM184B-modulates-endolysosomal-acidification-via?searchresult=1>) associated with this repository is the third to utilize function(s) from the
[TMEM](<https://github.com/eriklarsen4/TMEM>) package. The first publication has
its own dedicated repository, 
[eriklarsen4/Itch](<https://github.com/eriklarsen4/Itch>),
and the second publication's repository, `eriklarsen4/Hippo`, is under
development. Please see those repositories for each paper's
respective URL. Among other datasets and analyses, they contain multiple, functionally validated bulk RNAseq datasets.

Please see the [`FIREpHly` microscopy data vignette](<https://github.com/eriklarsen4/Endo/blob/dev/Rpkg/vignettes/FIREpHly.md>) and the [immunoprecipitation- label-free mass spec (`IPMS`) data vignette](<https://github.com/eriklarsen4/Endo/blob/dev/Rpkg/vignettes/IPMS.md>) for more detail on each dataset's data cleaning and analysis.

Additional computational biology work on the published data is outlined in broad strokes [here](<https://github.com/eriklarsen4/Endo/tree/dev/python#readme>), with other vignettes in progress. This work expands on the "gold standard" mixture modeling approach in affinity-precipitation mass spectrometry (`APMS` or `IPMS`) to infer proteins that truly interact with a target protein from a mixture of proteins.

Three packages are associated with this repository:

 1. A custom R package, [TMEM](<[https://github.com/eriklarsen4/TMEM](https://github.com/eriklarsen4/TMEM/pkgs/container/tmem)>), used for downstream gene set analyses and visualizations
 2. This repository's R package, [Endo](<https://github.com/eriklarsen4/Endo/pkgs/container/endo>), which gathers the data used for the publication
 3. A python package, [saint](<https://github.com/eriklarsen4/Endo/tree/dev/python>) (name change is forthcoming), which enables "plug-and-play" analysis of protein-protein interactions (`interactomes`)
    + this package contains two protein-interaction inference pipelines-- the `SAINT` algorithm, and a more complex iteration of it-- enabling benchmarking relative to the original gold-standard in interactomics
