---
title: "Spectral Clustering"
author: "Erik Larsen"
date: "`r Sys.Date()`"
always_allow_html: true
output:
  html_document:
    toc: true
    fig_width: 10
    fig_height: 10
  pdf_document:
    fig_width: 10
    fig_height: 10
    latex_engine: xelatex
  md_document:
    variant: gfm
---

# Overview

This markdown was developed to formally document an ad-hoc cluster analysis for 
TMEM184B immunoprecipitation tandem mass-spec experiments for 
**Dr. Martha Bhattacharya's lab** at the **University of Arizona** for the
endolysosomal paper under revision (2025-02-28) for the Journal of Cell Science

This vignette documents the `Spectral Clustering` of this IPMS data

```{r reset-shiny-mode, include=FALSE}
options(htmlwidgets.shinyMode = FALSE)
```

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.width = 10,
  fig.height = 12,
  dev = c("png", "pdf"),
  fig.path = "~/R/SpectralClustering/"
)
```

## Attach Packages

```{r Attach Packages, message=F,warning=F}
# data/wrangling
library(tidyverse)
library(Endo)
library(DT)

# clustering
library(pheatmap)
library(dynamicTreeCut)
  ## (kNN)
library(FNN)

# GSEA/viz
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(AnnotationDbi)
library(cowplot)

# progress bar
library(progressr)

# parallelization
library(future)
library(future.apply)

# network viz
library(visNetwork)

options(stringsAsFactors = FALSE)
```

## Load the Pre-Computed Data

```{r Data Load}
# Load precomputed objects from preprocessing
load("/R/precompute_objects.rda")
```

```{r Second Data Load, echo = FALSE}
rm(IPMS_mat_log)
rm(IPMS_mat_Z)
hSAINT <- hSAINT |> 
  dplyr::mutate(stat_class = dplyr::case_when(score > 0.5 & 
                                                lit_flag + CRAPome_flag <= 1 ~ 'High SAINT, Low Contamination',
                                              score < 0.5 &
                                                lit_flag + CRAPome_flag <= 1 ~ 'Low SAINT, Low Contamination',
                                              score > 0.5 &
                                                lit_flag + CRAPome_flag > 1 ~ 'High SAINT, High Contamination',
                                              score < 0.5 &
                                                lit_flag + CRAPome_flag > 1 ~ 'Low SAINT, High Contamination'))
```

## Transform the Data to Distance Matrix

```{r Compute Distance Matrix}
IPMS_dist <- 1 - stats::cor(t(IPMS_mat),
                            use = 'pairwise.complete.obs') |> 
  as.matrix()
```


```{r, results='asis', echo=FALSE}
fmt <- knitr::pandoc_to()

if (identical(fmt, "html")) {
  cat('# Spectral Clustering {.tabset .tabset-pills .tabset-fade}\n')
} else if (identical(fmt, "markdown")) {
  cat('<details><summary><strong>Spectral Clustering</strong></summary>\n')
} else {
  cat('# Spectral Clustering\n')
}
```

Spectral clustering *uses the same correlation-based distance representation as [hierarchical clustering](HierarchicalClustering.Rmd).*

We compute the pairwise distance matrix as 
**1 - Pearson correlation between proteins**, matching the above vignette, exactly.

This method was chosen because:

  + the wild variation in protein abundance detected by spectrometers in these experiments
  + contamination proteins often inflate abundances
  + immunoprecipitation (bait) biases (differential abundance for certain proteins, depending on the antibody + target/tag)
  + we want to emphasize the pattern in terms of signal, not scale

Correlation distance:

+ centers and scales each protein
+ removes/accounts for magnitude effects, focusing on the shape of the profile
  + is therefore robust to contaminants and bait biases

This is standard in proteomics clustering

```{r, results='asis', echo=FALSE}
fmt <- knitr::pandoc_to()

if (identical(fmt, "html")) {
  cat('# Construct Similarity Matrices {.tabset .tabset-pills .tabset-fade}\n')
} else if (identical(fmt, "markdown")) {
  cat('<details><summary><strong>Construct Similarity Matrices</strong></summary>\n')
} else {
  cat('# Construct Similarity Matrices\n')
}
```

Build a similarity matrix from the correlation matrix of raw counts

  -> will be used for clustering (closeness between data, e.g. proteins)
  
  Two options:
  
  1. **Gaussian similarity matrix** (using the Radial Basis Function kernel)
    + ideal for smooth, low-noise datasets
    + smooth, continuous calculations of distances without hard cutoffs
    
  2. **k-Nearest Neighbors adjacency graph**
    + ideal for sparse or data with heteregeneous density; noisy data (e.g. 
    IPMS datasets with lots of contaminants)
    
  We will compare and select the method that yields the most stable Laplacian
  spectrum
  
## Gaussian (Radial Basis Function) Kernel Similarity Matrix

*S_ij* = exp(-(*D_ij*^2)/2*sigma*^2)
  
```{r Sim Matrix Construction w Gaussian kernel aka radial basis kernel}
# choose a kernel width (sigma) based on the distance distribution
sigma <- stats::median(IPMS_dist)

# compute the Gaussian (RBF) similarity matrix
sim_mat <- exp(-(IPMS_dist^2)/(2*sigma^2))
```
  
## kNN Adjacency

```{r Sim Matrix Construction w kNN}
# compute k-nearest neighbors from the correlation-distance matrix
kNN_idx <- FNN::get.knn(IPMS_dist, k = 10)$nn.index

# convert to binary adjacency matrix
adjacency <- matrix(0, nrow = nrow(kNN_idx), ncol = nrow(kNN_idx))

for (i in seq_len(nrow(kNN_idx))) {
  adjacency[i, kNN_idx[i, ]] <- 1
}

# make symmetrical
adjacency_sym <- (adjacency + t(adjacency)) > 0
```

```{r, results='asis', echo=FALSE}
fmt <- knitr::pandoc_to()

if (identical(fmt, "html")) {
  cat('# Compute the Normalized Laplacian {.tabset .tabset-pills .tabset-fade}\n')
} else if (identical(fmt, "markdown")) {
  cat('<details><summary><strong>Compute the Normalized Laplacian</strong></summary>\n')
} else {
  cat('# Compute the Normalized Laplacian\n')
}
```

Compute the graph `Laplacian` of the similarity/adjacency matrix

+ this is the formal encoding of the graph structure based on distribution of 
space (similarity distances between nodes)

+ **The Laplacian**, in general, serves as a transformation of a matrix, 
replacing diagonal values with the row-sum of the matrix

  + the transformation replaces the off-diagonal values of the matrix with the 
  corresponding matrix values multiplied by -1 ( (off-diagonal matrix value) * (-1) ). 
  
  -> in other words, this captures how each protein's similarity value differs 
  from a weighted combination of all other proteins, where the weights are the 
  similarity matrix entries, themselves:
  
    - **L = D - W**

  + where:
  
    + *W* is either the **Gaussian similarity matrix** (*S*, `sim_mat`) or **kNN adjacency matrix**
    (*A*, `adjacency_sym`)
    + *D* is the degree matrix (diagonal, i.e. row-sum of *W*, or the total
    connection strength of a node):
    
      + D[*i*,*i*] = Sum_j*W*[*i*,*j*]
      + `Gaussian`-based degrees vary smoothly; `kNN`-based degrees are roughly
      constant because each row has the same number of non-zero entries
    
+ For spectral clustering, the standard approach is to use the **normalized Laplacian**,
which re-scales the contribution of each row/column by its degree to account for
uneven similarity magnitudes across proteins:

  + *L_norm* = *D*^(-1/2)(*D*-*W*)*D*^(-1/2)
  or *L_norm* = *I*-*D*^(-1/2)*WD*^(-1/2)
  
+ **The Laplacian** is the operator whose eigenvectors reveal:

  + connected components
  + cluster structure
  
+ **The normalized Laplacian** is essentially a function that transforms a
similarity matrix ("graph" or space) in the same way as the **unnormalized Laplacian**
but with an additional correction step

+ The key difference is that the **normalized Laplacian**
*penalizes proteins with large total similarity (large row-sums)*: the 
normalization re-scales each row and column by the square root of its degree

  + This makes proteins with very different total similarity strengths contribute
  on a comparable on the same scale, preventing high-degree proteins from 
  dominating transformation

## Normalized Laplacian for Gaussian Similarity Matrix

Derive the normalized Laplacian from the Gaussian kernel on the correlation
matrix

```{r Construct the Gaussian Normalized Laplacian}
# for the Gaussian RBF kernel similarity matrix
W <- sim_mat

# degree matrix
D <- diag(rowSums(W))

# normalized Laplacian
D_inv_sqrt <- diag(1 / sqrt(rowSums(W)))
gaussian_norm_L <- diag(nrow(W)) - D_inv_sqrt %*% W %*% D_inv_sqrt
```

## Normalized Laplacian for kNN Adjacency

Derive the normalized Laplacian from the k-Nearest Neighbors adjacency matrix

```{r Construct the kNN Normalized Laplacian}
## for the kNN adjacency matrix
W <- adjacency_sym

# degree matrix
D <- diag(rowSums(W))

# normalized Laplacian
D_inv_sqrt <- diag(1 / sqrt(rowSums(W)))
kNN_norm_L <- diag(nrow(W)) - D_inv_sqrt %*% W %*% D_inv_sqrt
```

```{r, results='asis', echo=FALSE}
fmt <- knitr::pandoc_to()

if (identical(fmt, "html")) {
  cat('# Normalized Laplacian Spectra {.tabset .tabset-pills .tabset-fade}\n')
} else if (identical(fmt, "markdown")) {
  cat('<details><summary><strong>Normalized Laplacian Spectra</strong></summary>\n')
} else {
  cat('# Normalized Laplacian Spectra\n')
}
```

Conduct **eigendecomposition** and extract **eigenvectors** for each similarity matrix 
(graph). Compare and select for **embeddings** upon which to cluster

  + **embeddings** are the numerical representations of each/every protein in 
  lower-dimensional space
  
    + an **embedding** is constructed from the **eigenvectors** of the 
    `graph Laplacian`
    
    + an **eigenvector** encodes how one protein relates to every other protein
    in the similarity graph
  
  + this transformation places proteins in a visualizable, geometric coordinate 
  system where proteins with similar interaction profiles are closer together 
  and dissimilar proteins are further apart
  
  + an embedding (represented by `k`) with multiple eigenvectors has multiple 
  dimensions-- multiple clusters-- by which the proteins relate: 
  **the number of eigenvectors = k = clusters within the graph Laplacian**

Extract the `eigenvectors` and `eigenvalues` of each similarity graph method

## Eigendecomposition of Gaussian Normalized Laplacian

```{r Extract Gaussian Eigencomponents, echo = F}
gaussian_eigenvalues <- gaussian_eig$values
gaussian_eigenvectors <- gaussian_eig$vectors
```

```{r Conduct Eigendecomposition on Gaussian RBF kernel similarity matrix, eval = F}
# compute eigendecomposition of the Gaussian RBF kernel-based normalized Laplacian
gaussian_eig <- eigen(gaussian_norm_L, symmetric = TRUE)

# extract eigenvalues and eigenvectors
gaussian_eigenvalues <- gaussian_eig$values
gaussian_eigenvectors <- gaussian_eig$vectors
```

## Eigendecomposition of kNN Normalized Laplacian

```{r Extract kNN Eigencomponents, echo = F}
kNN_eigenvalues <- kNN_eig$values
kNN_eigenvectors <- kNN_eig$vectors
```

```{r Conduct Eigendecomposition on kNN adjacency matrix, eval = F}
# compute eigendecomposition of the kNN adjacency matrix-based normalized Laplacian
kNN_eig <- eigen(kNN_norm_L, symmetric = TRUE)

# extract eigenvalues and eigenvectors
kNN_eigenvalues <- kNN_eig$values
kNN_eigenvectors <- kNN_eig$vectors
```

```{r, results='asis', echo=FALSE}
fmt <- knitr::pandoc_to()

if (identical(fmt, "html")) {
  cat('# Automatic Similarity Graph Choice {.tabset .tabset-pills .tabset-fade}\n')
} else if (identical(fmt, "markdown")) {
  cat('<details><summary><strong>Automatic Similarity Graph Choice</strong></summary>\n')
} else {
  cat('# Automatic Similarity Graph Choice\n')
}
```

Compute the **eigengap**-- the largest difference between 
consecutive `Laplacian` **eigenvalues** in ascending order-- for each graph

+ The index (**eigengap**) represents the most stable number of clusters for the 
graph

```{r Choose the Graph that Optimizes the Laplacian}
# compute eigengaps and identify the largest eigen gap for each graph
  # sort the eigenvalues
sorted_gaussian_eigenvalues <- sort(gaussian_eigenvalues)
  # find the biggest eigengap
gaussian_best_gap <- max(diff(sorted_gaussian_eigenvalues))
  # find k where this occurs
gaussian_k <- which(diff(sorted_gaussian_eigenvalues) == gaussian_best_gap)

  # sort the eigenvalues
sorted_kNN_eigenvalues <- sort(kNN_eigenvalues)
  # find the biggest eigengap
kNN_best_gap <- max(diff(sorted_kNN_eigenvalues))
  # find k where this occurs
kNN_k <- which(diff(sorted_kNN_eigenvalues) == kNN_best_gap)

# choose the graph with the larger eigengap
if (gaussian_best_gap > kNN_best_gap) {
  chosen_graph <- "gaussian" 
  } else {
  chosen_graph <- "kNN"
}

# store chosen Laplacian and eigenvectors
if (chosen_graph == "gaussian") {
  L_chosen <- gaussian_norm_L
  eigvals_chosen <- gaussian_eigenvalues
  eigvecs_chosen <- gaussian_eigenvectors
  W <- sim_mat
  D <- diag(rowSums(W))
  D_inv_sqrt <- diag(1 / sqrt(rowSums(W)))
  k <- gaussian_k
} else {
  L_chosen <- kNN_norm_L
  eigvals_chosen <- kNN_eigenvalues
  eigvecs_chosen <- kNN_eigenvectors
  k <- kNN_k
}

chosen_graph
```

```{r, results='asis', echo=FALSE}
fmt <- knitr::pandoc_to()

if (identical(fmt, "html")) {
  cat('# Spectral Embedding {.tabset .tabset-pills .tabset-fade}\n')
} else if (identical(fmt, "markdown")) {
  cat('<details><summary><strong>Spectral Embedding</strong></summary>\n')
} else {
  cat('# Spectral Embedding\n')
}
```
## Plot Embedding Dimensions (k)

Show the `eigengap` among `eigenvalues` for each graph (computed above)

```{r SortedEigenvalues, fig.height=8,fig.width=8}
sorted_kNN_eigenvalues |> 
  as.data.frame() |> 
  dplyr::mutate(sorted_gaussian_eigenvalues = 
    sorted_gaussian_eigenvalues
  ) |> 
  dplyr::mutate(k = dplyr::row_number(), .before = 1) |> 
  tidyr::pivot_longer(cols = c(dplyr::contains("eig")), names_to = 'graph', values_to = 'value') |> 
  dplyr::group_by(graph) |>
  ggplot2::ggplot(aes(x = k, y = value)) +
  ggplot2::geom_point(aes(group = graph, color = graph)) +
  ggplot2::labs(title = 'Sorted Eigenvalues\nby Similarity Graph Method',
                y = 'Eigenvalue',
                color = 'Similarity\nGraph\nMethod') +
  ggplot2::theme_bw()+
  ggplot2::theme(title = ggplot2::element_text(face = 'bold', size = 12),
                 axis.text = ggplot2::element_text(face = 'bold', size = 9),
                 axis.title = ggplot2::element_text(face = 'bold', size = 11),
                 legend.title = ggplot2::element_text(face = 'bold', size = 11),
                 legend.text = ggplot2::element_text(face = 'bold', size = 9),
                 strip.background = ggplot2::element_rect(color = 'black',
                                                          fill = 'white'),
                 strip.text = ggplot2::element_text(face = 'bold', size = 9)) +
  ggplot2::facet_wrap(~graph)

```

These spectra are not unexpected-- `kNN` is very smooth, and the `Gaussian` is 
very flat, both suggesting there aren't many clear structural boundaries in the
data, which isn't uncommon in IPMS interactomes due to the substantial overlap 
in protein networks (pathways, complexes, and differential functions)

```{r Remove waste, echo = FALSE}
if (chosen_graph == "gaussian") {
  rm(kNN_norm_L, kNN_eig, kNN_eigenvalues, kNN_eigenvectors,
     adjacency, adjacency_sim, kNN_graph, kNN_idx, kNN_k, sorted_kNN_eigenvalues)
} else {
  rm(sim_mat, gaussian_norm_L, gaussian_eig, gaussian_eigenvectors, 
     gaussian_eigenvalues, gaussian_best_gap, sigma, gaussian_k,
     sorted_gaussian_eigenvalues)
}
```

## Compute the Spectral Embedding

```{r Compute the Embedding by Subsetting the Eigenvectors through the Eigengap}
# spectral embedding
embedding <- eigvecs_chosen[, k, drop = FALSE]
```

## Results and Interpretation

We have determined `k` = **1** because the `Gaussian similarity matrix` showed the 
larger **eigengap** (`r round(gaussian_best_gap, 2) > round(kNN_best_gap, 2) `), 
and this occurred at `k` = **`r gaussian_k`** in the sorted spectrum

As shown in the [eigenvalue plot above](#plot-embedding-dimensions-(k)), the largest slope 
change (largest **eigengap**) occurs between the first two eigenvalues (i.e.
`k` = **`r gaussian_k`**): 

  + there is little structural difference within the data-- very smooth for both
  methods, which isn't uncommon in IPMS interactomes (see [section above](#plot-embedding-dimensions-(k))
  for rationale)

This indicates that the `Laplacian spectrum` contains **no evidence of more than one meaningful cluster**, 
and under the eigengap rule, this yields a *1-dimensional embedding*

Because `k` = **1**, the `spectral embedding` reduces to the first **eigenvector** of the
chosen graph `Laplacian` (`r stringr::str_to_sentence(chosen_graph)`). For 
completeness and downstream coherence analysis, we retain the full **eigenvector** 
matrix, but only the first column carries spectral meaning

```{r, results='asis', echo=FALSE}
fmt <- knitr::pandoc_to()

if (identical(fmt, "html")) {
  cat('# k-means Clustering {.tabset .tabset-pills .tabset-fade}\n')
} else if (identical(fmt, "markdown")) {
  cat('<details><summary><strong>k-means Clustering</strong></summary>\n')
} else {
  cat('# k-means Clustering\n')
}
```

Conduct `k-means clustering` on the `eigenvectors`/`embedding`

## Extract Clusters

(1 for all proteins)

```{r k means on embeddings}
set.seed(4)
labels <- stats::kmeans(embedding, centers = k)$cluster

clusters_df <- labels |> 
  as.data.frame() |> 
  dplyr::mutate(Protein = colnames(IPMS_dist), .before = 1) |> 
  dplyr::rename(cluster = 2)

head(clusters_df)
```

```{r, results='asis', echo=FALSE}
fmt <- knitr::pandoc_to()

if (identical(fmt, "html")) {
  cat('# Cluster Descriptions {.tabset .tabset-pills .tabset-fade}\n')
} else if (identical(fmt, "markdown")) {
  cat('<details><summary><strong>Cluster Descriptions</strong></summary>\n')
} else {
  cat('# Cluster Descriptions\n')
}
```

## Join Clusters to hSAINT

Include score summaries (mean and median `posterior probability scores`, `ranks`, and literature or CRAPome `contamination rates`)

```{r Summarize SAINT Scores and Contamination Rates}
spectral_cluster_scores <- labels |> 
  as.data.frame() |> 
  dplyr::mutate(Protein = colnames(IPMS_dist), .before = 1) |> 
  dplyr::select(Protein) |> 
  dplyr::mutate(cluster = labels) |> 
  dplyr::inner_join(
    hSAINT
  ) |> 
  dplyr::arrange(desc(score)) |> 
  dplyr::group_by(cluster) |> 
  dplyr::summarize(mean_score = mean(score),
                   median_score = stats::median(score),
                   score_mass = sum(score),
                   n = dplyr::n(),
                   lit_flag_perc = (sum(lit_flag)/n)*100,
                   CRAPome_flag_perc = (sum(CRAPome_flag)/n)*100) |> 
  dplyr::arrange(desc(mean_score)) |> 
  dplyr::mutate(mean_rank = dplyr::row_number()) |> 
  dplyr::arrange(desc(median_score)) |> 
  dplyr::mutate(median_rank = dplyr::row_number()) |> 
  dplyr::arrange(desc(score_mass)) |> 
  dplyr::mutate(total_rank = dplyr::row_number()) |> 
  dplyr::arrange(mean_rank)

spectral_cluster_scores
```

Join the cluster scores to the cluster (per protein) dataframe and the `hSAINT` results for classification (for posterity)

```{r Join the Summary dataframe to the Protein Level dataframes}
spectral_protein_cluster_df <- clusters_df |> 
  dplyr::inner_join(
    hSAINT
  ) |> 
  dplyr::full_join(
    spectral_cluster_scores
  ) |> 
  dplyr::arrange(desc(score)) |> 
  dplyr::mutate(
    SAINT_classification = dplyr::case_when(
      
      mean_rank <= 5 ~ 'High SAINT',
      mean_rank > 5 ~ 'Low SAINT'
      ),
    
    contams_classification = dplyr::case_when(
      
      lit_flag_perc > 5 ~'High Contaminants',
      TRUE ~ 'Low Contaminants')) |> 
  
  dplyr::mutate(classification = paste0(SAINT_classification,
                                        ", ",
                                        contams_classification)) |> 
  dplyr::select(-SAINT_classification,
                -contams_classification) |> 
  dplyr::relocate(classification, .after = cluster)

head(spectral_protein_cluster_df)
```

## Cluster GO GSEA

```{r Downstream Analysis}
EntrezIDs <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = colnames(IPMS_dist),
  columns = 'ENTREZID',
  keytype = 'SYMBOL'
) |> 
  dplyr::filter(!is.na(2))

EntrezMap <- stats::setNames(EntrezIDs$ENTREZID, EntrezIDs$SYMBOL)
```

Create enrichment functions

```{r Create Enrichment Functions}
background_genes <- EntrezIDs |> # use the EntrezIDs mapped from the modules_df
  dplyr::pull(ENTREZID)

get_GO <- function(genes, ontology, background_genes) {
  clusterProfiler::enrichGO(
    gene = genes,
    universe = background_genes,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = ontology,
    pAdjustMethod = 'BH',
    qvalueCutoff = 0.05,
    readable = T
  )
}

get_kegg <- function(genes, EntrezMap) {
  entrez <- EntrezMap[genes]
  entrez <- entrez[!is.na(entrez)]
  if (length(entrez)==0) {
    return(NULL)
  }
  clusterProfiler::enrichKEGG(
    gene = entrez,
    organism = 'hsa',
    pAdjustMethod = 'BH',
    qvalueCutoff = 0.05
  )
}

get_reactome <- function(genes, EntrezMap) {
  entrez <- EntrezMap[genes]
  entrez <- entrez[!is.na(entrez)]
  if (length(entrez)==0) {
    return(NULL)
  }
  ReactomePA::enrichPathway(
    gene = entrez,
    organism = 'human',
    pAdjustMethod = 'BH',
    qvalueCutoff = 0.05,
    readable = T
  )
}
```

Run the functions on the dataset in parallel (multiple CPUs)

```{r Parallelize Enrichment, eval = FALSE}
future::plan(future::multisession, workers = 6)
progressr::handlers(progressr::handler_txtprogressbar())

cluster_ids <- unique(labels)

progressr::with_progress({
  p <- progressr::progressor(steps = length(cluster_ids))

  enrichment_results <- future.apply::future_lapply(
    cluster_ids,
    function(mod) {

      # tick the progress bar
      p(sprintf("cluster %s", mod))

      genes <- labels |>
        as.data.frame() |> 
        dplyr::mutate(Protein = colnames(IPMS_dist), .before = 1) |> 
        dplyr::rename(cluster = labels) |> 
        dplyr::filter(cluster == mod) |>
        dplyr::inner_join(EntrezIDs |> dplyr::rename(Protein = 1)) |> 
        dplyr::distinct(ENTREZID) |>
        dplyr::pull(ENTREZID)

      list(
        cluster = mod,
        go_bp  = get_GO(genes, ontology = "BP", background_genes = background_genes),
        go_cc  = get_GO(genes, ontology = "CC", background_genes = background_genes),
        go_mf  = get_GO(genes, ontology = "MF", background_genes = background_genes)
        # add KEGG/Reactome later
      )
    }
  )
})

names(enrichment_results) <- cluster_ids
rm(p)
```

Gather enrichment results

```{r Consolidate Enrichment Results, eval = FALSE}
spectral_enrichment_results_df <- enrichment_results |> 
  purrr::pluck(1) |> 
  purrr::list_merge() %>%
  purrr::map_df(., .f = as.data.frame) |> 
  dplyr::rename(cluster = 1, GO_Term = Description, GO_ID = ID)

spectral_enrichment_results_df
```

```{r Load Module GO GSEA, eval = TRUE, echo = FALSE}
load("~/R/SpectralClustering/spectralCluster_GO_GSEA.rda")
```

**No GO Terms were enriched** in the only cluster of the spectral clustering method. 
As in [Endo/Rpkg/vignette/HierarchicalClustering.Rmd](/vignette/HierarchicalClustering.Rmd),
we would run GO enrichment on the clusters that did not return enrichment results,
using only each cluster's proteins as the background for enrichment

+ In this case, that remains the same background protein list, since it is the 
only cluster in the data

+ **This yields the same result**

## Cluster Barplot

For posterity, since there are no enriched GO Terms from the only cluster in the data

```{r ClusterBarPlot}
spectral_enrichment_results_df |> 
  dplyr::full_join(spectral_cluster_scores) |> 
  dplyr::filter(!is.na(GO_Term)) |> 
  dplyr::mutate(
    is_enriched = 'Not Enriched'
    ) |> 
  dplyr::group_by(cluster) |> 
  dplyr::arrange(p.adjust) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = -log10(p.adjust),
      y = forcats::fct_rev(forcats::fct_infreq(GO_Term, mean_score)),
      size = FoldEnrichment,
      alpha = mean_score,
      shape = is_enriched,
    )
  ) +
  ggplot2::geom_vline(xintercept = -log10(0.05),
                      linetype = 'dashed',
                      color = 'firebrick',
                      size = 1,
                      alpha = 0.5) +
  ggplot2::geom_point() +
  ggplot2::labs(title = 'GO Enrichment Across Clusters\nDerived from Eigendecomposition and k-means (Spectral Clustering)\nin TMEM184B IPMS',
                x = '-log10(Adj.P-val)',
                y = 'GO Term',
                alpha = 'Cluster\nMean\nHierarchical\nSAINT\nScore',
                size = 'FoldEnrichment',
                shape = 'IsEnrichedCluster') +
  ggplot2::theme_bw() +
  ggplot2::theme(title = ggplot2::element_text(face = 'bold', size = 12),
                 axis.text = ggplot2::element_text(face = 'bold', size = 10),
                 strip.text = ggplot2::element_text(face = 'bold', size = 7), 
                 strip.background = ggplot2::element_rect(color = 'black',
                                                          fill = 'white')) +
  ggplot2::guides(#shape = ggplot2::guide_legend(override.aes = list(size = 4, shape = c(19, 17))),
                  color = ggplot2::guide_legend(override.aes = list(size = 4,
                                                                    shape = 19,
                                                                    color = c('navy',
                                                                              'darkgoldenrod3',
                                                                              'firebrick'))),
                  alpha = ggplot2::guide_legend(override.aes = list(size = 5,
                                                                    shape = 19)),
                  size = ggplot2::guide_legend(override.aes = list(shape = 19))
                  # size = 'none'
                  ) +
  ggplot2::scale_fill_manual(values = c('navy',
                                        'darkgoldenrod3',
                                        'firebrick'))
```

## Downsample Heatmap of Similarity Matrix

Subset the heatmap of the similarity matrix (for rendering/interpretability)

```{r SimMatHeatmap}
idx <- sample(seq_len(nrow(sim_mat)), 1000)
pheatmap::pheatmap(
  sim_mat[idx, idx],
  color = grDevices::colorRampPalette(c('navy',
                                        'white',
                                        'darkgoldenrod'))((200)),
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = F,
  show_colnames = F,
  main = 'Gaussian RBF Kernel Similarity Matrix\nTMEM184B IPMS TSCs\n(subset for interpretability and render speed)'
)
```

The `Gaussian similarity matrix` is very smooth, with broad regions of near-1 
similarity and several thin, interleaving low-similarity (blue) stripes.

The diagonal is not visually distinct because many proteins have near-identical
similarity profiles, producing relatively large high-similarity (gold) regions.

This pattern is characteristic of a 1-dimensional spectral structure with no 
modularity and is not inconsistent with interactome profiles (see 
 [Compute the Spectral Embedding subsection](#compute-the-spectral-embedding))

## PCA on Embedding

Conduct `PCA` on the spectral embedding

+ since `k` = **1**, *there will only be one component* (**PC1**)

+ we, therefore, plot the distribution of this component across the protein index 
of the embedding

```{r PCA on SpectralEmbedding}
pca_embedding <- stats::prcomp(embedding, center = TRUE, scale = FALSE)

pc1 <- pca_embedding$x[, 1]

df_pc1 <- tibble::tibble(
  index = seq_along(pc1),
  PC1 = pc1
)
```

```{r PCA results from SpectralEmbedding}
pca_embedding$x <- pca_embedding$x[1:5, 1, drop = FALSE]

head(pca_embedding)
```


```{r SpectralEmbeddingPCA}
cowplot::plot_grid(
  df_pc1 |> 
  ggplot2::ggplot(ggplot2::aes(x = index, y = PC1)) +
  ggplot2::geom_point(size = 0.6, color = 'navy') +
  ggplot2::labs(title = 'Distribution of PC1 from\nPCA of the Gaussian Spectral Embedding\n(k = 1; PC1',
                x = 'Protein Index',
                y = 'PC1') +
  ggplot2::theme_bw() +
  ggplot2::theme(title = ggplot2::element_text(face = 'bold', size = 12),
                 axis.title = ggplot2::element_text(face = 'bold', size = 10),
                 axis.text = ggplot2::element_text(face = 'bold', size = 9)),
  
  df_pc1 |> 
  ggplot2::ggplot(ggplot2::aes(x = PC1)) +
  ggplot2::geom_density(fill = 'navy', alpha = 0.5)+
  ggplot2::labs(title = 'Density of Protein PC1 values\nfrom PCA of the Gaussian Spectral Embedding\n(k = 1; PC1)',
                x = 'PC1',
                y = 'density') +
  ggplot2::theme_bw() +
  ggplot2::theme(title = ggplot2::element_text(face = 'bold', size = 12),
                 axis.title = ggplot2::element_text(face = 'bold', size = 10),
                 axis.text = ggplot2::element_text(face = 'bold', size = 9))
  , 
  ncol = 2,
  align = 'h'
)
```


## Construct the Cluster GO (similarity) Network

For posterity, *since there are no enriched terms*, **no GO network can be constructed**

```{r Construct GO Network Object}
spectral_enrichment_results_df |> 
  dplyr::full_join(spectral_cluster_scores) |> 
  dplyr::mutate(node_to = GO_Term,
                node_from = cluster,
                weight = -log10(p.adjust)) |> 
  dplyr::filter(weight > -log10(0.05))
```

```{r, results='asis', echo=FALSE}
fmt <- knitr::pandoc_to()

if (identical(fmt, "html")) {
  cat('# Results Tables {.tabset .tabset-pills .tabset-fade}\n')
} else if (identical(fmt, "markdown")) {
  cat('<details><summary><strong>Results Tables</strong></summary>\n')
} else {
  cat('# Results Tables\n')
}
```

To conclude, in this dataset, spectral clustering on the `Gaussian kernel graph` 
produced *a single coherent cluster*. The **eigengap** analysis supported this 
outcome and the 1-dimensional spectral embedding showed **no evidence of separable substructure**

Because all proteins fell into this one cluster, GO enrichment analysis returned 
**no significant terms**-- there is no contrast set against which to detect 
over-representation. Some datasets do not contain multiple functional clusters 
at the resolution provided by their similarity structure

The spectral workflow still provides a reproducible, graph-based view of the 
data, even when the optimal partition is a single group
