IPMS vignette
================
Erik Larsen
2025-01-22

- [Overview](#overview)
- [Set up Environment](#set-up-environment)
- [Wrangle the TSC Data](#wrangle-the-tsc-data)
- [Filter Wrangled IPMS Counts Data](#filter-wrangled-ipms-counts-data)
  - [Apply filtering to TSCs](#apply-filtering-to-tscs)
- [Upload GO Results](#upload-go-results)
  - [Cell Component Results](#cell-component-results)
  - [Biological Process Results](#biological-process-results)
  - [Molecular Function Results](#molecular-function-results)
  - [Combine the Results](#combine-the-results)
- [Integrate with Custom GO Query
  Function](#integrate-with-custom-go-query-function)
  - [Query Bioconductor and GO dbs](#query-bioconductor-and-go-dbs)
  - [Integrate Query Results into Counts
    DF](#integrate-query-results-into-counts-df)
- [Plot Results](#plot-results)
  - [Plot Full Scatter](#plot-full-scatter)
  - [vATPase Scatter (zoom)](#vatpase-scatter-zoom)
  - [Regulation of Macroautophagy
    (zoom)](#regulation-of-macroautophagy-zoom)
  - [Protein Transport (zoom)](#protein-transport-zoom)
  - [GO Bar Plot](#go-bar-plot)
  - [GO CC Bar Plot](#go-cc-bar-plot)
  - [GO BP Bar Plot](#go-bp-bar-plot)

## Overview

This markdown was developed to formally document the data analysis
pipeline for TMEM184B immunoprecipitation tandem mass-spec experiments
for **Dr. Martha Bhattacharya’s lab** at the **University of Arizona**
for the [endolysosomal paper](<https://journals.biologists.com/jcs/article/doi/10.1242/jcs.263908/368440/Transmembrane-protein-184B-TMEM184B-modulates>)

Until a container is built, the `R version` used was **R 4.3.2** with
the packages installed below

## Set up Environment

Attach packages

``` r
packages <- c('Endo', 'ggplot2', 'rstatix', 'forcats', 'knitr', 'TMEM')

for (package in packages) {
  if (!require(package, character.only = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}
```

Load the IPMS TSC “raw” data (label-free spectral counts file of 8 samples)

``` r
data("IPMS_counts")
ALL_DF <- IPMS_counts
rm(IPMS_counts)
```

## Wrangle the TSC Data

Break the original wrangling call into commented chunks for clarity

``` r
ALL_DF <- ALL_DF |> 
  dplyr::select(-1) |> ## remove empty column
  dplyr::rename_with(
    
    ~ALL_DF |>
      dplyr::select(-1) |>
      dplyr::slice(1) |>
      unlist() |>
      as.character()
    ) |> # rename the columns using the first row
  
  dplyr::slice(-1) |> # remove the first row
  dplyr::slice_tail(n = -1) |> # remove the last row
  dplyr::select(-dplyr::contains("307-411")) |> # remove C terminus samples
  
  dplyr::rename(
    Proteins.Identified = dplyr::contains("Proteins")
    ) |>  # rename the protein ID column
  
  dplyr::rename(
    Accession.Number = dplyr::contains("Accession")
    ) |> # rename the protein accession number column
  
  dplyr::rename(
    MW = dplyr::contains("Molecular")
    ) |> # rename the molecular weight column
  
  dplyr::filter(
    !grepl(Accession.Number, pattern = "DECOY")
    ) |> # filter decoy proteins out
  
  dplyr::rename_with(
    ~gsub(.x, pattern = '[0-9]+_', replacement = '')
    ) |> # remove _'s from column names
  
  dplyr::rename_with(
    ~gsub(.x, pattern = '^m', replacement = '')
    ) |> # remove m's from column names (mTMEM)
  
  dplyr::rename(BirAV5_1 = `BirA-V1`, 
                BirAV5_2 = `BirA-V2`,
                TMEMV5_1 = `TMEM-V1`,
                TMEMV5_2 = `TMEM-V2`
                ) # rename the V5 columns by replacing hyphens with _'s
```

``` r
ALL_DF <- ALL_DF |>
  dplyr::mutate(
    HumanGeneID = stringr::str_extract(
      Proteins.Identified,
      pattern = 'GN=[:alnum:]+') |>
      gsub(Proteins.Identified, pattern = 'GN=', replacement = ''),
    .after = Proteins.Identified
    ) |> # create a column containing proteins' HUGO symbols
  
  dplyr::mutate(
    Proteins.Identified = gsub(Proteins.Identified,
                               pattern = ' O.+$',
                               replacement = '')
    ) |> # remove extra strings
  
  dplyr::mutate(
    MW = gsub(MW, pattern = ' kDa', replacement = '') |> as.numeric()
                ) |> # remove kDa's from MW column to enable math operations

  dplyr::mutate(
    dplyr::across(dplyr::contains("_"), ~dplyr::if_else(.x == '?', '0', .x)),
    dplyr::across(dplyr::contains("_"), ~dplyr::if_else(.x == '', '0', .x)),
    dplyr::across(dplyr::contains("_"), ~as.numeric(.x)
                  ), # replace ?'s, empty columns with 0's; convert to numeric
    dplyr::across(dplyr::contains("_"), ~dplyr::if_else(.x == 0, 0.1, .x)
                  ), # add 0.1 to 0's for FC calculations
    dplyr::across(dplyr::contains("_"), .fns = list("MW_norm" = ~.x/MW)
                  ) # add new columns that normalize TSC's by MW
    )
```

``` r
ALL_DF <- ALL_DF |> 
  dplyr::mutate(
    
    V5ControlAvg = (BirAV5_1_MW_norm + BirAV5_2_MW_norm)/2,
    mycControlAvg = (GFPmyc_1_MW_norm + GFPmyc_2_MW_norm)/2,
    V5TMEMAvg = (TMEMV5_1_MW_norm + TMEMV5_2_MW_norm)/2,
    mycTMEMAvg = (TMEMmyc_1_MW_norm + TMEMmyc_2_MW_norm)/2, # tag & bait avgs
    
    avgV5_FC = V5TMEMAvg/V5ControlAvg,
    avgmyc_FC = mycTMEMAvg/mycControlAvg, # compute the avg FC for each tag
    
    ControlAvg = (BirAV5_1_MW_norm+
                    BirAV5_2_MW_norm+
                    GFPmyc_1_MW_norm+
                    GFPmyc_2_MW_norm
                  )/4, # average MW-normalized average for control bait samples
    
    TMEMAvg = (TMEMV5_1_MW_norm+
                 TMEMV5_2_MW_norm+
                 TMEMmyc_1_MW_norm+
                 TMEMmyc_2_MW_norm
               )/4, # average MW-normalized average for target bait samples
    
    avgFC = TMEMAvg/ControlAvg # compute the grand FC
    
    ) |>
  dplyr::mutate(
    
    MouseGeneID = paste(substr(HumanGeneID, 1, 1),
                        tolower(substr(HumanGeneID, 2, 10)), sep = ''),
    
    MouseGeneID = dplyr::case_when(
      grepl(Accession.Number, pattern = 'p53_HUMAN') ~ "Trp53",
      HumanGeneID == 'birA' ~ HumanGeneID,
      Proteins.Identified == 'Green fluorescent protein' ~ HumanGeneID,
      TRUE ~ MouseGeneID)
    ) |> # create mouse EntrezID column
  
  dplyr::relocate(MouseGeneID, .after = HumanGeneID) |>
  dplyr::arrange(desc(avgFC)) # sort by grand FC from highest to lowest)
```

## Filter Wrangled IPMS Counts Data

Filter using aggressive heuristics; supported by Dr. Paul Langlais (University of
Arizona College of Molecular Medicine Proteomics; personal
communication)

``` r
Compiled_Candidate_List <- ALL_DF |>
  dplyr:::arrange(desc(avgFC)) |>
  dplyr::filter(
    avgFC > 2,
    avgV5_FC > 2,
    avgmyc_FC > 2, 
    # > one-fold increase of MW-normalized count avg for each tag, and across all
    
    (TMEMV5_1 & TMEMV5_2 > 5),
    
    (TMEMmyc_1 & TMEMmyc_2 > 5), # > 5 TSCs across duplicates for each tag AND
    
    (
      TMEMV5_1 > 10 |
        TMEMV5_2 > 10 |
        TMEMmyc_1 > 10 |
        TMEMmyc_2 > 10
      ) # at least 10 TSCs in at least one sample
    
    ) |>
  
  # Remove Mulcahy contams (2008 paper)
  dplyr::filter(
    
    !grepl(
      HumanGeneID,
      pattern = 'AHNAK|ALB|ANXA(1|2)|ASS1|ATP5|CFL1|CLTC|CSDA|EDARRAD|ENO|FARS(A|B)|FASN|FKSG30|GAPDH|GFAP|ILF(2|3)|LDHA|LGALS(1|3)|PCB(1|2)|PDIA6|PHB|PHB2|POTE2|PPIA|PRDX(1|2|3|4)|SERPIN(H1|A11)|SLC25A|TUFM|TXN|UBA52')
    
    ) |>
  
  # Remove Dunham contams (2012 review)
  dplyr::filter(
    
    !grepl(
      HumanGeneID,
      pattern = 'EEF1|EIF|HNRN|HSPA|KRT|TUB(A|B)|RPL|HIST2|H2A|RRP|(X|I)PO|NUP|DNAJ')
    
    ) |>
  
  dplyr::filter(HumanGeneID != "TMEM184B") |>
  dplyr::select(HumanGeneID) |>
  unlist() |> as.character()
```

Import CRAPome results (load from the `Endo` package)

``` r
data("CRAPome_results")
```

Wrangle CRAPome results

``` r
Crapome_results <- CRAPome_results |>
  dplyr::filter(!grepl(Mapped, pattern = '/|identifier')) |>
  dplyr::select(1:7) |>
  dplyr::rename(User_input = User,
                Mapped_Gene_Symbol = Input,
                Num_Exp_Found = Mapped,
                Num_Exp_Total = Symbol,
                Avg_SC = Num,
                Max_SC = of) |>
  dplyr::mutate(dplyr::across(c(3,5:7), ~as.numeric(.x))) |>
  dplyr::select(-Gene)
rm(CRAPome_results)
```

Identify proteins in **\> 5%** of CRAPome studies and remove them

``` r
Crapome_hits <- Crapome_results |>
  dplyr::filter(Num_Exp_Found <= 35.8) |> # fewer than 5% of experiments
  dplyr::select(Mapped_Gene_Symbol) |>
  unlist() |>
  as.character()
```

### Apply filtering to TSCs

Use CRAPome list to annotate the TSC data

``` r
ALL_DF_filtered <- ALL_DF |>
  dplyr::mutate(poi = 'Background', .after = MouseGeneID) |>
  dplyr::full_join(Crapome_hits |>
                     as.data.frame() |>
                     dplyr::mutate(HumanGeneID = Crapome_hits)
                   ) |> 
  dplyr::mutate(poi = dplyr::case_when(
    
    !is.na(Crapome_hits) ~ 'Candidates',
    HumanGeneID == 'TMEM184B' ~ 'TMEM184B',
    TRUE ~ poi)
    
    ) |>
  dplyr::select(-Crapome_hits)
```

## Upload GO Results

Upload the result `.txt` files from the `geneontology.org` query

### Cell Component Results

Load and wrangle cell components gene ontology results (in `Endo` package)

``` r
data("GO_CC_results")

GO_CC_results <- GO_CC_results |> 
  dplyr::rename_with(
    ~gsub(.x, pattern = '\\.\\.', replacement = '')
    ) |> # remove .'s from GO CC column
  
  dplyr::rename_with(
    ~gsub(.x, pattern = 'Homo.sapiens.|upload_1|\\.$', replacement = '')
    ) |> # remove excess strings
  
  dplyr::rename(
    GO_Term = GO.cellular.component.complete
                ) |> # rename the CC column
  
  dplyr::rename(
    GO_Term_Size = REFLIST20580
                ) |> # rename the GO Term Size column
  
  dplyr::rename(
    Overlap = `138`
                ) |> # rename the overlap column
  dplyr::mutate(
    dplyr::across(c(2:4,7,8), ~as.numeric(.x))
                ) |> # convert to numeric
  
  # dplyr::slice(c(1:925)) |> 
  tidyr::separate(
    GO_Term, into = c("GO_Term", "GO_Term_ID"), sep = '\\s\\(GO:'
                  ) |> # split the Term ID from the Term
  
  dplyr::mutate(
    GO_Term_ID = paste0('GO:',
                        gsub(GO_Term_ID, pattern = '\\)', replacement = ''))
                ) |>
  
  dplyr::mutate(GO_Type = "CC",
                .after = GO_Term_ID
                ) |> # create the category column
  
  dplyr::mutate(
    fold.Enrichment = stringr::str_extract(fold.Enrichment,
                                           pattern = '[:alnum:]+') |> 
      as.numeric(fold.Enrichment)
    ) |> # convert the FE column to numeric
  
  dplyr::arrange(
    FDR, desc(fold.Enrichment)
    ) # sort by the highest FE and lowest adjusted p-value
```

### Biological Process Results

Load and wrangle biological process gene ontology results (in `Endo`
package)

``` r
data("GO_BP_results")

GO_BP_results <- GO_BP_results |> 
  dplyr::rename_with(
    ~gsub(.x, pattern = '\\.\\.', replacement = '')
    ) |> # remove the .'s  from the GO BP column
  
  dplyr::rename_with(
    ~gsub(.x, pattern = 'Homo.sapiens.|upload_1|\\.$', replacement = '')
    ) |> # remove excess strings
  
  dplyr::rename(
    GO_Term = GO.biological.process.complete
                ) |> # rename the BP column
  
  dplyr::rename(
    GO_Term_Size = REFLIST20580
    ) |> # rename the GO Term Size column
  
  dplyr::rename(
    Overlap = `138`
    ) |> # rename the overlap column
  
  dplyr::mutate(
    dplyr::across(c(2:4,7,8), ~as.numeric(.x))
    ) |> # convert to numeric
  
  tidyr::separate(
    GO_Term, into = c("GO_Term", "GO_Term_ID"), sep = '\\s\\(GO:'
    ) |> # split the Term ID from the Term
  
  dplyr::mutate(
    GO_Term_ID = paste0('GO:',
                        gsub(GO_Term_ID, pattern = '\\)', replacement = ''))
    ) |>
  
  dplyr::mutate(
    GO_Type = "BP", .after = GO_Term_ID
    ) |> # create the category column
  
  dplyr::mutate(
    fold.Enrichment = stringr::str_extract(fold.Enrichment,
                                           pattern = '[:alnum:]+') |>
      as.numeric(fold.Enrichment)
    ) |> # convert the FE column to numeric
  
  dplyr::arrange(
    FDR, desc(fold.Enrichment)
    ) # sort by the highest FE and lowest adjusted p-value
```

### Molecular Function Results

Load and wrangle molecular function gene ontology results (in `Endo`
package)

``` r
data("GO_MF_results")

GO_MF_results <- GO_MF_results |> 
  dplyr::rename_with(
    ~gsub(.x, pattern = '\\.\\.', replacement = '')
    ) |> # remove the .'s from the GO BP column
  
  dplyr::rename_with(
    ~gsub(.x, pattern = 'Homo.sapiens.|upload_1|\\.$', replacement = '')
    ) |> # remove excess strings
  
  dplyr::rename(
    GO_Term = GO.molecular.function.complete
    ) |> # rename the MF column
  
  dplyr::rename(
    GO_Term_Size = REFLIST20580
    ) |> # rename the GO Term Size column
  
  dplyr::rename(
    Overlap = `138`
    ) |> # rename the overlap column
  
  dplyr::mutate(
    dplyr::across(c(2:4,7,8), ~as.numeric(.x))
    ) |> # convert to numeric
  
  tidyr::separate(
    GO_Term, into = c("GO_Term", "GO_Term_ID"), sep = '\\s\\(GO:'
    ) |> # split the Term ID from the Term
  
  dplyr::mutate(
    GO_Term_ID = paste0('GO:',
                        gsub(GO_Term_ID, pattern = '\\)', replacement = ''))
    ) |>
  dplyr::mutate(
    GO_Type = "MF", .after = GO_Term_ID
    ) |> # create the category column
  
  dplyr::mutate(
    fold.Enrichment = stringr::str_extract(fold.Enrichment,
                                           pattern = '[:alnum:]+') |> 
      as.numeric(fold.Enrichment)
    ) |> # convert the FE column to numeric
  
  dplyr::arrange(
    FDR, desc(fold.Enrichment)
    ) # sort by the highest FE and lowest adjusted p-value
```

### Combine the Results

``` r
GO_results <- GO_BP_results |>
  dplyr::full_join(GO_CC_results) |>
  dplyr::full_join(GO_MF_results) |>
  dplyr::filter(!grepl(GO_Term, pattern = 'unclassified', ignore.case = T))
```

## Integrate with Custom GO Query Function

Call a custom function (`TMEM` package) created to aggregate all GO information
for a provided list of gene/protein identifiers

- queries both `Bioconductor` and `GO` dbs

### Query Bioconductor and GO dbs

Use the `TMEM` package query function to acquire `GO` info for plotting

``` r
x <- TMEM::get_GO_info(list_of_interest = Crapome_hits, species = 'human')
```

### Integrate Query Results into Counts DF

``` r
GO_info_by_term_df_sig <- GO_results |>
  dplyr::mutate(db = 'geneontology.org') |>
  dplyr::full_join(x$GO_info_by_term_df |> 
                     dplyr::mutate(db = 'bioconductor')) |>
  dplyr::filter(!is.na(GO_Term)) |>
  
  # dplyr::filter(!is.na(GO_Type)) |> # for discrepancies between 
  # geneontology.org and bioconductor; default to geneontology.org
  
  dplyr::mutate(fold.Enrichment = as.numeric(fold.Enrichment)) |>
  dplyr::arrange(desc(fold.Enrichment), FDR) |>
  dplyr::filter(!grepl(GO_Term, pattern = "Unclassified", ignore.case = T)) |>
  dplyr::mutate(sig = dplyr::case_when(
    
    FDR <= 0.001 ~ "***",
    dplyr::between(FDR, 0.001, 0.0499) ~ "**",
    dplyr::between(FDR, 0.01, 0.05) ~ "*",
    TRUE ~ "NS"), .after = FDR
  ) # create a statistical significance column for filtering/annotation
```

## Plot Results

Create DFs for each plot

### Plot Full Scatter

- Full overview

``` r
## redefine protiens of interest for labeling (control baits)
ALL_DF_filtered_full <- ALL_DF_filtered |>
  dplyr::mutate(
    
    poi = dplyr::case_when(
      grepl(Proteins.Identified, pattern = 'Green fluorescent') ~ 'Control myc',
      grepl(Proteins.Identified, pattern = 'BirA') ~ 'Control V5',
      TRUE ~ poi)
    
    ) # full scatter
```

``` r
ggplot(data = ALL_DF_filtered_full) +
  geom_point(data = subset(ALL_DF_filtered_full,
                           poi == 'Background'),
             aes(x = ControlAvg,
                 y = TMEMAvg,
                 color = poi), alpha = 0.1, size = 1.5) +
  
  geom_point(data = subset(ALL_DF_filtered_full,
                           poi == 'Candidates'),
             aes(x = ControlAvg,
                 y = TMEMAvg,
                 color = poi), alpha = 0.3, size = 3) +
  
  geom_point(data = subset(ALL_DF_filtered_full,
                           poi == 'Control myc'),
             aes(x = ControlAvg,
                 y = TMEMAvg,
                 color = poi), alpha = 1.0, size = 3) +
  
  geom_point(data = subset(ALL_DF_filtered_full,
                           poi == 'Control V5'),
             aes(x = ControlAvg,
                 y = TMEMAvg,
                 color = poi), alpha = 1.0, size = 3) +
  
  geom_point(data = subset(ALL_DF_filtered_full,
                           poi == 'TMEM184B'),
             aes(x = ControlAvg,
                 y = TMEMAvg,
                 color = poi), alpha = 1.0, size = 3) +
  coord_cartesian(xlim = c(0,17), ylim = c(0,17)) +
  labs(title = 'TMEM184B IPMS Total Spectral Counts',
       x = 'Average Control IP Prey Spectra (TSC)',
       y = 'Average TMEM184B IP Prey Spectra (TSC)',
       col = 'Protein\nCategory') +
  #guides(col = 'none')+
  geom_abline(slope = 2,
              intercept = 0,
              linetype = 'dashed',
              color = 'black') +
  geom_text(x = 6.5,
            y = 15,
            color = 'black',
            label = 'FC \u2265 2.0',
            fontface = 'bold',
            size = 4.5) +
  scale_color_manual(
    values = c('gray48', 'navy', 'green1', 'deeppink3', 'darkgoldenrod3')
    ) +
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5, size = 20),
        axis.text = element_text(face = 'bold', size = 12),
        axis.title = element_text(face = 'bold', size = 15),
        legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.5)
        )
```

![](https://github.com/eriklarsen4/Endo/analysis/R/IPMS_plots/Full%20Scatter-1.jpeg)<!-- -->

### vATPase Scatter (zoom)

- vATPase hits

``` r
ALL_DF_filtered_vATPASE_acidification <- ALL_DF_filtered |> 
  dplyr::mutate(poi = dplyr::case_when(
    
    (grepl(HumanGeneID, pattern = 'ATP6V0A(1|2)') |
       grepl(HumanGeneID, pattern = 'DMXL|CLN6')) &
      grepl(HumanGeneID, pattern = paste0(Crapome_hits, collapse = '$|')) ~ 
      
      'Candidates in\nVacuolar Acidification',
    
    grepl(HumanGeneID, pattern = 'ATP6V') &
      !grepl(HumanGeneID, pattern = paste0(Crapome_hits, collapse = '$|')) ~ 
      
      'Other vATPase\nComplex Subunits',
    
    TRUE ~ 
      
      poi)) # vATPase scatter
```

``` r
ggplot(data = ALL_DF_filtered_vATPASE_acidification) +
  
  geom_point(data = subset(
    ALL_DF_filtered_vATPASE_acidification,
    poi == 'Background'),
                      aes(x = ControlAvg,
                          y = TMEMAvg,
                          color = poi), alpha = 0.1, size = 1) +
  
  geom_point(data = subset(
    ALL_DF_filtered_vATPASE_acidification,
    poi == 'Candidates'),
                      aes(x = ControlAvg,
                          y = TMEMAvg,
                          color = poi), alpha = 0.5, size = 3) +
  
  geom_point(data = subset(
    ALL_DF_filtered_vATPASE_acidification,
    poi == 'Other vATPase\nComplex Subunits'),
                      aes(x = ControlAvg,
                          y = TMEMAvg,
                          color = poi), alpha = 1.0, size = 3) +
  
  geom_point(data = subset(
    ALL_DF_filtered_vATPASE_acidification,
    poi == 'Candidates in\nVacuolar Acidification'),
                      aes(x = ControlAvg,
                          y = TMEMAvg,
                          color = poi), alpha = 1.0, size = 3) +
  
  geom_point(data = subset(
    ALL_DF_filtered_vATPASE_acidification,
    poi == 'TMEM184B'),
             aes(x = ControlAvg,
                 y = TMEMAvg,
                 color = poi), alpha = 1.0) +
  
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  labs(title = 'TMEM184B IPMS Total Spectral Counts',
       x = 'Average Control IP Prey Spectra (TSC)',
       y = 'Average TMEM184B IP Prey Spectra (TSC)',
       col = 'Protein\nCategory') +
  #guides(col = 'none')+
  geom_abline(slope = 2,
                       intercept = 0,
                       linetype = 'dashed',
                       color = 'black',
                       size = 1) +
  geom_text(x = 1.2,
                     y = 3,
                     color = 'black',
                     label = 'FC \u2265 2.0',
                     fontface = 'bold',
                     size = 4.5) +
  scale_color_manual(
    values = c('gray48', 'dodgerblue', 'navy', 'magenta', 'darkgoldenrod3')
    ) +
  theme_classic()+
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 20),
    axis.text = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 15),
    legend.text = element_text(face = 'bold', size = 14),
    legend.title = element_blank(),
    legend.position = c(0.84, 0.5)
  )
```

![](https://github.com/eriklarsen4/Endo/analysis/R/IPMS_plots/vATPase%20Scatter%20zoom-1.jpeg)<!-- -->

### Regulation of Macroautophagy (zoom)

- Macroautophagy hits

``` r
ALL_DF_filterd_reg_macro <- ALL_DF_filtered |>
  dplyr::mutate(poi = dplyr::case_when(
    
    grepl(HumanGeneID,
          pattern = 'ATP6V0A(1|2)|EXOC(1|8)|ATG2A|SPTLC2|TSC2|CALCOCO2|AMBRA1'
          ) &
      grepl(HumanGeneID,
            pattern = paste0(Crapome_hits,
                             collapse = '$|')
            ) ~ 
      
      'Candidates in\nRegulation of Macroautophagy',
    
    TRUE ~ 
      
      poi)
    ) # regulation of macroautophagy
```

``` r
ggplot(data = ALL_DF_filterd_reg_macro) +
  
  geom_point(data = subset(
    ALL_DF_filterd_reg_macro,
    poi == 'Background'),
             aes(x = ControlAvg,
                 y = TMEMAvg,
                 color = poi), alpha = 0.1, size = 1) +
  
  geom_point(data = subset(
    ALL_DF_filterd_reg_macro,
    poi == 'Candidates'),
             aes(x = ControlAvg,
                 y = TMEMAvg,
                 color = poi), alpha = 0.5, size = 3) +
  
  geom_point(data = subset(
    ALL_DF_filterd_reg_macro,
    poi == 'Candidates in\nRegulation of Macroautophagy'),
             aes(x = ControlAvg,
                 y = TMEMAvg,
                 color = poi), alpha = 1.0, size = 3) +
  
  geom_point(data = subset(
    ALL_DF_filterd_reg_macro,
    poi == 'TMEM184B'),
             aes(x = ControlAvg,
                 y = TMEMAvg,
                 color = poi), alpha = 1.0, size = 3) +
  
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  labs(title = 'TMEM184B IPMS Total Spectral Counts',
       x = 'Average Control IP Prey Spectra (TSC)',
       y = 'Average TMEM184B IP Prey Spectra (TSC)',
       col = 'Protein\nCategory') +
  #guides(col = 'none')+
  geom_abline(slope = 2,
              intercept = 0,
              linetype = 'dashed',
              color = 'black',
              size = 1) +
  geom_text(x = 1.2,
            y = 3,
            color = 'black',
            label = 'FC \u2265 2.0',
            fontface = 'bold',
            size = 4.5) +
  scale_color_manual(
    values = c('gray48', 'dodgerblue', 'navy', 'darkgoldenrod3')
    ) +
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5, size = 20),
        axis.text = element_text(face = 'bold', size = 12),
        axis.title = element_text(face = 'bold', size = 15),
        legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.3)
  )
```

![](https://github.com/eriklarsen4/Endo/analysis/R/IPMS_plots/Macro%20zoom-1.jpeg)<!-- -->

### Protein Transport (zoom)

- Intracellular protein transport hits

``` r
ALL_DF_filterd_IC_protein_transport <- ALL_DF_filtered |>
  dplyr::mutate(
    poi = dplyr::case_when(
      grepl(
        HumanGeneID,
        pattern = '^(CCDC22)$|^(GLE1)$|^(TRAPPC3)$|^(RBM15B)$|^(EHD3)$|^(COG3)$|^(AUP1)$|^(TSC2)$|^(EHD1)$|^(VPS4A)$|^(MTX1)$|^(NOP9)$|^(TANC2)$|^(NDC1)$|^(DERL1)$|^(PEX3)$|^(EHD4)$|^(SMG6)$|^(VPS50)$|^(VPS4A)$') ~ 
      
      'Candidates in\nIntracellular Transport',
    
    TRUE ~ 
      
      poi)) # protein transport
```

``` r
ggplot(data = ALL_DF_filterd_IC_protein_transport) +
  
  geom_point(data = subset(ALL_DF_filterd_IC_protein_transport,
                           poi == 'Background'),
             aes(x = ControlAvg,
                 y = TMEMAvg,
                 color = poi), alpha = 0.1, size = 1) +
  
  geom_point(data = subset(ALL_DF_filterd_IC_protein_transport,
                           poi == 'Candidates'),
             aes(x = ControlAvg,
                 y = TMEMAvg,
                 color = poi), alpha = 0.5, size = 3) +
  
  geom_point(data = subset(ALL_DF_filterd_IC_protein_transport,
                           poi == 'Candidates in\nIntracellular Transport'),
             aes(x = ControlAvg,
                 y = TMEMAvg,
                 color = poi), alpha = 1.0, size = 3) +
  
  geom_point(data = subset(ALL_DF_filterd_IC_protein_transport,
                           poi == 'TMEM184B'),
             aes(x = ControlAvg,
                 y = TMEMAvg,
                 color = poi), alpha = 1.0, size = 3) +
  
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  labs(title = 'TMEM184B IPMS Total Spectral Counts',
       x = 'Average Control IP Prey Spectra (TSC)',
       y = 'Average TMEM184B IP Prey Spectra (TSC)',
       col = 'Protein\nCategory') +
  #guides(col = 'none')+
  geom_abline(slope = 2,
              intercept = 0,
              linetype = 'dashed',
              color = 'black',
              size = 1) +
  geom_text(x = 1.2,
            y = 3,
            color = 'black',
            label = 'FC \u2265 2.0',
            fontface = 'bold',
            size = 4.5) +
  scale_color_manual(
    values = c('gray48', 'dodgerblue', 'navy', 'darkgoldenrod3')
    ) +
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5, size = 20),
        axis.text = element_text(face = 'bold', size = 12),
        axis.title = element_text(face = 'bold', size = 15),
        legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_blank(),
        legend.position = c(0.78, 0.4)
  )
```

![](https://github.com/eriklarsen4/Endo/analysis/R/IPMS_plots/IC%20transport%20zoom-1.jpeg)<!-- -->

### GO Bar Plot

``` r
GO_info_by_term_df_sig |> 
  dplyr::filter(sig != 'NS') |> 
  dplyr::arrange(desc(fold.Enrichment)) |> 
  dplyr::group_by(GO_Type) |> 
  dplyr::slice(c(1:10)) |> # top ten GO terms for each type
  
  # dplyr::mutate(GO_Term = dplyr::case_when(
  #   grepl(GO_Term,
  #         pattern = 'via multivesicular') ~
  # 
  #     gsub(GO_Term,
  #          pattern = '(^.+)(\\s)(via.+$)', replacement = '\\1\\2\n\\3'),
  # 
  #   TRUE ~
  # 
  #     GO_Term)) |>
  
  ggplot(aes(x = -log10(FDR),
             y = forcats::fct_rev(forcats::fct_infreq(GO_Term, -log10(FDR))),
             color = 'black',
             fill = GO_Type,
             alpha = fold.Enrichment)) +
  
  geom_bar(color = 'black', stat = 'identity') +
  labs(title = 'TMEM184B IPMS Candidate GO Analysis',
       x = bquote(~ -log[10] ~ '(Adj.P-value)'),
         # expression(paste0("-log",[10],'(Adj.P-value)')),
       y = 'GO Term',
       fill = 'GO\nCategory',
       size = '# Genes in List',
       alpha = 'Fold\nEnrichment') +
  geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'black') +
  geom_text(aes(label = 'FDR\n\u2264 0.05'),
            size = 5,
            color = 'black',
            x = 3,
            y = 4,
            show.legend = FALSE) +
  geom_text(aes(label = sig,
                fontface = 'bold'),
            size = 2.5,
            color = 'black',
            angle = 90,
            vjust = 1,
            alpha = 1) +
  # geom_text(aes(label = paste0(Overlap, ' / ', GO_Term_Size),
  #               fontface = 'bold'),
  #           hjust = 0,
  #           size = 5,
  #           color = 'black',
  #           alpha = 0.7) +
  theme_bw()+
  theme(plot.title = element_text(face = 'bold', size = 13, hjust = 0.5),
        axis.title = element_text(face = 'bold', size = 12),
        axis.text = element_text(face = 'bold', size = 10, vjust = 0.5),
        axis.ticks = element_line(linetype = 'solid'),
        legend.title = element_text(face = 'bold', size = 10),
        legend.text = element_text(face = 'bold', size = 9)) +
  guides(col = 'none',
         fill = guide_legend(
           override.aes = list(
             size = 2.5, 
             shape = 19,
             color = 'black',
             fill = c('navy', 'darkgoldenrod3', 'darkgray'))
           ),
         alpha = guide_legend(override.aes = list(size = 2.5, shape = 19))) +
  scale_color_manual(values = c('navy', 'darkgoldenrod3', 'darkgray')) +
  scale_fill_manual(values = c('navy', 'darkgoldenrod3', 'darkgray'))
```

![](https://github.com/eriklarsen4/Endo/analysis/R/IPMS_plots/GO%20Bar%20Plot-1.jpeg)<!-- -->


Interpretations:  
- in general, higher fold enrichment is common among GO terms with few annotated
members:  
  + explicitly, having a higher percentage of genes/proteins from the list of
  interest appear in a GO term (because of it having few gene/protein members)
  leads to that term having a high fold enrichment
  
  + conversely, large, vague and general GO terms tend to have weak fold
  enrichment  
  
  + therefore, large GO terms tend to have high significance and low fold
  enrichment, and small GO terms tend to have high fold enrichment and lower
  significance
  
- thus, depending on the experiment and the questions, the most informative
results for unbiased screens are those that balance these two features

- of course, the most informative and powerful results usually are combined with
other (possibly *a priori*) information

- in our case, we had other preliminary evidence that suggested ties to pH 
and/or ionic differences in mutant mouse neurons, with phenotypes resembling
those of autophagy perturbation

### GO CC Bar Plot

``` r
GO_info_by_term_df_sig |> 
  dplyr::filter(sig != "NS") |> 
  dplyr::arrange(desc(fold.Enrichment)) |> 
  dplyr::filter(GO_Type == 'CC') |> 
  dplyr::slice(c(1:10)) |> # top 10 GO terms in cellular components
  dplyr::mutate(GO_Term = dplyr::case_when(
    
    grepl(GO_Term, pattern = 'via multivesicular') ~ 
      
      gsub(GO_Term,
           pattern = '(^.+)(\\s)(via.+$)', replacement = '\\1\\2\n\\3'),
    
    TRUE ~ 
      
      GO_Term)
    
    ) |>
  ggplot(aes(x = -log10(FDR),
             y = forcats::fct_rev(forcats::fct_infreq(GO_Term, -log10(FDR))),
             color = 'navy',
             fill = 'navy',
             alpha = -log10(FDR))) +
  geom_bar(color = 'black',
           stat = 'identity') +
  labs(title = 'TMEM184B IPMS Candidate GO Analysis',
       x = bquote(~ -log[10] ~ '(Adj.P-value)'),
       y = 'GO Cell Component Term',
       alpha = bquote(~ -log[10] ~ '(Adj.P-value)')) +
  theme_bw()+
  theme(plot.title = element_text(face = 'bold', size = 13),
        axis.title = element_text(face = 'bold', size = 12),
        axis.text = element_text(face = 'bold', size = 12, vjust = 0.5),
        axis.ticks = element_line(linetype = 'solid'),
        legend.title = element_text(face = 'bold', size = 12),
        legend.text = element_text(face = 'bold', size = 12)) +
  scale_color_manual(values = c('navy')) +
  scale_fill_manual(values = c('navy')) +
  guides(col = 'none',
         fill = 'none',
         alpha = guide_legend(
           
           override.aes = list(size = 2.5,
                               shape = 19,
                               color = 'black',
                               fill = 'navy'))) +
  geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'black') +
  geom_text(aes(label = 'FDR\n\u2264 0.05'),
            size = 5,
            color = 'black',
            x = 3,
            y = 4,
            show.legend = FALSE) +
  geom_text(aes(label = sig,
                fontface = 'bold'),
            size = 5,
            color = 'black',
            angle = 90,
            vjust = 1,
            alpha = 1)
```

![](https://github.com/eriklarsen4/Endo/analysis/R/IPMS_plots/GO%20CC%20Bar%20Plot-1.jpeg)<!-- -->

### GO BP Bar Plot

``` r
GO_info_by_term_df_sig |> 
  dplyr::filter(sig != "NS") |> 
  dplyr::arrange(desc(fold.Enrichment)) |> 
  dplyr::filter(GO_Type == 'BP') |> 
  dplyr::slice(c(1:10)) |> # top ten GO terms in biological process
  dplyr::mutate(GO_Term = dplyr::case_when(
    
    grepl(GO_Term, pattern = 'via multivesicular') ~ 
      
      gsub(GO_Term,
           pattern = '(^.+)(\\s)(via.+$)', replacement = '\\1\\2\n\\3'),
    
    TRUE ~ 
      
      GO_Term)
    ) |>
  ggplot(aes(x = -log10(FDR),
             y = forcats::fct_rev(forcats::fct_infreq(GO_Term, -log10(FDR))),
             color = 'navy',
             fill = 'navy',
             alpha = -log10(FDR))) +
  geom_bar(color = 'black', stat = 'identity') +
  labs(title = 'TMEM184B IPMS Candidate GO Analysis',
       x = bquote(~ -log[10] ~ '(Adj.P-value)'),
       y = 'GO Biological Process Term',
       alpha = bquote(~ -log[10] ~ '(Adj.P-value)')) +
  theme_bw()+
  theme(plot.title = element_text(face = 'bold', size = 13),
        axis.title = element_text(face = 'bold', size = 12),
        axis.text = element_text(face = 'bold', size = 12, vjust = 0.5),
        axis.ticks = element_line(linetype = 'solid'),
        legend.text = element_text(face = 'bold', size = 12),
        legend.title = element_text(face = 'bold', size = 12)) +
  scale_color_manual(values = c('navy')) +
  scale_fill_manual(values = c('navy')) +
  guides(col = 'none',
         fill = 'none',
         alpha = guide_legend(override.aes = list(size = 2.5,
                                                  shape = 19,
                                                  color = 'black',
                                                  fill = 'navy'))) +
  geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'black') +
  geom_text(aes(label = 'FDR\n\u2264 0.05'),
            size = 5,
            color = 'black',
            x = 2.25,
            y = 4,
            show.legend = FALSE) +
  geom_text(aes(label = sig,
                fontface = 'bold'),
            size = 5,
            color = 'black',
            angle = 90,
            vjust = 1,
            alpha = 1)
```

![](https://github.com/eriklarsen4/Endo/analysis/R/IPMS_plots/GO%20BP%20Bar%20Plot-1.jpeg)<!-- -->
