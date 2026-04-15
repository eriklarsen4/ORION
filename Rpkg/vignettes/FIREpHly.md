FIREpHly
================
Erik Larsen
2025-01-22

### Overview

This markdown was developed to formally document the data analysis of
FIRE-pHly data for **Dr. Martha Bhattacharya’s lab** at the **University
of Arizona** for the [endolysosomal paper](<https://journals.biologists.com/jcs/article/doi/10.1242/jcs.263908/368440/Transmembrane-protein-184B-TMEM184B-modulates>)

Until an image is built, the `R version` used was **R 4.3.2** with the
packages installed below

### Set up Environment

#### Attach packages

``` r
packages <- c('Endo',
              'tidyverse',
              'ggplot2',
              'ggtext',
              # 'readr',
              'ggsignif',
              'emmeans',
              'broom'
              )
for ( package in packages ) {
  if (!require(package, character.only = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}
```

### Import FIRE-pHly Data

<span style="color:black;font-size: 18px">Directly load the
`tidied_puncta_data`</span>

``` r
data(tidied_puncta_data)
Puncta_data <- tidied_puncta_data
rm(tidied_puncta_data)
```

<span style="color:black;font-size: 18px">Code chunks that clean the raw
data from separate Excel sheets into the data object,
`tidied_puncta_data`, have been included in the [R Markdown
vignette](https://github.com/eriklarsen4/Endo/Rpkg/vignettes/FIREpHly.Rmd)</span>

### Visualize

<span style="color:black;font-size: 16px">Investigate the distributions
of measured variables</span>

#### Puncta Area Histogram

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Puncta%20Area%20Histogram%20full-1.jpeg)<!-- -->

#### Puncta Area Histogram (zoom)

<span style="color:black;font-size: 16px">Crop to better-discern
differences in the majority of the data</span>

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Puncta%20Area%20Histogram%20cropped-1.jpeg)<!-- -->

#### Puncta Area Density

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Puncta%20Area%20Density%20full-1.jpeg)<!-- -->

#### Puncta Area Density (zoom)

<span style="color:black;font-size: 16px">Crop to better-discern
differences in the majority of the data</span>

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Puncta%20Area%20Density%20cropped-1.jpeg)<!-- -->

#### Puncta Area Violin

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Puncta%20Area%20Violin-1.jpeg)<!-- -->

#### Puncta Acidity Histogram

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Puncta%20Acidity%20Histogram%20full-1.jpeg)<!-- -->

#### Puncta Acidity Histogram (zoom)

<span style="color:black;font-size: 16px">Crop to better-discern
differences in the majority of the data</span>

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Puncta%20Acidity%20Histogram%20cropped-1.jpeg)<!-- -->

#### Puncta Acidity Density

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Puncta%20Acidity%20Density%20full-1.jpeg)<!-- -->

#### Puncta Acidity Density (zoom)

<span style="color:black;font-size: 16px">Crop to better-discern
differences in the majority of the data</span>

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Puncta%20Acidity%20Density%20cropped-1.jpeg)<!-- -->

#### Puncta Area and Acidity Contour

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Puncta%20Area%20and%20Acidity%202d%20density%20full-1.jpeg)<!-- -->

#### Puncta Area and Acidity Contour (zoom)

<span style="color:black;font-size: 16px">Crop to better-discern
differences in the majority of the data</span>

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Puncta%20Area%20and%20Acidity%202d%20density%20cropped-1.jpeg)<!-- -->

### Statistics

<span style="color:black;font-size: 16px">Perform weighted linear
regressions on puncta data</span>

<span style="color:black;font-size: 18px">It appropriately accounts
for:</span>

- <span style="color:black;font-size: 16px">experimenter bias</span>  
- <span style="color:black;font-size: 16px">inter-mouse
  variability</span>  
- <span style="color:black;font-size: 16px">(genotypic) inter-puncta
  variability</span>

#### Puncta Acidity

    ## 
    ## Call:
    ## lm(formula = Cell_Avg_Puncta_GreenRedRatio ~ Genotype, data = Puncta_avg_acidity, 
    ##     weights = mouse_puncta_w)
    ## 
    ## Weighted Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.85269 -0.36344 -0.05489  0.14189  1.12221 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   1.4302     0.1643   8.703 2.08e-08 ***
    ## GenotypeWT   -0.7026     0.2324  -3.023  0.00647 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5197 on 21 degrees of freedom
    ## Multiple R-squared:  0.3032, Adjusted R-squared:   0.27 
    ## F-statistic: 9.139 on 1 and 21 DF,  p-value: 0.006472

<span style="color:navy;font-size: 18px">Strong significant difference
between genotypes (mutant more alkaline)</span>

<span style="color:black;font-size: 16px">Check shaprio-wilk</span>

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  as.numeric(summary(puncta_acidity_lm)$residuals)
    ## W = 0.95054, p-value = 0.3003

<span style="color:black;font-size: 18px">Not non-normal</span>

#### Puncta Area

    ## 
    ## Call:
    ## lm(formula = Cell_Avg_Puncta_Area ~ Genotype, data = Puncta_avg_area, 
    ##     weights = mouse_puncta_w)
    ## 
    ## Weighted Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.17678 -0.09906 -0.02857  0.08363  0.30133 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.55752    0.04490  12.418 3.86e-11 ***
    ## GenotypeWT   0.01227    0.06349   0.193    0.849    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.142 on 21 degrees of freedom
    ## Multiple R-squared:  0.001775,   Adjusted R-squared:  -0.04576 
    ## F-statistic: 0.03734 on 1 and 21 DF,  p-value: 0.8486

<span style="color:black;font-size: 18px">No significant difference
between genotypes in auto-detected puncta in terms of their areas</span>

<span style="color:black;font-size: 16px">Check shapiro-wilk</span>

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  as.numeric(summary(puncta_area_lm)$residuals)
    ## W = 0.92732, p-value = 0.09572

<span style="color:black;font-size: 18px">Not non-normal</span>

#### Puncta Number (wlm)

    ## 
    ## Call:
    ## lm(formula = puncta_n ~ Genotype, data = Puncta_n, weights = mouse_puncta_w)
    ## 
    ## Weighted Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -54.201 -20.327  -1.051  27.863  52.679 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   24.692      2.402  10.281 1.19e-09 ***
    ## GenotypeWT    -6.492      3.165  -2.051    0.053 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 34.28 on 21 degrees of freedom
    ## Multiple R-squared:  0.1669, Adjusted R-squared:  0.1272 
    ## F-statistic: 4.207 on 1 and 21 DF,  p-value: 0.05295

<span style="color:black;font-size: 18px">No significant difference
between genotypes</span>

#### Puncta Number (Welch’s t)

``` r
## more traditional Welch's t
t.test(
  x = c(Puncta_data |>
          dplyr::group_by(Genotype) |>
          dplyr::distinct(Puncta_Count) |>
          dplyr::ungroup() |>
          dplyr::filter(Genotype == "WT") |>
          dplyr::select(Puncta_Count) |>
          unlist() |> as.numeric()
),
y = c(Puncta_data |>
          dplyr::group_by(Genotype) |>
          dplyr::distinct(Puncta_Count) |>
          dplyr::ungroup() |>
          dplyr::filter(Genotype == "MUT") |>
          dplyr::select(Puncta_Count) |>
          unlist() |> as.numeric()
),
alternative = 'two.sided', var.equal = F
)
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  c(as.numeric(unlist(dplyr::select(dplyr::filter(dplyr::ungroup(dplyr::distinct(dplyr::group_by(Puncta_data, Genotype), Puncta_Count)), Genotype == "WT"), Puncta_Count)))) and c(as.numeric(unlist(dplyr::select(dplyr::filter(dplyr::ungroup(dplyr::distinct(dplyr::group_by(Puncta_data, Genotype), Puncta_Count)), Genotype == "MUT"), Puncta_Count))))
    ## t = -1.3095, df = 17.876, p-value = 0.207
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -13.683890   3.178839
    ## sample estimates:
    ## mean of x mean of y 
    ##  21.11111  26.36364

<span style="color:black;font-size: 18px">No significant difference
between genotype-dependent puncta number auto-detected by computer
software</span>

<span style="color:black;font-size: 16px">Check shapiro-wilk</span>

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  as.numeric(summary(puncta_n_lm)$residuals)
    ## W = 0.95369, p-value = 0.3484

<span style="color:black;font-size: 18px">Not non-normal</span>

#### Cell Body Area

``` r
  ## Welch's t test
t.test(x = c(Puncta_data |>
               dplyr::group_by(Genotype, cell_num) |>
               dplyr::distinct(Cell_Body_Area) |>
               dplyr::ungroup() |>
               dplyr::filter(Genotype == "WT") |>
               dplyr::select(Cell_Body_Area) |>
               unlist() |> as.numeric()
),
y = c(Puncta_data |>
        dplyr::group_by(Genotype, cell_num) |>
        dplyr::distinct(Cell_Body_Area) |>
        dplyr::ungroup() |>
        dplyr::filter(Genotype == "MUT") |>
        dplyr::select(Cell_Body_Area) |>
        unlist() |> as.numeric()
),
alternative = 'two.sided', var.equal = F)
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  c(as.numeric(unlist(dplyr::select(dplyr::filter(dplyr::ungroup(dplyr::distinct(dplyr::group_by(Puncta_data, Genotype, cell_num), Cell_Body_Area)), Genotype == "WT"), Cell_Body_Area)))) and c(as.numeric(unlist(dplyr::select(dplyr::filter(dplyr::ungroup(dplyr::distinct(dplyr::group_by(Puncta_data, Genotype, cell_num), Cell_Body_Area)), Genotype == "MUT"), Cell_Body_Area))))
    ## t = -2.2192, df = 20.999, p-value = 0.03762
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -174.804621   -5.676421
    ## sample estimates:
    ## mean of x mean of y 
    ##  124.1645  214.4051

<span style="color:darkred;font-size: 18px">Mutant cell bodies are
significantly larger than wild-type cell bodies; this is likely due to
sampling bias and likely does not relate to pH-dependent results since
those are weighted by the number of puncta, which is not significantly
different</span>

<span style="color:black;font-size: 16px">Check shapiro-wilk</span>

``` r
  ## Examine residuals
    ## extract from lm (~ same as Welch's t test)
cell_area_lm <- lm(data = Puncta_data |>
                     dplyr::group_by(Genotype, cell_num) |>
                     dplyr::distinct(Cell_Body_Area) |>
                     dplyr::ungroup(),
                   formula = Cell_Body_Area ~ Genotype)

  ## Shapiro.wilk
shapiro_cell_area <- shapiro.test(
  summary(cell_area_lm)$residuals |> as.numeric())

shapiro.test(summary(cell_area_lm)$residuals |> as.numeric())
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  as.numeric(summary(cell_area_lm)$residuals)
    ## W = 0.90911, p-value = 0.03913

<span style="color:darkred;font-size: 18px">Non-normal residuals,
suggesting cell body size difference is difficult to properly discern a
significant difference due to outliers or, more broadly,
differently-distributed cell body areas. This is obvious and, again, not
necessarily relevant</span>

<span style="color:black;font-size: 16px">Significant. Show the
residuals</span>

#### Cell Body Area residuals dist

``` r
  ## Cell Body area residuals distributions
Puncta_data |>
  dplyr::group_by(Genotype, cell_num) |>
  dplyr::distinct(Cell_Body_Area) |>
  dplyr::ungroup() |>
  dplyr::mutate(Cell_Body_Area = as.numeric(Cell_Body_Area),
                resids = lm(data = Puncta_data |>
                              dplyr::group_by(Genotype, cell_num) |>
                              dplyr::distinct(Cell_Body_Area) |>
                              dplyr::ungroup(),
                            formula = as.numeric(Cell_Body_Area) ~ Genotype
                            )$residuals |> as.numeric()
                ) |>
  ggplot(aes(x = resids, y = Genotype, color = Genotype, fill = Genotype)) +
  geom_violin(alpha = 0.3) +
  stat_summary(fun.y = 'median',
               geom = 'point',
               color = 'black',
               size = 5,
               show.legend = FALSE) +
  guides(alpha = 'none') +
  geom_jitter(data = Puncta_data |>
                dplyr::group_by(Genotype, cell_num) |>
                dplyr::distinct(Cell_Body_Area) |>
                dplyr::ungroup() |>
                dplyr::mutate(
                  Cell_Body_Area = as.numeric(Cell_Body_Area),
                  resids = lm(data = Puncta_data |>
                                dplyr::group_by(Genotype, cell_num) |>
                                dplyr::distinct(Cell_Body_Area) |>
                                dplyr::ungroup(),
                              formula = as.numeric(Cell_Body_Area) ~ Genotype
                              )$residuals |> as.numeric()
                ),
              aes(x = resids,
                  y = Genotype), alpha = 1, width = 0.01, height = 0.1) +
  labs(title = 'Hippocampal Neuron Cell Body\nlm fit residual distributions',
       x = 'lm fit residuals',
       y = 'Genotype') +
  theme_bw() +
  theme(title = element_text(face = 'bold', size = 13),
        axis.title = element_text(face = 'bold', size = 11),
        axis.text = element_text(face = 'bold', size = 11),
        panel.grid = element_blank(),
        legend.text = element_text(face = 'bold', size = 10)
  ) +
  guides(alpha = 'none') +
  scale_color_manual(values = c('navy', 'darkgoldenrod3')) +
  scale_fill_manual(values = c('navy', 'darkgoldenrod3'))
```

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Cell%20Body%20Area%20lm%20residuals%20dist-1.jpeg)<!-- -->

<span style="color:darkred;font-size: 18px">The residual average cell
body size for each mouse is not the same between genotypes</span>

<span style="color:black;font-size: 16px">What’s the shape of the dv,
taking both genotypes into account?</span>

#### Cell Body Area density

``` r
Puncta_data |> dplyr::distinct(Cell_Body_Area) |> 
  ggplot(aes(x = as.numeric(Cell_Body_Area))) +
  geom_density(fill = 'black', alpha = 0.5) +
  labs(x = expression(paste("Cell Body Area (", mu, m^2, ")")))
```

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Cell%20Body%20Area%20density-1.jpeg)<!-- -->

<span style="color:darkred;font-size: 18px">Negative binomial.
Unsurprising, particularly for small sample sizes; the borderline
difference suggests the sample sizes are nearly powerful enough</span>

<span style="color:black;font-size: 18px">Still, consider
alternatives:</span>

- <span style="color:black;font-size: 16px">a `glm` with an
  **inverse-Gaussian fit** (probably the appropriate
  alternative)</span>  
- <span style="color:black;font-size: 16px">`Wilcoxon-ranked sum (Mann-Whitney)`
  (probably the next appropriate alternative)</span>  
- <span style="color:black;font-size: 16px">`Komolgorov-Smirnov`</span>

<span style="color:black;font-size: 16px">First, the `glm`</span>

``` r
  ## GLM
cell_area_glm <- glm(data = Puncta_data |>
                       dplyr::group_by(Genotype, cell_num) |>
                       dplyr::distinct(Cell_Body_Area) |>
                       dplyr::ungroup(),
                     formula = as.numeric(Cell_Body_Area) ~ Genotype,
                     family = 'Gamma')

summary(cell_area_glm)
```

    ## 
    ## Call:
    ## glm(formula = as.numeric(Cell_Body_Area) ~ Genotype, family = "Gamma", 
    ##     data = dplyr::ungroup(dplyr::distinct(dplyr::group_by(Puncta_data, 
    ##         Genotype, cell_num), Cell_Body_Area)))
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 0.0046641  0.0007643   6.102 4.69e-06 ***
    ## GenotypeWT  0.0033898  0.0016878   2.008   0.0576 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.3491339)
    ## 
    ##     Null deviance: 8.1667  on 22  degrees of freedom
    ## Residual deviance: 6.5386  on 21  degrees of freedom
    ## AIC: 272.89
    ## 
    ## Number of Fisher Scoring iterations: 5

<span style="color:navy;font-size: 18px">No significant difference
between genotypes when using a `glm` with an **inverse-Gaussian** to
account for different (non-normal) distribution shapes</span>

<span style="color:black;font-size: 18px">What about the
`MannWhitney`/`Wilcoxon-ranked sum`?</span>

- <span style="color:black;font-size: 16px">Compares the medians of two
  groups because the data violate assumptions (such as normality, as in
  this case)</span>  
- <span style="color:black;font-size: 16px">often, “outliers” skew
  distributions and/or unequal variances between groups</span>

``` r
wilcox.test(data = Puncta_data |>
              dplyr::group_by(Genotype, cell_num) |>
              dplyr::distinct(Cell_Body_Area) |>
              dplyr::ungroup(),
            as.numeric(Cell_Body_Area) ~ Genotype)
```

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  as.numeric(Cell_Body_Area) by Genotype
    ## W = 103, p-value = 0.01778
    ## alternative hypothesis: true location shift is not equal to 0

<span style="color:darkred;font-size: 18px">Significant difference
between the genotypes. Not sure how this contributes to pH difference in
endolysosomes marked by the FIREpHly marker (Lamp1)</span>

<span style="color:black;font-size: 18px">What about the
`Kolmogorov-Smirnov`?</span>

- <span style="color:black;font-size: 16px">`KS test` does *not* assume
  the sample means of two groups come from the same distribution (a more
  conservative, non-parametric test)</span>

<span style="color:black;font-size: 16px">We assume this in genetic
studies (hence the `glm`)</span>

``` r
  ## KS
ks.test(x = c(Puncta_data |>
                dplyr::group_by(Genotype, cell_num) |>
                dplyr::distinct(Cell_Body_Area) |>
                dplyr::ungroup() |>
                dplyr::filter(Genotype == "WT") |>
                dplyr::select(Cell_Body_Area) |>
                unlist() |> as.numeric()
        ),
        y = c(Puncta_data |>
                dplyr::group_by(Genotype, cell_num) |>
                dplyr::distinct(Cell_Body_Area) |>
                dplyr::ungroup() |>
                dplyr::filter(Genotype == "MUT") |>
                dplyr::select(Cell_Body_Area) |>
                unlist() |> as.numeric()
        ),
        alternative = 'two.sided')
```

    ## 
    ##  Exact two-sample Kolmogorov-Smirnov test
    ## 
    ## data:  c(as.numeric(unlist(dplyr::select(dplyr::filter(dplyr::ungroup(dplyr::distinct(dplyr::group_by(Puncta_data, Genotype, cell_num), Cell_Body_Area)), Genotype == "WT"), Cell_Body_Area)))) and c(as.numeric(unlist(dplyr::select(dplyr::filter(dplyr::ungroup(dplyr::distinct(dplyr::group_by(Puncta_data, Genotype, cell_num), Cell_Body_Area)), Genotype == "MUT"), Cell_Body_Area))))
    ## D = 0.64615, p-value = 0.009936
    ## alternative hypothesis: two-sided

<span style="color:darkred;font-size: 18px">Significant difference
between the genotypes.</span>

<span style="color:black;font-size: 18px">Defer to GraphPad Prism
statistical results, though I am satisfied by the `glm` with the
inverse-Gaussian fit; again, this is a minor issue, though I will note
it here</span>

<span style="color:black;font-size: 16px">Is there a relationship
between cell body size and puncta size?</span>
<span style="color:black;font-size: 16px">Do the genotypes differ? (see
below)</span>

#### Cell Body Area vs Avg Puncta Area scatter w/fit

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Cell%20Body%20Area%20Scatter%20w%20lm%20fit-1.jpeg)<!-- -->

#### Cell Body Area vs Avg Puncta Area summary

``` r
summary(lm(data = Puncta_data |>
     as.data.frame() |>
     dplyr::mutate(
       Puncta_Area = as.numeric(Puncta_Area),
       Cell_Body_Area = as.numeric(Cell_Body_Area),
       cell_num = as.numeric(cell_num)) |>
     dplyr::group_by(Genotype, cell_num) |>
     dplyr::mutate(mean_Puncta_Area = mean(Puncta_Area)) |>
     dplyr::distinct(Cell_Body_Area, mean_Puncta_Area) |>
     dplyr::ungroup(),
   
   formula = mean_Puncta_Area ~ 
     Cell_Body_Area*Genotype
     ))
```

    ## 
    ## Call:
    ## lm(formula = mean_Puncta_Area ~ Cell_Body_Area * Genotype, data = dplyr::ungroup(dplyr::distinct(dplyr::mutate(dplyr::group_by(dplyr::mutate(as.data.frame(Puncta_data), 
    ##     Puncta_Area = as.numeric(Puncta_Area), Cell_Body_Area = as.numeric(Cell_Body_Area), 
    ##     cell_num = as.numeric(cell_num)), Genotype, cell_num), mean_Puncta_Area = mean(Puncta_Area)), 
    ##     Cell_Body_Area, mean_Puncta_Area)))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.20418 -0.10868 -0.04210  0.08739  0.32149 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                5.505e-01  9.855e-02   5.586 2.19e-05 ***
    ## Cell_Body_Area             6.032e-05  4.115e-04   0.147    0.885    
    ## GenotypeWT                 4.689e-02  1.354e-01   0.346    0.733    
    ## Cell_Body_Area:GenotypeWT -2.979e-04  7.519e-04  -0.396    0.696    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1584 on 19 degrees of freedom
    ## Multiple R-squared:  0.008792,   Adjusted R-squared:  -0.1477 
    ## F-statistic: 0.05617 on 3 and 19 DF,  p-value: 0.982

<span style="color:black;font-size: 18px">No clear relationship between
average puncta area and cell body area between the two genotypes
(though, mutants have larger cell bodies)</span>

<span style="color:black;font-size: 16px">Thus, mutant endolysosomal
cell bodies may be significantly larger; however, there is no difference
in the puncta size between the two genotypes</span>

### Stats Plots

<span style="color:black;font-size: 16px">Visualize the statistical test
results</span>

#### Puncta Acidity Bar

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Puncta%20Acidity%20Bar-1.jpeg)<!-- -->

#### Puncta Area Bar

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Puncta%20Area%20Bar-1.jpeg)<!-- -->

#### Puncta Number Bar

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Puncta%20Number%20Bar-1.jpeg)<!-- -->

#### Cell Body Area Bar (Welch’s t)

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Cell%20Body%20Area%20Welch%20Bar-1.jpeg)<!-- -->

#### Cell Body Area Bar (glm)

![](https://github.com/eriklarsen4/Endo/analysis/R/FIREpHly_plots/Cell%20Body%20Area%20Gamma%20Bar-1.jpeg)<!-- -->

### Results Tables

#### Tidied Puncta Data

    ##   Genotype `cell num` `Puncta ID` `Puncta Circularity` `Puncta Area`
    ## 1       WT          1           1               0.9913        0.4290
    ## 2       WT          1           2               0.9913        0.4290
    ## 3       WT          1           3               0.9913        0.4290
    ## 4       WT          1           4               0.9103        0.8023
    ## 5       WT          1           5               0.9913        0.4290
    ## 6       WT          1           6               0.9913        0.4290
    ##   `Intra Puncta Green Intensity` `Intra Puncta Red Intensity` `GreenRed Ratio`
    ## 1                        39.9078                      32.6997           1.2204
    ## 2                        39.3823                      43.0102           0.9156
    ## 3                        48.0068                     112.3959           0.4271
    ## 4                        41.3431                      56.3704           0.7334
    ## 5                        47.8498                      68.9454           0.6940
    ## 6                        45.9352                      79.5324           0.5776
    ##   `Cell Body Area` `#` `Cell Avg Puncta Green Intensity`
    ## 1         123.3492  10                           41.4696
    ## 2         123.3492  10                           41.4696
    ## 3         123.3492  10                           41.4696
    ## 4         123.3492  10                           41.4696
    ## 5         123.3492  10                           41.4696
    ## 6         123.3492  10                           41.4696
    ##   `Cell Avg Puncta Red Intensity` `# Normalized to Cell Area`
    ## 1                          63.452                      8.1071
    ## 2                          63.452                      8.1071
    ## 3                          63.452                      8.1071
    ## 4                          63.452                      8.1071
    ## 5                          63.452                      8.1071
    ## 6                          63.452                      8.1071
    ##   `Cell Avg Puncta GreenRedRatio`
    ## 1                          0.6536
    ## 2                          0.6536
    ## 3                          0.6536
    ## 4                          0.6536
    ## 5                          0.6536
    ## 6                          0.6536

#### Table of Stats

    ## # A tibble: 6 × 18
    ##   target_var  term  estimate std.error statistic p.value r.squared adj.r.squared
    ##   <chr>       <chr>    <dbl>     <dbl>     <dbl>   <dbl>     <dbl>         <dbl>
    ## 1 puncta_aci… (Int…   1.43      0.164     8.70   0        NA             NA     
    ## 2 puncta_aci… Geno…  -0.703     0.232    -3.02   0.00647  NA             NA     
    ## 3 puncta_aci… <NA>   NA        NA         9.14   0.00647   0.303          0.270 
    ## 4 puncta_area (Int…   0.558     0.0449   12.4    0        NA             NA     
    ## 5 puncta_area Geno…   0.0123    0.0635    0.193  0.849    NA             NA     
    ## 6 puncta_area <NA>   NA        NA         0.0373 0.849     0.00177       -0.0458
    ## # ℹ 10 more variables: sigma <dbl>, df <dbl>, logLik <dbl>, AIC <dbl>,
    ## #   BIC <dbl>, deviance <dbl>, df.residual <dbl>, nobs <dbl>,
    ## #   null.deviance <dbl>, df.null <dbl>

#### Table of Predicted Values

``` r
lm_results_preds_and_resids %>%
  dplyr::mutate(across(c(2:ncol(.),-9), ~ as.numeric(.x))) %>%
  dplyr::mutate(across(c(2:ncol(.),-9), ~ round(.x, digits = 4))) |> 
  head()
```

    ## # A tibble: 6 × 14
    ##   Genotype `(weights)` .fitted  .resid   .hat .sigma .cooksd .std.resid
    ##   <chr>          <dbl>   <dbl>   <dbl>  <dbl>  <dbl>   <dbl>      <dbl>
    ## 1 WT             0.550   0.728 -0.074  0.0549  0.532  0.0003     -0.109
    ## 2 WT             1.32    0.728 -0.419  0.132   0.520  0.0749     -0.993
    ## 3 WT             0.989   0.728  0.186  0.0989  0.531  0.0078      0.376
    ## 4 WT             1.15    0.728  0.0599 0.115   0.532  0.0011      0.132
    ## 5 WT             1.26    0.728  0.541  0.126   0.512  0.113       1.25 
    ## 6 WT             0.769   0.728 -0.490  0.0769  0.523  0.0308     -0.860
    ## # ℹ 6 more variables: target_var <chr>, Cell_Avg_Puncta_GreenRedRatio <dbl>,
    ## #   Cell_Avg_Puncta_Area <dbl>, puncta_n <dbl>, Cell_Body_Area_lm <dbl>,
    ## #   Cell_Body_Area_glm <dbl>

#### Shapiro Results Table

``` r
shapiro_test_results %>%
  dplyr::mutate(across(c(1:2), ~ round(.x, digits = 4))) |> 
  head()
```

    ## # A tibble: 4 × 4
    ##   statistic p.value method                      target_var    
    ##       <dbl>   <dbl> <chr>                       <chr>         
    ## 1     0.950  0.300  Shapiro-Wilk normality test puncta_acidity
    ## 2     0.909  0.0391 Shapiro-Wilk normality test cell_body_area
    ## 3     0.927  0.0957 Shapiro-Wilk normality test puncta_area   
    ## 4     0.954  0.348  Shapiro-Wilk normality test puncta_n
