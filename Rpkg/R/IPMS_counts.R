#' IPMS Spectral Counts from the JCS Study
#'
#' Raw immunoprecipitation mass spectrometry (IPMS) spectral count data 
#' originating from the JCS publication associated with the ORION project. 
#' This dataset represents the unprocessed spectral counts prior to any 
#' filtering, normalization, or IRIS-specific manipulation performed in the 
#' ORION analysis pipeline.
#'
#' The version included in this package is the *published form* of the 
#' dataset. Any downstream transformations (e.g., filtering, CRAPome 
#' adjustments, IRIS-ready formatting) are performed in the analysis scripts 
#' under \code{analysis/R/} and are not included in the package.
#'
#' @format A data frame or tibble containing raw spectral counts. The exact 
#' column structure reflects the published supplementary data.
#'
#' @usage data(IPMS_counts)
#'
#' @source JCS Publication: <https://journals.biologists.com/jcs/article/138/15/jcs263908/368852/TMEM184B-modulates-endolysosomal-acidification-via?searchresult=1>
#'
"IPMS_counts"
