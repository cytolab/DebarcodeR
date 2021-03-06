% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morphology_corr_knignenburg.R
\name{morphology_corr.knijnenburg}
\alias{morphology_corr.knijnenburg}
\title{Corrects morphology based on scatter or uptake control.}
\usage{
morphology_corr.knijnenburg(
  fcb,
  uptake,
  channel,
  fsc_ssc = c(fsc = "FSC-A", ssc = "SSC-A"),
  subsample = 10000,
  updateProgress = NULL,
  ret.model = FALSE
)
}
\arguments{
\item{fcb}{the barcoded dataframe, post compensation and preprocessing}

\item{uptake}{Optional: a dataframe consisting of all cells barcoded with a single level of the barcoding dye}

\item{channel}{The name (string) of the channel to be corrected, ie. the column name in 'fcb_df'}

\item{subsample}{Integer, number of cells to sample (with replacement) for the morphology correction, defaults to 10,000.}

\item{updateProgress}{used in reactive context (shiny) to return progress information to GUI#'}
}
\value{
a tibble/data.frame with the selected channel corrected for fsc and ssc
}
\description{
Corrects morphology based on scatter or uptake control.
}
\seealso{
\code{\link{selectDenseScatterArea}} \code{\link{doRegressConstrained}}
}
