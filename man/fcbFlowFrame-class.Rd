% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcbFlowFrame.R
\name{fcbFlowFrame-class}
\alias{fcbFlowFrame-class}
\alias{fcbFlowFrame}
\title{fcbFlowFrame
-----------------------------------------------------------------------------
A containiner for barcoded flow cytometry data, with slots for the barcoded
flowFrame, a single level 'uptake control', and a slot to contain the results
of the debarcoding functions contained within debarcoder.}
\usage{
fcbFlowFrame(x, barcodes = list())
}
\description{
fcbFlowFrame
-----------------------------------------------------------------------------
A containiner for barcoded flow cytometry data, with slots for the barcoded
flowFrame, a single level 'uptake control', and a slot to contain the results
of the debarcoding functions contained within debarcoder.
}
\section{Slots}{

\describe{
\item{\code{barcoded.ff}}{{Object of class\code{flowFrame} containing the barcoded
data, approriately compensated, transformed, and gated}}

\item{\code{uptake.ff}}{{object of class \code{flowFrame} containing cells which
have been barcoded with a single level of each barcoding channel as well
as stained with the approriate barcoding controls}}

\item{\code{barcodes}}{{a list of barcodes, each one named for the channel from which
it was derived}}

\item{\code{platemap}}{{a platemap for conditions per barcode level}}
}}

