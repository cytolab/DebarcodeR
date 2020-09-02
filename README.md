# DebarcodeR
## Benjamin Reisman, 2020
DebarcodeR is an R package for demultiplexing fluorescent cell barcoded (FCB) barcoded flow cytometer data. Fluorescent cell barcoding is a technique for implementing pooled sample multiplexing in fluoresence flow cytometery using amine reactive dyes to label each sample with unique barcodes by varying the concentration of one or more dyes. Sample can the be pooled, stained, and acquired in a single tube, which reduces instrument time and reagent consumption and improves data robustness. Following acquisition, the data must be 'debarcoded' (also known as devoncolution or demultiplexing) which traditionally has required manual biaxial gating. DebarcodeR implements an algorithm for automating this demultiplexing process which improves reproducibility and enables high throughput data processing. For more information about DebarcodeR, see our preprint (link) and for more information about FCB, see Krutzik et al. 2006.

# Installation
DebarcodeR reuses many of the classes and methods for working with FCS files implemented in the flowCore. Most workflows will also require the use of one or more of the complementary cytoverse pacakges such as CytoML, FlowWorkspace, and ggCyto. A Bioconductor version of DebarcodeR is in the works, but for now can be installed from github. 
```{r}
# install flowCore first
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("flowCore")

#Followed by DebarcodeR
remotes::install_github(bjreisman/DebarcodeR)
```
