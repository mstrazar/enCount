# Read user-defined r library location
.libPaths("/home/enuser/.R/")
sessionInfo()

# Install via CRAN
install.packages("http://hartleys.github.io/QoRTs/QoRTs_LATEST.tar.gz",
    repos=NULL,type="source");

# Install via Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("JunctionSeq", ask=FALSE)
biocLite("DEXSeq", ask=FALSE)