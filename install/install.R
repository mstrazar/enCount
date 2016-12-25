# Read user-defined r library location
.libPaths("/home/enuser/.R/")
sessionInfo()

# Install via Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("ShortRead"))
biocLite(c("DESeq2"))
biocLite(c("DEXSeq"))
biocLite(c("JunctionSeq"))

# Install manually via source
# Assume QoRTs_1.1.8.tar.gz is downloaded
install.packages("QoRTs_1.1.8.tar.gz", repos=NULL,type="source");
