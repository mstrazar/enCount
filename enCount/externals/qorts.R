# Estimate size factors with QoRTs
library("QoRTs")

in_dir     = args[1]    # directory contains raw count files
in_decoder = args[2]    # decoder is the metadata file indicating sample IDs
out_file   = args[3]    # size factor file to be generated

res = read.qc.results.data(in_dir, decoder = in_decoder,
    calc.DESeq2 = TRUE, calc.edgeR = TRUE);
get.size.factors(res, outfile=out_file);
