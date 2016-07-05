# Estimate size factors with QoRTs

args = commandArgs(TRUE)
in_dir     = args[1]    # directory contains raw count files
in_decoder = args[2]    # decoder is the metadata file indicating sample IDs
out_file   = args[3]    # size factor file to be generated


print(in_dir)
print(in_decoder)
print(out_file)

library("QoRTs")

directory = in_dir
decoder.data = read.table(in_decoder, header=T, stringsAsFactors=F);

# calc.DESeq2 = TRUE results in error for zero gene counts
res = read.qc.results.data(directory, decoder=decoder.data,
    calc.DESeq2 = FALSE, calc.edgeR = TRUE,);

get.size.factors(res, outfile=out_file, sf.method=c("edgeR"));
