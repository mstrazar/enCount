# Perform JunctionSeq differential Usage analysis

args = commandArgs(TRUE)
in_decoder   = args[1]   # decoder is the metadata file indicating sample IDs
in_gff       = args[2]   # input gff file with novel junctions
in_count_dir = args[3]   # input count files directory
out_dir      = args[4]   # Output directory for jscs analysis


library("JunctionSeq")


decoder.file <- in_decoder
decoder <- read.table(in_decoder, header=TRUE)

gff.file <- in_gff

countFiles <- paste0(in_count_dir, decoder$sample.ID,
    "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")


jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
            sample.names = decoder$sample.ID,
            condition=factor(decoder$group.ID),
            flat.gff.file = gff.file, nCores = 1,
            analysis.type = "junctionsAndExons");

prefix = paste0(out_dir, "test/")
dir.create(prefix, recursive=TRUE)
writeCompleteResults(jscs, outfile.prefix=prefix, save.jscs=TRUE)

prefix = paste0(out_dir, "plots/")
dir.create(prefix, recursive=TRUE)
buildAllPlots(jscs=jscs, outfile.prefix=prefix,
    use.plotting.device = "png", FDR.threshold=0.01);
