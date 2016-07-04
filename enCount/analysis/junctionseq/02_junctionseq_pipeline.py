"""
Prototype JunctionSeq analysis pipeline.
"""
import enCount.externals.junctionseq as jseq
import enCount.config as cfg
import os
import csv
import sys
import multiprocessing as mp


try:
    n_workers = int(sys.argv[1])
except:
    print("Run JunctionSeq pipeline.")
    print("Usage: %s <n_workers>" % sys.argv[0])
    quit(1)

# Input
in_dir_bam      = os.path.join(cfg.data_root, "bam")
in_dir_control  = os.path.join(cfg.data_root, "controls")
in_metafile     = os.path.join(in_dir_bam, "metadata.tsv")
in_gtf          = cfg.start_gtf

# Output
data_output     = os.path.join(cfg.data_root, "junctionseq/")
out_dir_raw_cts = os.path.join(data_output,   "rawCts/")
out_dir_cts     = os.path.join(data_output,   "cts/")
out_gtf_dir     = os.path.join(data_output,   "gtf/")
out_jscs_dir    = os.path.join(data_output,   "jscs/")
out_size_factors = os.path.join(out_jscs_dir, "size_factors.txt")

# Read experiment ids
experiments = sorted(list(set([row[jseq.EXPERIMENT_ACCESSION] for row in
               csv.DictReader(open(in_metafile), delimiter="\t")])))

# Open worker pool
pool = mp.Pool(n_workers)

# 1. Generate decoders
# r = jseq.generate_decoders(in_metafile, in_dir_control, out_jscs_dir)
# assert r == 0

# 2.
# TODO: call counting in separate processes per sample


# 3. Call merge counts in parallel per experiment
# Use multiprocessing ; make sure arguments are in correct order
in_decoder_sample = os.path.join(out_jscs_dir, "decoder.bySample.txt")
inputs = []
for exp_dir in experiments:
    in_decoder_uid = os.path.join(out_jscs_dir, exp_dir, "decoder.byUID.txt")
    inputs.append((out_dir_raw_cts, in_decoder_uid, out_dir_cts,))

results = pool.starmap(jseq.run_QoRTs_merge_counts, inputs)
assert all(map(lambda r: r == 0, results))


# 4. Size factors are required for correctly merging novel splice sites
in_decoder_uid = os.path.join(out_jscs_dir, "decoder.byUID.txt")
jseq.run_QoRTs_size_factors(in_dir=out_dir_raw_cts,
                       in_decoder=in_decoder_uid,
                       out_file=out_size_factors)


# 5. Discover novel splice junctions for all experiments to generate a common
# .gff file
r = jseq.run_QoRTs_novel_splices(in_dir=out_dir_cts,
                                in_gtf=in_gtf,
                                in_size_factors=out_size_factors,
                                out_dir=out_gtf_dir)
assert r == 0

# 6. Run a JunctionSeq analysis in parallel
in_gff = os.path.join(out_gtf_dir, "withNovel.forJunctionSeq.gff.gz")
inputs = []
for exp_dir in experiments:
        in_decoder_sample = os.path.join(out_jscs_dir, exp_dir,
                                         "decoder.bySample.txt")
        out_dir = os.path.join(out_jscs_dir, exp_dir)
        inputs.append((in_decoder_sample, in_gff, out_gtf_dir, out_dir))

results = pool.starmap(jseq.run_JunctionSeq_analysis, inputs)
assert all(map(lambda r: r == 0, results))

