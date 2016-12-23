from subprocess import call as sp_call

def gtf_to_gff(in_gtf, out_gff):
    """
    Prepare a gff (non-overlapping bins) based on a gtf.
    Calls provided (python2 only) dexseq_prepare_annotation.py script.

    :param in_gtf:
        Input .gtf file.
    :param out_gff:
        Output .gff file.
    :return:
        External process call status.
    """
    return sp_call(["/usr/bin/python",
                    "/home/enuser/.R/DEXSeq/python_scripts/dexseq_prepare_annotation.py",
                    in_gtf, out_gff])