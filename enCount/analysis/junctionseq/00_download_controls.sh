#!/bin/bash
# Generating metadata files with experimental design in separate files.
# Apparently, the only way to get the correct control is to request the specific
# constrol as a http request.

# We start with input METADATA file. Correct metadata files will be generated
# accordingly.

CONTROLID="Non-specific target control-human"
BATCH="Homo_sapiens.GRCh37.75.altevent.unique"

RNASEQDIR="/n/users/martins/Dev/data/encode/rnaseq/bam/"
RNASEQDATA="$RNASEQDIR"/metadata.tsv
DEXSEQDIR="/n/users/martins/Dev/data/encode/rnaseq/junctionseq/$BATCH/"
CONTROLDIR="/n/users/martins/Dev/data/encode/rnaseq/controls"


ACCESSIONS=`grep "$CONTROLID" $RNASEQDATA | awk -F "\t" '{ print $4 }' | sort | uniq `

# For each control experiment, download the target (knockdown) experiment
for ACC in $ACCESSIONS ; do
    echo Downloading accession $ACC
    PARAMETERS="type=Experiment&assay_title=shRNA+RNA-seq&assembly=hg19&replicates.library.biosample.life_stage=adult&replicates.library.nucleic_acid_term_name=polyadenylated+mRNA&files.file_type=bam&&possible_controls.accession=$ACC"
    wget "https://www.encodeproject.org/metadata/$PARAMETERS/metadata.tsv" -O $CONTROLDIR/$ACC.controls.tsv
done

TARGETS=`grep -v "$CONTROLID" $RNASEQDATA | grep -v "File accession" | awk -F "\t" '{ print $4 }' | sort | uniq ` ;

# For each target (knockdown), append the appropriate controls to file determining experimental design
for TGT in $TARGETS ; do
    target=`grep "$TGT" $RNASEQDATA | awk -F "\t" '{ print $16 }' | head -n 1`
    echo Finding controls for target $TGT $target

    mkdir -p $DEXSEQDIR/$target/
    METAFILE=$DEXSEQDIR/$target/metadata.tsv

    head -n 1 $RNASEQDATA > $METAFILE
    grep "$TGT" $RNASEQDATA >> $METAFILE

    for ACC in $ACCESSIONS ; do
        test=`grep "$TGT" $CONTROLDIR/$ACC.controls.tsv | wc -l`

        if [[ $test > 0 ]] ; then
            # Found proper controls - put into metafile
            grep "$ACC" $RNASEQDATA >> $METAFILE
        fi
    done
done

