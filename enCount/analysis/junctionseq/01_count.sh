#!/bin/bash
IN_METADATA=$1
IN_BAM=$2
IN_GTF=$3
OUTDIR=$4
N=$5
QORTS_JAR=../../externals/libs/QoRTs.jar

# Kill child processes on SIGSTOP
# Remove remaining temporary files
trap ctrl_c INT
function ctrl_c() {
        echo "** Trapped CTRL-C"

        # Kill child processes
        pkill --signal 9 -P $$
        exit
}

if ! [ -e $IN_GTF ] ; then
    echo "Usage: $0 <input metadata> <input bam directory> <input gtf> <output dir> <no. parallel processes>"
    exit 1
fi

count = 0
while read line ; do
    acc=`echo $line | cut -f1 -d" "`
    f=$IN_BAM/$acc.bam

    if ! [ -e $f ]; then
        echo "File $f does not exist, skipping ..."
        continue
    fi

    # Run QoRTs count
    OUT_DIR_ACC=$OUTDIR/$acc/
    mkdir -p $OUT_DIR_ACC

    echo java -jar $QORTS_JAR QC --stranded $f $IN_GTF $OUT_DIR_ACC
    java -jar $QORTS_JAR QC --stranded $f $IN_GTF $OUT_DIR_ACC 2>>$0.stderr 1>/dev/null &
    count=$(($count+1))

    if [[ $(($count % $N)) == 0 ]] ; then
            echo "Waiting ..."
            wait
    fi

done < $IN_METADATA


