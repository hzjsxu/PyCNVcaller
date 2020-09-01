#!/bin/sh
#---------------------------------------------------------------------------
# File: Individual_process.sh
# Created Date: 2020-08-30
# Author: jsxu
# Contact: <hzaujsxu@163.com>
# Last Modified: Sunday August 30th 2020 12:48:30 pm
# modified based on CNVcaller
#---------------------------------------------------------------------------

# samtools executables must be on path.
which samtools > /dev/null || exit 1
# CNVcaller is the installation directory for CNVcaller - set it as an enviromental variable.
export CNVcaller=/store/jsxu/biosoft/CNVcaller   # your CNVcaller path
echo "CNVcaller install directory $CNVcaller."
if [ ! -f "$CNVcaller/bin/1.1.CNVprocess.py" ]; then
    echo "You should set CNVcaller installation directory."
    exit 1
fi

#############################
usage() {
    echo "Usage: $0 -b <BAM> -h <header> -d <dup> -s <sex_chromosome>"
    echo "required options:"
        echo "-b|--bam      alignment file in BAM format using BWA" 
        echo "-h|--header   header of BAM file, the prefix of the output file (same with SM tag of the input BAM file)"
        echo "-d|--dup      duplicated window record file used for copy number correction"
        echo "-s|--sex      the name of sex chromosome"
    1>&2
    exit 1
}

OPTS=$(getopt -o b:h:d:s: --long bam:,header:,dup:,sex: -- "$@")
eval set -- "$OPTS"
while true ; do
    case "$1" in
        -b|--bam )
            bam=$2 ; shift 2
            ;;
        -h|--header )
            header=$2 ; shift 2
            ;;
        -d|--dup )
            dup=$2; shift 2
            ;;
        -s|--sex )
            sex=$2; shift 2
            ;;
        --)
            shift ; break
            ;;
        *)
            echo "Option error!"
            exit 1 ;;
    esac
done

if [ -z "$bam" ] || [ -z "$header" ] || [ -z "$dup" ] || [ -z "$sex" ]; then usage; exit 1; fi
refdb=$(ls `pwd`/referenceDB*)
windowsize=$(echo $refdb | grep -oP "referenceDB.\d+" | grep -oP "\d+")
stepsize=$(echo $refdb | grep -oP "_\d+" | grep -oP "\d+")
echo "refdb $refdb"
echo "bam $bam"
echo "header $header"
echo "dup $dup"
echo "sex $sex"
echo "window size $windowsize"
echo "step size $stepsize"
mkdir -p RD_raw
cd RD_raw || exit 1
python $CNVcaller/bin/1.1.CNVprocess.py -b $bam -r $refdb -w $windowsize -s $stepsize
echo "raw reads count finished."
cd ..
mkdir -p RD_absolute
python $CNVcaller/bin/1.2.CNVcollect.py -raw RD_raw/${header}_raw -link $dup -o RD_absolute/$header
echo "absolute correct finished."
mkdir -p RD_normalized
cd RD_normalized || exit 1
python $CNVcaller/bin/1.3.CNVnormalize.py -i ../RD_absolute/$header -w $windowsize -s $sex 
cd ..
echo "normalization finished!"