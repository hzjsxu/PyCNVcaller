#!/bin/bash
#---------------------------------------------------------------------------
# File: CNV.discovery.sh
# Created Date: 2020-08-30
# Author: jsxu
# Contact: <hzaujsxu@163.com>
# Last Modified: Sunday August 30th 2020 7:39:47 pm
# modified based on CNVcaller
#---------------------------------------------------------------------------

# CNVcaller installation directory.
export CNVcaller=/store/jsxu/biosoft/CNVcaller
echo "CNVcaller install directory $CNVcaller."

if [ ! -f "$CNVcaller/bin/2.1.CNVDiscoveryMerge.py" ]; then
    echo "You should set CNVcaller installation directory."
    exit 1
fi
#####################################
usage() {
    echo "Usage:  $0 -l <normalizedFileList> -e <excludedFileList> -f <frequency> -h <homozygous> -r <pearsonCorrelation> -p <primaryCNVR> -m <mergedCNVR>"
    echo "required options:"
        
        echo "-l|--normalizedFileList  individual normalized copy number files list, with absolute path"

        echo "-e|--excludedFileList    the samples in this list will be excluede from CNVR detection,"
        echo "                         their genotyping are reported based on the CNVR boundaries defined by other samples."
        echo "                         This option is applicable to the outgroup individual"

        echo "-f|--frequency           minimum frequency of gain/loss individuals when define a candidate CNV window"
        echo "                         [recommend 0.1]"

        echo "-h|--homozygous          minimum number of homozygous gain/loss individuals when define a candidate CNV window"
        echo "                         [recommend 3]"

        echo "-r|--pearsonCorrelation  minimum of pearson correlation coefficient between the two adjacent non-overlapping windows"
        echo "                         during CNVR discovery"
        echo "                         [recommend:]"
        echo "                         0.5 for sample size (0, 30]"
        echo "                         0.4 for sample size (30, 50]"
        echo "                         0.3 for sample size (50, 100]"
        echo "                         0.2 for sample size (100, 200]"
        echo "                         0.15 for sample size (200, 500]"
        echo "                         0.1 for sample size (500,+∞)"
        echo "-p|--primaryCNVR         primary CNVR result"

        echo "-m|--mergedCNVR          merged CNVR result"
    1>&2
    exit 1
}

##################
OPTS=$(getopt -o l:e:f:h:r:p:m: --long normalizedFileList:,excludedFileList:,frequency:,homozygous:,pearsonCorrelation:,primaryCNVR:,mergedCNVR: -- "$@")
if [ $? != 0 ]; then usage; exit 1; fi
eval set -- "$OPTS"
while true ; do
    case "$1" in
        -l|--normalizedFileList )
            normalizedFileList=$2 ; shift 2
            ;;
        -e|--excludedFileList )
            excludedFileList=$2 ; shift 2
            ;;
        -f|--frequency )
            frequency=$2 ; shift 2
            ;;
        -h|--homozygous )
            homozygous=$2 ; shift 2
            ;;
        -r|pearsonCorrelation )
            pearsonCorrelation=$2 ; shift 2
            ;;
        -p|--primaryCNVR )      
            primaryCNVR=$2 ; shift 2
            ;;
        -m|--mergedCNVR )      
            mergedCNVR=$2 ; shift 2
            ;;
        --) 
            shift ; break
            ;;
        *) 
            echo "$1 Internal error!" 
            exit 1
            ;;
    esac
done

if [ -z $normalizedFileList ] || [ -z $excludedFileList ] || [ -z $frequency ] || [ -z $homozygous ] || [ -z $pearsonCorrelation ] || [ -z $primaryCNVR ] || [ -z $mergedCNVR ]
    then usage
    exit 1
fi

refdb=$(ls `pwd`/referenceDB*)
windowsize=$(echo $refdb | grep -oP "referenceDB.\d+" | grep -oP "\d+")
stepsize=$(echo $refdb | grep -oP "_\d+" | grep -oP "\d+")
echo "referenceDB $refdb";
echo "normalizedFileList $normalizedFileList";
echo "excludedFileList $excludedFileList";
echo "windowsize $windowsize";
echo "stepsize $stepsize"
echo "frequency $frequency";
echo "homozygous $homozygous";
echo "pearsonCorrelation $pearsonCorrelation";
echo "primaryCNVR $primaryCNVR";
echo "mergedCNVR $mergedCNVR";
python $CNVcaller/bin/2.1.CNVDiscoveryMerge.py $refdb $normalizedFileList $excludedFileList $primaryCNVR -w $windowsize -f $frequency -h $homozygous -r $pearsonCorr    elation
echo "primary CNVR define finished!"

python $CNVcaller/bin/2.2.CNVRRedundancy.py -i $primaryCNVR -cor $pearsonCorrelation -o $mergedCNVR
echo "merge CNVR finished!"