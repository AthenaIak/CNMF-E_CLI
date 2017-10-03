#!/bin/sh
echo ${BASH_VERSION}
# system specific locations
FIJIDIR="/home/athina/Downloads/Fiji.app"
MATLABDIR="/usr/local/MATLAB/R2015b/bin"
CURRDIR=$PWD

# command line arguments
if [ $# -ne 6 ]; then
	echo "Illegal number of parameters!"
	echo "Use bash ./cnmf_script.sh [input directory] [mouseID] [SessionID] [Trial number] [trial/rest] [analysis tag]"
	exit 1
fi

INDIR=$1
if [ $5 == "trial" ]; then 
	isTrial=1
elif [ $5 == "rest" ]; then 
	isTrial=0
else
	echo "Wrong input for the 5th argument. It should be trial or rest"
	exit 1
fi
RECID=`(grep "$2,$3,$4,$isTrial" $INDIR/rec_ids.dat) | awk -F"," '{print $5}'`
TAG=$6

# find out how many files comprise this movie
cmdOutput=`find $INDIR -name recording_$RECID*.tif`
numFiles=0
for entry in $cmdOutput
do
	numFiles=`expr $numFiles + 1`
done

# construct the string necessary to create the corresponding matlab cell of movie files
unCompressedFiles=\'$INDIR/recording_$RECID.tif\'
for i in `seq 2 $numFiles`
do
	unCompressedFiles="$unCompressedFiles, '$INDIR/recording_$RECID-00`expr $i - 1`.tif'"
done
unCompressedFiles="{$unCompressedFiles}"
#echo $unCompressedFiles

# temporally downsample the movies using matlab (generates the pp files)
cd $MATLABDIR
echo matlab -nosplash -nodesktop -r "cd('${CURRDIR}');cnmfe_setup;movieFiles=${unCompressedFiles};downsampling(movieFiles,4);exit();"

# find out how many preprocessed (pp) files were generated
cmdOutput=`find $INDIR -name pp_recording_$RECID*.tif`
numPPfiles=0
for entry in $cmdOutput
do
	numPPfiles=`expr $numPPfiles + 1`
done

# motion correct using ImageJ
#~/Downloads/Fiji.app$ ./ImageJ-linux64 -macro ~/Data/Macro.ijm ~/Data+test_TRP+2
cd $FIJIDIR
# ./ImageJ-linux64 -macro ${CURRDIR}/Macro.ijm ${INDIR}+${RECID}+${numPPfiles}
numMCfiles=$numPPfiles

# construct the string necessary to create the corresponding matlab cell of motion corrected movie files
mcFiles=\'$INDIR/mc_recording_$RECID-1.tif\'
for i in `seq 2 $numMCfiles`
do
	mcFiles="$mcFiles, '$INDIR/recording_$RECID-$i.tif'"
done
mcFiles="{$mcFiles}"
#echo $mcFiles

# use matlab to run cnmf-e
echo matlab -nosplash -nodesktop -r "cd('${CURRDIR}');movieFiles=${mcFiles};run_analysis(movieFiles,'${TAG}');exit();"

# check if results file was successfully created. if so, delete all tiff (and mat?) files with this RECID

#cd $CURRDIR # not required at the end of program (automatically returns to previous environment location)

