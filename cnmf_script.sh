#!/bin/sh
echo ${BASH_VERSION}
# system specific locations
FIJIDIR="/home/athina/Downloads/Fiji.app"
MATLABDIR="/usr/local/MATLAB/R2015b/bin"
CURRDIR=$PWD

# command line arguments
if [ $# -ne 6 -a $# -ne 5 ]; then
	echo "Illegal number of parameters!"
	echo "USE:"
	echo "bash ./cnmf_script.sh [input directory] [mouseID] [SessionID] [Trial number] [trial/rest] [analysis tag]"
	echo "or"
	echo "bash ./cnmf_script.sh [input directory] [mouseID] [SessionID] [Day (1/2/3/4/5)] [analysis tag]"
	exit 1
fi

INDIR=$1
if [ $# -eq 6 ]; then
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
	cmdOutput=`find $INDIR -name recording_$RECID*.tif`
	# find out how many files comprise this movie
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
else
	let a="($4 - 1) * 5"
	RECID=`(grep "$2,$3,$(expr $a + 1),1" $INDIR/rec_ids.dat) | awk -F"," '{print $5}'`
	RECIDtmp=`(grep "$2,$3,$(expr $a + 1),0" $INDIR/rec_ids.dat) | awk -F"," '{print $5}'`
	TAG=$5
		
	cmdOutput=`find $INDIR -name recording_$RECID*.tif`
	numFiles=0
	for entry in $cmdOutput
		do
			numFiles=`expr $numFiles + 1`
		done
	if [ $numFiles -ne 0 ]; then
			unCompressedFiles=\'$INDIR/recording_$RECID.tif\'
	fi
	#normally we would exit if $numFiles -e 0
	for i in `seq 2 $numFiles`
		do
			unCompressedFiles="$unCompressedFiles, '$INDIR/recording_$RECID-00`expr $i - 1`.tif'"
		done
		
	cmdOutput=`find $INDIR -name recording_$RECIDtmp*.tif`
	numFiles=0
	for entry in $cmdOutput
		do
			numFiles=`expr $numFiles + 1`
		done
	if [ $numFiles -ne 0 ]; then
		unCompressedFiles="$unCompressedFiles, '$INDIR/recording_$RECIDtmp.tif'"
	fi
	# normally we would exit if $numFiles -e 0
	for i in `seq 2 $numFiles`
		do
			unCompressedFiles="$unCompressedFiles, '$INDIR/recording_$RECIDtmp-00`expr $i - 1`.tif'"
		done
	if [ $4 -le 4 ]; then 
		for n in `seq 2 5`
			do
				for isTrial in 1 0
					do
						RECIDtmp=`(grep "$2,$3,$(expr $a + $n),$isTrial" $INDIR/rec_ids.dat) | awk -F"," '{print $5}'`

						cmdOutput=`find $INDIR -name recording_$RECIDtmp*.tif`
						numFiles=0
						for entry in $cmdOutput
						do
							numFiles=`expr $numFiles + 1`
						done
						if [ $numFiles -ne 0 ]; then
							unCompressedFiles="$unCompressedFiles, '$INDIR/recording_$RECIDtmp.tif'"
						fi

						# normally we would exit if $numFiles -e 0
						for i in `seq 2 $numFiles`
						do
							unCompressedFiles="$unCompressedFiles, '$INDIR/recording_$RECIDtmp-00`expr $i - 1`.tif'"
						done
				done
		done
		
	fi
fi

unCompressedFiles="{$unCompressedFiles}"
#echo $unCompressedFiles

# temporally downsample the movies using matlab (generates the pp files)
cd $MATLABDIR
matlab -nosplash -nodesktop -r "cd('${CURRDIR}');cnmfe_setup;movieFiles=${unCompressedFiles};downsampling(movieFiles,4);exit();"

# find out how many preprocessed (pp) files were generated
cmdOutput=`find $INDIR -name pp_recording_$RECID*.tif  | grep "\-[1-9]"`
numPPfiles=0
for entry in $cmdOutput
do
	numPPfiles=`expr $numPPfiles + 1`
done

# motion correct using ImageJ
cd $FIJIDIR
./ImageJ-linux64 -macro ${CURRDIR}/Macro.ijm ${INDIR}+${RECID}+${numPPfiles}
numMCfiles=$numPPfiles

# construct the string necessary to create the corresponding matlab cell of motion corrected movie files
mcFiles=\'$INDIR/mc_recording_$RECID-1.tif\'
for i in `seq 2 $numMCfiles`
do
	mcFiles="$mcFiles, '$INDIR/mc_recording_$RECID-$i.tif'"
done
mcFiles="{$mcFiles}"
#echo $mcFiles

# use matlab to run cnmf-e
matlab -nosplash -nodesktop -r "cd('${CURRDIR}');movieFiles=${mcFiles};run_analysis(movieFiles,'${TAG}');exit();"

# check if results file was successfully created. if so, delete all tiff (and mat?) files with this RECID

#cd $CURRDIR # not required at the end of program (automatically returns to previous environment location)

#matlab -nosplash -nodesktop -r "cd('${CURRDIR}');cnmfe_setup;path='${INDIR}/$RECID-${TAG}';view_analysis(path,'neurons');view_analysis(path,'ncomp');view_analysis(path,'nhollow');exit();"

