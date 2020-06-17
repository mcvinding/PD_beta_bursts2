# Get Freesurfer summary statistics of ROI.
# @mc_vinding

startdir=$(pwd)
source /home/mikkel/PD_long/scripts/setUpFreesurfer.sh
cd $SUBJECTS_DIR
# Define dirs
OUTDIR=/home/mikkel/PD_longrest/meg_data

#######################################################################################
# Define subjects
SUBJECTS=( $(awk -F "\"*;\"*" '{print $2}' ~/PD_long/subj_data/subjects_and_dates.csv) ); unset 'SUBJECTS[0]' 		# REMOVE header
INCLUDE=( $(awk -F "\"*;\"*" '{print $4}' ~/PD_long/subj_data/subjects_and_dates.csv) ); unset 'INCLUDE[0]' 		# REMOVE header

# Remove non-included
DEL=XXXX

for i in "${!INCLUDE[@]}"; do
   if [ ${INCLUDE[$i]} == "0" ]; then
   	SUBJECTS[$i]=$DEL
	fi
done

SUBJECTS=( "${SUBJECTS[@]/$DEL}" )

#######################################################################################
# RUN
for SUB in ${SUBJECTS[@]}; do
	echo "Processing subject $SUB"
	SUBDIR=$OUTDIR/$SUB
	OUTPUTLH=$SUBDIR/$SUB.lh.sensmotor.stats

	if [ -f $OUTPUTLH ];	then
	echo "Output $OUTPUTLH exists. Del to run again!"
	continue
	fi
	
	# Left hemi
	echo Label file: $SUBDIR/label/lh.sensmotor.label
	echo Output file: $OUTPUTLH
	mris_anatomical_stats \
		-l $SUB/label/lh.sensmotor.label \
		-f $OUTPUTLH \
		$SUB \
		lh

	# Right hemi
	OUTPUTRH=$SUBDIR/$SUB.rh.sensmotor.stats
	
	echo Label file: $SUB/label/rh.sensmotor.label
	echo Output file: $OUTPUTRH
	mris_anatomical_stats \
		-l $SUB/label/rh.sensmotor.label \
		-f $OUTPUTRH \
		$SUB \
		rh
		
done

cd $startdir

#END
