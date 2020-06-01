# Get cortical thickness of ROI (created in petmeg_get_betatrace.py script).

startdir=$(pwd)
source /home/mikkel/PD_long/scripts/setUpFreesurfer.sh
cd $SUBJECTS_DIR

SUBJECTS=('0406' '0525' '0551' '0552' '0560' '0565' '0569' '0584' '0586' '0588')

OUTDIR=/home/mikkel/PD_longrest/meg_data

for SUB in ${SUBJECTS[@]}; do
	echo "Processing subject $SUB"
	SUBDIR=$OUTDIR/$SUB
	OUTPUTLH=$SUBDIR/$SUB.lh.sensmotor.stats

	echo Label file: $SUB/label/lh.sensmotor.label
	echo Output file: $OUTPUTLH
	mris_anatomical_stats \
		-l $SUB/label/lh.sensmotor.label \
		-f $OUTPUTLH \
		$SUB \
		lh

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
