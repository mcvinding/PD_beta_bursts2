# Get cortical thickness of ROI (created in petmeg_get_betatrace.py script).

startdir=$(pwd)
source /home/mikkel/PD_long/scripts/setUpFreesurfer.sh
cd $SUBJECTS_DIR

SUBJECTS=('0525' '0551' '0552' '0560' '0565' '0569' '0584' '0586' '0588')

for SUB in ${SUBJECTS[@]}; do
	echo "Processing subject $SUB"

	echo Label file: $SUB/label/lh.sensmotor.label
	echo Output file:  /PETMEG/fs_stats/${SUB:0:4}.lh.sensmotor.stats
	mris_anatomical_stats \
		-l $SUB/label/lh.sensmotor.label \
		-f ~/PETMEG/fs_stats/${SUB:0:4}.lh.sensmotor.stats \
		$SUB \
		lh

	echo Label file: $SUB/label/rh.sensmotor.label
	echo Output file:  /PETMEG/fs_stats/${SUB:0:4}.rh.sensmotor.stats
	mris_anatomical_stats \
		-l $SUB/label/rh.sensmotor.label \
		-f ~/PETMEG/fs_stats/${SUB:0:4}.rh.sensmotor.stats \
		$SUB \
		rh
		
done

cd $startdir
