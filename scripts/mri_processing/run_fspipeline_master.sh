#!/bin/bash
# Run Freesurfer recon-all on subjects.

# Values
megID="$1"
T1seq="$2"

## Warnings
if [[ -z "$2" ]]; then
	echo "Warning: missing sequence number"
	exit 1
fi

## Setup Freesurfer
export FREESURFER_HOME=/opt/freesurfer
source /home/mikkel/PD_long/scripts/mri_processing/setUpFreesurfer.sh

## Paths
export SUBJECTS_DIR=/home/mikkel/PD_long/fs_subjects_dir
export TOOL_DIR=/home/mikkel/mri_scripts

# Run
$TOOL_DIR/NatMEG_recon $megID $T1seq $SUBJECTS_DIR

#END
