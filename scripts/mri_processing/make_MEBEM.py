#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Sort multiecho DICOM for all subjects. This use the multiecho BEM from 
https://github.com/mcvinding/mri_stuff/tree/master/multiecho_bem
"""
# Import
import os.path as op
import sys
import subprocess
sys.path.append('/home/mikkel/mri_scripts/multiecho_bem')   # Path to folder with scripts
from multiechoBEM_funs import sort_MEdicom, run_MEBEM, copyBEM2folder

#%% Options
# Where to put files
subjects_dir = '/home/mikkel/PD_long/fs_subjects_dir'       # Path to Freesurfer SUBJECTS_DIR
dicom_folder = '/home/mikkel/PD_long/MRI'                   # Intermediate path for sorted DICOM. Will create a subfolder with name <subj>
raw_folder = '/archive/20079_parkinsons_longitudinal/MRI'

# Input filenames (subject name and name of folder with raw DICOM)
subs_and_folders = {
                    '0522':'00000004',
                    '0523':'00000004',
                    '0524':'00000004',
                    '0525':'00000004',
                    '0528':'00000004',
                    '0529':'00000004',
                    '0530':'00000004',                    
                    '0531':'00000004',
                    '0532':'00000004',
                    '0533':'00000004',
                    '0534':'00000004',                    
                    '0535':'00000004',
                    '0536':'00000004',
                    '0537':'00000004',
                    '0540':'00000004',
                    '0541':'00000004',
                    '0542':'00000004',
                    '0543':'00000006',
                    '0544':'00000004',
                    '0545':'00000004',
                    '0546':'00000004',
                    '0547':'00000004',
                    '0548':'00000004',
                    '0549':'00000004',
                    '0550':'00000004',
                    '0552':'00000004',
                    '0553':'00000004',
                    '0554':'00000004',
                    '0556':'00000004',
                    '0557':'00000004',
                    '0559':'00000004',
                    '0561':'00000004',
                    '0562':'00000004',
                    '0563':'00000004',
                    '0565':'00000004',
                    '0566':'00000004',      
                    '0567':'00000004',                
                    '0568':'00000004',
                    '0570':'00000004',
                    '0571':'00000004',
                    '0572':'00000004',
                    '0573':'00000004',
                    '0574':'00000004',
                    '0575':'00000004',
                    '0576':'00000004',
                    '0578':'00000004',                    
                    '0579':'00000004',
                    '0580':'00000004',
                    '0581':'00000004',
                    '0582':'00000004',
                    '0583':'00000004',
                    '0584':'00000004',
                    '0585':'00000004',                    
                    '0587':'00000004',
                    '0589':'00000004',
                    '0590':'00000004',
                    '0591':'00000004',
                    '0592':'00000004',
                    '0593':'00000004',
                    '0594':'00000004',                
                    '0596':'00000004',
                    '0597':'00000004',
                    '0598':'00000004',
                    '0599':'00000004',
                    '0600':'00000004',
                    '0601':'00000004',
                    '0602':'00000004',
#                    '0603_12917':'00000004',
                    '0603':'00000004',                    
                    '0604':'00000004',  # Inspect
                    '0605':'00000004',
                    '0606':'00000004',
                    '0607':'00000004',
                    '0608':'00000004',
                    '0610':'00000004',
                    '0611':'00000004',
                    '0612':'00000004',
                    '0614':'00000004',
                    '0615':'00000005',
                    '0616':'00000004',
                    '0617':'00000004',
                    '0618':'00000004',
                    '0619':'00000400',
                    '0620':'00000400',
                    '0621':'00000004',
                    '0622':'00000004',
                    '0623':'00000004',
                    '0624':'00000004',
                    '0625':'00000004',
                    '0626':'00000004',
                    # '0627':'00000004',  # Folde rname on Archive is wrong. Do manually.
                    '0628':'00000003',  # ! Correct number?
                    '0629':'00000004',
                    '0630':'00000400',
                    '0632':'00000004',
                    '0633':'00000004',
                    '0634':'00000004',
                    '0635':'00000004',
                    '0636':'00000004',
                    '0637':'00000400',
                    '0638':'00000400',
                    '0640':'00000400',
                    '0641':'00000400',
                    '0642':'00000400',
                    '0643':'00000400',
                    '0645':'00000400',
                    '0647':'00000400',
                    '0648':'00000400',
#                    '0649_13297':'00000400',
                    '0649':'00000400',
                    '0650':'00000004'
                    }

# Append full filepath
subs_and_folders = {k:v for (k, v) in zip(subs_and_folders.keys(), (op.join(raw_folder,'NatMEG_'+f,subs_and_folders[f]) for f in subs_and_folders.keys()))}

# Generic output filenames (Intermediate path with sorted DICOM)
subs_and_dicompaths = {k:v for (k, v) in zip(subs_and_folders.keys(), (op.join(dicom_folder,f) for f in subs_and_folders.keys()))}

#%% RUN
replace = False

for sub in subs_and_folders:
    print('Processing subj: '+sub)
    
    if op.exists(op.join(subjects_dir, sub, 'bem', 'inner_skull.surf')):
        print("Files appear to already exist for sub "+sub+". Delete to run again")
        continue
        
#    #Sort DICOMS
#    print('##### SORTING ME DICOMS #####')
#    sort_MEdicom(subs_and_folders[sub], subs_and_dicompaths[sub], replace=replace)
    
    # Run ME BEM
    print('##### MAKING BEM #####')
    run_MEBEM(sub, subs_and_dicompaths[sub], subjects_dir)
    
    # Copy BEM surfaces to folder for further processing in MNE-PY
    print('##### COPY BEM FILES #####')
    copyBEM2folder(sub, subjects_dir, replace=True)
    
    # Add a high-resulution headmodel for better coalignment
    print('##### MAKING SCALP SURF #####')
    cmd = 'mne make_scalp_surfaces -s '+sub+' -d '+subjects_dir+' --force'        #Command to execute
    val = subprocess.call(cmd, shell=True)
        
    print('Done subj: '+sub)

#END