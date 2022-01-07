#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Setup and metavariables for PDbb2.

Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M., Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease. medRxiv.org. https://doi.org/10.1101/2021.06.27.21259592

@author: mcvinding
"""
import csv
import os.path as op
import sys

#%% Paths
if sys.platform == 'linux':
    raw_path        = '/archive/20079_parkinsons_longitudinal/MEG/'
    old_raw_path    = '/archive/20055_parkinson_motor/MEG'
    meg_path        = '/home/mikkel/PD_longrest/meg_data'
#    trans_path      = '/home/mikkel/PD_long/trans_files'           # !! Files are in subj_folder
    # old_trans_path  = '/home/mikkel/PD_motor/tap/trans_files'
    fs_subjects_dir = '/home/mikkel/PD_long/fs_subjects_dir'
    subj_data_path  = '/home/mikkel/PD_long/subj_data/'
    grp_data_path   = '/home/mikkel/PD_longrest/groupanalysis'
else:
#    raw_path        = '/archive/20079_parkinsons_longitudinal/MEG/'
#    old_raw_path    = '/archive/20055_parkinson_motor/MEG'
    meg_path        = 'X:\\PD_longrest\\script\\meg_data'
#    trans_path      = 'X:\\PD_longrest\\script\\trans_files'
    # old_trans_path  = '/home/mikkel/PD_motor/tap/trans_files'
    fs_subjects_dir = 'X:\\PD_long\\fs_subjects_dir'
    subj_data_path  = 'X:\\PD_long\\subj_data'

#%% Read subjects
subj_file = op.join(subj_data_path, 'subjects_and_dates.csv')

with open(subj_file, newline='') as csvfile:
    tmp = csv.reader(csvfile, delimiter=';', quotechar='"')
    subjid = []
    date = []
    is_bb = []
    for ii, row in enumerate(tmp):
        subjid.append(row[1])
        date.append(row[2])
        is_bb.append(row[3])
        
is_bb = [x == '1' for x in is_bb]
tt = [x == '' for x in date]
idxer = [a and not b for a, b in zip(is_bb, tt)]
    
subjid_flt = [i for (i, v) in zip(subjid, idxer) if v]
subjid_flt = ['0'+s for s in subjid_flt]

date_flt = [i for (i, v) in zip(date, idxer) if v]

subjects_and_dates = [op.join('NatMEG_'+s, d) for (s, d) in zip(subjid_flt, date_flt)]
subjects = subjid_flt.copy()

old_subjs = ['0313',
             '0314',
             '0319',
             '0320',
             '0322',
             '0325',
             '0328',
             '0332',
             '0333',
             '0339',
             '0340',
             '0342',
             '0343',
             '0352',
             '0353',
             '0355',
             '0366',
             '0377',
             '0392',
             '0397',
             '0398',
             '0406']

#%% Filenames
filestring          = 'rest_ec_mc_avgtrans_tsss_corr95'     # Unique string for finding raw files
old_filestring      = 'rest_ec_2_mc_trans_tsss'             # Unique string for finding raw files

empt_filestring     = 'empty_room_before_tsss'
old_empt_filestring = 'empty_room2_before_tsss'


exceptions = {
    '0582':'rest_ec_tsss.fif',              # No cHPI
    '0604':'rest_ec_tsss_mc.fif',           # Error due to too many 'autobad' channels. Manual MaxFilter.
    '0322':'rest_eo2_mc_trans_tsss',        # Wrong filename during recording
    '0352':'rest_ec_2_mc_trans_tsss_max95'  # Multiple MaxFilter files
    }


#%% Grouping
data_file = op.join(grp_data_path, 'alldata_subj2.csv')

with open(data_file, newline='') as csvfile:
    tmp = csv.reader(csvfile, delimiter=',', quotechar='"')
    group = []
    for ii, row in enumerate(tmp):
        group.append(row[10])

group = group[1:]

#%% Settings
spacing = 'ico4'

#END