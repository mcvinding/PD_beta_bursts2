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
import pandas as pd

#%% Paths
if sys.platform == 'linux':
    raw_path        = '/archive/20079_parkinsons_longitudinal/MEG/'
    old_raw_path    = '/archive/20055_parkinson_motor/MEG'
    meg_path        = '/home/mikkel/PD_longrest/meg_data'
    fs_subjects_dir = '/home/mikkel/PD_long/fs_subjects_dir'
    subj_data_path  = '/home/mikkel/PD_long/subj_data/'
    grp_data_path   = '/home/mikkel/PD_longrest/groupanalysis'
else:
    meg_path        = 'X:\\PD_longrest\\script\\meg_data'
    fs_subjects_dir = 'X:\\PD_long\\fs_subjects_dir'
    subj_data_path  = 'X:\\PD_long\\subj_data'

#%% Read subjects
subj_file = op.join(subj_data_path, 'subjects_and_dates.csv')

tmp = pd.read_csv(subj_file)      
is_bb = [x == 1 for x in tmp['rest_ec']]
tmp['id'] = ['0'+str(s) for s in tmp['id']]
subjects_and_dates = [op.join('NatMEG_'+s, str(d)) for (s, d) in zip(tmp['id'][is_bb], tmp['date'][is_bb])]
subjects = list(tmp['id'][is_bb])

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
tmp2 = pd.read_csv(data_file)
group = tmp2['group']

#%% Settings
spacing = 'ico4'

#END