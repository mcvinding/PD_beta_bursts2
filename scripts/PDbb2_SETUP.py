#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Setup and metavariables for PDbb2. Created on Thu May 14 12:02:50 2020
@author: mikkel
"""
import csv
import os.path as op

#%% Paths
raw_path = '/archive/20079_parkinsons_longitudinal/MEG/'
meg_path = '/home/mikkel/PD_longrest/meg_data'


#%% Read subjects
subj_file = '/home/mikkel/PD_long/subj_data/subjects_and_dates.csv'

with open(subj_file, newline='') as csvfile:
    tmp = csv.reader(csvfile, delimiter=',', quotechar='"')
    subjid = []
    date = []
    is_bb = []
    for ii, row in enumerate(tmp):
        subjid.append(row[2])
        date.append(row[3])
        is_bb.append(row[4])
        
is_bb = [x == '1' for x in is_bb]
tt = [x == '' for x in date]
idxer = [a and not b for a, b in zip(is_bb, tt)]
    
subjid_flt = [i for (i, v) in zip(subjid, idxer) if v]
date_flt = [i for (i, v) in zip(date, idxer) if v]

subjects_and_dates = ['NatMEG_'+s+'/'+d for (s, d) in zip(subjid_flt, date_flt)]

#%% Filenames
filestring          = 'rest_ec_mc_avgtrans_tsss_corr95'     # Unique string for finding raw files

exceptions = {
    '0582':'rest_ec_tsss.fif',              # No cHPI
    '0604':'rest_ec_tsss_mc.fif'            # Error due to too many 'autobad' channels. Manual MaxFilter.
    }



#END