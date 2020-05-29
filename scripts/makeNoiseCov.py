#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Make noise convariance matrix from empty room data
* Need to only use first two minutes of raw files
@author: mikkel
"""
import os.path as op
from os import listdir
from mne.io import Raw
from mne import compute_raw_covariance, write_cov
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/mikkel/PD_longrest/scripts')
from PDbb2_SETUP import subjects_and_dates, meg_path, raw_path, old_raw_path, old_subjs, empt_filestring, old_empt_filestring

#%% Run options
overwrite = False      # Wheter files should be overwritten if already exist

# Missing list for debugging/diagnostics
missing_list = []

# %% Run through files
for subj_date in subjects_and_dates:   
    
    subj = subj_date.split('/')[0][-4:]
    print('Now processing subject '+subj)
    
    if subj in old_subjs:
        raw_fpath = op.join(old_raw_path, subj_date)
    else:
        raw_fpath = op.join(raw_path, subj_date)

    fig_path = op.join(meg_path, subj, 'ica')
    outFname = op.join(meg_path, subj, subj+'-cov.fif')         #Covfefe
        
    if op.exists(outFname) and not overwrite:
        print('Do not overwrite cov file')
        continue

    # Find files
    file_list = listdir(raw_fpath)

    if subj in old_subjs:
        inFiles = [op.join(raw_fpath,f) for f in file_list if old_empt_filestring in f]
    else:
        inFiles = [op.join(raw_fpath,f) for f in file_list if empt_filestring in f]
    inFiles.sort()
    
    if not inFiles:
        print('WARNING: NO FILE FOR SUBJ '+subj)
        missing_list += [subj]
        continue
    
    raw_temp = Raw(inFiles[0], preload=True)
    
    # raw_temp.crop(tmin=0, tmax=120)

    # Estimate cov
    noise_cov = compute_raw_covariance(raw_temp, tmin=0, tmax=120)
    
    # Plot for inspection
    cov_fig = noise_cov.plot(raw_temp.info)
    [fig.savefig(op.join(fig_path,'cov'+str(i)+'.png')) for i, fig in enumerate(cov_fig)]
    plt.close('all')
    
    # Save
    write_cov(outFname, noise_cov)
    
    print("Done with "+subj)