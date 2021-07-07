#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Make noise convariance matrix from empty room data.

Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M., Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease. medRxiv.org. https://doi.org/10.1101/2021.06.27.21259592

@author: mcvinding

"""
import os.path as op
from os import listdir
from mne.io import Raw
from mne import compute_raw_covariance, write_cov, pick_types
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/mikkel/PD_longrest/scripts')
from PDbb2_SETUP import subjects_and_dates, meg_path, raw_path, old_raw_path, old_subjs, empt_filestring, old_empt_filestring

#%% Run options
overwrite = False      # Wheter files should be overwritten if already exist

# Filter options
bandpass_freqs = [None, 48]
notch_freqs = [50, 100, 150]

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

    fig_path = op.join(meg_path, subj, 'plots')
    outFname = op.join(meg_path, subj, subj+'-cov.fif')
    outTemp  = op.join(meg_path, subj, subj+'-empt-raw.fif')
        
    if op.exists(outFname) and not overwrite:
        print('Do not overwrite cov file')
        continue

    # Find files
    file_list = listdir(raw_fpath)

    # Errors in raw filenames
    if subj == '0581': 
        inFiles = [op.join(raw_fpath,f) for f in file_list if 'empty_room_after_tsss' in f][0] 
    elif subj in ['0319', '0320', '0398']:
        inFiles = [op.join(raw_fpath,f) for f in file_list if 'empty_room2_after_tsss' in f][0]
    elif subj == '0328':
        inFiles = '/archive/20055_parkinson_motor/MEG/NatMEG_0327/160901/empty_room1_after_tsss.fif'   # Same day
    elif subj in ['0342','0343']:        
        inFiles = [op.join(raw_fpath,f) for f in file_list if 'empty_room_12_after_tsss' in f][0] 
    elif subj in old_subjs:
        inFiles = [op.join(raw_fpath,f) for f in file_list if old_empt_filestring in f][0]
    else:
        inFiles = [op.join(raw_fpath,f) for f in file_list if empt_filestring in f][0]
    
    if not inFiles:
        print('WARNING: NO FILE FOR SUBJ '+subj)
        missing_list += [subj]
        continue
    
    raw_temp = Raw(inFiles, preload=True)
    
    # Filter data
    print('Filtering....')
    picks_meg = pick_types(raw_temp.info, meg=True, eeg=False, eog=False, emg=False, misc=False, 
                           stim=False, exclude='bads')
    raw_temp.notch_filter(notch_freqs, n_jobs=3, picks=picks_meg)
    raw_temp.filter(bandpass_freqs[0], bandpass_freqs[1], n_jobs=3, picks=picks_meg)

    # Estimate cov
    noise_cov = compute_raw_covariance(raw_temp, tmin=0, tmax=120, rank='info')

    # Plot for inspection
    cov_fig = noise_cov.plot(raw_temp.info)
    [fig.savefig(op.join(fig_path,'cov'+str(i)+'.png')) for i, fig in enumerate(cov_fig)]
    plt.close('all')
    
    # Save
    write_cov(outFname, noise_cov)
    raw_temp.save(outTemp, overwrite=overwrite)
    
    print("Done with "+subj)
