#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Run ICA and remove components correlted with EOG and ECG. Save figures in subj 
folder for inspection.
@author: mikkel
"""
import matplotlib.pyplot as plt
from mne.preprocessing import create_ecg_epochs, create_eog_epochs, ICA
from mne.io import read_raw_fif
from mne import pick_types
from os import listdir, mkdir
import os.path as op
import numpy as np
import sys
sys.path.append('/home/mikkel/PD_longrest/scripts')
from PDbb2_SETUP import subjects_and_dates, raw_path, meg_path, exceptions, filestring

# %% run options
overwrite_old_files = False                                     # Wheter files should be overwritten if already exist

# bandpass filter
Filter = [None, 48]

# Number of components to reject
n_max_ecg = 3
n_max_eog = 2

# extreme jump rejection
reject = dict(grad=4000e-13,    # T / m (gradiometers)
              mag=5e-12,        # T (magnetometers)
              )  
    
#%% RUN
missing_list = []
for subj_date in subjects_and_dates:   

    # Find all files for same subject to run ICA on. The first step is to accomodate
    # files that have been split in two or more files.
    subj = subj_date.split('/')[0][-4:]
    print('Now processing subject '+subj)
    
    # Define paths (this outght to be done in the config file)
    raw_fpath = op.join(raw_path, subj_date)
    file_list = listdir(raw_fpath)
    if subj in exceptions:
        inFiles = [op.join(raw_fpath,f) for f in file_list if exceptions[subj] in f]
        inFiles.sort()
    else:
        inFiles = [op.join(raw_fpath,f) for f in file_list if filestring in f]
        inFiles.sort()
    
    if not inFiles:
        print('WARNING: NO FILE FOR SUBJ '+subj)
        missing_list += [subj]
        continue
    
    sub_path        = op.join(meg_path, subj)
    ica_path        = op.join(meg_path, subj, 'ica')
    out_icaFname    = op.join(ica_path, 'comp-ica.fif')            # Name of ICA component (saved for bookkeeping)
    outfname        = op.join(sub_path, subj+'-ica-raw.fif')
    
    if op.exists(outfname) and not overwrite_old_files:
        print('Output '+outfname+' already exists')
        continue
    
    # Make dirs if they do not exist
    if not op.exists(sub_path):
        mkdir(sub_path)
    if not op.exists(ica_path):
        mkdir(ica_path)      
    
    # Load data (This loop will make sure that split files are read toghether)
    for i, fname in enumerate(inFiles):
        print('loading '+str([f for f in inFiles]))
        if i < 1:
            raw = read_raw_fif(fname, preload=True)
        else:
            raw.append(read_raw_fif(fname, preload=True))
            
    # mark channels as bad if necessary - Maybe bad channels can be written to a dict and loaded.
#        %raw.info['bads'] = ['MEG2332']                         # change for subjects after summer break (08/19)
#        print (raw.info['bads'])
     
    # Filter data
    print('Filtering....')
    picks_meg = pick_types(raw.info, meg=True, eeg=False, eog=False, emg=False, misc=False, 
                           stim=False, exclude='bads')
    raw_filt = raw.copy().filter(Filter[0], Filter[1], n_jobs=3, picks=picks_meg)
    raw.filter(None, Filter[1], n_jobs=3, picks=picks_meg)
    
    # RUN ICA
    ica = ICA(n_components=0.95, method='fastica', random_state=0)

    ica.fit(raw_filt, picks=picks_meg, decim=3, reject=reject, verbose=True)
#    ica.labels_ = dict()        
    
    # Plot and save
    ica_fig = ica.plot_components()
    [fig.savefig(op.join(ica_path,'ICA_allComp'+str(i)+'.png')) for i, fig in enumerate(ica_fig)]
    
    ica.save(out_icaFname)
    print('ICA comp saved as '+out_icaFname)

    # REMOVE COMPONENTS
    picks_eXg = pick_types(raw.info, meg=False, eeg=False, eog=True, ecg = True, emg=False, misc=False, stim=False, exclude='bads')
    raw.filter(Filter[0], Filter[1], n_jobs=3, picks=picks_eXg)
    raw.notch_filter(50, n_jobs=3, picks=picks_eXg)             # Remove residual 50Hz line noise
    
    # Find ECG artifacts
    ecg_epochs = create_ecg_epochs(raw, ch_name='ECG003', tmin=-.5, tmax=.5)    #, picks=picks)
    ecg_inds, ecg_scores = ica.find_bads_ecg(ecg_epochs, method='ctps', verbose=False)
    
    # Update reject info
    ica.exclude += ecg_inds[:n_max_ecg]

    # Plot ECG ICs for inspection
    ecg_scores_fig = ica.plot_scores(ecg_scores, exclude=ecg_inds, title='Component score (ecg)', show=True)
    ecg_scores_fig.savefig(op.join(ica_path,'ICA_ecg_comp_score.png'))
    plt.close()
    
    if ecg_inds:
        show_picks = np.abs(ecg_scores).argsort()[::-1][:5]
        
        ecg_comp_fig = ica.plot_components(ecg_inds, title='ecg comp', colorbar=True, show=False)
        ecg_comp_fig.savefig(op.join(ica_path,'ICA_ecg_comp_topo.png'))
        plt.close()
        
        # estimate average artifact
        ecg_evoked = ecg_epochs.average()
        
        # plot ECG sources + selection
        ecg_evo_fig1 = ica.plot_overlay(ecg_evoked, exclude=ecg_inds, show=False)    
        ecg_evo_fig1.savefig(op.join(ica_path, 'ICA_ecg_overlay.png'))
        plt.close()
    
    # Find EOG artifacts 
    eog_epochs = create_eog_epochs(raw, reject=reject, tmin=-.5, tmax=.5)  # get single EOG trials
    eog_inds, eog_scores = ica.find_bads_eog(raw)
    
    # Update reject info
    ica.exclude += eog_inds[:n_max_eog]

    # Plot EOG ICs for inspection
    eog_scores_fig = ica.plot_scores(eog_scores, exclude=eog_inds, title='Component score (ecg)', show=True)
    eog_scores_fig.savefig(op.join(ica_path, 'ICA_eog_comp_score.png'), show=False)
    plt.close()

    if eog_inds:
        eog_comp_fig = ica.plot_components(eog_inds, title="Sources related to EOG artifacts", colorbar=True,show=False)
        eog_comp_fig.savefig(op.join(ica_path, 'ICA_eog_comp.png'))
        plt.close()
                                  
        # plot EOG sources + selection
        eog_evo_fig = ica.plot_overlay(eog_epochs.average(), exclude=eog_inds, show=False)  # plot EOG cleaning
        eog_evo_fig.savefig(op.join(ica_path, 'ICA_eog_overlay.png'))
        plt.close()
            
    # Apply the solution to Raw, Epochs or Evoked like this:
    raw_ica = ica.apply(raw)
    raw_ica.save(outfname, overwrite=overwrite_old_files)

    print('----------- FINISHED '+subj+' -----------------')
    plt.close('all')

#END