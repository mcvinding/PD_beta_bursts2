#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Run ICA and remove components correlted with EOG and ECG. Save figures in subj 
folder for inspection.
@author: mikkel
"""
import matplotlib.pyplot as plt
#import mne
from mne.preprocessing import create_ecg_epochs, create_eog_epochs, ICA, annotate_muscle_zscore
from mne.io import read_raw_fif, RawArray
from mne import pick_types, find_events, make_fixed_length_events, Epochs
from os import listdir, mkdir
import os.path as op
import numpy as np
import pandas as pd
import sys
sys.path.append('/home/mikkel/PD_longrest/scripts')
from PDbb2_SETUP import subjects_and_dates, raw_path, meg_path, exceptions, filestring, old_subjs, old_raw_path, old_filestring

# %% run options
overwrite_old_files = False     # Wheter files should be overwritten if already exist

# Muscle artefact detection
threshold_muscle = 5  # z-score

# bandpass filter
bandpass_freqs = [None, 48]
notch_freqs = [50, 100, 150]

# Number of components to reject
n_max_ecg = 3
n_max_eog = 2

# extreme jump rejection
reject = dict(grad=4000e-13,    # T / m (gradiometers)
              mag=5e-12,        # T (magnetometers)
              )  

startTrigger = 1
stopTrigger  = 64

# Logs for data cleaning reporting (merge to DataFrame later)
log_subj = []   # Subj id
log_len = []    # Length of cleanded segment
log_drop = []   # MNE drop log
log_ica = []    # N ICA components removed

#%% RUN
for subj_date in subjects_and_dates:   

    subj = subj_date.split('/')[0][-4:]
    log_subj += [subj]
    print('Now processing subject '+subj)
    
    # Define paths and filenames (this outght to be done in the config file)
    if subj in old_subjs:
        raw_fpath   = op.join(old_raw_path, subj_date)
    else:
        raw_fpath   = op.join(raw_path, subj_date)

    sub_path        = op.join(meg_path, subj)
    ica_path        = op.join(sub_path, 'ica')
    fig_path        = op.join(sub_path, 'plots')
    out_icaFname    = op.join(ica_path, 'comp-ica2.fif')         # Name of ICA component (saved for bookkeeping)
    outfname        = op.join(sub_path, subj+'-ica-raw2.fif')
    
    if op.exists(outfname) and not overwrite_old_files:
        print('Output '+outfname+' already exists')
        continue
    
    # Make dirs if they do not exist
    if not op.exists(sub_path):
        mkdir(sub_path)
    if not op.exists(ica_path):
        mkdir(ica_path)      
    if not op.exists(fig_path):
        mkdir(fig_path)      

    # Find all files for same subject to run ICA on. The first step is to accomodate
    # files that have been split in two or more files.
    file_list = listdir(raw_fpath)
    if subj in exceptions:
        inFiles = [op.join(raw_fpath,f) for f in file_list if exceptions[subj] in f]
    elif subj in old_subjs:
        inFiles = [op.join(raw_fpath,f) for f in file_list if old_filestring in f]
    else:
        inFiles = [op.join(raw_fpath,f) for f in file_list if filestring in f]
  
    # Load data
    fname = inFiles[0]
    raw = read_raw_fif(fname, preload=True)
                    
    # Initial rejection of segments with muscle movements
    annot_muscle, scores_muscle = annotate_muscle_zscore(raw, ch_type="mag", threshold=threshold_muscle, min_length_good=0.2)
    
    # Plot as save for inspection
    fig, ax = plt.subplots()
    ax.plot(raw.times, scores_muscle)
    ax.axhline(y=threshold_muscle, color='r')
    ax.set(xlabel='time, (s)', ylabel='zscore', title='Muscle activity (Z-score)')
    fig.savefig(op.join(fig_path,'muscle_artefact.png'))
    plt.close()
    
    raw.set_annotations(annot_muscle)
    
    # Inspect
    # raw.plot()
        
    # Filter data
    print('Filtering....')
    picks_meg = pick_types(raw.info, meg=True, eeg=False, eog=False, emg=False, misc=False, 
                           stim=False, exclude='bads')
    raw.notch_filter(notch_freqs, n_jobs=3, picks=picks_meg)                    # Remove residual linenoise
    raw.filter(bandpass_freqs[0], bandpass_freqs[1], n_jobs=3, picks=picks_meg)

    # Find events and crop data
    eve = find_events(raw, stim_channel='STI101')
    
    # Inspect
    # raw.plot(eve)
    
    # Trigger exceptions
    if subj == '0333':                    # Missing triggers
        startSam = 20000+raw.first_samp
        stopSam = 180102+startSam       
    elif subj == '0523':
        startSam = eve[eve[:,2] == startTrigger,0][0]
        stopSam = startSam+180102
    elif subj == '0524':
        startSam = eve[eve[:,2] == startTrigger,0][0]
        stopSam = startSam+180102        
    elif subj == '0529':                 # Missing stop trigger
        startSam = eve[eve[:,2] == startTrigger,0][0]
        stopSam = 180102+startSam    
    elif subj == '0548':                 # No start trigger. Start trigger val is stop trigger.
        stopSam = eve[eve[:,2] == startTrigger,0][0]
        startSam = stopSam-180102               
    elif subj == '0583':                 # No start trigger. Start trigger val is stop trigger.
        stopSam = eve[eve[:,2] == startTrigger,0][0]
        startSam = stopSam-180102                       
    elif subj == '0590':                 # No start trigger.
        startSam = raw.first_samp+1000             
        stopSam = startSam+180102
    elif subj == '0605':                 # Triggers numbers
        startSam = eve[eve[:,2] == 14593,0][0]
        stopSam = eve[eve[:,2] == 14656,0][0]
    elif subj == '0615':                 # No triggers.
        startSam = raw.first_samp+10000             
        stopSam = startSam+180102        
    else:
        # if not len(eve) == 2:
            # trigger_err += [subj]            
            # continue
        startSam = eve[eve[:,2] == startTrigger,0][0]
        stopSam = eve[eve[:,2] == stopTrigger,0][0]
    
    raw.crop(tmin=(startSam - raw.first_samp ) / raw.info['sfreq'], tmax=(stopSam - raw.first_samp ) / raw.info['sfreq'])

    # PSD for diagnostics
    fig = raw.plot_psd(tmax=np.inf, fmax=55, dB=True)
    fig.savefig(op.join(fig_path,'PSD_raw.png'))
    plt.close()
    
    # Reject bad segments
    dummyeve = make_fixed_length_events(raw, id=99, duration=1)
    pseudoepo = Epochs(raw, dummyeve, tmin=0., tmax=0.999, baseline=None, reject=reject, reject_by_annotation=True)
    pseudoepo.drop_bad()
    n_chans = len(pseudoepo.info['ch_names'])
    data = pseudoepo.get_data().transpose(1,0,2).reshape(n_chans,-1)
    raw_cln = RawArray(data, raw.info)
    log_len += [raw_cln.last_samp/raw_cln.info['sfreq']]
    log_drop += [pseudoepo.drop_log_stats()]

    # RUN ICA
    ica = ICA(n_components=0.95, method='fastica', random_state=0)
    ica.fit(raw_cln, picks=picks_meg, decim=3, reject=reject, verbose=True, reject_by_annotation=True)    

    # REMOVE COMPONENTS
    picks_eXg = pick_types(raw_cln.info, meg=False, eeg=False, eog=True, ecg = True, emg=False, misc=False, stim=False, exclude='bads')
    raw_cln.filter(bandpass_freqs[0], bandpass_freqs[1], n_jobs=3, picks=picks_eXg)
    raw_cln.notch_filter(50, n_jobs=3, picks=picks_eXg)             # Remove residual 50Hz line noise
    
    # Find ECG artifacts
    ecg_epochs = create_ecg_epochs(raw_cln, ch_name='ECG003', tmin=-.5, tmax=.5)    #, picks=picks)
    ecg_inds, ecg_scores = ica.find_bads_ecg(ecg_epochs, method='ctps', verbose=False)
    
    # Update reject info
    ica.exclude += ecg_inds[:n_max_ecg]

    # Plot ECG ICs for inspection
    ecg_scores_fig = ica.plot_scores(ecg_scores, exclude=ecg_inds, title='Component score (ECG)', show=True)
    ecg_scores_fig.savefig(op.join(fig_path,'ICA_ecg_comp_score.png'))
    plt.close()
    
    if ecg_inds:
        # show_picks = np.abs(ecg_scores).argsort()[::-1][:5]
        
        # ecg_comp_fig = ica.plot_components(ecg_inds, title='ecg comp', colorbar=True, show=False)
        # ecg_comp_fig.savefig(op.join(fig_path,'ICA_ecg_comp_topo.png'))
        # plt.close()
                
        # plot ECG sources + selection
        ecg_evoked = ecg_epochs.average()
        ecg_evo_fig1 = ica.plot_overlay(ecg_evoked, exclude=ecg_inds, show=False)    
        ecg_evo_fig1.savefig(op.join(fig_path, 'ICA_ecg_overlay.png'))
        plt.close()
    
    # Find EOG artifacts 
    eog_epochs = create_eog_epochs(raw_cln, reject=reject, tmin=-.5, tmax=.5)  # get single EOG trials
    eog_inds, eog_scores = ica.find_bads_eog(eog_epochs)
    
    # Update reject info
    ica.exclude += eog_inds[:n_max_eog]

    # Plot EOG ICs for inspection
    eog_scores_fig = ica.plot_scores(eog_scores, exclude=eog_inds, title='Component score (EOG)', show=True)
    eog_scores_fig.savefig(op.join(fig_path, 'ICA_eog_comp_score.png'), show=False)
    plt.close()

    if eog_inds:
        # eog_comp_fig = ica.plot_components(eog_inds, title="Sources related to EOG artifacts", colorbar=True,show=False)
        # eog_comp_fig.savefig(op.join(fig_path, 'ICA_eog_comp.png'))
        # plt.close()
                                  
        # plot EOG sources + selection
        eog_evo_fig = ica.plot_overlay(eog_epochs.average(), exclude=eog_inds, show=False)  # plot EOG cleaning
        eog_evo_fig.savefig(op.join(fig_path, 'ICA_eog_overlay.png'))
        plt.close()
        
    # Plot components
    ica_fig = ica.plot_components()
    [fig.savefig(op.join(fig_path,'ICA_allComp'+str(i)+'.png')) for i, fig in enumerate(ica_fig)]
    log_ica.append(len(ica.exclude))

    # Apply  ICA to Raw
    raw_ica = ica.apply(raw_cln)
    
#    # Inspect
#    raw_ica.plot(eve)

    # Save
    raw_ica.save(outfname, overwrite=overwrite_old_files)
    ica.save(out_icaFname)
    
    # PSD for diagnostics
    fig = raw_ica.plot_psd(tmax=np.inf, fmax=55, dB=True)
    fig.savefig(op.join(fig_path,'PSD_cln.png'))
    plt.close()

    print('----------- FINISHED '+subj+' -----------------')
    plt.close('all')
    
# %% Save data log
dd = list(zip(log_subj, log_len, log_drop, log_ica))
df = pd.DataFrame(dd, columns=['subj','data_length','drop_pct','ica_remove'])
df.to_csv('/home/mikkel/PD_longrest/groupanalysis/data_log.csv', index = False, header=True)

print('----------- FINISHED ALL -----------------')
#END