#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mikkel
"""
import numpy as np
import os.path as op
import mne
from mne.time_frequency import psd_array_welch
import sys
sys.path.append('/home/mikkel/PD_longrest/scripts')
from PDbb2_SETUP import subjects, meg_path
import scipy.io as sio
import matplotlib.pyplot as plt

#%% Settings
overwrite = False

#%% RUN

for subj in subjects:
        
    hemi = 'lh'
    subj_path   = op.join(meg_path, subj)
    tcfname     = op.join(subj_path, subj+'-ts-rawtc-'+hemi+'.mat')            # Raw time-series
    psdfname    = op.join(subj_path, subj+'-ts-psd'+hemi+'.???')
    
    # Load data
    dat = sio.loadmat(tcfname)['label_tc'][0]    
    
    psd, freqs = psd_array_welch(dat, sfreq=1000, fmin=0.1, fmax=40, n_fft=4096, n_overlap=2048)
    
    plt.plot(freqs,np.log(psd))

#END