#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mikkel
"""
#import numpy as np
import os.path as op
#import mne
from mne.time_frequency import psd_array_welch
import sys
sys.path.append('/home/mikkel/PD_longrest/scripts')
from PDbb2_SETUP import subjects, meg_path
import scipy.io as sio
import matplotlib.pyplot as plt

#%% Settings
overwrite = False

freq_range = [0.5, 45]
t_win = 3072
              
#%% RUN
all_psd = []

for subj in subjects:
        
    hemi = 'lh'
    subj_path   = op.join(meg_path, subj)
    tcfname     = op.join(subj_path, subj+'-ts-rawtc-'+hemi+'.mat')            # Raw time-series
    psdfname    = op.join(subj_path, subj+'-ts-psd'+hemi+'.mat')
    
    # Load data
    dat = sio.loadmat(tcfname)['label_tc'][0]    
    
    psdx = dict()
    psdx['psd'], psdx['freqs'] = psd_array_welch(dat, sfreq=1000, fmin=freq_range[0], fmax=freq_range[1], n_fft=t_win, n_overlap=t_win/2)
    
#    plt.plot(freqs, np.log10(psd))
    sio.savemat(psdfname, psdx)

#END