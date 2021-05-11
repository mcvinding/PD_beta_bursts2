#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot PSD split by group. @author: mikkel
"""
import numpy as np
import matplotlib.pyplot as plt
import os.path as op
import pandas as pd
import scipy.io as sio
import sys
sys.path.append('/home/mikkel/PD_longrest/scripts')
from PDbb2_SETUP import subjects, meg_path

#%% Read group data
datafile = '/home/mikkel/PD_longrest/groupanalysis/alldata_subj2.csv'
alldata = pd.read_csv(datafile)

alldata['subj'] = '0'+alldata['subj'].astype(str)

#%% Read PSD data
all_psd = np.zeros((len(subjects), 137))
all_freqs = np.zeros((len(subjects), 137))

for ii, subj in enumerate(alldata.subj):

    hemi = 'lh'
    subj_path   = op.join(meg_path, subj)
    psdfname    = op.join(subj_path, subj+'-ts-psd2'+hemi+'.mat')
    
    # Load data
    all_psd[ii,:] = sio.loadmat(psdfname)['psd'][0]
    all_freqs[ii,:] = sio.loadmat(psdfname)['freqs'][0]
    
    
#%% Split by group
freqs = all_freqs[1,:]
psd_ctrl = all_psd[alldata['group']=='control']
psd_ptns = all_psd[alldata['group']=='patient']
psd_ctrl_avg = np.mean(psd_ctrl, axis=0)
psd_ptns_avg = np.mean(psd_ptns, axis=0)

psd_ctrl_sem = 1.96*np.std(psd_ctrl, axis=0)/np.sqrt(len(psd_ctrl))
psd_ptns_sem = 1.96*np.std(psd_ptns, axis=0)/np.sqrt(len(psd_ptns))

plt.plot(freqs, np.log(psd_ptns_avg), 'b')
plt.fill_between(freqs, np.log(psd_ptns_avg-psd_ptns_sem), np.log(psd_ptns_avg+psd_ptns_sem), color='b', alpha=.1)
plt.plot(freqs, np.log(psd_ctrl_avg), 'r')
plt.fill_between(freqs, np.log(psd_ctrl_avg-psd_ctrl_sem), np.log(psd_ctrl_avg+psd_ctrl_sem), color='r', alpha=.1)
plt.legend(['PD','HC'], fontsize = 12)
plt.xlabel("Frequency (Hz)", size = 16)
plt.ylabel("Log-power", size = 16)

plt.savefig('/home/mikkel/PD_longrest/figures/psd_by group.jpg', dpi=600)
