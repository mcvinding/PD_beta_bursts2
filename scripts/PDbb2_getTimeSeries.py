#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Get ROI time course from ROI
@author: mikkel
"""
import numpy as np
import os.path as op
import mne
import sys
sys.path.append('/home/mikkel/PD_longrest/scripts')
from PDbb2_SETUP import subjects, meg_path, fs_subjects_dir, spacing
from sensorymotorROI import make_sensorymotorROI
# from scipy.signal import hilbert
import scipy.io as sio

#%% Run settings
overwrite = False

no_stc = []

#%% Run
for subj in subjects:
    print(subj)
    subj_path   = op.join(meg_path, subj)
    rawfile     = op.join(subj_path, subj+'-ica-raw.fif')
    covfile     = op.join(subj_path, subj+'-cov.fif') 
    fwdfile     = op.join(subj_path, subj+'-'+spacing+'-fwd.fif')
    srcfile     = op.join(subj_path, subj+'-'+spacing+'-src.fif')
    stcfile     = op.join(subj_path, subj+'-dspm-lh.stc')

    # outhilbt = op.join(subj_path, subj+'-ts-hilbt')      # Hilbert envelope
    outrawtc = op.join(subj_path, subj+'-ts-rawtc')      # Raw time-series
    # outrawft = op.join(subj_path, subj+'-ts-rawft')      # Band-pass filtered time-series

    if op.exists(outrawtc) and not overwrite:
        print('File '+outrawtc+' exists. Continue!')
        
    # hilb = dict()
    rawtc = dict()
    # filtc = dict()
    
    if op.exists(srcfile):
        stc = mne.read_source_estimate(stcfile)
        src = mne.read_source_spaces(srcfile)
        # print('yes')
    else:
        no_stc += [subj]
        continue
        
    for hemi in ['lh','rh']:

        lab = make_sensorymotorROI(subj, fs_subjects_dir, hemi=hemi)

        # Save label
        lab.save(op.join(fs_subjects_dir, subj, 'label',hemi+'.sensmotor.label'))

        # Extract label time-series
        label_tc = stc.extract_label_time_course(lab, src, mode='pca_flip')[0,:]
        label_tc  = np.float64(label_tc)
        
        # filtc = mne.filter.filter_data(label_tc, 1000, 13, 30, method='fir',n_jobs=3)

        # analytic = hilbert(filtc)
        # envelope = np.abs(analytic)         
            
    # if not op.exists(outhilbt) or overwrite:    
        # sio.savemat(outhilbt+'-'+hemi+'.mat', dict(envelope=envelope))
    if not op.exists(outrawtc) or overwrite:
        sio.savemat(outrawtc+'-'+hemi+'.mat', dict(label_tc=label_tc))
    # if not op.exists(outrawft) or overwrite:
        # sio.savemat(outrawft+'-'+hemi+'.mat', dict(filtc=filtc))
        
#END