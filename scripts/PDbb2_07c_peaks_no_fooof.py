#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Additional analysis of the power spectral density of the ROI time series to estimate alpha and beta peak power and peak frequency without removal of the 1/f-component
Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M., Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease. medRxiv.org. https://doi.org/10.1101/2021.06.27.21259592

@author: joswal

"""
# Import fooof functions for creating spectra and managing parameters
import os
import numpy as np
import scipy.io as sio
import os.path as op
import csv
from scipy.signal import find_peaks

# load all files with PSD filestring
study = 'parkinsons_longitudinal'#'pd_longitudinal'#'parkinson_motor' # 'pd_longitudinal' #
archive_nr = '20079' # '20055' # '20079'
psd_path = '/archive/'+archive_nr + '_' + study + '/analysis/meg_data'
subj_folders = [name for name in os.listdir(psd_path) if os.path.isdir(os.path.join(psd_path, name))]
             
summary_all = []

for subj in subj_folders:
    psdfile    = op.join(psd_path,subj, subj+'-ts-psd2lh.mat')
    power = sio.loadmat(psdfile)['psd'][0]
    frequency = sio.loadmat(psdfile)['freqs'][0]
     
    beta_frequency = frequency[frequency> 13]
    beta_power = power[frequency>13]
    beta_power = beta_power[beta_frequency< 30]
    beta_frequency = beta_frequency[beta_frequency<30]
     
    alpha_frequency = frequency[frequency> 8]
    alpha_power = power[frequency> 8]
    alpha_power = alpha_power[alpha_frequency< 13]
    alpha_frequency = alpha_frequency[alpha_frequency<13]

     
    # Identify peaks within the beta frequency range
    peak_indices_beta, _ = find_peaks(beta_power, height=beta_power.max() * 0.1)
    # Find the peak frequency within the beta range
    raw_beta_cf = beta_frequency[peak_indices_beta[np.argmax(beta_power[peak_indices_beta])]]
    # Cget power at the peak frequency
    raw_beta_pw = beta_power[peak_indices_beta[np.argmax(beta_power[peak_indices_beta])]]
   
    # Identify peaks within the beta frequency range
    peak_indices_alpha, _ = find_peaks(alpha_power, height=alpha_power.max() * 0.1)
    # Find the peak frequency within the beta range
    raw_alpha_cf = alpha_frequency[peak_indices_alpha[np.argmax(alpha_power[peak_indices_alpha])]]
    # Cget power at the peak frequency
    raw_alpha_pw= alpha_power[peak_indices_alpha[np.argmax(alpha_power[peak_indices_alpha])]]
   
    summary_subj = [subj,raw_alpha_cf, raw_alpha_pw, raw_beta_cf, raw_beta_pw]
    summary_all.append(summary_subj)

output = '/home/joswal/pd_longitudinal/subj_data/peaks_no_fooof.csv'
with open(output, 'w', ) as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        for word in summary_all:
            wr.writerow([word])
