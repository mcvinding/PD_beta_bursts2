#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of the power spectral density of the ROI time series: use FOOOF /https://fooof-tools.github.io/fooof/index.html) 
to estimate alpha and beta power and 1/f broadband intercept and slope.

Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M., Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease. medRxiv.org. https://doi.org/10.1101/2021.06.27.21259592

@author: mcvinding
"""
# Import fooof functions for creating spectra and managing parameters
import numpy as np
import pandas
import scipy.io as sio
import matplotlib.pyplot as plt
import os.path as op
from os import mkdir
import sys
sys.path.append('/home/mikkel/PD_longrest/scripts')
from PDbb2_SETUP import subjects, meg_path

# Import fooof functions
from fooof import FOOOF
from fooof import FOOOFGroup
from fooof.bands import Bands
from fooof.analysis import get_band_peak_fg #get_band_peak_fm
#from fooof.plts.spectra import plot_spectrum, plot_spectra
#from fooof.plts.spectra import plot_spectrum_shading, plot_spectra_shading
from fooof.plts.periodic import plot_peak_fits, plot_peak_params
from fooof.plts.aperiodic import plot_aperiodic_params, plot_aperiodic_fits

#%% SETTINGS
# Define frequency bands of interest
bands = Bands({'delta' : [1, 4],
               'theta' : [4, 7],
               'alpha' : [7, 13],
               'beta'  : [13, 30],
               'gamma' : [30, 45]})

# Constrain analysis
max_n_peaks     = 8
peak_wlim       = [0.75, 12]
freq_range      = [1, 45]
peak_threshold  = 2
min_peak_height = 0.05

#%% Load and collect data
all_psd = np.zeros((len(subjects), 137))
all_freqs = np.zeros((len(subjects), 137))

for ii, subj in enumerate(subjects):

    hemi = 'lh'
    subj_path   = op.join(meg_path, subj)
    psdfname    = op.join(subj_path, subj+'-ts-psd2'+hemi+'.mat')
    figdir      = op.join(subj_path, 'plots')
    if not op.exists(figdir):
        mkdir(figdir)
    
    # Load data
    all_psd[ii,:] = sio.loadmat(psdfname)['psd'][0]
    all_freqs[ii,:] = sio.loadmat(psdfname)['freqs'][0]

    # Individual FOOOF analysis  (for inspection)
    fm = FOOOF(max_n_peaks=max_n_peaks, peak_threshold=peak_threshold, peak_width_limits=peak_wlim, 
                aperiodic_mode='fixed', min_peak_height=min_peak_height)

    # Report: fit the model, print the resulting parameters, and plot the reconstruction
    fm.report(all_freqs[ii,], all_psd[ii,], freq_range)
    plt.title(subj)
    # plt.savefig(op.join(figdir, 'fooof_psd2.jpg'))
    # plt.close()

#%% FOOOF analysis
fg = FOOOFGroup(max_n_peaks=max_n_peaks, peak_threshold=peak_threshold, peak_width_limits=peak_wlim, 
                aperiodic_mode='fixed', min_peak_height=min_peak_height)
fg.fit(all_freqs[0], all_psd)

# Get parameters
aper  = fg.get_params('aperiodic_params')
delta = get_band_peak_fg(fg, bands.delta)
theta = get_band_peak_fg(fg, bands.theta)
alpha = get_band_peak_fg(fg, bands.alpha)
beta  = get_band_peak_fg(fg, bands.beta)
gamma = get_band_peak_fg(fg, bands.gamma)

fg.plot()

plot_peak_params(alpha, freq_range=bands.alpha)
plot_peak_params(beta, freq_range=bands.beta)
plot_aperiodic_params(aper)

#%% Get model error
fooof_err = [None]*len(subjects)
for ii, subj in enumerate(subjects):
    fooof_err[ii] = fg.get_fooof(ind=ii).r_squared_
    
fooof_err = np.array(fooof_err)

#%% Export
df_dct = {'subj':subjects,
      'a_intercept': aper[:,0],
      'a_slope': aper[:,1],
      'delta_cf': delta[:,0],
      'delta_pw': delta[:,1],
      'delta_bw': delta[:,2],
      'theta_cf': theta[:,0],
      'theta_pw': theta[:,1],
      'theta_bw': theta[:,2],
      'alpha_cf': alpha[:,0],
      'alpha_pw': alpha[:,1],
      'alpha_bw': alpha[:,2],
      'beta_cf': beta[:,0],
      'beta_pw': beta[:,1],
      'beta_bw': beta[:,2],
      'gamma_cf': gamma[:,0],
      'gamma_pw': gamma[:,1],
      'gamma_bw': gamma[:,2],
      'rsqrd': fooof_err}

df = pandas.DataFrame(df_dct)

# save to csv    
df.to_csv('/home/mikkel/PD_longrest/groupanalysis/fooof_df2.csv', index=False, sep=';')


#%% Plot example
fm.report(all_freqs[0], all_psd[1])
plt.savefig('/home/mikkel/PD_longrest/figures/psd_example.jpg', dpi=800)

#END
