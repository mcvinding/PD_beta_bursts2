#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create BEM solution, source model and, forward model
@author: mikkel
"""
import mne
import sys
import os.path as op
sys.path.append('/home/mikkel/PD_longrest/scripts')
from PDbb2_SETUP import subjects, meg_path, fs_subjects_dir, spacing

#%% Run settings
overwrite = False

no_trans = []
no_bem = []

#%% RUN
for subj in subjects:
    
    # Filenames
    sub_megdir = op.join(meg_path, subj)
    sub_mridir = op.join(fs_subjects_dir, subj)
    
    # Input
    trans_fname     = op.join(sub_megdir, subj+'-trans.fif')
    bem_surf_fname  = op.join(sub_mridir, 'bem', 'inner_skull.surf')
    raw_fname       = op.join(sub_megdir, subj+'-ica-raw.fif')
    
    # output
    bem_sol_fname   = op.join(sub_megdir, subj+'-bem-sol.fif')
    src_fname       = op.join(sub_megdir, subj+'-'+spacing+'-src.fif')
    fwd_fname       = op.join(sub_megdir, subj+'-'+spacing+'-fwd.fif')
    
    print('Processing subject:', subj)
    
    # Check if ok for RUN
    if op.exists(fwd_fname) and not overwrite:
        print(fwd_fname+' already exists. Continue!')
        continue
    if not op.exists(trans_fname):
        print(trans_fname+' doen not exists for subj '+subj+'. Continue!')
        no_trans += [subj]
        continue
    if not op.exists(bem_surf_fname):
        print(bem_surf_fname+' doen not exists for subj '+subj+'. Continue!')
        no_bem += [subj]
        continue
    if not op.exists(raw_fname):          
        print(raw_fname+' doen not exists for subj '+subj+'. Continue!')
        continue    
    
    # #%Make dirs
    # if not op.exists(sub_megdir):
    #     os.mkdir(sub_megdir)
    
    # Create source space in individual subject
    # The is-else is to avoid overwiritng if the script is run several times
    if not op.exists(src_fname) or overwrite:
        src = mne.setup_source_space(subj, 
                                     spacing=spacing, 
                                     subjects_dir=fs_subjects_dir, 
                                     n_jobs=3, 
                                     add_dist=False)
        
        mne.write_source_spaces(src_fname, src, overwrite=overwrite)
    else:
        src = mne.read_source_spaces(src_fname)
        
    # Read BEM surfaces and make the BEM solution
    # Plot for inspection
    # mne.viz.plot_bem(subject=subj, subjects_dir=subjects_dir, brain_surfaces='white', orientation='coronal')
    
    if not op.exists(bem_sol_fname) or overwrite:
        bem_mod = mne.make_bem_model(subj,
                                     ico=4, 
                                     subjects_dir=fs_subjects_dir, 
                                     conductivity=(0.3,))
        
        bem = mne.make_bem_solution(bem_mod)
        mne.write_bem_solution(bem_sol_fname, bem)
    else:
        bem = mne.read_bem_solution(bem_sol_fname)
                     
    # Read transformation (i.e. the coregistration)    
    trans = mne.read_trans(trans_fname)


    # compute the forward operator, commonly referred to as leadfield matrix
    fwd = mne.make_forward_solution(raw_fname,
                                    trans=trans,
                                    src=src,
                                    bem=bem,
                                    meg=True, eeg=False,
                                    mindist=1.0,
                                    n_jobs=3)
    
    mne.write_forward_solution(fwd_fname, fwd, overwrite=overwrite)

#END