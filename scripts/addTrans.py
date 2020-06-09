# -*- coding: utf-8 -*-
"""
Get transformation MRI <-> MEG with MNE coregistration
@author: mikkel
"""
import mne
import sys
import os.path as op
sys.path.append('X:\PD_longrest\scripts')
from PDbb2_SETUP import meg_path, fs_subjects_dir, subjects

#%% Coregistration
# Make sure to give consistent naming when saving -trans files, e.g. "0523-trans.fif"

# Specify subject id
subj = '0650'

#%% Do coregistration
print('------------------ Sub: '+subj+' ------------------')
if sys.platform == 'linux':
    inst = op.join(meg_path, subj, subj+'-ica-raw.fif')
else:
    inst = 'X:/PD_longrest/meg_data/'+subj+'/'+subj+'-ica-raw.fif'
    
mne.gui.coregistration(subject=subj, subjects_dir=fs_subjects_dir, head_high_res=True, inst=inst)

#%% Inspect
# The transformation file obtained by coregistration
#trans = op.join(meg_path, subj, 'trans.fif')
#
#info = mne.io.read_info(raw_fname)
## Here we look at the dense head, which isn't used for BEM computations but
## is useful for coregistration.
#mne.viz.plot_alignment(info, trans, subject=subject, dig=True,
#                       meg=['helmet', 'sensors'], subjects_dir=subjects_dir,
#                       surfaces='head-dense')

#END

#END