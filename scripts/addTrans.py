# -*- coding: utf-8 -*-
"""
Get transformation MRI <-> MEG with MNE coregistration
@author: mikkel
"""
import mne
import sys
import os.path as op
sys.path.append('/home/mikkel/PETMEG/scripts')
from PETMEG_config import subjects, meg_path, subjects_dir

#%% Coregistration
# Make sure to give consistent naming when saving -trans files, e.g. "0523-trans.fif"

# Specify subject id
subj = '0588'

# Do coregistration
print('------------------ Sub: '+subj+' ------------------')
mne.gui.coregistration(subject=subj, 
               subjects_dir=subjects_dir, 
               head_high_res=True)

#%% Inspect
# The transformation file obtained by coregistration
trans = op.join(meg_path, subj, 'trans.fif')

info = mne.io.read_info(raw_fname)
# Here we look at the dense head, which isn't used for BEM computations but
# is useful for coregistration.
mne.viz.plot_alignment(info, trans, subject=subject, dig=True,
                       meg=['helmet', 'sensors'], subjects_dir=subjects_dir,
                       surfaces='head-dense')

#END

#END