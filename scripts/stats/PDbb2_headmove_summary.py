#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get summaries of head movement in scanner per group.
@author: mikkel

For information  please refer to the paper :
    Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M., Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease. medRxiv.org. https://doi.org/10.1101/2021.06.27.21259592
    
"""
import mne
import os.path as op
import numpy
import pandas as pd
import os
import numpy as np
import sys
sys.path.append('/home/mikkel/PD_longrest/scripts')
from PDbb2_SETUP import subjects_and_dates, raw_path, exceptions, filestring, old_subjs, old_raw_path, old_filestring, group

#%% Init
# Triggers
startTrigger = 1
stopTrigger  = 64

# Group index
ptnsidx = [i for i,x in enumerate(group) if x=='patient']
ctrlidx = [i for i,x in enumerate(group) if x=='control']
grp = np.ones(len(group))
grp[ctrlidx] = 2

#%% Read and collect headmovement from cHPI
total_move = [0]*len(subjects_and_dates)
avg_move = [0]*len(subjects_and_dates)

for ii, subj_date  in enumerate(subjects_and_dates):
    subj = subj_date.split('/')[0][-4:]
    print('Now processing subject '+subj)
    
    if subj == '0582':
        total_move[ii] = float('nan')
        avg_move[ii] = float('nan')
        continue
    
    if subj in old_subjs:
        raw_fpath   = op.join(old_raw_path, subj_date)
    else:
        raw_fpath   = op.join(raw_path, subj_date)
    
    file_list = os.listdir(raw_fpath)
    if subj in exceptions:
        inFiles = [op.join(raw_fpath,f) for f in file_list if exceptions[subj] in f]
    elif subj in old_subjs:
        inFiles = [op.join(raw_fpath,f) for f in file_list if old_filestring in f]
    else:
        inFiles = [op.join(raw_fpath,f) for f in file_list if filestring in f]
    inFile = inFiles[0]
    
    # Read MEG data
    raw = mne.io.read_raw(inFile, preload=True)
    eve = mne.find_events(raw, stim_channel='STI101', initial_event=True)
    
    # Trigger exceptions
    if subj == '0333':                                                  # Missing triggers
        startSam = 20000+raw.first_samp
        stopSam = 180102+startSam       
    elif subj in ['0322','0523','0524','0525','0529']:
        startSam = eve[eve[:,2] == startTrigger,0][0]
        stopSam = startSam+180102     
    elif subj in ['0548','0583']:                                       # No start trigger. Start trigger val is stop trigger.
        stopSam = eve[eve[:,2] == startTrigger,0][0]
        startSam = stopSam-180102                                     
    elif subj == '0590':                                                # No start trigger.
        startSam = raw.first_samp+1000             
        stopSam = startSam+180102
    elif subj == '0605':                                                # Triggers numbers
        startSam = eve[eve[:,2] == 14593,0][0]
        stopSam = eve[eve[:,2] == 14656,0][0]
    elif subj == '0615':                                                # No triggers.
        startSam = raw.first_samp+10000             
        stopSam = startSam+180102
    else:
        startSam = eve[eve[:,2] == startTrigger,0][0]
        stopSam = eve[eve[:,2] == stopTrigger,0][0]
        
    # Crop and find movement in window
    raw.crop(tmin=(startSam - raw.first_samp ) / raw.info['sfreq'], tmax=(stopSam - raw.first_samp ) / raw.info['sfreq'])
    
    quat, times = raw.pick_types(chpi=True).get_data(return_times=True)
    
    qT = numpy.vstack((times, quat)).transpose()
    # fig = mne.viz.plot_head_positions(qT, mode='traces', show=False)#, info=info)
    
    d1 = quat[0:3,0:-2]
    d2 = quat[0:3,1:-1]
        
    dnorm = numpy.sqrt(numpy.sum((d1-d2)**2, axis=0))
    dcum = numpy.cumsum(dnorm)
    
    avg = max(dcum)/(max(times)-min(times))
        
    text = 'Moved a cummulative total of %.2f cm during the session\nTotal recording time: %.2f min. (%.1f s)\nAverage movement: %.2f mm/s\n'
    print(text % (max(dcum)*100, max(times)/60, max(times), avg*1000))
    
    avg_move[ii] = avg*1000         # mm/s
    total_move[ii] = max(dcum)*100  # cm

#%% stat
import scipy.stats as stats

df = pd.DataFrame({'Group':grp, 'move':total_move, 'avgmove':avg_move})
df = df.dropna()

med =  df.groupby(['Group']).apply(lambda x: np.median(x['avgmove']))
men =  df.groupby(['Group']).apply(lambda x: np.mean(x['avgmove']))
std =  df.groupby(['Group']).apply(lambda x: np.std(x['avgmove']))

U, p = stats.mannwhitneyu(df['avgmove'][df['Group']==2], df['avgmove'][df['Group']==1])
t, p = stats.ttest_ind((df['avgmove'][df['Group']==1]),(df['avgmove'][df['Group']==2]), equal_var=False)
r, p = stats.ranksums(df['avgmove'][df['Group']==2], df['avgmove'][df['Group']==1])


#%% Make plot
import seaborn as sns
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(4, 5))
custom_params = {"axes.spines.right": False, "axes.spines.top": False, "figure.figsize":(5, 5)}
sns.set_theme(style='white', rc=custom_params)
# plt.ylabel('Average movement (mm/s)', weight='bold')
fig = sns.swarmplot(x='Group', y='avgmove', data=df, dodge=False, size=6, hue='Group', legend=False)
fig.set_ylabel('Average movement (mm/s)', weight='bold')
fig.set_xlabel('Group', weight='bold')
zz = fig.get_figure()
zz.tight_layout()
zz.savefig("/home/mikkel/PD_longrest/figures/move.png", dpi=600) 

#END