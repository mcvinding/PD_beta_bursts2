%% Setup PDbb2 in MATLAB

function [subjects, meg_path] = PDbb2_SETUP()
%% Paths
% raw_path        = '/archive/20079_parkinsons_longitudinal/MEG/'
% old_raw_path    = '/archive/20055_parkinson_motor/MEG'
meg_path        = '/home/mikkel/PD_longrest/meg_data';
% trans_path      = '/home/mikkel/PD_long/trans_files'           # !! Files are in subj_folder
% old_trans_path  = '/home/mikkel/PD_motor/tap/trans_files'
% fs_subjects_dir = '/home/mikkel/PD_long/fs_subjects_dir'
subj_data_path  = '/home/mikkel/PD_long/subj_data/';

%% Read subjects
subj_file = fullfile(subj_data_path, 'subjects_and_dates.csv');

subj_dat = readtable(subj_file);
subj_names = table2cell(subj_dat(:,1));
include = subj_dat(:,4);

for ss = 1:length(subj_names)
    subj_names{ss} = ['0', num2str(subj_names{ss})];
end

idxer = table2array(include)==1;
subjects = subj_names(idxer);

end
%END