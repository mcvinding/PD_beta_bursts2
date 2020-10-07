function [subjects, dirs, subj_date] = PDbb2_SETUP()
% Setup PDbb2 in MATLAB

    % Paths
    dirs.raw_path        = '/archive/20079_parkinsons_longitudinal/MEG/';
    dirs.old_raw_path    = '/archive/20055_parkinson_motor/MEG';
    dirs.meg_path        = '/home/mikkel/PD_longrest/meg_data';
    % old_trans_path  = '/home/mikkel/PD_motor/tap/trans_files'
    % fs_subjects_dir = '/home/mikkel/PD_long/fs_subjects_dir'
    dirs.subj_data_path  = '/home/mikkel/PD_long/subj_data/';
    dirs.group_path      = '/home/mikkel/PD_longrest/groupanalysis';

    % Read subjects
    subj_file = fullfile(dirs.subj_data_path, 'subjects_and_dates.csv');

    subj_dat = readtable(subj_file);
    subj_names = table2cell(subj_dat(:,1));
    include = subj_dat(:,4);
    idxer = table2array(include)==1;

    % Arrange subject id
    for ss = 1:length(subj_names)
        subj_names{ss} = ['0', num2str(subj_names{ss})];
    end

    subjects = subj_names(idxer);

    % Arrange subject and adate folders
    dd = (table2cell(subj_dat(:,3)));
    for ss = 1:length(dd)
        subject_date{ss} = fullfile(subj_names{ss}, num2str(dd{ss}));
    end
    subj_date = subject_date(idxer);

end
%END