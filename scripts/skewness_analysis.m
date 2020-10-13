% Analyse the skewness of beta and mu+beta time series for each subject
clear all; close all
addpath('/home/mikkel/fieldtrip/fieldtrip/')
ft_defaults
[subjects, dirs] = PDbb2_SETUP();

%% Run
skw_b = nan(length(subjects),1);
krt_b = nan(length(subjects),1);
skw_u = nan(length(subjects),1);
krt_u = nan(length(subjects),1);

for ss = 1:length(subjects)
    subj = subjects{ss};
    fprintf('Reading subj %s.\n', subj)
    
    % Load hilbert data
    infile_b = fullfile(dirs.meg_path, subj, 'roidata_beta_hlbt.mat');
    infile_m = fullfile(dirs.meg_path, subj, 'roidata_mube_hlbt.mat');
    
    % Get data
    load(infile_b); % Hilbert beta
    dat_b = roidata_beta_hlb.trial{1}(1,:);
    skw_b(ss) = skewness(dat_b);
    krt_b(ss) = kurtosis(dat_b);
    
    load(infile_m); % Hilbert mu+beta
    dat_u = roidata_mube_hlb.trial{1}(1,:);
    skw_u(ss) = skewness(dat_u);
    krt_u(ss) = kurtosis(dat_u);
    
    clear dat*
end

% SAVE
save('/home/mikkel/PD_longrest/groupanalysis/skw.mat', 'skw_b','skw_u')
save('/home/mikkel/PD_longrest/groupanalysis/krt.mat', 'krt_b','krt_u')
disp('done')

%END