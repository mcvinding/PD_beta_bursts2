% First find common threshold for defining beta bursts. The napply threshold 
% to individual data.
% set paths
clear all; close all;
addpath('/home/mikkel/fieldtrip/fieldtrip/')
ft_defaults
addpath('/home/mikkel/beta_bursts/functions')
addpath('/home/mikkel/PD_longrest/scripts/')
[subjects, dirs] = PDbb2_SETUP();

%% Settings
overwrite = 0;   % Overwirte old files 0=false or 1=true

steps   = 0:0.1:4;
labels  = {'lh_roi'};
fsample = 1000; %Hz

%% Run across thrershold to find optimum
rhomat_beta = zeros(length(subjects),length(steps), 2); % For beta band filtered data
rhomat_mube = zeros(length(subjects),length(steps), 2); % For mu+beta band filtered data

for ss = 1:length(subjects)
    subj = subjects{ss};
    fprintf('Now loading subj %s...\n', subj)
    % Input
    infile_rh = fullfile(dirs.meg_path, subj,[subj,'-ts-rawtc2-rh.mat']);
    infile_lh = fullfile(dirs.meg_path, subj,[subj,'-ts-rawtc2-lh.mat']);
    
    % Output
    outfname_raw = fullfile(dirs.meg_path, subj, 'roidata2.mat');
    outfname_beta_hlb = fullfile(dirs.meg_path, subj, 'roidata_beta_hlbt2.mat');
    outfname_mube_hlb = fullfile(dirs.meg_path, subj, 'roidata_mube_hlbt2.mat');

    % Load data
%     load(infile_rh);
%     rh_dat = label_tc;
    load(infile_lh);
    lh_dat = label_tc;
    clear label_tc
    
    % Make pseudo data
    roidata = [];
    roidata.trial      = {lh_dat};
    roidata.time       = {(1:length(lh_dat))/fsample};
    roidata.label      = labels;
    roidata.fsample    = fsample;
    roidata.dimord     = 'chan_time';
    roidata.sampleinfo = [1, length(lh_dat)];
    
    % Save
    if ~exist(outfname_raw, 'file') || overwrite
        fprintf('saving %s...', outfname_raw); save(outfname_raw, 'roidata'); disp('done')
    end
    
    if ~exist(outfname_beta_hlb, 'file') || overwrite

        % Band-pass to beta band
        cfg = [];
        cfg.bpfilter    = 'yes';
        cfg.bpfreq      = [13 30];
        cfg.hilbert     ='abs';
        roidata_beta_hlb = ft_preprocessing(cfg, roidata);

        % Band-pass to mu+beta band
        cfg = [];
        cfg.bpfilter    = 'yes';
        cfg.bpfreq      = [8 30];
        cfg.hilbert     ='abs';
        roidata_mube_hlb = ft_preprocessing(cfg, roidata);
        
        % Save
        fprintf('saving %s...', outfname_beta_hlb);
        save(outfname_beta_hlb, 'roidata_beta_hlb'); disp('done')
        fprintf('saving %s...', outfname_mube_hlb);
        save(outfname_mube_hlb, 'roidata_mube_hlb'); disp('done')
        
    elseif exist(outfname_beta_hlb, 'file') && ~overwrite
        load(outfname_beta_hlb);
        load(outfname_mube_hlb);
    end
          
    % Find correlations across thresholds
    for ii = 1:length(labels)
        cfg = [];
        cfg.length      = 3;
        cfg.overlap     = 0;
        cfg.steps       = steps;       
        cfg.cutofftype  = 'med';
        cfg.corrtype    = 'amp';
        cfg.channel     = labels{ii};
        [~, rhomat_beta(ss,:,ii)] = find_betaevents(cfg, roidata_beta_hlb);
        
        cfg = [];
        cfg.length      = 3;
        cfg.overlap     = 0;
        cfg.steps       = steps;       
        cfg.cutofftype  = 'med';
        cfg.corrtype    = 'amp';
        cfg.channel     = labels{ii};
        [~, rhomat_mube(ss,:,ii)] = find_betaevents(cfg, roidata_mube_hlb);
        disp('done')
    end
end

save('/home/mikkel/PD_longrest/groupanalysis/rhomats2.mat','rhomat_beta', 'rhomat_mube')
disp('done')

%% Find cutoff
% load('/home/mikkel/PD_longrest/groupanalysis/rhomats')
rhomat_beta = rhomat_beta(:,:,1)
rhomat_mube = rhomat_mube(:,:,1)
threshold = find_threshold(rhomat_beta, steps, 1); title('Threshold')
threshold = find_threshold(rhomat_mube, steps, 1); title('Threshold')

%END