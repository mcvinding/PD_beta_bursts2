% First find common threshold for defining beta bursts. The napply threshold 
% to individual data.
% set paths
clear all
close all
addpath('/home/mikkel/fieldtrip/fieldtrip/')
ft_defaults
addpath('/home/mikkel/beta_bursts/functions')
addpath('/home/mikkel/PD_longrest/scripts/')
[subjects, megdir] = PDbb2_SETUP();

%% Settings
overwrite = 1;   % Overwirte old files 0=false or 1=true

steps   = 0:0.1:4;
labels  = {'lh_roi','rh_roi'};
fsample = 1000; %Hz

%% Find peaks
rhomat = zeros(length(subjects),length(steps), 2);
for ss = 1:length(subjects)
    subj = subjects{ss};
    fprintf('Now loading subj %s...\n', subj)
    % Input
    infile_rh = fullfile(megdir, subj,[subj,'-ts-rawtc-rh.mat']);
    infile_lh = fullfile(megdir, subj,[subj,'-ts-rawtc-lh.mat']);
    
    % Output
    outfname_raw = fullfile(megdir, subj, 'roidata.mat');
    outfname_hlb = fullfile(megdir, subj, 'roidata_hlbt.mat');
        
    % Load data
    load(infile_rh);
    rh_dat = label_tc;
    load(infile_lh);
    lh_dat = label_tc;
    clear label_tc
    
    % Make pseudo data
    roidata = [];
    roidata.trial      = {[lh_dat; rh_dat]};
    roidata.time       = {(1:length(lh_dat))/fsample};
    roidata.label      = labels;
    roidata.fsample    = fsample;
    roidata.dimord     = 'chan_time';
    roidata.sampleinfo = [1, length(lh_dat)];
    
    % Save
    if ~exist(outfname_raw, 'file') || overwrite
        fprintf('saving %s...', outfname_raw); save(outfname_raw, 'roidata'); disp('done')
    end
    
    % Band-pass to beta band
    cfg = [];
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [13 30];
    cfg.hilbert     ='abs';
    
    roidata_hlb = ft_preprocessing(cfg, roidata);
    
    % Save
    if ~exist(outfname_hlb, 'file') || overwrite
        fprintf('saving %s...', outfname_hlb); save(outfname_hlb, 'roidata_hlb'); disp('done')
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
        [~, rhomat(ss,:,ii)] = find_betaevents(cfg, roidata_hlb);
        disp('done')
    end
end

save('/home/mikkel/PD_longrest/groupanalysis/rhomats.mat','rhomat')
disp('done')

%% Find cutoff

cutoff_mdamp = find_threshold(rhomat, steps, 1); title('med amp')

%END