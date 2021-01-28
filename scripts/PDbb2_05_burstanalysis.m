% Find beta bursts for each subject
% To do: clean up in this script
clear all; close all
addpath('/home/mikkel/fieldtrip/fieldtrip/')
ft_defaults
addpath('/home/mikkel/beta_bursts/functions')
addpath('/home/mikkel/PD_longrest/scripts/')
[subjects, dirs] = PDbb2_SETUP();

%% Settings
overwrite = 1;   % Overwirte old files 0=false or 1=true

fsample = 1000;
labels  = {'lh_roi'};

%% Get beta summary
% General settings
cfg = [];
cfg.halfmax     = 'mixed';
cfg.makeplot    = 'no';
cfg.channel     = 'lh_roi';

% Run
for ss = 1:length(subjects)
    subj = subjects{ss};
    fprintf('Reading subj %s (%i of %i)\n', subj, ss, length(subjects))

    % Input
    infile_lh = fullfile(dirs.meg_path, subj,[subj,'-ts-rawtc2-lh.mat']);

    % Output
    outfname_raw = fullfile(dirs.meg_path, subj, 'roidata2.mat');
    outfname_beta_hlb = fullfile(dirs.meg_path, subj, 'roidata_beta_hlbt2.mat');
    outfname_mube_hlb = fullfile(dirs.meg_path, subj, 'roidata_mube_hlbt2.mat');

    % Load
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
    
    % Different outputs mu (u) and beta (b), avg thresholds (m1), 2x median (m2),
    % and percentile (pc)
    outfname_b_m2 = fullfile(dirs.meg_path, subj,[subj,'-b_m2-burst2.mat']);
    outfname_b_pc = fullfile(dirs.meg_path, subj,[subj,'-b_pc-burst2.mat']);
    outfname_u_m2 = fullfile(dirs.meg_path, subj,[subj,'-u_m2-burst2.mat']);
    outfname_u_pc = fullfile(dirs.meg_path, subj,[subj,'-u_pc-burst2.mat']);
 
    % Init
    burstsummary_b_m2 = [];
    burstsummary_b_pc = [];
    burstsummary_u_m2 = [];   
    burstsummary_u_pc = [];

    if exist(outfname_b_m2,'file') && ~overwrite
        warning('File %s exists. Continue!', outfname_b_m2);
        continue
    end
    
    % BURST ANALSYIS
    % 2x median beta
    cfg.steps       = 2;
    cfg.cutofftype  = 'med';
    cfg.mindist     = 1/13;
    burstsummary_b_m2 = find_betaevents(cfg, roidata_beta_hlb);
        
    % 80% beta
    cfg.cutofftype  = 'pct';
    cfg.steps       = 80;
    cfg.mindist     = 1/13;
    burstsummary_b_pc = find_betaevents(cfg, roidata_beta_hlb);

    % 2x median mu+beta
    cfg.steps       = 2;
    cfg.cutofftype  = 'med';
    cfg.mindist     = 1/13;
    burstsummary_u_m2 = find_betaevents(cfg, roidata_mube_hlb);

    % 80% mu+beta
    cfg.cutofftype  = 'pct';
    cfg.steps       = 80;
    cfg.mindist     = 1/13;
    burstsummary_u_pc = find_betaevents(cfg, roidata_mube_hlb);

%     end

    % Save
    fprintf('saving... ')
    save(outfname_b_m2, 'burstsummary_b_m2')
    save(outfname_b_pc, 'burstsummary_b_pc')
    save(outfname_u_m2, 'burstsummary_u_m2')
    save(outfname_u_pc, 'burstsummary_u_pc')
    disp('done')
end

%END