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

labels  = {'lh_roi','rh_roi'};

%% Choose threshold
load('/home/mikkel/PD_longrest/groupanalysis/rhomats.mat')
cutoff_beta = find_threshold(rhomat_beta, 0:0.1:4, 1); title('Threshold')
cutoff_mube = find_threshold(rhomat_mube, 0:0.1:4, 1); title('Threshold')

%% Get beta summary
% General settings
cfg = [];
cfg.halfmax     = 'mixed';
cfg.makeplot    = 'no';
cfg.channel     = 'lh_roi';

% neve_lh = nan(size(subjects));
% neve_rh = nan(size(subjects));
% cut_lh = nan(size(subjects));
% cut_rh = nan(size(subjects));

for ss = 41:length(subjects)
    subj = subjects{ss};
    fprintf('Reading subj %s (%i of %i)\n', subj, ss, length(subjects))

    % Load hilbert data
    infile_b = fullfile(dirs.meg_path, subj, 'roidata_beta_hlbt.mat');
    infile_u = fullfile(dirs.meg_path, subj, 'roidata_mube_hlbt.mat');

    % Different outputs mu (u) and beta (b), avg thresholds (m1), 2x median (m2),
    % and percentile (pc)
    outfname_b_m2 = fullfile(dirs.meg_path, subj,[subj,'-b_m2-burst.mat']);
    outfname_b_m1 = fullfile(dirs.meg_path, subj,[subj,'-b_m1-burst.mat']);
    outfname_b_pc = fullfile(dirs.meg_path, subj,[subj,'-b_pc-burst.mat']);
    outfname_u_m2 = fullfile(dirs.meg_path, subj,[subj,'-u_m2-burst.mat']);
    outfname_u_m1 = fullfile(dirs.meg_path, subj,[subj,'-u_m1-burst.mat']);
    outfname_u_pc = fullfile(dirs.meg_path, subj,[subj,'-u_pc-burst.mat']);
    
    % Init
    burstsummary_b_m1 = [];
    burstsummary_b_m2 = [];
    burstsummary_b_pc = [];
    burstsummary_u_m1 = [];
    burstsummary_u_m2 = [];   
    burstsummary_u_pc = [];

%     if exist(outfname,'file') && ~overwrite
%         warning('File %s exists. Continue!', outfname);
%         continue
%     end
%     
    % Load data
    load(infile_b); % Hilbert beta
    load(infile_u); % Hilbert mu+beta

    % Find beta events
%     for ii = 1:length(labels)
%         cfg.channel = labels{ii};
        ii=1;
        
        % Standard beta
        cfg.cutofftype  = 'med';
        cfg.steps       = cutoff_beta;
        cfg.mindist     = 1/13;
        [burstsummary_b_m1{ii}] = find_betaevents(cfg, roidata_beta_hlb);
               
        % 2x median beta
        cfg.steps       = 2;
        cfg.cutofftype  = 'med';
        cfg.mindist     = 1/13;
        [burstsummary_b_m2{ii}] = find_betaevents(cfg, roidata_beta_hlb);
        
        % 80% beta
        cfg.cutofftype  = 'pct';
        cfg.steps       = 80;
        cfg.mindist     = 1/13;
        [burstsummary_b_pc{ii}] = find_betaevents(cfg, roidata_beta_hlb);
        
        % Standard mu+beta
        cfg.cutofftype  = 'med';
        cfg.steps       = cutoff_mube;
        cfg.mindist     = 1/13;
        [burstsummary_u_m1{ii}] = find_betaevents(cfg, roidata_mube_hlb);
               
        % 2x median mu+beta
        cfg.steps       = 2;
        cfg.cutofftype  = 'med';
        cfg.mindist     = 1/13;
        [burstsummary_u_m2{ii}] = find_betaevents(cfg, roidata_mube_hlb);
        
        % 80% mu+beta
        cfg.cutofftype  = 'pct';
        cfg.steps       = 80;
        cfg.mindist     = 1/13;
        [burstsummary_u_pc{ii}] = find_betaevents(cfg, roidata_mube_hlb);

%     end

    % Save
    fprintf('saving... ')
    save(outfname_b_m1, 'burstsummary_b_m1')
    save(outfname_b_m2, 'burstsummary_b_m2')
    save(outfname_b_pc, 'burstsummary_b_pc')
    save(outfname_u_m1, 'burstsummary_u_m1')
    save(outfname_u_m2, 'burstsummary_u_m2')
    save(outfname_u_pc, 'burstsummary_u_pc')
    disp('done')
    
%     % Intermediate summaries for inspection
%     neve_lh(ss) = burstsummary{1}.n_events;
%     neve_rh(ss) = burstsummary{2}.n_events;
%     cut_lh(ss) = burstsummary{1}.cutoff;
%     cut_rh(ss) = burstsummary{2}.cutoff;
end

%END