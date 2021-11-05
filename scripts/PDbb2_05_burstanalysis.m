%% Find beta bursts for each subject
% Import ROI time series. Filter and find bursts.
% The tool for finding bursts are found herE: https://github.com/mcvinding/beta_bursts
%
% Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M., Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease. medRxiv.org. https://doi.org/10.1101/2021.06.27.21259592
%
%@author: mcvinding

clear all; close all
addpath('/home/mikkel/fieldtrip/fieldtrip/')
ft_defaults
addpath('/home/mikkel/beta_bursts/functions')
addpath('/home/mikkel/PD_longrest/scripts/')
[subjects, dirs] = PDbb2_SETUP();

%% Settings
overwrite = 0;   % Overwirte old files 0=false or 1=true

fsample = 1000;
labels  = {'lh_roi'};

%% Get beta summary
% General settings
cfg = [];
cfg.halfmax     = 'mixed';
cfg.makeplot    = 'no';
cfg.channel     = 'lh_roi';

%% Run
for ss = 1:length(subjects)
    subj = subjects{ss};
    fprintf('Reading subj %s (%i of %i)\n', subj, ss, length(subjects))

    % Input
    infile_lh = fullfile(dirs.meg_path, subj,[subj,'-ts-rawtc2-lh.mat']);

    % Output
    outfname_raw = fullfile(dirs.meg_path, subj, 'roidata2.mat');
    outfname_mube_hlb = fullfile(dirs.meg_path, subj, 'roidata_mube_hlbt2.mat');

    % Load
    load(infile_lh);
    lh_dat = label_tc;
    clear label_tc
    
    % Make FT data
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
    
    if ~exist(outfname_mube_hlb, 'file') || overwrite
         
        % Band-pass to mu band
        cfg = [];
        cfg.bpfilter    = 'yes';
        cfg.bpfreq      = [8 30];
        cfg.hilbert     ='abs';
        roidata_mube_hlb = ft_preprocessing(cfg, roidata);
        
        % Save
        fprintf('saving %s...', outfname_mube_hlb);
        save(outfname_mube_hlb, 'roidata_mube_hlb'); disp('done')
        
    elseif exist(outfname_beta_hlb, 'file') && ~overwrite
        load(outfname_mube_hlb);
     end
    
    % Outputs mu (u) 2x median (m2)
    outfname = fullfile(dirs.meg_path, subj,[subj,'-u_m2-burst2.mat']);
 
    % Init
    burstsummary_u_m2 = [];   

    if exist(outfname,'file') && ~overwrite
        warning('File %s exists. Continue!', outfname);
        continue
    end
    
    % BURST ANALSYIS
    % Median + 2x median mu
    cfg.steps       = 2;
    cfg.cutofftype  = 'med';
    cfg.mindist     = 1/13;
    burstsummary_u_m2 = find_betaevents(cfg, roidata_mube_hlb);

    % Save
    fprintf('saving... ')
    save(outfname, 'burstsummary_u_m2')
    disp('done')
end

%END
