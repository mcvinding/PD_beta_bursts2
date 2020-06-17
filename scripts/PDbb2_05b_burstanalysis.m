% Find beta bursts for each subject
clear all
close all
addpath('/home/mikkel/fieldtrip/fieldtrip/')
ft_defaults
addpath('/home/mikkel/beta_bursts/functions')
addpath('/home/mikkel/PD_longrest/scripts/')
[subjects, dirs] = PDbb2_SETUP();

%% Settings
overwrite = 0;   % Overwirte old files 0=false or 1=true

labels  = {'lh_roi','rh_roi'};

%% Choose threshold
load('/home/mikkel/PD_longrest/groupanalysis/rhomats.mat')
cutoff = find_threshold(rhomat, 0:0.1:4, 1); title('Threshold')

%% Get beta summary
cfg = [];
cfg.length      = 3;
cfg.overlap     = 0;
cfg.steps       = cutoff;
cfg.corrtype    = 'amp';
cfg.cutofftype  = 'med';
cfg.halfmax     = 'yes';
cfg.makeplot    = 'no';

neve_lh = nan(size(subjects));
neve_rh = nan(size(subjects));
cut_lh = nan(size(subjects));
cut_rh = nan(size(subjects));

for ss = 1:length(subjects)
    burstsummary = [];
    subj = subjects{ss};
    fprintf('Reading subj %s.\n', subj)
    infile = fullfile(dirs.megdir, subj, 'roidata_hlbt.mat');
    outfname = fullfile(dirs.megdir, subj,[subj,'-burstsummary.mat']);
    
    if exist(outfname,'file') && ~overwrite
        warning('File %s exists. Continue!', outfname);
        continue
    end
    
    % Load data
    load(infile);
    
    % Find beta events
    for ii = 1:length(labels)
        cfg.channel = labels{ii};
        
        [burstsummary{ii}] = find_betaevents(cfg, roidata_hlb);
    end

    % Save
    fprintf('saving %s... ', outfname)
    save(outfname, 'burstsummary')
    disp('done')
    
    % Intermediate summaries for inspection
    neve_lh(ss) = burstsummary{1}.n_events;
    neve_rh(ss) = burstsummary{2}.n_events;
    cut_lh(ss) = burstsummary{1}.cutoff;
    cut_rh(ss) = burstsummary{2}.cutoff;
end

%END