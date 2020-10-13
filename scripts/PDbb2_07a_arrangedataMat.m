% Get various summaries of beta burst (explorative)
clear all
close all
addpath('/home/mikkel/fieldtrip/fieldtrip/')
ft_defaults
addpath('/home/mikkel/beta_bursts/functions')
addpath('/home/mikkel/PD_longrest/scripts/')
[subjects, dirs] = PDbb2_SETUP();

%% Init arrays
% Different outputs mu (u) and beta (b), avg thresholds (m1), 2x median (m2),
% and percentile (pc)
    
% N events
nevent_b_m1 = nan(length(subjects), 1);
nevent_b_m2 = nan(length(subjects), 1);
nevent_b_pc = nan(length(subjects), 1);
nevent_u_m1 = nan(length(subjects), 1);
nevent_u_m2 = nan(length(subjects), 1);
nevent_u_pc = nan(length(subjects), 1);
% nevent_rh = nan(length(subjects), 1);

% Event length
len_b_m1 = [];
len_b_m2 = [];
len_b_pc = [];
len_u_m1 = [];
len_u_m2 = [];
len_u_pc = [];
% len_rh = [];

% Mak event peak
max_b_m1 = [];
max_b_m2 = [];
max_b_pc = [];
max_u_m1 = [];
max_u_m2 = [];
max_u_pc = [];
% max_rh = [];

% Time until event
tue_b_m1 = [];
tue_b_m2 = [];
tue_b_pc = [];
tue_u_m1 = [];
tue_u_m2 = [];
tue_u_pc = [];
% toe_rh = [];

% Subject identifier
sub_b_m1 = [];
sub_b_m2 = [];
sub_b_pc = [];
sub_u_m1 = [];
sub_u_m2 = [];
sub_u_pc = [];
% sub_rh = [];

% cut = nan(length(subjects),2);

%% Run
for ii = 1:length(subjects)
    subj = subjects{ii};
    fprintf('Reading subj %s (%i of %i)\n', subj, ii, length(subjects))
    
    % Load data
    infname_b_m2 = fullfile(dirs.meg_path, subj,[subj,'-b_m2-burst.mat']);
    infname_b_m1 = fullfile(dirs.meg_path, subj,[subj,'-b_m1-burst.mat']);
    infname_b_pc = fullfile(dirs.meg_path, subj,[subj,'-b_pc-burst.mat']);
    infname_u_m2 = fullfile(dirs.meg_path, subj,[subj,'-u_m2-burst.mat']);
    infname_u_m1 = fullfile(dirs.meg_path, subj,[subj,'-u_m1-burst.mat']);
    infname_u_pc = fullfile(dirs.meg_path, subj,[subj,'-u_pc-burst.mat']);
    
    fprintf('Load data... ')
    load(infname_b_m1)
    load(infname_b_m2)
    load(infname_b_pc)
    load(infname_u_m1)
    load(infname_u_m2)
    load(infname_u_pc)
    disp('done')
    
%     % # Cut off #
%     cut(ii,1) = burstsummary{1}.cutoff;
%     cut(ii,2) = burstsummary{2}.cutoff;

    % # N events #
    nevent_b_m1(ii) = burstsummary_b_m1{1}.n_events;
    nevent_b_m2(ii) = burstsummary_b_m2{1}.n_events;
    nevent_b_pc(ii) = burstsummary_b_pc{1}.n_events;
    nevent_u_m1(ii) = burstsummary_u_m1{1}.n_events;
    nevent_u_m2(ii) = burstsummary_u_m2{1}.n_events;
    nevent_u_pc(ii) = burstsummary_u_pc{1}.n_events;

%     nevent_rh(ii) = burstsummary{2}.n_events;
    
    % # Event duration #
    len_b_m1 = [len_b_m1; burstsummary_b_m1{1}.bdat.evelen(1:end-1)];
    len_b_m2 = [len_b_m2; burstsummary_b_m2{1}.bdat.evelen(1:end-1)];
    len_b_pc = [len_b_pc; burstsummary_b_pc{1}.bdat.evelen(1:end-1)];
    len_u_m1 = [len_u_m1; burstsummary_u_m1{1}.bdat.evelen(1:end-1)];
    len_u_m2 = [len_u_m2; burstsummary_u_m2{1}.bdat.evelen(1:end-1)];
    len_u_pc = [len_u_pc; burstsummary_u_pc{1}.bdat.evelen(1:end-1)];

%     len_rh = [len_rh; burstsummary{2}.bdat.evelen(1:end-1)];
  
    % # Max peak #
    max_b_m1 = [max_b_m1; burstsummary_b_m1{1}.bdat.maxpk(1:end-1)];
    max_b_m2 = [max_b_m2; burstsummary_b_m2{1}.bdat.maxpk(1:end-1)];
    max_b_pc = [max_b_pc; burstsummary_b_pc{1}.bdat.maxpk(1:end-1)];
    max_u_m1 = [max_u_m1; burstsummary_u_m1{1}.bdat.maxpk(1:end-1)];
    max_u_m2 = [max_u_m2; burstsummary_u_m2{1}.bdat.maxpk(1:end-1)];
    max_u_pc = [max_u_pc; burstsummary_u_pc{1}.bdat.maxpk(1:end-1)];

%     max_rh = [max_rh; burstsummary{2}.bdat.maxpk(1:end-1)];


    % # Time between events #
    tt_b_m1 = zeros(length(burstsummary_b_m1{1}.bdat.event-1)-1,1);  
    tt_b_m2 = zeros(length(burstsummary_b_m2{1}.bdat.event-1)-1,1);  
    tt_b_pc = zeros(length(burstsummary_b_pc{1}.bdat.event-1)-1,1);  
    tt_u_m1 = zeros(length(burstsummary_u_m1{1}.bdat.event-1)-1,1);  
    tt_u_m2 = zeros(length(burstsummary_u_m2{1}.bdat.event-1)-1,1);  
    tt_u_pc = zeros(length(burstsummary_u_pc{1}.bdat.event-1)-1,1);  

    for k = 1:length(tt_b_m1)
        tt_b_m1(k) = (burstsummary_b_m1{1}.bdat.event(k+1,1)-burstsummary_b_m1{1}.bdat.event(k,2))/1000;
    end
    for k = 1:length(tt_b_m2)
        tt_b_m2(k) = (burstsummary_b_m2{1}.bdat.event(k+1,1)-burstsummary_b_m2{1}.bdat.event(k,2))/1000;
    end
    for k = 1:length(tt_b_pc)
        tt_b_pc(k) = (burstsummary_b_pc{1}.bdat.event(k+1,1)-burstsummary_b_pc{1}.bdat.event(k,2))/1000;
    end
    for k = 1:length(tt_u_m1)
        tt_u_m1(k) = (burstsummary_u_m1{1}.bdat.event(k+1,1)-burstsummary_u_m1{1}.bdat.event(k,2))/1000;
    end
    for k = 1:length(tt_u_m2)
        tt_u_m2(k) = (burstsummary_u_m2{1}.bdat.event(k+1,1)-burstsummary_u_m2{1}.bdat.event(k,2))/1000;
    end
    for k = 1:length(tt_u_pc)
        tt_u_pc(k) = (burstsummary_u_pc{1}.bdat.event(k+1,1)-burstsummary_u_pc{1}.bdat.event(k,2))/1000;
    end    
    
%     tt_rh = zeros(length(burstsummary{2}.bdat.event-1)-1,1);
%     for k = 1:length(burstsummary{2}.bdat.event)-1
%         tt_rh(k) = (burstsummary{2}.bdat.event(k+1,1)-burstsummary{2}.bdat.event(k,2))/1000;
%     end
    
%     toe_lh = [toe_lh; tt_lh];
%     toe_rh = [toe_rh; tt_rh];
    tue_b_m1 = [tue_b_m1; tt_b_m1];
    tue_b_m2 = [tue_b_m2; tt_b_m2];
    tue_b_pc = [tue_b_pc; tt_b_pc];
    tue_u_m1 = [tue_u_m1; tt_u_m1];
    tue_u_m2 = [tue_u_m2; tt_u_m2];
    tue_u_pc = [tue_u_pc; tt_u_pc];
%     sub_lhT = [sub_lhT; repmat(subj, length(tt_lh),1)];
%     sub_rhT = [sub_rhT; repmat(subj, length(tt_rh),1)];

    
    % # Subject identifier variable #
    sub_b_m1 = [sub_b_m1; repmat(subj,length(burstsummary_b_m1{1}.bdat.evelen)-1, 1)];
    sub_b_m2 = [sub_b_m2; repmat(subj,length(burstsummary_b_m2{1}.bdat.evelen)-1, 1)];
    sub_b_pc = [sub_b_pc; repmat(subj,length(burstsummary_b_pc{1}.bdat.evelen)-1, 1)];
    sub_u_m1 = [sub_u_m1; repmat(subj,length(burstsummary_u_m1{1}.bdat.evelen)-1, 1)];
    sub_u_m2 = [sub_u_m2; repmat(subj,length(burstsummary_u_m2{1}.bdat.evelen)-1, 1)];
    sub_u_pc = [sub_u_pc; repmat(subj,length(burstsummary_u_pc{1}.bdat.evelen)-1, 1)];

%     sub_rh = [sub_rh; repmat(subj,length(burstsummary{2}.bdat.evelen)-1, 1)];

end

%% Arrange data
% neve    = [nevent_lh; nevent_rh];
% hemiN   = cat(1, repmat({'lh'}, size(subjects)), repmat({'rh'}, size(subjects)));
% subjsN  = [subjects; subjects];
% leneve  = [len_lh; len_rh];
% hemi    = cat(1, repmat({'lh'}, size(len_lh)), repmat({'rh'}, size(len_rh)));
% subjs   = [sub_lh; sub_rh];
% toeeve  = [toe_lh; toe_rh];
% % hemiT   = cat(1, repmat({'lh'}, size(toe_lh)), repmat({'rh'}, size(toe_rh)));
% % subjsT  = [sub_lhT; sub_rhT];
% maxeve  = [max_lh; max_rh];

%% Save
fprintf('Save data... ');
% save(fullfile(dirs.group_path, 'neve_data.mat'), 'neve', 'hemiN','subjsN');
% save(fullfile(dirs.group_path, 'leneve_data.mat'), 'leneve', 'hemi','subjs');
% save(fullfile(dirs.group_path, 'toeeve_data.mat'), 'toeeve', 'hemi','subjs');
% save(fullfile(dirs.group_path, 'maxeve_data.mat'), 'maxeve', 'hemi','subjs');

% Save N data
save(fullfile(dirs.group_path, 'neve_b_m1_data.mat'), 'nevent_b_m1','subjects');
save(fullfile(dirs.group_path, 'neve_b_m2_data.mat'), 'nevent_b_m2','subjects');
save(fullfile(dirs.group_path, 'neve_b_pc_data.mat'), 'nevent_b_pc','subjects');
save(fullfile(dirs.group_path, 'neve_u_m1_data.mat'), 'nevent_u_m1','subjects');
save(fullfile(dirs.group_path, 'neve_u_m2_data.mat'), 'nevent_u_m2','subjects');
save(fullfile(dirs.group_path, 'neve_u_pc_data.mat'), 'nevent_u_pc','subjects');

% Save length data
save(fullfile(dirs.group_path, 'leneve_b_m1.mat'), 'len_b_m1','sub_b_m1');
save(fullfile(dirs.group_path, 'leneve_b_m2.mat'), 'len_b_m2','sub_b_m2');
save(fullfile(dirs.group_path, 'leneve_b_pc.mat'), 'len_b_pc','sub_b_pc');
save(fullfile(dirs.group_path, 'leneve_u_m1.mat'), 'len_u_m1','sub_u_m1');
save(fullfile(dirs.group_path, 'leneve_u_m2.mat'), 'len_u_m2','sub_u_m2');
save(fullfile(dirs.group_path, 'leneve_u_pc.mat'), 'len_u_pc','sub_u_pc');

% Save max data
save(fullfile(dirs.group_path, 'maxeve_b_m1.mat'), 'max_b_m1','sub_b_m1');
save(fullfile(dirs.group_path, 'maxeve_b_m2.mat'), 'max_b_m2','sub_b_m2');
save(fullfile(dirs.group_path, 'maxeve_b_pc.mat'), 'max_b_pc','sub_b_pc');
save(fullfile(dirs.group_path, 'maxeve_u_m1.mat'), 'max_u_m1','sub_u_m1');
save(fullfile(dirs.group_path, 'maxeve_u_m2.mat'), 'max_u_m2','sub_u_m2');
save(fullfile(dirs.group_path, 'maxeve_u_pc.mat'), 'max_u_pc','sub_u_pc');

% Save tue data
save(fullfile(dirs.group_path, 'tueeve_b_m1.mat'), 'tue_b_m1','sub_b_m1');
save(fullfile(dirs.group_path, 'tueeve_b_m2.mat'), 'tue_b_m2','sub_b_m2');
save(fullfile(dirs.group_path, 'tueeve_b_pc.mat'), 'tue_b_pc','sub_b_pc');
save(fullfile(dirs.group_path, 'tueeve_u_m1.mat'), 'tue_u_m1','sub_u_m1');
save(fullfile(dirs.group_path, 'tueeve_u_m2.mat'), 'tue_u_m2','sub_u_m2');
save(fullfile(dirs.group_path, 'tueeve_u_pc.mat'), 'tue_u_pc','sub_u_pc');

disp('done')

%END