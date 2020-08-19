% Get various summaries of beta burst (explorative)
clear all
close all
addpath('/home/mikkel/fieldtrip/fieldtrip/')
ft_defaults
addpath('/home/mikkel/beta_bursts/functions')
addpath('/home/mikkel/PD_longrest/scripts/')
[subjects, dirs] = PDbb2_SETUP();

%% N events
nevent_lh = nan(length(subjects), 1);
nevent_rh = nan(length(subjects), 1);
len_lh = [];
len_rh = [];
sub_lh = [];
sub_rh = [];
toe_lh = [];
toe_rh = [];
% sub_lhT = [];
% sub_rhT = [];
max_lh = [];
max_rh = [];

for ii = 1:length(subjects)
    subj = subjects{ii};
    fprintf('Reading subj %s (%i of %i)\n', subj, ii, length(subjects))
    load(fullfile(dirs.meg_path, subj, [subj,'-burstsummary.mat']))
      
    % N events
    nevent_lh(ii) = burstsummary{1}.n_events;
    nevent_rh(ii) = burstsummary{2}.n_events;
    
    % Event duration
        % Event duration
    len_lh = [len_lh; burstsummary{1}.bdat.evelen(1:end-1)];
    len_rh = [len_rh; burstsummary{2}.bdat.evelen(1:end-1)];
    
    sub_lh = [sub_lh; repmat(subj,length(burstsummary{1}.bdat.evelen)-1, 1)];
    sub_rh = [sub_rh; repmat(subj,length(burstsummary{2}.bdat.evelen)-1, 1)];
    
    % Time between events
    tt_lh = zeros(length(burstsummary{1}.bdat.event-1)-1,1);
    tt_rh = zeros(length(burstsummary{2}.bdat.event-1)-1,1);

    for k = 1:length(burstsummary{1}.bdat.event)-1
        tt_lh(k) = (burstsummary{1}.bdat.event(k+1,1)-burstsummary{1}.bdat.event(k,2))/1000;
    end
    
    for k = 1:length(burstsummary{2}.bdat.event)-1
        tt_rh(k) = (burstsummary{2}.bdat.event(k+1,1)-burstsummary{2}.bdat.event(k,2))/1000;
    end
    
    toe_lh = [toe_lh; tt_lh];
    toe_rh = [toe_rh; tt_rh];
    
%     sub_lhT = [sub_lhT; repmat(subj, length(tt_lh),1)];
%     sub_rhT = [sub_rhT; repmat(subj, length(tt_rh),1)];
    
    % Max peak
    max_lh = [max_lh; burstsummary{1}.bdat.maxpk(1:end-1)];
    max_rh = [max_rh; burstsummary{2}.bdat.maxpk(1:end-1)];

end

%% Arrange data
neve    = [nevent_lh; nevent_rh];
hemiN   = cat(1, repmat({'lh'}, size(subjects)), repmat({'rh'}, size(subjects)));
subjsN  = [subjects; subjects];
leneve  = [len_lh; len_rh];
hemi    = cat(1, repmat({'lh'}, size(len_lh)), repmat({'rh'}, size(len_rh)));
subjs   = [sub_lh; sub_rh];
toeeve  = [toe_lh; toe_rh];
% hemiT   = cat(1, repmat({'lh'}, size(toe_lh)), repmat({'rh'}, size(toe_rh)));
% subjsT  = [sub_lhT; sub_rhT];
maxeve  = [max_lh; max_rh];

% Save
save(fullfile(dirs.group_path, 'neve_data.mat'), 'neve', 'hemiN','subjsN');
save(fullfile(dirs.group_path, 'leneve_data.mat'), 'leneve', 'hemi','subjs');
save(fullfile(dirs.group_path, 'toeeve_data.mat'), 'toeeve', 'hemi','subjs');
save(fullfile(dirs.group_path, 'maxeve_data.mat'), 'maxeve', 'hemi','subjs');
disp('done')

%END