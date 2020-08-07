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

for ii = 1:length(subjects)
    subj = subjects{ii};
    disp(subj);
    load(fullfile(dirs.meg_path, subj, [subj,'-burstsummary.mat']))
      
    % N events
    nevent_lh(ii) = burstsummary{1}.n_events;
    nevent_rh(ii) = burstsummary{2}.n_events;
    
    % Event duration
end

% Arrange data
neve    = [nevent_lh; nevent_rh];
hemi    = cat(1, repmat({'lh'}, size(subjects)), repmat({'rh'}, size(subjects)));
subjs   = [subjects; subjects];
% neve_table = table(neve, hemi, subjs);

% Save
% save(fullfile(dirs.group_path, 'neve_table.mat'), 'neve_table');
save(fullfile(dirs.group_path, 'neve_data.mat'), 'neve', 'hemi','subjs');

%% Event duration
clear hemi subjs
len_lh = [];
len_rh = [];
sub_lh = [];
sub_rh = [];

for ii = 1:length(subjects)
    subj = subjects{ii};
    disp(subj);
    load(fullfile(dirs.meg_path, subj, [subj,'-burstsummary.mat']))
    
    len_lh = [len_lh; burstsummary{1}.bdat.evelen];
    len_rh = [len_rh; burstsummary{2}.bdat.evelen];
    
    sub_lh = [sub_lh; repmat(subj,length(burstsummary{1}.bdat.evelen),1)];
    sub_rh = [sub_rh; repmat(subj,length(burstsummary{2}.bdat.evelen),1)];

end

% Arrange data
leneve  = [len_lh; len_rh];
hemi    = cat(1, repmat({'lh'}, size(len_lh)), repmat({'rh'}, size(len_rh)));
subjs   = [sub_lh; sub_rh];

% Save for export
save(fullfile(dirs.group_path, 'leneve_data.mat'), 'leneve', 'hemi','subjs');

%% Time to next event
clear hemi subjs
toe_lh = [];
toe_rh = [];
sub_lh = [];
sub_rh = [];

for ii = 1:length(subjects)
    subj = subjects{ii};
    disp(subj);
    load (fullfile(dirs.meg_path, subj, [subj,'-burstsummary.mat']))
    
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
    
    sub_lh = [sub_lh; repmat(subj, length(tt_lh),1)];
    sub_rh = [sub_rh; repmat(subj, length(tt_rh),1)];
end

% Arrange data
toeeve  = [toe_lh; toe_rh];
hemi    = cat(1, repmat({'lh'}, size(toe_lh)), repmat({'rh'}, size(toe_rh)));
subjs   = [sub_lh; sub_rh];

% Save for export
save(fullfile(dirs.group_path, 'toeeve_data.mat'), 'toeeve', 'hemi','subjs');

%% Max peak in events
clear hemi subjs
max_lh = [];
max_rh = [];
sub_lh = [];
sub_rh = [];

for ii = 1:length(subjects)
    subj = subjects{ii};
    disp(subj);
    load (fullfile(dirs.meg_path, subj, [subj,'-burstsummary.mat']))
    
    max_lh = [max_lh; burstsummary{1}.bdat.maxpk];
    max_rh = [max_rh; burstsummary{2}.bdat.maxpk];
    
    sub_lh = [sub_lh; repmat(subj,length(burstsummary{1}.bdat.maxpk),1)];
    sub_rh = [sub_rh; repmat(subj,length(burstsummary{2}.bdat.maxpk),1)];

    % Save pkidx for reading and plotting
    maxidx_lh = burstsummary{1}.bdat.maxidx;
    maxidx_rh = burstsummary{2}.bdat.maxidx;
    save(fullfile(dirs.meg_path, subj, 'pkidx_lh.mat'), 'maxidx_lh')
    save(fullfile(dirs.meg_path, subj, 'pkidx_rh.mat'), 'maxidx_rh')
end

maxeve  = [max_lh; max_rh];
hemi    = cat(1, repmat({'lh'}, size(max_lh)), repmat({'rh'}, size(max_rh)));
subjs   = [sub_lh; sub_rh];

% Save for export
save(fullfile(dirs.group_path, 'maxeve_data.mat'), 'maxeve', 'hemi','subjs');

% END