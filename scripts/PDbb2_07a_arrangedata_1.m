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
b1_lh = nan(length(subjects), 1);
b2_lh = nan(length(subjects), 1);
b3_lh = nan(length(subjects), 1);
b1_rh = nan(length(subjects), 1);
b2_rh = nan(length(subjects), 1); 
b3_rh = nan(length(subjects), 1);

for ii = 1:length(subjects)
    subj = subjects{ii};
    disp(subj);
    load(fullfile(dirs.meg_path, subj, [subj,'-burstsummary.mat']))
        
    nevent_lh(ii) = burstsummary{1}.n_events;
    nevent_rh(ii) = burstsummary{2}.n_events;
    
    b1_lh(ii) = sum(burstsummary{1}.bdat.maxidx < 60000);
    b2_lh(ii) = sum(burstsummary{1}.bdat.maxidx > 60000 & burstsummary{1}.bdat.maxidx < 120000);
    b3_lh(ii) = sum(burstsummary{1}.bdat.maxidx > 120000);
    b1_rh(ii) = sum(burstsummary{2}.bdat.maxidx < 60000);
    b2_rh(ii) = sum(burstsummary{2}.bdat.maxidx > 60000 & burstsummary{2}.bdat.maxidx < 120000);
    b3_rh(ii) = sum(burstsummary{2}.bdat.maxidx > 120000);
end

neve    = [nevent_lh; nevent_rh];
hemi    = cat(1, repmat({'lh'}, size(subjects)), repmat({'rh'}, size(subjects)));
subjs   = [subjects; subjects];
neve_table = table(neve, hemi, subjs);

% Save
save(fullfile(dirs.group_path, 'neve_table.mat'), 'neve_table');

%% Event duration
lenmean1 = zeros(length(subs), 1);   % subjects x number of steps
lenmedn1 = zeros(length(subs), 1);   % subjects x number of steps
lenmean2 = zeros(length(subs), 1);   % subjects x number of steps
lenmedn2 = zeros(length(subs), 1);   % subjects x number of steps

len1 = [];
len2 = [];
sub1 = [];
sub2 = [];

for ii = 1:length(subs)
    load(fullfile(dirs.megDir,subs{ii},'subvals.mat'))
    
    len1 = [len1; burstsummary{1}.bdat.evelen];
    len2 = [len2; burstsummary{2}.bdat.evelen];
    
    sub1 = [sub1; repmat(subs{ii},length(burstsummary{1}.bdat.evelen),1)];
    sub2 = [sub2; repmat(subs{ii},length(burstsummary{2}.bdat.evelen),1)];

    lenmedn1(ii) = median(burstsummary{1}.bdat.evelen);
    lenmean1(ii) = mean(burstsummary{1}.bdat.evelen);
    lenmedn2(ii) = median(burstsummary{2}.bdat.evelen);
    lenmean2(ii) = mean(burstsummary{2}.bdat.evelen);
end

PDlenmn1 = lenmedn1(PDidx);
ctrllenmn1 = lenmedn1(ctrlidx);
PDlenmn2 = lenmedn2(PDidx);
ctrllenmn2 = lenmedn2(ctrlidx);

PDlenavg1 = mean(PDlenmn1);
ctrllenavg1 = mean(ctrllenmn1);
PDlenavg2 = mean(PDlenmn2);
ctrllenavg2 = mean(ctrllenmn2);

PDlensd1 = std(PDlenmn1);
ctrllensd1 = std(ctrllenmn1);
PDlensd2 = std(PDlenmn2);
ctrllensd2 = std(ctrllenmn2);

figure; 
subplot(1,2,1); histogram(PDlenmn1,20); hold on
subplot(1,2,1); histogram(ctrllenmn1,20); hold off
subplot(1,2,2); histogram(PDlenmn2,20); hold on
subplot(1,2,2); histogram(ctrllenmn2,20); hold off

% Save for export
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/lenevent.mat', ...
    'len1', 'len2', 'sub1', 'sub2');

%% Time to next event
toemedn1 = zeros(length(subs), 1);   % subjects x number of steps
toemedn2 = zeros(length(subs), 1);   % subjects x number of steps

alltoe1 = [];
alltoe2 = [];
sub1 = [];
sub2 = [];

for ii = 1:length(subs)
    load(fullfile(dirs.megDir,subs{ii},'subvals.mat'))
    
    toe1 = zeros(length(burstsummary{1}.bdat.event-1)-1,1);
    toe2 = zeros(length(burstsummary{2}.bdat.event-1)-1,1);

    for k = 1:length(burstsummary{1}.bdat.event)-1
        toe1(k) = (burstsummary{1}.bdat.event(k+1,1)-burstsummary{1}.bdat.event(k,2))/1000;
    end
    
    for k = 1:length(burstsummary{2}.bdat.event)-1
        toe2(k) = (burstsummary{2}.bdat.event(k+1,1)-burstsummary{2}.bdat.event(k,2))/1000;
    end
    
    toemedn1(ii) = median(toe1);
    toemedn2(ii) = median(toe2);
    
    alltoe1 = [alltoe1; toe1];
    alltoe2 = [alltoe2; toe2];
    
    sub1 = [sub1; repmat(subs{ii},length(toe1),1)];
    sub2 = [sub2; repmat(subs{ii},length(toe2),1)];
end

PDtoemn1 = toemedn1(PDidx);
ctrltoemn1 = toemedn1(ctrlidx);
PDtoemn2 = toemedn2(PDidx);
ctrltoemn2 = toemedn2(ctrlidx);

figure
subplot(1,2,1); histogram(PDtoemn1,10); hold on
subplot(1,2,1); histogram(ctrltoemn1,10); hold off
subplot(1,2,2); histogram(PDtoemn2,10); hold on
subplot(1,2,2); histogram(ctrltoemn2,10); hold off

PDtoeavg1 = mean(PDtoemn1);
ctrltoeavg1 = mean(ctrltoemn1);
PDtoesd1 = std(PDtoemn1);
ctrltoesd1 = std(ctrltoemn1);

% Save for export
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/toevent.mat', ...
    'alltoe1', 'alltoe2', 'sub1', 'sub2');

%% Max peak in events
maxmedn1 = zeros(length(subs), 1);   % subjects x number of steps
maxmedn2 = zeros(length(subs), 1);   % subjects x number of steps

max1 = [];
max2 = [];
sub1 = [];
sub2 = [];

for ii = 1:length(subs)
    load(fullfile(dirs.megDir,subs{ii},'subvals.mat'))
    
    max1 = [max1; burstsummary{1}.bdat.maxpk];
    max2 = [max2; burstsummary{2}.bdat.maxpk];
    
    sub1 = [sub1; repmat(subs{ii},length(burstsummary{1}.bdat.maxpk),1)];
    sub2 = [sub2; repmat(subs{ii},length(burstsummary{2}.bdat.maxpk),1)];

    maxmedn1(ii) = median(burstsummary{1}.bdat.maxpk);
    maxmedn2(ii) = median(burstsummary{2}.bdat.maxpk);
    
    % Save pkidx for reading and plotting in MNE-Py
    maxidx1 = burstsummary{1}.bdat.maxidx;
    maxidx2 = burstsummary{2}.bdat.maxidx;
    save(fullfile(dirs.megDir,subs{ii},'pkidx1.mat'),'maxidx1')
    save(fullfile(dirs.megDir,subs{ii},'pkidx2.mat'),'maxidx2')
end

PDmaxmn1 = maxmedn1(PDidx);
ctrlmaxmn1 = maxmedn1(ctrlidx);
PDmaxmn2 = maxmedn2(PDidx);
ctrlmaxmn2 = maxmedn2(ctrlidx);

PDlenavg1 = mean(PDmaxmn1);
ctrllenavg1 = mean(ctrlmaxmn1);
PDlenavg2 = mean(PDmaxmn2);
ctrllenavg2 = mean(ctrlmaxmn2);

PDlensd1 = std(PDmaxmn1);
ctrllensd1 = std(ctrlmaxmn1);
PDlensd2 = std(PDmaxmn2);
ctrllensd2 = std(ctrlmaxmn2);

figure; 
subplot(1,2,1); histogram(PDmaxmn1,10); hold on
subplot(1,2,1); histogram(ctrlmaxmn1,10); hold off
subplot(1,2,2); histogram(PDmaxmn2,10); hold on
subplot(1,2,2); histogram(ctrlmaxmn2,10); hold off

save('/home/mikkel/PD_motor/rest_ec/groupanalysis/pkmaxevent.mat', ...
    'max1', 'max2', 'sub1', 'sub2');

% END