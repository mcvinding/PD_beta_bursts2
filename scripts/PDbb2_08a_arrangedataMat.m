%% Get  summaries of burst and arrange data
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

%% Init arrays
% N events
nevent_u_m2 = nan(length(subjects), 1);

% Event length
len_u_m2 = [];
% Mak event peak
max_u_m2 = [];
% Time until event
tue_u_m2 = [];
% Subject identifier
sub_u_m2 = [];

%% Run
for ii = 1:length(subjects)
    subj = subjects{ii};
    fprintf('Reading subj %s (%i of %i)\n', subj, ii, length(subjects))
    
    % Load data
    infname_u_m2 = fullfile(dirs.meg_path, subj,[subj,'-u_m2-burst2.mat']);
    fprintf('Load data... ')
    load(infname_u_m2)
    disp('done')

    % # N events #
    nevent_u_m2(ii) = burstsummary_u_m2.n_events;
    
    % # Event duration #
    len_u_m2 = [len_u_m2; burstsummary_u_m2.bdat.evelen(1:end-1)];
  
    % # Max peak #
    max_u_m2 = [max_u_m2; burstsummary_u_m2.bdat.maxpk(1:end-1)];

    % # Time between events #
    tt_u_m2 = zeros(length(burstsummary_u_m2.bdat.event-1)-1,1);  

    for k = 1:length(tt_u_m2)
        tt_u_m2(k) = (burstsummary_u_m2.bdat.event(k+1,1)-burstsummary_u_m2.bdat.event(k,2))/1000;
    end
   
    tue_u_m2 = [tue_u_m2; tt_u_m2];
    
    % # Subject identifier variable #
    sub_u_m2 = [sub_u_m2; repmat(subj,length(burstsummary_u_m2.bdat.evelen)-1, 1)];
end

%% Save
fprintf('Save data... ');
% Save N data
save(fullfile(dirs.group_path, 'neve_u_m2_data2.mat'), 'nevent_u_m2','subjects');
% Save length data
save(fullfile(dirs.group_path, 'leneve_u_m22.mat'), 'len_u_m2','sub_u_m2');
% Save max data
save(fullfile(dirs.group_path, 'maxeve_u_m22.mat'), 'max_u_m2','sub_u_m2');
% Save tue data
save(fullfile(dirs.group_path, 'tueeve_u_m22.mat'), 'tue_u_m2','sub_u_m2');
disp('done')
%END
