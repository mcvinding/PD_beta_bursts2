%% Example plots
clear all; close all
addpath('/home/mikkel/fieldtrip/fieldtrip/')
ft_defaults
addpath('/home/mikkel/PD_longrest/scripts/')
addpath('/home/mikkel/matlab/export_fig/')

[subs, dirs] = PDbb2_SETUP();
cd(dirs.meg_path);

%% "raw" ROI bp-filterd data with envelope
subID = subs{1};
sub_dir = fullfile(dirs.meg_path, subID);
infiles = find_files(sub_dir, 'roidata_mube_hlbt2.mat');
load(fullfile(sub_dir,infiles{1}))
infiles = find_files(sub_dir, 'roidata2.mat');
load(fullfile(sub_dir,infiles{1}))

% Band-pass to mu band
cfg = [];
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [8 30];
roidata_flt = ft_preprocessing(cfg, roidata);
        
% Plot vars
raw_dat = roidata.trial{1};
flt_dat = roidata_flt.trial{1};
hlb_dat = roidata_mube_hlb.trial{1};
tim = 0:0.001:5;

%% Plot filtered data and enevelope
fig = figure; hold on
set(fig,'Position', [0 0 800 300], 'color','w');
plot(tim, flt_dat(7000:12000), 'b', 'lineWidth',1)
plot(tim, hlb_dat(7000:12000), 'r', 'lineWidth',1.5)
line([0, 5], [3, 3], 'linestyle','--', 'color', 'k')
set(gca, 'LineWidth', 1,'fontweight','bold','fontsize',14, ...
    'YTick', [], 'yticklabel','', 'YColor','none', ...
    'XTick', [0,1,2,3,4,5], 'xticklabel',{0,1,2,3,4,5});
xlabel('Time (s)','fontsize',16);
ylabel('','fontsize',16)

export_fig(fullfile(dirs.figures,'hilbt.png'), '-r500', '-p0.05', '-CMYK')

%% Plot raw
fig = figure; hold on
set(fig,'Position', [0 0 800 300], 'color','w');
plot(tim, raw_dat(7000:12000), 'b', 'lineWidth',1)
set(gca, 'LineWidth', 1,'fontweight','bold','fontsize',14, ...
    'YTick', [], 'yticklabel','', 'YColor','none', ...
    'XTick', [0,1,2,3,4,5], 'xticklabel',{0,1,2,3,4,5});
xlabel('Time (s)','fontsize',16);
ylabel('','fontsize',16)

export_fig(fullfile(dirs.figures,'rawbt.png'), '-r500', '-p0.05', '-CMYK')

%END