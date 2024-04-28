%%% This script plots spectrogram
%%% Run 'ERPD_preprocessing.m' before this script
%%% Necessary Functions: TDTbin2mat; locdetrend (Chronux_2_12), eegfilt and epoch (EEGLAB) 

%% load the data and define parameters 
close all
clear all
clc
tic
load('ERPDpath') % path of the output of 'ERPD_preprocessing.m'
channel = 30; % channel number 
trial = 1; % trial number
ToShow = [-20 150]; % the secs to show in the figure, before and after stim onset

%%
data = ERPD(channel,:,trial);
cycles = 7; % using the wavelet method: the number of cycles in each Morlet wavelet, held constant across frequencies.
epochlim = timelim*1000; % timelim in miliseconds 
frames = numel(data(1,:));% number of frames (samples) in each trial
[ersp,itc,powbase,timess,freqs] = newtimef(data(:,:), frames, epochlim, srate, cycles,...
    'baseline',[-50000 -2000],...
    'basenorm','off',...
    'trialbase','off',...
    'freqs',[0.1 125],...
    'timesout',500,...
    'padratio',2,...
    'plotersp','off',...
    'plotitc','off');

%% plot 
figure()
set(gca, 'Position', [0.02,0.1,1,0.82])
t = tiledlayout(10,1,'TileSpacing','none');
set(gcf,'Position',[50 50 1000 200])
ax1 = nexttile([8 1]);
time = timess(:,:)./1000; % in second
contourf(time,freqs,(ersp),20,'linestyle','none','LineWidth',0.01,'edgecolor','none')% SPECIFY THE # ~20-100 for resolution
name = sprintf('%s   Channel #%d   Trial #%d', Block, channel, trial);
title(name,'FontSize',10)
hold on 

line([0 0], [0.5 100],'LineStyle','--','Color','black') % stim onset
line([19 19], [0.5 100],'LineStyle','--','Color','black') % stim offset 
line([5 5], [0.5 100],'LineStyle','--','Color','black') %  move starts
line([10 10], [0.5 100],'LineStyle','--','Color','black') %  aRORR starts
line([23 23], [0.5 100],'LineStyle','--','Color','black') %  RORR starts
line([110 110], [0.5 100],'LineStyle','--','Color','black') % LORR starts  
line(ToShow, [4.5 4.5],'LineStyle','--','Color','black')
line(ToShow, [30 30],'LineStyle','--','Color','black')

ylabel('Frequency (Hz)','FontSize',10)%,'interpreter','latex'
set(gca,'Ydir','normal','xticklabel',[])%'Ydir','reverse',,'FontSize',10'TickLength',[0.03 0.03],'Linewidth',1
ylim ([0.5 125])
xlim (ToShow)
set(gca,'Yscale','log')
caxis ([-25 25])
colormap 'jet'
box off
set(get(gca,'XAxis'),'visible','off')
ax = gca;
ax.YAxis.Exponent = 0;
yticks([0.5 4.5 30 125])

%%
ax2 = nexttile([1 1]);
EMGw = abs(ERPD(35,:,trial)).*1000;

plot(times,EMGw')
xlim (ToShow)
ylim ([0 40000])
set(gca,'xticklabel',[],'yticklabel',[])%'Ydir','reverse',,'FontSize',10'TickLength',[0.03 0.03],'Linewidth',1
set(get(gca,'XAxis'),'visible','off')
box off

ax3 = nexttile([1 1]);

EMGn = abs(ERPD(36,:,trial)).*1000;
plot(times,EMGn')
xlim (ToShow)
ylim ([0 200000])
ylabel('EMGs (mV)','FontSize',10)%,'interpreter','latex'
box off

xlabel('Time (sec)','FontSize',10,'FontName','Times New Roman')% 'interpreter','latex',
%%
print(gcf,[save_path, '[170sec] ', name],'-dtiffn','-r300');
toc