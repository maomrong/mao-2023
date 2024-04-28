%%% This script preprocesses data for plotting spectrogram
%%% Run this script before 'Spectrogram.m'
%%% Necessary Functions: TDTbin2mat; locdetrend (Chronux_2_12), eegfilt and epoch (EEGLAB) %%%

%% load the data and define parameters 
close all
clear all
clc
tic

DSamp = 1; % N downsample
HP = 0.1; % high pass filter
LP = 125; % low pass filter
NCha = 36; % LFP1 = 16; LFP2 = 16; EEGs = 2; EMGs = 2;
timelim = [-80 250]; % the secs to preprocess, before and after stim onset
Block = 'blocknumber'; % enter the block number of the data

TDT_times = load('TDT_times.txt'); % stim TDT time in sec, error ±2sec
save_path = ('savepath'); % locate the path to save preprocessed data
data_path = ('datapath'); % locate the path to load data
%%

for i = 1:numel(TDT_times)
    Ti = TDT_times(i)-2;
    Tf = TDT_times(i)+2;
    data = TDTbin2mat(char(data_path), 'T1', Ti, 'T2', Tf,'TYPE',{'epocs'});
    stim = data.epocs.yLSR.onset(1,1);
    Ti = stim+timelim(1);
    Tf = stim+timelim(2);
    data = TDTbin2mat(char(data_path), 'T1', Ti, 'T2', Tf,'TYPE',{'streams'});
    fs = data.streams.EEGs.fs; % 1, sampling frequency rate
    LFPs = data.streams.LFPs.data(1:32,:);
    EEGs = data.streams.EEGs.data(1:2,:);
    EMGw = data.streams.EMGw.data(1,:);
    EMGn = data.streams.EMGs.data(1,:);

    data1 = double(cat(1,LFPs,EEGs,EMGw,EMGn));

    w2=2*LP/fs;
    w=w2;
    [B,A]=butter(3,w,'low');
%     w2=2*HP/fs;
%     w=[w2];
%     [C,D]=butter(3,w,'high');
    
    trial = nan(size(data1));
    MW=[2 0.1]; % detrend parameters 
    for n = 1:size(data1,1) 
       TT = data1(n,:);% select the stim + the channel
       TT = locdetrend(TT,fs,MW); % detrend the data (linear regresion)
       TT = filtfilt(B,A,TT); % 
       %TT = filtfilt(C,D,TT); % 
       TT = downsample(TT,DSamp);
       trial(n,:) = TT;
    end
    srate = fs/DSamp;

    ERPD(:,:,i) = trial;
    stim_times(i) = stim;
end

% times = (timelim(1):1/fs:timelim(2)); 
times =(timelim(1):1/(srate):timelim(2)-(1/(srate))); 

%% Mean ERP and/or zero mean  plot NREM
figure('position',[300 300 600 800])
imagesc(times,(1:1:36),mean(ERPD,3))%,'linestyle','none'
title ('All ERP','FontSize',14)%'interpreter','latex',
ylabel('Channels','FontSize',10)%,'interpreter','latex'
xlabel('Time (sec)','FontSize',10)% 'interpreter','latex',
h = colorbar;
set(get(h,'label'),'string','V');%,'$uV$','interpreter','latex'
xlim (timelim)
caxis ([-1000 1000])
colormap 'jet'

%% Save data and parameters
save ([save_path,'ERPD_',num2str(Block),'.mat'],'Block','srate','ERPD','HP','LP','timelim','stim_times','times','save_path','-v7.3');%
toc