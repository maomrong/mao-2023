%%% This script is used to generate the bar graph
%%% It outputs a table of results of gamma power during and post stim, SWA power during and post stim, in normed as well as asb values
%%% Necessary Functions: TDTbin2mat; locdetrend (Chronux_2_12) 

%% load the data and define parameters 
clear all
close all
clc
tic
% To specify for each exp 
save_path = ('savepath'); % locate the path to save preprocessed data
data_path = ('datapath'); % locate the path to load data
first_stim_time = ['firststim']; % first stim TDT time in sec, error ±2sec
stim_dur = [0 19]; % stim duration in secs
post_stim = [19 39]; % post-stim time window for SWA quantification, in sec from stim at 0
arrange = [4 6 5 2]; % pick and rearrange channels for power quantification, based on its surgery setup: local LFP-> mid EEG-> distal EEG-> distal LFP. 1=LFP f d, 2=v, 3=LFP p d, 4=v, 5=EEG f, 6=EEG p

LFPs_list = [3 14 19 30]; % TAKING 4 SPECIFIC CHANNELS
timelim = [-200 200];
baseline = [-50 -1];
delta = [0.5 4];
gamma = [70 100];
stim_1sec = [0 1]; % to explore the power change in the 1st sec of stim 
stim_8sec = [0 8]; % to explore the power change in the first 8sec of stim 
name_spect = ('spects');
name_result = ('power_results');

MW = [5 0.2]; % detrend parameters
DSamp = 1; % N downsample
LP = 100; % low pass filter

%% load data
% load the accurate 1st stim onset time
Ti = first_stim_time - 2;
Tf = first_stim_time + 2;
data = TDTbin2mat(char(data_path), 'T1', Ti, 'T2', Tf,'TYPE',{'epocs'}); % load the first pulse
stim = data.epocs.yLSR.onset(1,1);

% load data around it
Ti = stim + timelim(1);
Tf = stim + timelim(2);
data = TDTbin2mat(char(data_path), 'T1', Ti, 'T2', Tf,'TYPE',{'streams'},'STORE',{'LFPs','EEGs','EMGw','EMGs'}); % load the data surranding the first pulse
fs = data.streams.LFPs.fs; % 1, sampling frequency rate
LFPs = data.streams.LFPs.data(LFPs_list,:); % LFPs 
EEGs = data.streams.EEGs.data(1:2,:); % EEGs
EMGw = data.streams.EMGw.data(1,:); % EMG whiskers
EMGs = data.streams.EMGs.data(1,:); % EMG neck
data1 = double(cat(1,LFPs,EEGs,EMGw,EMGs)); % cancatenate all the data

% extract all the stim times 
Ti = stim+timelim(1);
Tf = stim+timelim(2);
data = TDTbin2mat(char(data_path), 'T1', Ti, 'T2', Tf,'TYPE',{'epocs'});
stim_onset = data.epocs.yLSR.onset;
stim_onset = stim_onset(1:2:numel(stim_onset)); % deal with TDT data
stim_offset = data.epocs.yLSR.offset;
stim_offset = stim_offset(1:2:numel(stim_offset)); % deal with TDT data

%% preprocess data
signals = nan(size(data1));
w2 = 2*LP/fs;
w = w2;
[B,A]=butter(3,w,'low');
for n = 1:size(data1,1) 
   TT = data1(n,:);% select the stim + the channel
   TT = locdetrend(TT,fs,MW); % detrend the data (linear reggresion)
   TT = filtfilt(B,A,TT); % 
   TT = downsample(TT,DSamp);
   signals(n,:) = TT;
end
srate = fs/DSamp;
time = timelim(1) + 1/srate:1/srate:timelim(2); % lay out time points of each data sample around stim, of the sampling rate

%% compute and plot raw spectrograms
figure()
t = tiledlayout(8,1);
set(gcf,'Position',[20 20 1000 1200]) %,'BackgroundColor',[0, 0, 0]
ylabel(t, 'Frequency (Hz)','FontName', 'Times New Roman','FontSize',10)

% the first channel (frontal dorsal)
nexttile(1)
[wt1,f] = cwt(signals(1,:),srate); % compute spectrogram
surface(time,f,abs(wt1)) % plot spectrogram
axis tight
shading flat
%xlabel('Time (s)')
ylabel('Frontal dor')
set(gca,'yscale','log','xtick',[])
set(gca, 'FontName', 'Times New Roman','FontSize',10)
h = colorbar;
set(get(h,'label'),'string','Magnitude');%,'$uV$','interpreter','latex'
ylim([0.5 125])
xlim([-10 10])
yticks([0.5 4.5 30 125])

nexttile(2)
[wt2,f] = cwt(signals(2,:),srate);
surface(time,f,abs(wt2))
axis tight
shading flat
%xlabel('Time (s)')
ylabel('Frontal ven')
set(gca,'yscale','log','xtick',[])
set(gca, 'FontName', 'Times New Roman','FontSize',10)
h = colorbar;
set(get(h,'label'),'string','Magnitude');%,'$uV$','interpreter','latex'
ylim([0.5 110])
xlim([-5 5])
yticks([0.5 4.5 30 125])

nexttile(3)
[wt3,f] = cwt(signals(3,:),srate);
surface(time,f,abs(wt3))
axis tight
shading flat
%xlabel('Time (s)')
ylabel('Parietal dor')
set(gca,'yscale','log','xtick',[])
set(gca, 'FontName', 'Times New Roman','FontSize',10)
h = colorbar;
set(get(h,'label'),'string','Magnitude');%,'$uV$','interpreter','latex'
ylim([0.5 125])
xlim([-10 10])
yticks([0.5 4.5 30 125])

nexttile(4)
[wt4,f] = cwt(signals(4,:),srate);
surface(time,f,abs(wt4))
axis tight
shading flat
%xlabel('Time (s)')
ylabel('Parietal ven')
set(gca,'yscale','log','xtick',[])
set(gca, 'FontName', 'Times New Roman','FontSize',10)
h = colorbar;
set(get(h,'label'),'string','Magnitude');%,'$uV$','interpreter','latex'
ylim([0.5 125])
xlim([-10 10])
yticks([0.5 4.5 30 125])

nexttile(5)
[wt5,f] = cwt(signals(5,:),srate);
surface(time,f,abs(wt5))
axis tight
shading flat
%xlabel('Time (s)')
ylabel('EEG frontal')
set(gca,'yscale','log','xtick',[])
set(gca, 'FontName', 'Times New Roman','FontSize',10)
h = colorbar;
set(get(h,'label'),'string','Magnitude');%,'$uV$','interpreter','latex'
ylim([0.5 125])
xlim([-10 10])
yticks([0.5 4.5 30 125])

nexttile(6)
[wt6,f] = cwt(signals(6,:),srate);
surface(time,f,abs(wt6))
axis tight
shading flat
%xlabel('Time (s)')
ylabel('EEG parietal')
set(gca,'yscale','log','xtick',[])
set(gca, 'FontName', 'Times New Roman','FontSize',10)
h = colorbar;
set(get(h,'label'),'string','Magnitude');%,'$uV$','interpreter','latex'
ylim([0.5 125])
xlim([-10 10])
yticks([0.5 4.5 30 125])

nexttile(7)
[wt7,f] = cwt(signals(7,:),srate);
surface(time,f,abs(wt7))
axis tight
shading flat
%xlabel('Time (s)')
ylabel('EMG wisk')
set(gca,'yscale','log','xtick',[])
set(gca, 'FontName', 'Times New Roman','FontSize',10)
h = colorbar;
set(get(h,'label'),'string','Magnitude');%,'$uV$','interpreter','latex'
ylim([0.5 125])
xlim([-10 10])
yticks([0.5 4.5 30 125])

nexttile(8)
[wt8,f] = cwt(signals(8,:),srate);
surface(time,f,abs(wt8))
axis tight
shading flat
xlabel('Time (s)')
ylabel('EMG neck')
set(gca,'yscale','log')
set(gca, 'FontName', 'Times New Roman','FontSize',10)
h = colorbar;
set(get(h,'label'),'string','Magnitude');%,'$uV$','interpreter','latex'
ylim([0.5 125])
xlim([-10 10])
yticks([0.5 4.5 30 125])

%% extract index 
idx_post = time>post_stim(1) & time<post_stim(2);
idx_1sec = time>stim_1sec(1) & time<stim_1sec(2);
idx_8sec = time>stim_8sec(1) & time<stim_8sec(2);
idx_stim = (time>stim_dur(1) & time<stim_dur(2));
idx_base = (time>baseline(1) & time<baseline(2)); % baseline index

idx_on = logical(zeros(size(time)));
for i = 1:size(stim_onset,1)
    %idx=(time>stim_onset(i,:)-stim & time<stim_offset(i,:)-stim); % baseline index
    idx=(time>(stim_onset(i,:)-stim)+((stim_offset(i,:)-stim_onset(i,:))/4)) & (time<(stim_offset(i,:)-stim)-((stim_offset(i,:)-stim_onset(i,:))/4)); % taking only the centeral half 
    idx_on = (idx_on | idx); % laser on index
end

idx_off = logical(zeros(size(time)));
for i = 1:size(stim_offset,1)-1
    idx=(time>(stim_offset(i,:)-stim)+((stim_onset(i+1,:)-stim_offset(i,:))/4)) & (time<(stim_onset(i+1,:)-stim)-((stim_onset(i+1,:)-stim_offset(i,:))/4)); % taking only the centeral half 
    idx_off = (idx_off | idx); % laser off index
end
% index of post-stim, excluding immed after stim termination
idx_no = logical(zeros(size(time)));
idx = (time>(stim_offset(numel(stim_offset),:)-stim)+((stim_onset(i+1,:)-stim_offset(i,:))/2)); % EXCLUDING THE LENGTH OF HALF-OFF 
idx_no = (idx_no | idx);

idx_delta = (f>=delta(1) & f<=delta(2)); % delta index
idx_gamma = (f>=gamma(1) & f<=gamma(2)); % gamma index

%% MAKE AND CHECK INDEX
figure(100)
scatter(stim_onset-stim, 1,'r')
hold on
scatter(stim_offset-stim,1,'k')
figure(2)
plot(idx_on)
hold on
plot(idx_off)
hold on
plot(idx_no)
xlim([50000 65000])

% stim-ON during the 1st sec
on_1sec = idx_1sec & idx_on;
on_8sec = idx_8sec & idx_on;
post_no = idx_post & idx_no;

% TO CHECK THE condition index
figure(5)
plot(idx_base)
hold on
plot(idx_on)
hold on
plot(on_1sec)
hold on
plot(on_8sec)
hold on 
plot(post_no)
xlim([35000 85000])

%% baseline extraction
base_power = cat(2,mean(abs(wt1(:,idx_base)),2),...
    mean(abs(wt2(:,idx_base)),2),...
    mean(abs(wt3(:,idx_base)),2),...
    mean(abs(wt4(:,idx_base)),2),...
    mean(abs(wt5(:,idx_base)),2),...
    mean(abs(wt6(:,idx_base)),2),...
    mean(abs(wt7(:,idx_base)),2),...
    mean(abs(wt8(:,idx_base)),2));

base_delta = mean(base_power(idx_delta,:),1);
base_gamma = mean(base_power(idx_gamma,:),1);
result_delta = base_delta;
result_gamma = base_gamma;

%% stim-ON across all beh states
stim_on_power = cat(2,mean(abs(wt1(:,idx_on)),2),...
    mean(abs(wt2(:,idx_on)),2),...
    mean(abs(wt3(:,idx_on)),2),...
    mean(abs(wt4(:,idx_on)),2),...
    mean(abs(wt5(:,idx_on)),2),...
    mean(abs(wt6(:,idx_on)),2),...
    mean(abs(wt7(:,idx_on)),2),...
    mean(abs(wt8(:,idx_on)),2));
stim_on_gamma = mean(stim_on_power(idx_gamma,:),1);
result_gamma = cat(2,result_gamma,stim_on_gamma);
stim_on_delta = mean(stim_on_power(idx_delta,:),1); 
result_delta = cat(2,result_delta,stim_on_delta);

%% stim-ON during the first 1 sec 
stim_1sec_power = cat(2,mean(abs(wt1(:,on_1sec)),2),...
    mean(abs(wt2(:,on_1sec)),2),...
    mean(abs(wt3(:,on_1sec)),2),...
    mean(abs(wt4(:,on_1sec)),2),...
    mean(abs(wt5(:,on_1sec)),2),...
    mean(abs(wt6(:,on_1sec)),2),...
    mean(abs(wt7(:,on_1sec)),2),...
    mean(abs(wt8(:,on_1sec)),2));
stim_1sec_gamma = mean(stim_1sec_power(idx_gamma,:),1);
result_gamma = cat(2,result_gamma,stim_1sec_gamma);
stim_1sec_delta = mean(stim_1sec_power(idx_delta,:),1); 
result_delta = cat(2,result_delta,stim_1sec_delta); 

%% stim-ON during the first 8 sec 
stim_8sec_power = cat(2,mean(abs(wt1(:,on_8sec)),2),...
    mean(abs(wt2(:,on_8sec)),2),...
    mean(abs(wt3(:,on_8sec)),2),...
    mean(abs(wt4(:,on_8sec)),2),...
    mean(abs(wt5(:,on_8sec)),2),...
    mean(abs(wt6(:,on_8sec)),2),...
    mean(abs(wt7(:,on_8sec)),2),...
    mean(abs(wt8(:,on_8sec)),2));
stim_8sec_gamma = mean(stim_8sec_power(idx_gamma,:),1);
result_gamma = cat(2,result_gamma,stim_8sec_gamma);
stim_8sec_delta = mean(stim_8sec_power(idx_delta,:),1);
result_delta = cat(2,result_delta,stim_8sec_delta); 

%% GAMMA and SWA of the post stim window
post_no_power = cat(2,mean(abs(wt1(:,post_no)),2),...
    mean(abs(wt2(:,post_no)),2),...
    mean(abs(wt3(:,post_no)),2),...
    mean(abs(wt4(:,post_no)),2),...
    mean(abs(wt5(:,post_no)),2),...
    mean(abs(wt6(:,post_no)),2),...
    mean(abs(wt7(:,post_no)),2),...
    mean(abs(wt8(:,post_no)),2));
post_no_gamma = mean(post_no_power(idx_gamma,:),1); 
result_gamma = cat(2,result_gamma,post_no_gamma); 
post_no_delta = mean(post_no_power(idx_delta,:),1);
result_delta = cat(2,result_delta,post_no_delta);

%% process gamma results
for i = 1:6 % channel index. 1=LFP f d, 2=v, 3=LFP p d, 4=v, 5=EEG f, 6=EEG p. had 'i = [2 4 5 6]'
    idx = [i:8:numel(result_gamma)]; % '8' b/c 'result_gamma' is laid out as 8 channels per condition
    abs_gamma(i,:) = result_gamma(1,idx);  % the power of channel i across all conditions
    lognorm_gamma(i,:) = log10(abs_gamma(i,:)./abs_gamma(i,1)); % take its lognorm
end
% process delta results
for i = 1:6 % channel index. 1=LFP f d, 2=v, 3=LFP p d, 4=v, 5=EEG f, 6=EEG p]. had 'i = [2 4 5 6]'
    idx = [i:8:numel(result_delta)]; % '8' b/c 'result_gamma' is laid out as 8 channels per condition
    abs_delta(i,:) = result_delta(1,idx);  % the power of channel i across all conditions
    lognorm_delta(i,:) = log10(abs_delta(i,:)./abs_delta(i,1)); % take its lognorm.
end
table1 = [abs_gamma abs_delta lognorm_gamma lognorm_delta];% 2024, combine the 4 tables
table2 = table1(:,[1:10 12:15 17:20]); % 2024 NEW, remove baseline zeros
table = table2(arrange,:); % PICK and rearrange channels: local LFP-> mid EEG-> distal EEG-> distal LFP. 1=LFP f d, 2=v, 3=LFP p d, 4=v, 5=EEG f, 6=EEG p
% to visualize the bar graph
figure
bar(table)
print(gcf,[save_path, 'bar', '.jpg'],'-djpeg','-r300'); 
save ([save_path 'bar'],'table')
toc