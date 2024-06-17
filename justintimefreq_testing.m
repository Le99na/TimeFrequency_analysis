
clear; clc; close all;
% Paths to toolboxes
TOOLBOXES = ['C:\Users\fpupillo\Data\BLIPS_EEG_Analysis\Toolboxes\'];

% Add the EEGLAB toolbox to the path
EEGLAB_TOOLBOX = [TOOLBOXES 'eeglab2022.1/'];
addpath(EEGLAB_TOOLBOX);

% Load EEGLAB
eeglab
%[ALLEEG,EEG,CURRENTSET] = eeglab('rebuild');

% Jusitn's toolbox
RIDDLER_TOOLBOX = [TOOLBOXES 'Riddler_Toolbox\'];

tmpDir = ['C:\Users\fpupillo\Data\BLIPS_EEG_Analysis\TIMEFREQ_Analyses\TIMEFREQ_Results\tmp\'];

%% Load EGI channel location information
chanlocFile = [RIDDLER_TOOLBOX 'EEG/EGI_128_channelLocations.mat'];
chanlocsStruct = load(chanlocFile);
chanlocs = chanlocsStruct.chanlocs;
dataChannelsFile = [RIDDLER_TOOLBOX 'EEG/EGI_128_dataChannels_idxs.mat'];
dataChannelsStruct = load(dataChannelsFile);
dataChannels = dataChannelsStruct.dataChannels;
numDataChannels = length(dataChannels);
data_chanlocs = chanlocs(dataChannels);


                        
EPOCH_RANGE = [-1.25 1.945];
EPOCH_BASELINE = [-0.5 -0.2];

%% Time Frequency Information
lowerBoundFreq = 2;
upperBoundFreq = 40;
% 1/f scaled time frequency estimate
numFreqs = 150;
% Exponent 0.05 based on Voytek, Gazzaley 2015
exponent = 0.05;
% Scaled frequencies to use to estimate power
[freq] = pinkFreqScale(lowerBoundFreq,upperBoundFreq,numFreqs,exponent);
% Wavelet size 5
numCycles = 3;
% Down sampled rate
sampleRate = 200;

% Load EEG file
eegFile = 'C:\Users\fpupillo\Data\BLIPS_EEG_Analysis\TIMEFREQ_Analyses\PRE_TIMEFREQ\GOOD_NEWS\goodnews_negative_events_BLIPS_P003_task.set';

EEG = pop_loadset('filename', eegFile);
time = EEG.times; 
eegFiles = {EEG.data};

fileNamePrefix = ['JR_timefreq_BLIPS_GN_P003'];

% Adjusted baseline period
baselinePeriod = EPOCH_BASELINE - EPOCH_RANGE(1);
% % Required fields for running the time frequency analysis
tfInfo = struct(...
    'sampleRate',sampleRate,...
    'numCycles',numCycles,...
    'fileName',fileNamePrefix,...
    'frequencies',freq,...
    'dataChannels',dataChannels,...
    'baselinePeriod',baselinePeriod,...
    'mirrorFlag',1);
outputDirectory = tmpDir;
[tfOutputFiles,~] = timeFreq_JR_testing(outputDirectory,eegFiles,tfInfo);

dataStruct = load(tfOutputFiles{1});
ersp = dataStruct.tf.data.ersp;
itpc = dataStruct.tf.data.itpc;
            % for dataChanIdx = 1:numDataChannels
            %     erspChan = ersp(:,:,dataChanIdx);
            %     itpcChan = itpc(:,:,dataChanIdx);
            % 
            %     if size(erspChan,1) == 600
            %         erspChan(1:20,:) = [];
            %         itpcChan(1:20,:) = [];
            %     end
            % 
            %     % Prepare output
            %     allCond_ERSP(:,:,dataChanIdx,condIdx)=erspChan';
            %     allCond_ITPC(:,:,dataChanIdx,condIdx)=itpcChan';
            % end
 % output file for time frequency results
timeFreqFile = ['JR_timefreq_BLIPS_GN_P003_results.mat'];
    % Save data out
    timeFreqStruct = struct(...
        'ERSP',ersp,...
        'ITPC',itpc);
    save(timeFreqFile,'-struct','timeFreqStruct');


% Define the time range you want to plot (in milliseconds)
startTime = -500; % Start time in milliseconds
endTime = 1000; % End time in milliseconds

% Find the corresponding indices in the time vector
startIndex = find(time >= startTime, 1);
endIndex = find(time <= endTime, 1, 'last');

% Extract the time vector for the specified time window
time_window = time(startIndex:endIndex);


% Select the time window from the ERSP and ITPC matrices
erspTimeWindow = ersp(startIndex:endIndex, :, :);
itpcTimeWindow = itpc(startIndex:endIndex, :, :);

% Select the channel index
channelIndex = 9; % Change this to the desired channel index (1 to 90)

% Plot the spectrogram for the selected time window and channel
figure;
subplot(2, 1, 1);
imagesc(squeeze(erspTimeWindow(:,:,channelIndex))); % Plot ERSP for the selected time window and channel
title(['ERSP Spectrogram - Channel ' num2str(channelIndex)]);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
colorbar;

subplot(2, 1, 2);
imagesc(squeeze(itpcTimeWindow(:,:,channelIndex))); % Plot ITPC for the selected time window and channel
title(['ITPC Spectrogram - Channel ' num2str(channelIndex)]);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
colorbar;



% Define time window of interest
time_window_start = -500; % in ms
time_window_end = 1000; % in ms

% Find the indices corresponding to the time window
start_idx = find(time >= time_window_start, 1);
end_idx = find(time <= time_window_end, 1, 'last');

% Extract the time vector for the specified time window
time_window = time(start_idx:end_idx);

% Plot the spectrogram for a specific channel
channel_idx = 9; % Replace with the index of the channel you want to plot

%esrp_avgchannel = mean(ersp(151:451, :, :), 3);

% Extract the data for the specified channel and time window
ersp_channel_window = ersp(start_idx:end_idx, :, channel_idx);

% Plot the spectrogram using pcolor
figure;
pcolor(time_window,freq,ersp_channel_window');
shading interp; % Interpolate colors between grid points for smoother appearance
colorbar; % Add color bar for reference
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title(['Spectrogram for Channel ' num2str(channel_idx) ' - Time Window: ' num2str(time_window_start) ' to ' num2str(time_window_end) ' ms']);
