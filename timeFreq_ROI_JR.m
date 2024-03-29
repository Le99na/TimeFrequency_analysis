function [tfOutputFiles] = timeFreq_ROI_JR(outputDirectory,tfInfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CREDITS                                                               %
%                                                                         %
%       By Justin Riddle, PhD                                             %
%       Flavio Frohlich's Lab                                             %
%       University of North Carolina at Chapel Hill                       %
%       Last updated on April 29th, 2021                                  %
%       Email for contact: justin_riddle@med.unc.edu                      %
%                                                                         %
%% Time Frequency Analysis                                                %
% INPUTS:                                                                 %
%   outputDirectory: sting, path to directory for output                  %
%   tfInfo: matlab structure with required and optional inputs            %
%       frequencies - frequencies for analysis                            %
%       mirrorFlag - mirror the data for better accuracy                  %
%       numCycles - double, or vector [1 x numFrequencies]                %
%           specifies how many cycles to be used for time freq analysis   %
%       sampleRate - acquistiion rate of data                             %
%       fileName - customizable name for data output                      %
%       baseline_timeIdxs - timeIdxs for baseline normalization           %
%       COND_NAMES: cellstring, names of each condition                   % 
%       cond_eegFiles: cellstring, files to run time frequency            %
%       ROI_NAMES: cellstring, names of each region of interst            %
%       roi_chanIdxs: cellstring, index to EEG files, check units         %
%                                                                         %
% OUTPUT:                                                                 %
%   tfOutputFiles: time frequency file for each region of interest        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Required fields for running the time frequency analysis
sampleRate = 200;
numCycles = 3; %high temporal resolution: 3-5, high frequency resolution: 5-8
fileNameOut = tfInfo.fileName;
% must choose a baseline period for normalization
baseline_timeIdxs = tfInfo.baseline_timeIdxs;
% condition names and data
COND_NAMES = tfInfo.COND_NAMES;
cond_eegFiles = tfInfo.cond_eegFiles;
% roi names and channel indexes
ROI_NAMES = tfInfo.ROI_NAMES;
roi_chanIdxs = tfInfo.roi_chanIdxs;
% frequencies to estimate time frequency
FREQ = [4,5,6,7,8,9,10,11,12]
numFreq = length(FREQ);

% Number of cycles to estimate tf can change for each frequency
if length(numCycles) ~= numFreq
    numCycles = ones(1,numFreq) * numCycles(1);
end

% Option to mirror the data when computing time-freq
% Defualt is to mirror the data
mirrorData_FLAG = 1;
if isfield(tfInfo,'mirrorFlag')
    mirrorData_FLAG = tfInfo.mirrorFlag;
end

% Assumes that if there are more than one EEG file that these are conditions
% If not by condition, then only run this function on a single EEG file
numConditions = length(cond_eegFiles);
assert(numConditions == length(COND_NAMES));
numROIs = length(roi_chanIdxs);
assert(numROIs == length(ROI_NAMES));

% Load EEG data into workspace
cond_eegData = cell(1,numConditions);
for condIdx = 1:numConditions
    eegFile = cond_eegFiles{condIdx};
    if ischar(eegFile)
        EEG = pop_loadset('filename',eegFile);
        condEEG = EEG.data;
    else
        % Option to just input a matrix of data
        condEEG = cond_eegFiles{condIdx};
    end
    cond_eegData{condIdx} = condEEG;
end

% initialize output variable
tfOutputFiles = cell(1,numROIs);
for roiIdx = 1:numROIs
    roiName = ROI_NAMES{roiIdx};
    chanIdxs = roi_chanIdxs{roiIdx};
    % Use channel number for the specific setup used
    tfOutFile = sprintf('%s%s_%s.mat',...
        outputDirectory,roiName,fileNameOut);
    tfOutputFiles{roiIdx} = tfOutFile;
    if exist(tfOutFile,'file')~=2

        % Print progress
        fprintf('ROI %s\n',roiName);
        
        % Structure to save time frequency results
        clear tfOut;
        
        % loop through the provided conditions
        for condIdx = 1:numConditions
            cond = COND_NAMES{condIdx};
            
            % load data for this condition
            condData = cond_eegData{condIdx};
            
            roi_condData = squeeze(mean(condData(chanIdxs,:,:),1));
            
            % find the save of this condition's data
            [numTimePts,numTrials] = size(roi_condData);
            
            % Run ersp for each channel separately
            amp_condition = NaN(numTimePts,numFreq,numTrials);
            phase_condition = NaN(numTimePts,numFreq,numTrials);
            
            % loop through each trial for this subject
            for trialIdx = 1:numTrials
                % data for only this trial
                trialData = roi_condData(:,trialIdx)';
                
                % If the data is being mirror flipped
                if mirrorData_FLAG
                    % mirror flip the trial data
                    mirrorData = [fliplr(trialData) trialData fliplr(trialData)];
                    mirrorIdxs = (numTimePts+1):(numTimePts*2);
                else
                    % Or just use the original data
                    mirrorData = trialData;
                    mirrorIdxs = 1:numTimePts;
                end
                
                % run time frequency analysis for each frequency
                for freqIdx = 1:numFreq
                    freq = FREQ(freqIdx);
                    % run wavelet based time frequency estimate
                    [phase,amp] = phaseAmp_calculator(sampleRate,freq,mirrorData,numCycles(freqIdx));
                    
                    % Prepare output
                    phase_condition(:,freqIdx,trialIdx) = phase(mirrorIdxs);
                    amp_condition(:,freqIdx,trialIdx) = amp(mirrorIdxs);
                end % freq loop
            end % trial loop

            % Average across trials
            amp_trialAvg = mean(amp_condition,3);
            tfOut.tf.(cond).ersp = amp_trialAvg;
            
            % Calculate inter-trial phase coherence
            itpc = NaN(numTimePts,numFreq);
            for timeIdx = 1:numTimePts
                for freqIdx = 1:numFreq
                    % vector of phase for each trial
                    phaseTrials = squeeze(phase_condition(timeIdx,freqIdx,:));
                    itpc(timeIdx,freqIdx) = abs(mean(exp(1i.*phaseTrials)));
                end % loop freq
            end % loop time
            tfOut.tf.(cond).itpc = itpc;
            
            % Baseline correction
            amp_baselineCorrected = NaN(numTimePts,numFreq);
            itpc_baselineCorrected = NaN(numTimePts,numFreq);

            % Calculate baseline
            amp_baselineAvg = mean(amp_trialAvg(baseline_timeIdxs,:),1);
            itpc_baselineAvg = mean(itpc(baseline_timeIdxs,:),1);

            % Baseline correct all data
            for timeIdx = 1:numTimePts
                timepoint_amp = amp_trialAvg(timeIdx,:);
                timepoint_amp = (timepoint_amp - amp_baselineAvg) ./ amp_baselineAvg;
                amp_baselineCorrected(timeIdx,:) = timepoint_amp;
                
                timepoint_itpc = itpc(timeIdx,:);
                timepoint_itpc = (timepoint_itpc - itpc_baselineAvg) ./ itpc_baselineAvg;
                itpc_baselineCorrected(timeIdx,:) = timepoint_itpc;
            end
            
            tfOut.tf.(cond).ersp_based = amp_baselineCorrected;
            tfOut.tf.(cond).itpc_based = itpc_baselineCorrected;
        end % condition
        
        tfOut.overview = struct(...
            'frequencies',FREQ,...
            'numCycles',numCycles,...
            'sampleRate',sampleRate);
        % save out time frequency results
        save(tfOutFile,'-struct','tfOut');
    end % time frequency file exists
end % data channels loop
end


% Plotting the spectrogram

% Plot the spectrogram
figure;
imagesc(1/200:1/200:300/200, FREQ, amp_trialAvg.');
set(gca, 'YDir', 'normal');
title('Time-Frequency Spectrogram');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;