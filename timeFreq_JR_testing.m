function [tfOutputFiles,dataChannels] = timeFreq_JR(outputDirectory,eegFiles,tfInfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CREDITS                                                               %
%                                                                         %
%       By Justin Riddle, PhD                                             %
%       Flavio Frohlich's Lab                                             %
%       University of North Carolina at Chapel Hill                       %
%       Last updated on February 11th, 2019                               %
%       Email for contact: justin_riddle@med.unc.edu                      %
%                                                                         %
%% Time Frequency Analysis                                                %
% INPUTS:                                                                 %
%   outputDirectory: sting, path to directory for output                  %
%   eegFiles: cellstring, the files on which time frequency will be run   %
%   tfInfo: matlab structure with required and optional inputs            %
%       frequencies - frequencies for analysis                            %
%       mirrorFlag - mirror the data for better accuracy                  %
%       conditions - labels for the eegFiles                              %
%       dataChannels - channels to which analysis is restricted           %
%       numCycles - double, or vector [1 x numFrequencies]                %
%           specifies how many cycles to be used for time freq analysis   %
%       sampleRate - acquistiion rate of data                             %
%       fileName - customizable name for data output                      %
%                                                                         %
% OUTPUT:                                                                 %
%   tfOutputFiles: final time frequency file for each eegFile provided    %
%   dataChannels: this is passed out for use by other functions           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Required fields for running the time frequency analysis
sampleRate = tfInfo.sampleRate;
numCycles = tfInfo.numCycles;
fileNameOut = tfInfo.fileName;
% frequencies to estimate time frequency
FREQ = tfInfo.frequencies;
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
numConditions = length(eegFiles);
if numConditions > 1
    if isfield(tfInfo,'conditions')
        COND = tfInfo.conditions;
    else
        % generate dummy filler condition names
        COND = cell(1,numConditions);
        for condIdx = 1:numConditions
            COND{condIdx} = sprintf('condition%02d',condIdx);
        end
    end
else
    COND = {'data'};
end % if multiple conditions

if isfield(tfInfo,'smallFiles')
    smallFilesFLAG = tfInfo.smallFiles;
else
    smallFilesFLAG = 0;
end

% Load EEG data into workspace
cond_eegData = cell(1,numConditions);
for condIdx = 1:numConditions
    eegFile = eegFiles{condIdx};
    if ischar(eegFile)
        EEG = pop_loadset('filename',eegFile);
        condEEG = EEG.data;
    else
        % Option to just input a matrix of data
        condEEG = eegFiles{condIdx};
    end
    [numChannels,~,~] = size(condEEG);
    cond_eegData{condIdx} = condEEG;
end

% need to establish the number of channels for analysis
if isfield(tfInfo,'dataChannels')
    dataChannels = tfInfo.dataChannels;
else
    dataChannels = 1:numChannels;
end
numDataChannels = length(dataChannels);

% initialize output variable
tfOutputFiles = cell(numDataChannels,1);

if isfield(tfInfo,'baselinePeriod')
    baselinePeriod = tfInfo.baselinePeriod;
    avgAndBaseline = 1;
else
    avgAndBaseline = 0;
end

if avgAndBaseline
    
    % Use channel number for the specific setup used
    tfOutFile = sprintf('%s%s_%03dchannels_%ifreqs_%icycles.mat',...
        outputDirectory,fileNameOut,numDataChannels,numFreq,numCycles(1));
    
    if exist(tfOutFile,'file')==2
        tfOutputFiles = {tfOutFile};
        return;
    end
end


% Loop through data channels
for dataChanIdx = 1:numDataChannels
    channel_name = dataChannels(dataChanIdx);
    
    % Use channel number for the specific setup used
    tfOutFile = sprintf('%s%s_%03dchan_%ifreq_%icycles.mat',...
        outputDirectory,fileNameOut,channel_name,numFreq,numCycles(1));
    tfOutputFiles{dataChanIdx} = tfOutFile;
    
    if exist(tfOutFile,'file')~=2
        
        % Print progress
        fprintf('Channel %03d of %i\n',dataChanIdx,numDataChannels);
        
        % Structure to save time frequency results
        clear tfOut;
        
        % loop through the provided conditions
        for condIdx = 1:numConditions
            cond = COND{condIdx};
            
            % load data for this condition
            condData = cond_eegData{condIdx};
            
            % find the save of this condition's data
            [~,numTimePts,numTrials] = size(condData);
            
            % Run ersp for each channel separately
            ersp_condition = NaN(numTimePts,numFreq,numTrials,2);
            
            % Size of matrix: numTimePoints x numTrials
            epochedChannelData = squeeze(condData(dataChanIdx,:,:));
            
            % loop through each trial for this subject
            for trialIdx = 1:numTrials
                % data for only this trial
                trialData = epochedChannelData(:,trialIdx)';
                
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
                    ersp_condition(:,freqIdx,trialIdx,1) = phase(mirrorIdxs);
                    ersp_condition(:,freqIdx,trialIdx,2) = amp(mirrorIdxs);
                end % freq loop
            end % trial loop
           
            if avgAndBaseline
                ampOnly = squeeze(ersp_condition(:,:,:,2));
                
                % Average across trials
                trialAvg = mean(ampOnly,3);
                
                % Calculate baseline
                baselineStart = max(round(baselinePeriod(1) * sampleRate),1);
                baselineStop = round(baselinePeriod(2) * sampleRate) - 1;
                baselineData = trialAvg(baselineStart:baselineStop,:);
                baselineAvg = mean(baselineData,1);
                
                % Baseline correct all data
                for timeIdx = 1:size(trialAvg,1)
                    timepoint = trialAvg(timeIdx,:);
                    timepoint = (timepoint - baselineAvg) ./ baselineAvg;
                    trialAvg(timeIdx,:) = timepoint;
                end
                
                tfOut.tf.(cond).ersp = trialAvg;
                
                % Calculate inter-trial phase coherence
                phaseOnly = squeeze(ersp_condition(:,:,:,1));
                itpc = NaN(numTimePts,numFreq);
                for timeIdx = 1:numTimePts
                    for freqIdx = 1:numFreq
                        % vector of phase for each trial
                        phaseTrials = squeeze(phaseOnly(timeIdx,freqIdx,:));
                        itpc(timeIdx,freqIdx) = abs(mean(exp(1i.*phaseTrials)));
                    end % loop freq
                end % loop time
                
                tfOut.tf.(cond).itpc = itpc;
                
            else
                tfOut.tf.(cond).amp = ersp_condition(:,:,:,1);
                tfOut.tf.(cond).phase = ersp_condition(:,:,:,2);
            end
        end % condition
        
        tfOut.overview = struct(...
            'frequencies',FREQ,...
            'numCycles',numCycles,...
            'sampleRate',sampleRate);
        % save out time frequency results
        save(tfOutFile,'-struct','tfOut');
    end % time frequency file exists
end % data channels loop

if avgAndBaseline && (~smallFilesFLAG)
    
        
    % Use channel number for the specific setup used
    tfOutFile = sprintf('%s%s_%03dchannels_%ifreqs_%icycles.mat',...
        outputDirectory,fileNameOut,numDataChannels,numFreq,numCycles(1));
    
    clear tfOut
    for condIdx = 1:numConditions
        cond = COND{condIdx};
        
        [~,numTimePts,~] = size(cond_eegData{condIdx});
        
        allERSP = NaN(numTimePts,numFreq,numDataChannels);
        allITPC = NaN(numTimePts,numFreq,numDataChannels);
        for dataChanIdx = 1:numDataChannels
            chanOutFile = tfOutputFiles{dataChanIdx};
            tfOut_chan = load(chanOutFile);
            ersp = tfOut_chan.tf.(cond).ersp;
            itpc = tfOut_chan.tf.(cond).itpc;
            allERSP(:,:,dataChanIdx) = ersp;
            allITPC(:,:,dataChanIdx) = itpc;
        end
        tfOut.tf.(cond).ersp = allERSP;
        tfOut.tf.(cond).itpc = allITPC;
    end
    tfOut.overview = tfOut_chan.overview;
    save(tfOutFile,'-struct','tfOut');
    
    for dataChanIdx = 1:numDataChannels
        chanOutFile = tfOutputFiles{dataChanIdx};
        delete(chanOutFile);
    end
    
    tfOutputFiles = {tfOutFile};
    
end

end % function ends