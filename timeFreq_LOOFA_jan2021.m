% Goal of the script - run single subject analysis
%   1. Run time frequency analysis
%   2. Run phase-amplitude coupling

clc; clear all; close all;

%% Variables that are liable to change
preprocVersion = 'june2019';
timeFreqVersion = 'jan2021';

%% Paths to Toolboxes
% Root is specific to the computer you are working on
COMPUTER = '/Users/justinri/';
ROOT = [COMPUTER 'Dropbox (Frohlich Lab)/'];
CODE = [ROOT 'Frohlich Lab Team Folder/Codebase/CodeJustin/'];
SCRIPTS = [CODE 'SSDE/'];
TOOLBOXES = [CODE 'Toolboxes/'];
% Toolboxes
EEGLAB_TOOLBOX = [TOOLBOXES 'eeglab14_0_0b/'];
addpath(EEGLAB_TOOLBOX);
[ALLEEG,EEG,CURRENTSET] = eeglab('rebuild');
pop_editoptions( 'option_storedisk', 1,...
    'option_savetwofiles', 0,  'option_computeica', 1, 'option_scaleicarms', 0);
RIDDLER_TOOLBOX = [CODE 'Riddler_Toolbox/'];
addpath(genpath(RIDDLER_TOOLBOX));

%% Experiment paths
DATA = [ROOT 'HumanStudies/SSDE_EEG/'];
PROC_EEG = [DATA 'Preproc_' preprocVersion '/'];
TF_EEG = [DATA 'TimeFreq/'];
mkdir_JR(TF_EEG);
SINGLE_TF = [TF_EEG 'Individual_' timeFreqVersion '/'];
mkdir_JR(SINGLE_TF);
%RAW_TASK = 'Z:/ClinicalTrials/Data/SSDE/Matlab Files/';
RAW_TASK = [DATA 'RawBehavior/'];

%% Subject Information
% Groups: depressed or healthy control
GROUPS = {'Depressed','Control'};
numGroups = length(GROUPS);

% Depressed subject and their age/sex matched healthy control
subjectInfoFile = 'subjectInfo.mat';
SUBJECT_INFO_struct = load(subjectInfoFile);
SUBJECT_INFO = SUBJECT_INFO_struct.SUBJECT_INFO;
numSubPerGroup = size(SUBJECT_INFO,1);
% Subject identification number
SUBJECTS = {SUBJECT_INFO{:,1} SUBJECT_INFO{:,4}};
% Depressed or control
GROUP_LABEL = [ones(1,numSubPerGroup) ones(1,numSubPerGroup)*2];
GROUP_LABEL = GROUPS(GROUP_LABEL);
% Number of subjects for the analysis
numSub = length(SUBJECTS);

%% Load EGI channel location information
chanlocFile = [RIDDLER_TOOLBOX 'EEG/EGI_128_channelLocations.mat'];
chanlocsStruct = load(chanlocFile);
chanlocs = chanlocsStruct.chanlocs;
dataChannelsFile = [RIDDLER_TOOLBOX 'EEG/EGI_128_dataChannels_idxs.mat'];
dataChannelsStruct = load(dataChannelsFile);
dataChannels = dataChannelsStruct.dataChannels;
numDataChannels = length(dataChannels);
data_chanlocs = chanlocs(dataChannels);

%% Tasks
% Emotional valence
%   Valence: pos, neg, neu
%   Future possibility to look at valence based on subject rating
EMO_COND = {...
    'tertile_pre','valence',{'positive','negative','neutral'};...
    'tertile_post','valence',{'positive','negative','neutral'}};
% Nested cell-string format for task conditions
%   Generate task averages
                        
EPOCH_RANGE = [-0.9 2.0];
EPOCH_BASELINE = [-0.8 -0.3];

%% Time Frequency Information
lowerBoundFreq = 2;
upperBoundFreq = 58;
% 1/f scaled time frequency estimate
numFreqs = 150;
% Exponent 0.05 based on Voytek, Gazzaley 2015
exponent = 0.05;
% Scaled frequencies to use to estimate power
[freq] = pinkFreqScale(lowerBoundFreq,upperBoundFreq,numFreqs,exponent);
% Wavelet size 5
numCycles = 5;
% Down sampled rate
sampleRate = 200;

%% Single Subject Time Frequency analysis
task = 'EMO';
SINGLE_TF_TASK = [SINGLE_TF task '/'];
mkdir_JR(SINGLE_TF_TASK);

% There are many ways to specify the task conditions
TASK_COND = EMO_COND;
numCondSpecs = size(TASK_COND,1);

% Analysis window varies based on task
TASK_RANGE = EPOCH_RANGE;
startEpoch = TASK_RANGE(1);
stopEpoch  = TASK_RANGE(2);
analysisWindow = stopEpoch - startEpoch;
ms = 1000;
timePoints = (startEpoch*ms):(ms/sampleRate):((stopEpoch*ms)-(ms/sampleRate));
numTimePoints = length(timePoints);

% Loop through specifications for this task
for specsIdx = 1:numCondSpecs
    COND_SPEC = TASK_COND(specsIdx,:);

    specName = COND_SPEC{1};
    epochName = COND_SPEC{2};
    tfName = [epochName '_' specName];

    SPEC_DIR = [SINGLE_TF_TASK tfName '/'];
    mkdir_JR(SPEC_DIR);

    % Loop through each subject
    timeFreqFiles = cell(1,numSub);
    for subIdx = 1:numSub

        subject = SUBJECTS{subIdx};
        SUB_SPEC_DIR = [SPEC_DIR subject '/'];
        mkdir_JR(SUB_SPEC_DIR);
        % Frequently used for naming sysstem
        subTask = [subject '_' task];

        % Directory for preprocessed data for this subject
        SUB_PROC = [PROC_EEG task '/' subject '/'];
        % Load preprocessed epoched data
        preprocFilename = ['mcealf_SSDE_' subTask '.set'];
        preprocFile = [SUB_PROC preprocFilename];
        if exist(preprocFile,'file')~=2
            % Alternative name
            preprocFilename = ['mceaXrelf_SSDE_' subTask '.set'];
            preprocFile = [SUB_PROC preprocFilename];
            if exist(preprocFile,'file')~=2
                fprintf('Missing preprocessing file for %s %s\n',subject,task);
                continue;
            end
        end
        epochRejectionFile = [SUB_PROC 'epochRejection.mat'];
        if exist(epochRejectionFile,'file')==2
            epochRejection = load(epochRejectionFile);
            rmTrialIdxs = epochRejection.rejectedEpochs;
        else
            rmTrialIdxs = [];
        end

        % output file for time frequency results
        timeFreqFile = [SUB_SPEC_DIR 'tf_' tfName '_' subTask '.mat'];
        timeFreqFiles{subIdx} = timeFreqFile;

        % Check if time freq file already exists
        if exist(timeFreqFile,'file')~=2

            fprintf('\n\nRunning time frequency for %s %s %s\n',subject,task,tfName);

            % Get the conditions for this task specification
            COND = COND_SPEC{3};
            numCond = length(COND);

            % Initialize output variables
            allCond_ERSP = NaN(numFreqs,numTimePoints,numDataChannels,numCond);
            allCond_ITPC = NaN(numFreqs,numTimePoints,numDataChannels,numCond);

            % load dataset
            EEG = pop_loadset('filename',preprocFile);
            [ALLEEG,EEG,~]=eeg_store(ALLEEG,EEG,0);
            data = EEG.data;
            actual_numEpochs = size(data,3);

            % Load all the task information
            specNameParts = regexp(specName,'_','split');
            % Either image valence is defined by IAPS or is
            % defined by subjective rating
            valenceDefinition = specNameParts{1};
            % Timing is either pre or post stimulation
            timing = specNameParts{2};
            pattern_taskFile = ['SSDE_EV_' timing '_*'];
            prefix_rating_taskFile = 'SSDE_EVrate_*';

            % Directory for raw task data
            RAW_TASK_SUB = [RAW_TASK subject '/'];
            taskFileFndr = dir([RAW_TASK_SUB pattern_taskFile]);
            if isempty(taskFileFndr)
                fprintf('Failled to find task file!\n');
                continue;
            end
            taskFile = [RAW_TASK_SUB taskFileFndr(1).name];
            taskData = load(taskFile);
            taskResults = taskData.subResults;

            % Figure out which trials are included in each condition
            allCond_trialIdxs = cell(1,numCond);

            numBlocks = 2;
            numTrials = 48;
            numTotalTrials = numBlocks * numTrials;

            dim = [1 numTotalTrials];

            % Valence - must be defined by IAPS or subjectively
            if strcmp(valenceDefinition,'iaps')
                cond_valence = [{taskData.subResults.blockInfo(1).trialInfo.valence}...
                    {taskData.subResults.blockInfo(2).trialInfo.valence}];
            elseif strcmp(valenceDefinition,'subjective') || strcmp(valenceDefinition,'tertile')
                
                % Use the ratings
                prefix_rating_taskFile = 'SSDE_EVrate_*';
                ratingFileFndr = dir([RAW_TASK_SUB prefix_rating_taskFile]);

                % Rating file needs to exist
                if isempty(ratingFileFndr)
                    fprintf('Cannot use subjective rating for %s %s\n',subject,task);
                    continue;
                end

                % Load the subjective ratings file
                ratingFile = [RAW_TASK_SUB ratingFileFndr(1).name];
                ratings = load(ratingFile);

                % Actual image presented
                imageNames = [{taskData.subResults.blockInfo(1).trialInfo.imageName}...
                    {taskData.subResults.blockInfo(2).trialInfo.imageName}];

                % Subjective valence ratings for each image
                rating_imageNames = {ratings.subResults.imageInfo.imageName};
                rating_responses = [ratings.subResults.imageInfo.response];

                
                if strcmp(valenceDefinition,'subjective')
                    negNeuThreshold = 3.5;
                    neuPosThreshold = 7.5;
                else
                    ratingTertiles = quantile(rating_responses,[0.333 0.666]);
                    negNeuThreshold = ratingTertiles(1);
                    neuPosThreshold = ratingTertiles(2);
                    
                    if ~isempty(find(strcmp({'P116','P135'},subject),1))
                        neuPosThreshold = neuPosThreshold - 1;
                    elseif ~isempty(find(strcmp({'P029','P043','P187','P107'},subject),1))
                        negNeuThreshold = negNeuThreshold + 1;
                    end
                end
                
                numNeg_sub = length(find(rating_responses < negNeuThreshold));
                numPos_sub = length(find(rating_responses > neuPosThreshold));
                numNeu_sub = length(rating_responses) - (numNeg_sub + numPos_sub);
                fprintf('%i neg, %i neu, %i pos\n',numNeg_sub,numNeu_sub,numPos_sub);
                valenceDefFile = [SUB_SPEC_DIR subject '_' valenceDefinition '_info.mat'];
                if exist(valenceDefFile,'file')==2
                    delete(valenceDefFile);
                end
                valenceStruct = struct(...
                    'negNeuThreshold',negNeuThreshold,...
                    'neuPosThreshold',neuPosThreshold,...
                    'numNeg_sub',numNeg_sub,...
                    'numPos_sub',numPos_sub,...
                    'numNeu_sub',numNeu_sub);
                save(valenceDefFile,'-struct','valenceStruct');
                
                % Determine valence based on subjective rating
                cond_valence = cell(1,numTotalTrials);
                for trialIdx = 1:numTotalTrials
                    trial_imageName = imageNames{trialIdx};
                    ratingIdx = find(strcmpi(rating_imageNames,trial_imageName),1);
                    rating_response = rating_responses(ratingIdx);
                    
                    % 1-3 negative or tertile
                    if rating_response < negNeuThreshold
                        cond_valence{trialIdx} = 'negative';
                    % 8-10 positive or tertile
                    elseif rating_response > neuPosThreshold
                        cond_valence{trialIdx} = 'positive';
                    else
                        % 4-7 neutral or tertile
                        cond_valence{trialIdx} = 'neutral';
                    end
                end % loop trials

            end % switch valence definition

            % Loop through conditions which are valence
            for condIdx = 1:numCond
                cond = COND{condIdx};

                % initialize trial indexes
                cond_trialIdxs = ones(dim);

                % Update based on the valence
                cond_trialIdxs = cond_trialIdxs & strcmpi(cond,cond_valence);

                doubleNumTrials = numTotalTrials*2;
                final_cond_trialIdxs = zeros(1,doubleNumTrials);
                if strcmp(timing,'pre')
                    final_cond_trialIdxs(1:numTotalTrials) = cond_trialIdxs;
                elseif strcmp(timing,'post')
                     final_cond_trialIdxs((numTotalTrials+1):doubleNumTrials) = cond_trialIdxs;
                end

                % Remove trials for rejected epochs
                if ~isempty(rmTrialIdxs)
                    final_cond_trialIdxs(rmTrialIdxs)=[];
                end
                doubleNumTrials = length(final_cond_trialIdxs);

                % For some participants, the beginning of data
                % collection is missing
                if actual_numEpochs ~= doubleNumTrials
                    error('This should work but double check each one');
                    diff_trialCount = (doubleNumTrials - actual_numEpochs);
                    final_cond_trialIdxs(1:diff_trialCount) = [];
                end

                % store conditions
                allCond_trialIdxs{condIdx} = final_cond_trialIdxs;
            end % loop conditions

            tmpDir = [SUB_SPEC_DIR 'tmp/'];
            mkdir_JR(tmpDir);
            
            % Loop through each in condition
            for condIdx = 1:numCond
                cond = COND{condIdx};

                % Figure out which trials to include in this condition
                cond_trialIdxs = allCond_trialIdxs{condIdx};

                if length(find(cond_trialIdxs)) < 5
                    continue;
                end
                % Load previous data and mirror flip
                cond_data = data(:,:,logical(cond_trialIdxs));

                % Instead of EEGLAB files, input brick of data
                eegFiles = {cond_data};

                % file name for saving out temporary files
                fileNamePrefix = ['SSDE_' subTask '_' cond];
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
                [tfOutputFiles,~] = timeFreq_JR(outputDirectory,eegFiles,tfInfo);
                
                dataStruct = load(tfOutputFiles{1});
                ersp = dataStruct.tf.data.ersp;
                itpc = dataStruct.tf.data.itpc;
                for dataChanIdx = 1:numDataChannels
                    erspChan = ersp(:,:,dataChanIdx);
                    itpcChan = itpc(:,:,dataChanIdx);
                    
                    if size(erspChan,1) == 600
                        erspChan(1:20,:) = [];
                        itpcChan(1:20,:) = [];
                    end
                    % Prepare output
                    allCond_ERSP(:,:,dataChanIdx,condIdx)=erspChan';
                    allCond_ITPC(:,:,dataChanIdx,condIdx)=itpcChan';
                end
                
                % Delete the file after extracting data
                delete(tfOutputFiles{1});
                
            end % loop conditions

            % Save data out
            timeFreqStruct = struct(...
                'ERSP',allCond_ERSP,...
                'ITPC',allCond_ITPC);
            save(timeFreqFile,'-struct','timeFreqStruct');

            % Delete the temporary directory
            rmdir(tmpDir,'s');
        end % time frequency files exists
    end % loop subjects
end % loop task specifications
