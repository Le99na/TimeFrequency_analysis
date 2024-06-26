% wipe out all previous data
clear all; clc; close all;

% Preprocessing version
preprocVersion = '02052023';

% Use ICLabel or not
useIClabel = 1;
skipManualIC = 0;

%% Step 0: get organized
% Paths to data
root = '/home/mc/Documents/code_matlab/GAD_TF/'
addpath(genpath(root))
ClassObj = NNC_path();

root_raw = ClassObj.DataRoot_GAD_TF_orig
root_proc = ClassObj.DataRoot_GAD_TF_processed 

%generate cell array of available participants
grps_strct=dir(root_raw)
grp_nams = {grps_strct.name}


longerThan3 = cellfun(@(x) length(x) > 3 && contains(lower(x),'task'), grp_nams)
GROUPS = grp_nams(find(longerThan3))

numGroups = length(GROUPS);

%read-in relevant file
fls_all = {dir(root_raw).name}
fl_rel_idx = find(contains(fls_all,'ERN_vals.set'))
fl_rel = fls_all(fl_rel_idx)

% Add the EEGLAB toolbox to the path
addpath(genpath(ClassObj.EEGLAB_TOOLBOX));

%add csc toolbox
addpath(genpath(ClassObj.CSC_TOOLBOX));

%add epi-eeg toolbox
addpath(genpath(ClassObj.epi_TOOLBOX));

%add power scripts toolbox
addpath(genpath(ClassObj.pwr_TOOLBOX));

% Add Mengsen butter function + rainbow maps

addpath(genpath(ClassObj.rainbow1)); 

addpath(genpath(ClassObj.rainbow2));

sampleRate = 200;
numCycles = 3; %high temporal resolution: 3-5, high frequency resolution: 5-8
FREQ = [1:100]
numFreq = length(FREQ);
baselinePeriod = 1:40

% Number of cycles to estimate tf can change for each frequency
if length(numCycles) ~= numFreq
    numCycles = ones(1,numFreq) * numCycles(1);
end

% Defualt is to mirror the data
mirrorData_FLAG = 1;


% Start EEGlab
eeglab

%% Go through each participant
for groupIdx = 1:length(fl_rel)
    % Name of the group which is used in folder organization

    path_out = [root_proc fl_rel{groupIdx}(1:10) '_nc3/']

    if exist(path_out,'dir')~=7
        mkdir(path_out);
    end
    
    EEG = pop_loadset([root_raw fl_rel{groupIdx}])

    answ = {EEG.event.type}
    idx_ep_cor = []
    idx_ep_icor = []
    count=1 %this is the count for the number of epochs
    for i = 1:numel(answ)
        if all(contains(answ{i},'PR_correct'))
            i
            count
            idx_ep_cor = [idx_ep_cor,count];
            count=count+1;
            answ{i}
        elseif all(contains(answ{i},'PR_incorrect'))
            i
            count
            idx_ep_icor = [idx_ep_icor,count];
            count=count+1;
            answ{i}
        end
    end

    EEG_corr = EEG.data([6,11,129],:,idx_ep_cor); %6:FCz, 11:Fz, 129:Cz
    EEG_icorr = EEG.data([6,11,129],:,idx_ep_icor);

    EEG.idx_cor = idx_ep_cor
    EEG.idx_icor = idx_ep_icor

    EEG = pop_saveset(EEG, 'filename', [fl_rel{groupIdx}, '_idx_cor_icor.set'], 'filepath', path_out);

    roi_Data_cor = squeeze(mean(EEG_corr(:,:,:),1));

    roi_Data_icor = squeeze(mean(EEG_icorr(:,:,:),1));

    % find the save of this condition's data
    [numTimePts_cor,numTrials_cor] = size(roi_Data_cor);
    [numTimePts_icor,numTrials_icor] = size(roi_Data_icor);

    % Run ersp for each channel separately
    amp_cor = NaN(numTimePts_cor,numFreq,numTrials_cor);
    amp_icor = NaN(numTimePts_icor,numFreq,numTrials_icor);
    phase_cor = NaN(numTimePts_cor,numFreq,numTrials_cor);
    phase_icor = NaN(numTimePts_icor,numFreq,numTrials_icor);
  
    % loop through each trial for this subject
    for trialIdx = 1:numTrials_cor
        % data for only this trial
        trialData = roi_Data_cor(:,trialIdx)';
        
        % If the data is being mirror flipped
        if mirrorData_FLAG
            % mirror flip the trial data
            mirrorData = [fliplr(trialData) trialData fliplr(trialData)];
            mirrorIdxs = (numTimePts_cor+1):(numTimePts_cor*2);
        else
            % Or just use the original data
            mirrorData = trialData;
            mirrorIdxs = 1:numTimePts_cor;
        end
        
        % run time frequency analysis for each frequency
        for freqIdx = 1:numFreq
            freq = FREQ(freqIdx);
            % run wavelet based time frequency estimate
            [phase,amp] = phaseAmp_calculator(sampleRate,freq,mirrorData,numCycles(freqIdx));
            
            % Prepare output
            phase_cor(:,freqIdx,trialIdx) = phase(mirrorIdxs);
            amp_cor(:,freqIdx,trialIdx) = amp(mirrorIdxs);
        end % freq loop
    end % trial loop

    % loop through each trial for this subject
    for trialIdx = 1:numTrials_icor
        % data for only this trial
        trialData = roi_Data_icor(:,trialIdx)';
        
        % If the data is being mirror flipped
        if mirrorData_FLAG
            % mirror flip the trial data
            mirrorData = [fliplr(trialData) trialData fliplr(trialData)];
            mirrorIdxs = (numTimePts_icor+1):(numTimePts_icor*2);
        else
            % Or just use the original data
            mirrorData = trialData;
            mirrorIdxs = 1:numTimePts_icor;
        end
        
        % run time frequency analysis for each frequency
        for freqIdx = 1:numFreq
            freq = FREQ(freqIdx);
            % run wavelet based time frequency estimate
            [phase,amp] = phaseAmp_calculator(sampleRate,freq,mirrorData,numCycles(freqIdx));
            
            % Prepare output
            phase_icor(:,freqIdx,trialIdx) = phase(mirrorIdxs);
            amp_icor(:,freqIdx,trialIdx) = amp(mirrorIdxs);
        end % freq loop
    end % trial loop

    % Average across trials
    amp_trialAvg_cor = mean(amp_cor,3);
    amp_trialAvg_icor = mean(amp_icor,3);
    
    % Baseline correction
    amp_baselineCorrected_cor = NaN(numTimePts_cor,numFreq);
    amp_baselineCorrected_icor = NaN(numTimePts_icor,numFreq);

    % Calculate baseline
    amp_baselineAvg_cor = mean(amp_trialAvg_cor(baselinePeriod,:),1);
    amp_baselineAvg_icor = mean(amp_trialAvg_icor(baselinePeriod,:),1);

    % Baseline correct all data
    for timeIdx = 1:numTimePts_cor
        timepoint_amp = amp_trialAvg_cor(timeIdx,:);
        timepoint_amp = (timepoint_amp - amp_baselineAvg_cor) ./ amp_baselineAvg_cor;
        amp_baselineCorrected_cor(timeIdx,:) = timepoint_amp;

        timepoint_amp = amp_trialAvg_icor(timeIdx,:);
        timepoint_amp = (timepoint_amp - amp_baselineAvg_icor) ./ amp_baselineAvg_icor;
        amp_baselineCorrected_icor(timeIdx,:) = timepoint_amp;

    end

    EEG.amp_avg_cor = amp_baselineCorrected_cor
    EEG.amp_avg_icor = amp_baselineCorrected_icor

    EEG = pop_saveset(EEG, 'filename', [EEG.filename(1:end-4), '_avg_amp.set'], 'filepath', path_out);

    % Plot the spectrogram
    figure;
    imagesc(1/200:1/200:300/200, FREQ(1:40), amp_baselineCorrected_cor.');
    set(gca, 'YDir', 'normal');
    title([EEG.filename(1:6) ' time-frequency spectrogram cor']);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar;
    % Save the figure as a JPEG
    output_filename = [path_out EEG.filename(1:6) '_time_frequency_spectrogram_cor.jpg'];  % Specify the output filename
    saveas(gcf, output_filename, 'jpeg');
    close all

    % Plot the spectrogram
    figure;
    imagesc(1/200:1/200:300/200, FREQ(1:40), amp_baselineCorrected_icor.');
    set(gca, 'YDir', 'normal');
    title([EEG.filename(1:6) ' time-frequency spectrogram icor']);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar;
    % Save the figure as a JPEG
    output_filename = [path_out EEG.filename(1:6) '_time_frequency_spectrogram_icor.jpg'];  % Specify the output filename
    saveas(gcf, output_filename, 'jpeg');
    close all
    
end

flds = dir(root_proc)

grp_nams = {flds.name}


longerThan3 = cellfun(@(x) length(x) > 3, grp_nams)
flds_rel = grp_nams(find(longerThan3))

lst_ps_icor = {}
lst_ps_cor = {}
for groupIdx = 1:length(flds_rel)
    pat = dir([root_proc flds_rel{groupIdx} '/'])
    nms = {pat.name}

    longerThan3 = cellfun(@(x) contains(x,'avg_amp.set'), nms)
    nam_rel = nms(find(longerThan3))

    EEG = pop_loadset([root_proc flds_rel{groupIdx} '/' nam_rel{1}])

    lst_ps_icor{groupIdx} = EEG.amp_avg_icor
    lst_ps_cor{groupIdx} = EEG.amp_avg_cor
end


diff_86_icor = lst_ps_icor{1}-lst_ps_icor{2}
diff_86_cor = lst_ps_cor{1}-lst_ps_cor{2}

diff_110_icor = lst_ps_icor{4}-lst_ps_icor{3}
diff_110_cor = lst_ps_cor{4}-lst_ps_cor{3}

 % Plot the spectrogram
figure;
imagesc(1/200:1/200:300/200, FREQ(1:40), diff_86_icor.');
set(gca, 'YDir', 'normal');
title('ptp86 difference time-frequency spectrogram icor');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
% Save the figure as a JPEG
output_filename = [root_proc 'ptp86_difference_time_frequency_spectrogram_icor.jpg'];  % Specify the output filename
saveas(gcf, output_filename, 'jpeg');
close all
figure;
imagesc(1/200:1/200:300/200, FREQ(1:40), diff_110_cor.');
set(gca, 'YDir', 'normal');
title('ptp110 difference time-frequency spectrogram cor');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
% Save the figure as a JPEG
output_filename = [root_proc 'ptp110_difference_time_frequency_spectrogram_cor.jpg'];  % Specify the output filename
saveas(gcf, output_filename, 'jpeg');

figure;
imagesc(1/200:1/200:300/200, FREQ(1:40), diff_86_cor.');
set(gca, 'YDir', 'normal');
title('ptp86 difference time-frequency spectrogram cor');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
% Save the figure as a JPEG
output_filename = [root_proc 'ptp86_difference_time_frequency_spectrogram_cor.jpg'];  % Specify the output filenam
figure;
imagesc(1/200:1/200:300/200, FREQ(1:40), diff_110_cor.');
set(gca, 'YDir', 'normal');e
saveas(gcf, output_filename, 'jpeg');
close all

figure;
imagesc(1/200:1/200:300/200, FREQ(1:40), diff_110_icor.');
set(gca, 'YDir', 'normal');
title('ptp110 difference time-frequency spectrogram icor');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
% Save the figure as a JPEG
output_filename = [root_proc 'ptp110_difference_time_frequency_spectrogram_icor.jpg'];  % Specify the output filename
saveas(gcf, output_filename, 'jpeg');
close all

figure;
imagesc(1/200:1/200:300/200, FREQ(1:40), diff_110_cor.');
set(gca, 'YDir', 'normal');
title('ptp110 difference time-frequency spectrogram cor');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
% Save the figure as a JPEG
output_filename = [root_proc 'ptp110_difference_time_frequency_spectrogram_cor.jpg'];  % Specify the output filename
saveas(gcf, output_filename, 'jpeg');
close all