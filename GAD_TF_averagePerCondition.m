path = "E:\TF_analysis\calculate_average\"

grps_strct=dir(path)
grp_nams = {grps_strct.name}

BS_icor = nan(300,150)
BS_cor = nan(300,150)
FU_icor = nan(300,150)
FU_cor = nan(300,150)

for groupIdx = 1:length(grp_nams)
    nam = grp_nams{groupIdx}

    if length(nam)>2 &&  contains(nam,'.set')

        fullPath = strcat(path, nam); 
    
        EEG = pop_loadset(char(fullPath))

        cor = EEG.amp_avg_cor
        icor = EEG.amp_avg_icor

        if contains(nam,'FU')
            FU_cor = cat(3,FU_cor,cor)
            FU_icor = cat(3,FU_icor,icor)
        else
            BS_cor = cat(3,BS_cor,cor)
            BS_icor = cat(3,BS_icor,icor)
        end

    end

end


mean_BS_icor = nanmean(BS_icor, 3);
mean_BS_cor = nanmean(BS_cor, 3);
mean_FU_cor = nanmean(FU_cor, 3);
mean_FU_icor = nanmean(FU_icor, 3);

sampleRate = 200;
lowerBoundFreq = 2;
upperBoundFreq = 40;
% 1/f scaled time frequency estimate
numFreqs = 150;
% Exponent 0.05 based on Voytek, Gazzaley 2015
exponent = 0.05;
% Scaled frequencies to use to estimate power
FREQ = pinkFreqScale(lowerBoundFreq,upperBoundFreq,numFreqs,exponent);


figure;
imagesc(1/200:1/200:300/200, FREQ, mean_BS_icor.');
set(gca, 'YDir', 'normal');
title('mean baseline incorrect');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
% Save the figure as a JPEG
output_filename = [char(path) 'mean_bs_icor.jpg'];  % Specify the output filename
saveas(gcf, output_filename, 'jpeg');
close all

figure;
imagesc(1/200:1/200:300/200, FREQ, mean_BS_cor.');
set(gca, 'YDir', 'normal');
title('mean baseline correct');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
% Save the figure as a JPEG
output_filename = [char(path) 'mean_bs_cor.jpg'];  % Specify the output filename
saveas(gcf, output_filename, 'jpeg');
close all

figure;
imagesc(1/200:1/200:300/200, FREQ, mean_FU_icor.');
set(gca, 'YDir', 'normal');
title('mean FollowUp incorrect');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
% Save the figure as a JPEG
output_filename = [char(path) 'mean_FU_icor.jpg'];  % Specify the output filename
saveas(gcf, output_filename, 'jpeg');
close all

figure;
imagesc(1/200:1/200:300/200, FREQ, mean_FU_cor.');
set(gca, 'YDir', 'normal');
title('mean FollowUp correct');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
% Save the figure as a JPEG
output_filename = [char(path) 'mean_FU_cor.jpg'];  % Specify the output filename
saveas(gcf, output_filename, 'jpeg');
close all