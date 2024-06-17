path = "E:\TF_analysis\calculate_average\"

grps_strct=dir(path)
grp_nams = {grps_strct.name}



idx = find(cellfun(@(x) contains(x, 'FU'), grp_nams));

nams_FU = grp_nams(idx)


diff_icor_all = nan(300,150)
diff_cor_all = nan(300,150)
for groupIdx = 1:length(nams_FU)
    nam = nams_FU{groupIdx}
    if length(nam)>2 &&  contains(nam,'.set')     
        prts = split(nam,'_')
        num = prts{1}(4:end)
        idx_bs = find(cellfun(@(x) contains(x, num) && ~contains(x, 'FU')  && ~contains(x, '.fdt'), grp_nams));
        if length(idx_bs)>1
            for idx_s = 1:length(idx_bs)
                prts_2 = split(grp_nams{idx_bs(idx_s)},'_')
                num_2 = prts_2{1}(4:end)
                if str2num(num_2) == str2num(num)
                    idx_bs = idx_bs(idx_s)
                end
            end
        else
             fullPath_bs = strcat(path, grp_nams{idx_bs});
        end
        fullPath_fu = strcat(path, nam);
        eeg_FU = pop_loadset(char(fullPath_fu))
        eeg_BS = pop_loadset(char(fullPath_bs))
        icor_FU = eeg_FU.amp_avg_icor
        cor_FU = eeg_FU.amp_avg_cor
        icor_BS = eeg_BS.amp_avg_icor
        cor_BS = eeg_BS.amp_avg_cor

        diff_icor = icor_FU - icor_BS
        diff_cor = cor_FU - cor_BS

        diff_icor_all = cat(3,diff_icor_all,diff_icor)
        diff_cor_all = cat(3,diff_cor_all,diff_cor)
    end
end

mean_diff_icor = nanmean(diff_icor_all, 3);
mean_diff_cor = nanmean(diff_cor_all, 3);

figure;
imagesc(1/200:1/200:300/200, FREQ, mean_diff_icor.');
set(gca, 'YDir', 'normal');
title('mean FU-BS incorrect');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
% Save the figure as a JPEG
output_filename = [char(path) 'mean_FU-BS_icor.jpg'];  % Specify the output filename
saveas(gcf, output_filename, 'jpeg');
close all

figure;
imagesc(1/200:1/200:300/200, FREQ, mean_diff_cor.');
set(gca, 'YDir', 'normal');
title('mean FU-BS correct');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
% Save the figure as a JPEG
output_filename = [char(path) 'mean_FU-BS_cor.jpg'];  % Specify the output filename
saveas(gcf, output_filename, 'jpeg');
close all