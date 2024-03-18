classdef NNC_path
    %UTILS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties  (Constant)

        %% Data Paths
        pathToParentFolder = 'C:\Users\magcam\Documents\code_matlab\GAD\'
        pathToParentParentFolder = 'C:\Users\magcam\Documents\code_matlab\'
        
        pathToToolbox = 'C:\Users\magcam\Documents\toolboxes\'
        EEGLAB_TOOLBOX = [NNC_path.pathToToolbox 'eeglab2022.1\'];
        CSC_TOOLBOX = [NNC_path.pathToToolbox 'csc-eeg-tools-develop/']
        epi_TOOLBOX = [NNC_path.pathToToolbox 'epi-eeg-tools-master/'];
        pwr_TOOLBOX = [NNC_path.pathToToolbox 'power_all/'];
        rainbow1 = [NNC_path.pathToToolbox 'rainbowspectra\dsputils\'];
        rainbow2 = [NNC_path.pathToToolbox 'rainbowspectra\plotutils\']; 
        
        DataRoot_GAD_healthy_original = 'C:\Users\magcam\Documents\GAD\healthy\DATA_Raw\'
        DataRoot_GAD_healthy_processed = 'C:\Users\magcam\Documents\GAD\healthy\DATA_preprocessed\'
        DataRoot_GAD_anxious_original = 'C:\Users\magcam\Documents\GAD\patients\DATA_raw\'
        DataRoot_GAD_anxious_processed = 'C:\Users\magcam\Documents\GAD\patients\DATA_preprocessed\'

        DataRoot_MBA_original = 'D:\emotion_regulation\raw\movie\Movie_EEG\'
        DataRoot_MBA_processed = 'D:\emotion_regulation\raw\movie\Movie_EEG\processed\'
        DataRoot_MBA_start_trigger = 'D:\emotion_regulation\raw\movie\Movie_EEG\start_trigger\'

        DataRoot_LUPUS_original = 'C:\Users\magcam\Documents\LUPUS\LUPUS_DATA_raw\'
        DataRoot_LUPUS_processed = 'C:\Users\magcam\Documents\LUPUS\LUPUS_DATA_preprocessed\'

        DataRoot_OUD_original = 'C:\Users\magcam\Documents\OUD\OUD_DATA_raw\'
        DataRoot_OUD_processed = 'C:\Users\magcam\Documents\OUD\ProcessedData\'
        
        pathToRawData_old = [NNC_path.DataRoot_GAD_healthy_original 'Data_Study_old\']
        pathToRawData_new = [NNC_path.DataRoot_GAD_healthy_original 'Data_Study_new\']
       
        pathToDataOutput =  [NNC_path.DataRoot_GAD_healthy_processed 'Output_Simulink\'];
        pathToDataOutput_Matlab =  [NNC_path.DataRoot_GAD_healthy_processed 'Output_Matlab\'];
        pathToLinearMappings = "LinearMappings\*.json"
        pathToEEWraperOutput = [NNC_path.DataRoot_GAD_healthy_processed 'Output_EEWrapper\'];       
       
        %% Configs Paths
        ActivityMappingConfigPath= "*\Config_ActivityMapping.json"
        MetRangesPath = "*\Config_MetRanges.json"

        CutsPath= "*\EndCuts.csv"

        %% Figure Paths
        Figures_PS_ = [NNC_path.DataRoot_GAD_healthy_processed 'EE_SONEE_Output_Figures\Figures_MeanMagSquare\'];
        Figures_LinearMappings = [NNC_path.DataRoot_GAD_healthy_processed 'EE_SONEE_Output_Figures\Figures_LinearMappings\'];
        Figures_EE_Performance= [NNC_path.pathToEEWraperOutput 'Figures\']

    end
    
    methods  (Static)
        %% Get SONEE TaskNames
        % Reads the SONEE task names from the activity Mapping File
        function [Tasks,Tasks_clean] = getSONEETaskNames()
            Mappingpath= dir(Utils_SONEE.ActivityMappingConfigPath);
            ActivityMappingConf= Utils.readJson(fullfile (Mappingpath.folder,Mappingpath.name));
            Mapping=ActivityMappingConf.ActivityIndexes.SONEE;
            Tasks = string(fieldnames(Mapping));
            Tasks_clean= (fieldnames(Mapping));
            for ind= 1: numel(Tasks_clean)
                Tasks_clean{ind} = Tasks_clean{ind}(2:end);
            end
            Tasks_clean= string(Tasks_clean);
        end
        %% get EEG Recordings
        function recs_ = getRecs(pth)
            all_fls = dir(pth)
            all_fls_nam = {all_fls.name}
            longerThan3 = cellfun(@(x) length(x) > 3, all_fls_nam)
            fls_eeg = all_fls_nam(find(longerThan3))
            lst_grp = {}
            for i = 1:length(fls_eeg)
                if contains(lower(fls_eeg{i}),'movieeeg')
                    lst_grp{i} = fls_eeg{i}
                end
            end
            recs_ = lst_grp(~cellfun('isempty',lst_grp))
        end
        %% getOverviewTable
        function ResultTable_= getOverviewTable(frequency)
            arguments
                frequency {mustBeMember(frequency,["","_0.1","_0.5","_1" ])} = ""                
            end
            ResultTableAll = load("SaveData\DataOverview"+frequency +".mat");                        % load Data sorted SONEE data
            ResultTable_=ResultTableAll.ResultTable;            
        end
        %% get plot 
        function ch= getCharacteristics(path)
            ch = readtable([path 'Characteristics_Summary_SONEE_V2.xlsx'],...
                'Range', 'B5:AG65', 'VariableNamingRule', 'modify');                                % read in participant characteristics table
        end
        %% get Photodiode signal of missing triggers
        function idx_miss = getPhotodiode(b,EEG)
            highs = abs(b.photodiode)>80000
    
            lst_on = []
            lst_off = []
            for i = 1:length(highs)
                if i==1
                    continue
                elseif highs(i)==1 && highs(i-1)~=1
                    lst_on = [lst_on,i]
                elseif highs(i)==0 && highs(i-1)~=0
                    lst_off = [lst_off,i-1]
                end
            end
            
            lats={EEG.event.latency}
            lats_ms =[cellfun(@(x) (x/200)*1000, lats)]
            
            lst_miss_on = []
            lst_miss_off = []
            for i = 1:length(lst_on)
                logic_on = (lats_ms > (lst_on(i)-50)) & (lats_ms < (lst_on(i)+50));
                logic_off = (lats_ms > (lst_off(i)-50)) & (lats_ms < (lst_off(i)+50));
                if sum(logic_on)==0
                    lst_miss_on = [lst_miss_on,i]
                elseif sum(logic_off)==0
                    lst_miss_off = [lst_miss_off,i]
                end
            end
            
            mss_vals_on = lst_on(lst_miss_on)
            mss_vals_off = lst_off(lst_miss_off)
            
            new_lats = [lats_ms,mss_vals_on,mss_vals_off]
            
            [sorted_values, sorted_indices] = sort(new_lats);
            idx_miss = find(sorted_indices>length(lats_ms))
        end
        %% get TestTrainSet
        function [TestSet, TrainSet] = getTestTrainSet()
            SettingPath= dir(Utils_SONEE.TestTrainSetPath);
            SplitData= Utils.readJson(fullfile(SettingPath.folder,SettingPath.name));
            TestSet= SplitData.Test;
            TrainSet= SplitData.Train;           
        end
        %% Get binary lists
        function idx_opt= getBinaries(lst_lat_ms,lst_xlsx)
            lst={}
            for i=1:length(lst_lat_ms)
                if i ==1
                    lst{i} = 1
                elseif lst_lat_ms{i}-lst_lat_ms{i-1} < 1000
                    lst{i} = 2
                elseif lst_lat_ms{i}-lst_lat_ms{i-1} > 1000
                    lst{i} = 1
                end
            end
            
            lst_manual_trggrs = []
            for i = 1:length(lst_xlsx)
                if all(contains(lst_xlsx{i},'PR_')) 
                    lst_manual_trggrs=[lst_manual_trggrs,2]
                elseif all(contains(lst_xlsx{i},'stimulusstart')) 
                    lst_manual_trggrs=[lst_manual_trggrs,1]
                elseif all(contains(lst_xlsx{i}, 'stimulusend'))
                    lst_manual_trggrs=[lst_manual_trggrs,2]
                elseif all(contains(lst_xlsx{i}, 'blockstart'))
                    lst_manual_trggrs=[lst_manual_trggrs,1]
                elseif all(contains(lst_xlsx{i}, 'blockend'))
                    lst_manual_trggrs=[lst_manual_trggrs,2]
                end
            end  

            c=1
            count = 0
            idx_opt=[]
            for i=1:length(lst_manual_trggrs)
                if lst_manual_trggrs(i) ~= lst{c} %&& all(contains(lst_xlsx_def{i},'blockend'))
                    idx_opt = [idx_opt,i-count];
                    c=c-1
                    count=count+1
                end
                c=c+1;
            end

        end
        %% getCuts 
        function CutTime= getCuts(table_,id_, activityString)
            cutTable= table_;
            activityTableID = cellfun(@(x) str2num(x(end-1:end)), cutTable.Properties.VariableNames(1:end-1))  ;
            activityID= str2num(activityString(1:2));
            Index= find(activityTableID==activityID);
            if isempty(Index)
                CutTime=0;
            else
                relevantEntry= cutTable(cutTable.ID ==id_ ,Index);
                CutTime=relevantEntry{1,1};
                if isnan(CutTime)
                    CutTime=0;                
                end
            end
        end

        %% replace Empty with NAN
        function out= replaceEmptyWithNan(cellarray)
            mycellarray= cellarray;
            emptyIndex = cellfun('isempty', mycellarray);     % Find indices of empty cells
            mycellarray(emptyIndex) = {nan};                    % Fill empty cells with 0
            out=mycellarray;        
        end
        %% Split Participants
        function [Test, Train, seed_,useSubset]=SplitParticipants (option_,seed_)
            arguments
                option_ {mustBeMember(option_,["","Seed"])} = ""
                seed_ = 0
            end

            switch option_
                case ""
                    Train = [1:60]'; 
                    Test = [1:60]'; 
                case "Seed"
                    rng(seed_) % set the seed for reproducibility
                    Train = randperm(60, 48)';                                       % Select Random 48 participants from the 60 for Training
                    Test = setdiff([1:60],Train)';
            end
            if strcmp(option_,"")
                useSubset = false;
            else
                useSubset= true;
            end
        end
        
    %% Return the Mapping
    function [expectedClusterMap, activityMap]= getActivityMapping(MappingConfig,MappingVersion)   
        ActivityMappingConf= MappingConfig;             % load the ActivityMappingConfig and show the Info
        if strcmp(MappingVersion(1),'_')
            mapV =MappingVersion(2:end);   
        else
            mapV=MappingVersion;
        end
        fNames= fieldnames(ActivityMappingConf.ActivityIndexes);
        fieldIndex= -1;
        for fI = 1:numel(fNames)   
            if contains(mapV,fNames{fI})
                fieldIndex=fI;
                break
            else
                continue
            end
        end
        if fieldIndex ==-1
            warning("Mapping Not Found")
            expectedClusterMap= struct();
            return
        end
        expectedClusterMap=ActivityMappingConf.ActivityIndexes.(fNames{fieldIndex});
        activityMap  =ActivityMappingConf.ActivityMapping.(fNames{fieldIndex});
    end        
        %% 
%%%%%%%%%%%%%%%%%%% Functions from Nick's Performance Script %%%%%%%%%%%%%%
%%    classifyIntensityByMET
    function T = classifyIntensityByMET(T, inVarName, outVarName)
    T.(outVarName)(lt(T.(inVarName), 1.6)) = "vLow";
    T.(outVarName)(le(T.(inVarName), 3.1) & ge(T.(inVarName), 1.6)) = "Low";
    T.(outVarName)(le(T.(inVarName), 4.7) & gt(T.(inVarName), 3.1)) = "Mod";
    T.(outVarName)(le(T.(inVarName), 6.7) & gt(T.(inVarName), 4.7)) = "Vig";
    T.(outVarName)(gt(T.(inVarName), 6.8)) = "Vig"; % should be vVig
    end      
    %% BA Plot
    function [m, lLOA, uLOA, err, t_mean] = ba_plot(t, e, produce_plot, method, plot_bfl, varargin)
        if ~isempty(varargin)
            varargin = varargin{:};
            groups = unique(varargin);
            g = length(groups);
        else
            varargin = cell(length(t), 1);
            varargin(:) = {'default'};
            groups = {'default'};
            g = 1;
        end
        idx = ~isnan(t) & ~isnan(e);
        t = t(idx);
        e = e(idx);
        varargin = varargin(idx);
        y = e - t;
        m = mean(y);
        lLOA = m - 1.96*std(y);
        uLOA = m + 1.96*std(y);
        err = 1.96*std(y);
        t_mean = mean(t);
        if produce_plot
            switch method
                case 'mean'
                    x = mean([e,t],2);
                case 'true'
                    x = t;
            end
            for i = 1 : g
                idx = strcmp(varargin, groups{i});
                scatter(x(idx), y(idx), 16, 'filled');
                hold on
            end
            yline(0, ':k', 'LineWidth', 2)
            yline(m,'-r', 'LineWidth', 2)
            yline(lLOA, '--r')
            yline(uLOA, '--r')
            l = get(gca, 'XLim');
            text(l(2), m, num2str(round(m,2)), 'Color', 'red');
        
            if plot_bfl
                c = polyfit(x, y, 1);
                xFit = linspace(l(1), l(2), 1000);
                yFit = polyval(c , xFit);
                plot(xFit, yFit, 'g', 'LineWidth', 2);
            end
        
            ylabel('difference between measures')
            xlabel([method, ' value'])
            box on
            if g > 1
                L = legend(replace(groups,'_',' '));
%                 L.String = L.String(1:2);
            L.Location  = 'bestOutside';
            end
        end
    end

    %%
    end
end

