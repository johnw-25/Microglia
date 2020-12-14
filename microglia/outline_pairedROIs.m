% % % Use this script to outline how to look at paired data efficiently
% % % A mix of pseudocode and just explaining some ideas w/ words

%% Load Suite2p Data
clc
main_suite2p = uigetdir();
% Get all the subfolders
videos = dir(main_suite2p);
SubFold = videos([videos.isdir]); % keep only the directories
SubFold = SubFold(3:end,:);
filelist = dir(fullfile(main_suite2p, '**\*.mat'));  % get list of mat files in any subfolder

% Loop on each folder
MyData = [];
for i = 1:length(filelist) % 1 to number of Fall.mat files
  filetoread = fullfile(filelist(i).folder,'Fall.mat'); % extract suite2p data from video subfolders
  % store each Fall.mat file into MyData
  fileID = fopen(filetoread);
  MyData{end+1} = load(filetoread); % the format depends of your files
  fclose(fileID);
end

data_names = cell(size(MyData,2),2); %create 2xn cell array to concatenate video names and corresponding traces
data_names(:,2) = MyData';
% create a table b/c matlab is weird and this is the best way i can think
% to access folder names
folder_data = struct2table(SubFold);
folder_data.name = str2double(folder_data.name);

if length(folder_data.name) ~= length(data_names(:,1))
    for ii = 1:length(folder_data.name)
        video_folder = fullfile(SubFold(ii).folder, SubFold(ii).name);
        data_check = ~isempty(dir(fullfile(video_folder, '**\*.mat')));
        if data_check == 0 % no Fall.mat file
            folder_data(ii,:) = [];
        end
    end
end
% turn folder names into cell array for input into data array
% % names = cellstr(folder_data.name(3:end));
% % % 2 column array with video name and corresponding video data from suite2p
data_names(:,1) = num2cell(folder_data.name(1:end));

%% get all rows with paired_exclude0_include1 == 1
[fileName, pathName] = uigetfile('*.xlsx');
test_data = readtable(fullfile(pathName,fileName),'Sheet',2);
% convert Final_calcium_process column  to matlab indeces
index = test_data.PAIRSexclude0_include1;
index(isnan(index)) = 0;
pairsIndex = logical(index);
pairsTable = test_data(pairsIndex,:);
suite2p_names = folder_data.name(1:end);
%convert test_data to struct - easier to index into
for jj = 1:length(suite2p_names)
    % Find which rows which match to current suite2p_name
    subTable = pairsTable(suite2p_names(jj) == pairsTable.Image,:);
    somaRois = subTable.Final_calcium_soma;
    processRois = subTable.Final_calcium_process;
    full_F = data_names{jj,2}.F;
    tempFieldName = strcat('T',num2str(suite2p_names(jj)));
    if isempty(subTable) % check if findtest is empty, cont if it is
        continue
    end
    
    for k = 1:size(subTable,1)
        currentTrace = subTable{k,1};
        currentTrace = currentTrace{1};
        TRACES.(tempFieldName).(currentTrace).timehitR = subTable.TimeHitRadius1(k);
        if ~isnan(somaRois(k))
            TRACES.(tempFieldName).(currentTrace).Soma.ROI = somaRois(k);
            TRACES.(tempFieldName).(currentTrace).Soma.traceF = full_F(somaRois(k)+1,:);
%         elseif ~isempty(subTable.Best_caROI(k)) && isnumeric(subTable.Best_caROI(k)) && strcmp(subTable.Somma_process{k},'s')
%             TRACES.(tempFieldName).(currentTrace).Soma.ROI = subTable.Best_caROI(k);
%             TRACES.(tempFieldName).(currentTrace).Soma.traceF = full_F(subTable.Best_caROI(k)+1,:);
        else
            continue % no pair
        end
        
        if ~isnan(processRois)
            TRACES.(tempFieldName).(currentTrace).Process.ROI = processRois(k);
            TRACES.(tempFieldName).(currentTrace).Process.traceF = full_F(processRois(k)+1,:);
%         elseif ~isempty(subTable.Best_caROI(k)) && isnumeric(subTable.Best_caROI(k)) && strcmp(subTable.Somma_process{k},'p')
%             TRACES.(tempFieldName).(currentTrace).Process.ROI = subTable.Best_caROI(k);
%             TRACES.(tempFieldName).(currentTrace).Process.traceF = full_F(subTable.Best_caROI(k)+1,:);
        else
            continue % no pair
        end
    end
end
%% So data is now easily accessible
