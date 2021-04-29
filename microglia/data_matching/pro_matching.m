function [process_data,excel_table] = pro_matching(fileName, pathName, main_suite2p)
% This function takes a filename and pathname for an excel file (.xlsx) and
% matches the data in a specified sheet to the Matlab output data from
% Suite2p. Written by John Wagner
% sheet = input('Please enter the number of the sheet you want to load. Must be an integer. \n');
test_data = readtable(fullfile(pathName,fileName),'Sheet',2);
% convert Final_calcium_process column  to matlab indeces
test_data.Final_calcium_process = test_data.Final_calcium_process+1;

%% Load Suite2p Data
clc
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
%% Match suite2p data to excel data
% put video names in their own vector
suite2p_names = folder_data.name(1:end);
% create new structure and phase out test_data
excel_structure = table2struct(test_data);
% remove rows that are marked for exclusion
index = test_data.Pro_exclude0_include1;
excel_structure = excel_structure(index == 1);
% Remove all rows with NaN -- we know these won't match
excel_structure = excel_structure(~isnan([excel_structure.Final_calcium_process]));
% let's test some logic here
for jj = 1:length(suite2p_names)
    % Find which rows which match to current suite2p_name
    findtest = find(suite2p_names(jj) == [excel_structure.Image]);
    if isempty(findtest) == 1 % check if findtest is empty, cont if it is
        continue
    end
    % Extract ROI identifiers from rows found above
    testroi = [excel_structure(findtest).Final_calcium_process];
    % Extract full F array for current video the loop is on
    full_F = data_names{jj,2}.F;
    % Subset of full_F -- only rows for each ROI
    match_F = full_F(testroi,:);
    for kk = 1:size(match_F,1)
        temp_F = match_F(kk,:);
        excel_structure(findtest(kk)).suite2pF = temp_F;
    end
end


% delete the rows with empty suite2pF field because there is no matched
% data
% convert excel_struct back to a table format
excel_table = struct2table(excel_structure);
% convert table to cell
excel_cells = table2cell(excel_table);
% create cell array with same #of rows as excel_cells/table
% final_matched data will have all the matched data -- map_ca sorts
final_matched_data = cell(size(excel_table,1),11);
% % % % Insert excel and suite2p data into columns of final_matched
% % ADD MEAN_velocity OR TIME_HITS_R_1/2 (1 of 2 columns)
% % AND CELL ID
maxD = excel_table.MaxD2Rum - excel_table.radius_1_from_area;
radius = excel_table.radius_1_from_area;
ROIs = excel_table.Final_calcium_process;
ROIs = num2cell(ROIs);
mouseNums = num2cell(excel_table.Mouse);
cellID = num2cell(excel_table.CELL_ID);
images = num2cell(excel_table.Image);

final_matched_data(:,1) = excel_table.TrackID;
final_matched_data(:,2) = ROIs;
final_matched_data(:,3) = excel_table.suite2pF;
final_matched_data(:,4) = excel_table.Groups;
final_matched_data(:,5) = num2cell(excel_table.DPI);
final_matched_data(:,6) = num2cell(excel_table.TimeHitRadius1);
final_matched_data(:,7) = mouseNums;
final_matched_data(:,8) = cellID;
final_matched_data(:,9) = images;
final_matched_data(:,10) = num2cell(maxD);
final_matched_data(:,11) = num2cell(radius);
% % % % done inserting excel data -- now delete rows w/ empty suite2pF
empty_cells = any(cellfun('isempty',final_matched_data),2);
final_matched_data(empty_cells,:) = [];

% % start with empty cell array size of final_matched_data
sorted_data = cell(size(final_matched_data));

% col. 4 has groups
TMEV_logicals = contains(final_matched_data(:,4),'TMEV'); 
% rows where TMEV occurs
TMEV_rows = find(TMEV_logicals);
% extract main table data
TMEV_extracted = final_matched_data(TMEV_logicals,:);
% repeat for PBS group
PBS_logicals = contains(final_matched_data(:,4),'PBS'); 
PBS_rows = find(PBS_logicals);
PBS_extracted = final_matched_data(PBS_logicals,:);

sorted_data(1:size(PBS_extracted,1),:) = PBS_extracted;
sorted_data(size(PBS_extracted,1)+1:end,:) = TMEV_extracted;

process_data = sortrows(sorted_data,5); % finished
end
