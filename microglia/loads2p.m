function [data_names, MyData, iscell_list] = loads2p()
% Path of the main folder : Rotation
main_suite2p = uigetdir;
% Get all the subfolders
videos = dir(main_suite2p);
SubFold = videos([videos.isdir]); % keep only the directories
filelist = dir(fullfile(main_suite2p, '**\*.mat'));  % get list of mat files in any subfolder
iscell_files = dir(fullfile(main_suite2p, '**\iscell.npy')); % extract all iscell files 
% Loop on each folder
MyData = cell(length(filelist),1);
iscell_list = cell(length(filelist),1);
for i = 1:length(filelist) % 1 to number of Fall.mat files
  filetoread = fullfile(filelist(i).folder,'Fall.mat'); % extract suite2p data from video subfolders
  iscell_toread = fullfile(filelist(i).folder,'iscell.npy');
  % store each Fall.mat file into MyData
  fileID = fopen(filetoread);
  % really only need iscellID if there are not-cell ROIs in suite2p
  iscellID = fopen(iscell_toread);
  MyData{end+1} = load(filetoread); % the format depends of your files
  iscell_list{end+1} = readNPY(iscell_toread);
  fclose(fileID);
  fclose(iscellID);
end
data_names = cell(size(MyData,2),4); %create 2xn cell array to concatenate video names and corresponding traces
data_names(:,2) = MyData';
% create a table b/c matlab is weird and this is the best way i can think
% to access folder names
folder_data = struct2table(SubFold);
folder_data.name = str2double(folder_data.name);
% turn folder names into cell array for input into data array
% % names = cellstr(folder_data.name(3:end));
% % % 2 column array with video name and corresponding video data from suite2p
data_names(:,1) = num2cell(folder_data.name(3:end));
end