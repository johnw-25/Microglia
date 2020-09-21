%% use this script to call excel_matching %%
% obtain file and path for excel workbook to pull data from
[fileName, pathName] = uigetfile('*.xlsx');

% first query is an integer for excel sheet number
% second is the directory for Fall.mat files
[sorted_data,excel_table1] = excel_matching(fileName, pathName);

[sorted_data2,excel_table2] = second_best_matching(fileName, pathName);

[sorted_data3,excel_table3] = third_best_matching(fileName, pathName);



%% shifted
map_ca_rev4(sorted_data2);
% unshifted_map_rev3(sorted_data);
%% unshifted
