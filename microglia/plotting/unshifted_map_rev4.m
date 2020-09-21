function [soma_data, process_data] = unshifted_map_rev4()
% this function takes the output from the "excel_matching" function and
% creates a heatmap proportional to the number of cell arrays within
% excel_matching. Nested sorting built into function. Written by John Wagner. I do not claim that this is an
% optimal solution for our data, but I do ensure that it is interesting and
% a bit convoluted.

[fileName, pathName] = uigetfile('*.xlsx');
% Path of the main folder : Rotation
main_suite2p = uigetdir;
[soma_data, soma_table] = soma_matching(fileName, pathName, main_suite2p);
[process_data, pro_table] = pro_matching(fileName, pathName, main_suite2p);
%% SOMA shifting: start by lining up Ca traces based on soma_timeHitRadius1
soma_frames = max(cellfun('size',soma_data(:,3),2));
fs = 1.04; % frames per second (sampling rate)
max_soma_time = soma_frames./fs;
% initialize sorted_F and insert each vector w/ for loop
sorted_F = zeros(size(soma_data,1),soma_frames);
% reason for this is to accomodate traces of different lengths, shorter
% traces will be padded with zeros
for ii = 1:size(sorted_F,1)
temp_F = cell2mat(soma_data(ii,3));
sorted_F(ii,1:length(temp_F)) = temp_F;
end
clear temp_F % conserve some memory
% % create dF/F using first 50 seconds
soma_dF_sorted = zeros(size(sorted_F));
soma_timeHitRadius1 = cell2mat(soma_data(:,6)); % soma_time in seconds
soma_time_shifts = round(soma_timeHitRadius1.*fs);
n = 0; % how much to move signal back in soma_time
% creating soma_time vectors because we'll use these for labeling horiz. axis
soma_time = [-soma_frames:soma_frames]./fs;
soma_time = round(soma_time-n-1,1);
soma_time = soma_time(rem(soma_time, 100) == 0);
soma_time_label = string([soma_time]);
for jj = 1:size(sorted_F,1)
% % calculate dF/F first
    F0 = mean(sorted_F(jj,1:50));
    soma_dF_sorted(jj,:) = (sorted_F(jj,:)-F0)./F0;
end

%% PROCESS shifting
process_frames = max(cellfun('size',process_data(:,3),2));
fs = 1.04; % frames per second (sampling rate)
max_process_time = process_frames./fs;
% initialize sorted_F and insert each vector w/ for loop
sorted_F = zeros(size(process_data,1),process_frames);
% reason for this is to accomodate traces of different lengths, shorter
% traces will be padded with zeros
for ii = 1:size(sorted_F,1)
temp_F = cell2mat(process_data(ii,3));
sorted_F(ii,1:length(temp_F)) = temp_F;
end
clear temp_F % conserve some memory (marginal, but good practice)
% % create dF/F using first 50 seconds
process_dF_sorted = zeros(size(sorted_F));
process_timeHitRadius1 = cell2mat(process_data(:,6)); % process_time in seconds
process_time_shifts = round(process_timeHitRadius1.*fs);
n = 0; % how much to move signal back in process_time
% creating process_time vectors because we'll use these for labeling horiz. axis
process_time = [-process_frames:process_frames]./fs;
process_time = round(process_time-n-1,1);
process_time = process_time(rem(process_time, 100) == 0);
process_time_label = string([process_time]);
for jj = 1:size(sorted_F,1)
% % calculate dF/F first
    F0 = mean(sorted_F(jj,1:50));
    process_dF_sorted(jj,:) = (sorted_F(jj,:)-F0)./F0;
end
%% this section just extracts the sorted traces into separate arrays based on soma/process and pathology
% % % soma PBS/TMEV
soma_TMEV_logicals = contains(soma_data(:,4),'TMEV'); 
% rows where TMEV occurs
TMEV_rows = find(soma_TMEV_logicals);
soma_TMEV_dF = soma_dF_sorted(soma_TMEV_logicals,:);
% extract main table data
soma_TMEV_extracted = soma_data(soma_TMEV_logicals,:);
for k = 1:size(soma_TMEV_dF,1)
soma_TMEV_extracted{k,3} = soma_TMEV_dF(k,:);
end
% repeat for PBS group
soma_PBS_logicals = contains(soma_data(:,4),'PBS'); 
PBS_rows = find(soma_PBS_logicals);
soma_PBS_dF = soma_dF_sorted(soma_PBS_logicals,:);
soma_PBS_extracted = soma_data(soma_PBS_logicals,:);
for k = 1:size(soma_PBS_dF,1)
soma_PBS_extracted{k,3} = soma_PBS_dF(k,:);
end

soma_PBS_videos = soma_PBS_extracted(:,1);
soma_PBS_dpi = num2str(cell2mat(soma_PBS_extracted(:,5)));
soma_PBS_ROI = num2str(cell2mat(soma_PBS_extracted(:,2)));
soma_PBS_labels = cell(size(soma_PBS_videos,1),1);
for i = 1:length(soma_PBS_labels)
soma_PBS_labels(i) = strcat(soma_PBS_videos(i),' , ', soma_PBS_dpi(i,:),' , ', soma_PBS_ROI(i,:));
end

soma_TMEV_videos = soma_TMEV_extracted(:,1);
soma_TMEV_dpi = num2str(cell2mat(soma_TMEV_extracted(:,5)));
soma_TMEV_ROI = num2str(cell2mat(soma_TMEV_extracted(:,2)));
soma_TMEV_labels = cell(size(soma_TMEV_videos,1),1);
for i = 1:length(soma_TMEV_labels)
soma_TMEV_labels(i) = strcat(soma_TMEV_videos(i),' , ', soma_TMEV_dpi(i,:),' , ', soma_TMEV_ROI(i,:));
end

p_TMEV_logicals = contains(process_data(:,4),'TMEV'); 
% rows where TMEV occurs
TMEV_rows = find(p_TMEV_logicals);
p_TMEV_dF = process_dF_sorted(p_TMEV_logicals,:);
% extract main table data
p_TMEV_extracted = process_data(p_TMEV_logicals,:);
for k = 1:size(p_TMEV_dF,1)
p_TMEV_extracted{k,3} = p_TMEV_dF(k,:);
end
% repeat for PBS group
p_PBS_logicals = contains(process_data(:,4),'PBS'); 
PBS_rows = find(p_PBS_logicals);
p_PBS_dF = process_dF_sorted(p_PBS_logicals,:);
p_PBS_extracted = process_data(p_PBS_logicals,:);
for k = 1:size(p_PBS_dF,1)
p_PBS_extracted{k,3} = p_PBS_dF(k,:);
end


p_PBS_videos = p_PBS_extracted(:,1);
p_PBS_dpi = num2str(cell2mat(p_PBS_extracted(:,5)));
p_PBS_ROI = num2str(cell2mat(p_PBS_extracted(:,2)));
p_PBS_labels = cell(size(p_PBS_videos,1),1);
for i = 1:length(p_PBS_labels)
p_PBS_labels(i) = strcat(p_PBS_videos(i),' , ', p_PBS_dpi(i,:),' , ', p_PBS_ROI(i,:));
end

p_TMEV_videos = p_TMEV_extracted(:,1);
p_TMEV_dpi = num2str(cell2mat(p_TMEV_extracted(:,5)));
p_TMEV_ROI = num2str(cell2mat(p_TMEV_extracted(:,2)));
p_TMEV_labels = cell(size(p_TMEV_videos,1),size(p_TMEV_videos,2));
for i = 1:length(p_TMEV_labels)
p_TMEV_labels(i) = strcat(p_TMEV_videos(i),' , ', p_TMEV_dpi(i,:),' , ', p_TMEV_ROI(i,:));
end

%% let the heatmapping begin!
% % % SOMA PBS % % %
figure()
imagesc([-max_soma_time max_soma_time], [1 size(soma_PBS_dF,1)],soma_PBS_dF);
yticks([1:size(soma_PBS_videos,1)]);
yticklabels(soma_PBS_labels);
% Make new tick marks.
xlim([-round(max_soma_time) round(max_soma_time)]);
xl = xlim(); % Get existing range.
numberOfXTickMarks = length(soma_time);
xTickLocations = linspace(xl(1), xl(2), numberOfXTickMarks);
set(gca,'XTick', xTickLocations);
% Make new labels for the new tick marks
set(gca,'XTickLabel', soma_time_label);
set(gca,'YDir','reverse');
title('PBS Soma'); xlabel('Time (sec)'); ylabel('Video, DPI, ROI #');
set(gca,'TickLabelInterpreter','none');
lowestValue = min(soma_PBS_dF(:));
highestValue = max(soma_PBS_dF(:));
myColorMap = jet(256);
caxis(gca,[lowestValue*1.5, highestValue*1.5]);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar
hold on
line([0, 0], [0 size(soma_PBS_dF,2)], 'Color', 'red','LineStyle','--');
saveas(gcf,fullfile(pathName,'U_Soma_PBS'));
% % % RROCESS PBS % % %
% ax(2) = subplot(2,2,2); % top right plot will be process pbs
figure()
imagesc([-max_process_time max_process_time], [1 size(p_PBS_dF,1)], p_PBS_dF);
yticks([1:size(p_PBS_videos,1)]);
yticklabels(p_PBS_labels);
% Make new tick marks.
xlim([-round(max_process_time) round(max_process_time)]);
xl = xlim(); % Get existing range.
xTickLocations = linspace(xl(1), xl(2), numberOfXTickMarks);
set(gca,'XTick', xTickLocations);
% Make new labels for the new tick marks
set(gca,'XTickLabel', process_time_label);
set(gca,'YDir','reverse');
title('PBS Process'); xlabel('Time (sec)'); ylabel('Video, DPI, ROI #');
set(gca,'TickLabelInterpreter','none');
myColorMap = jet(256);
caxis(gca,[lowestValue*1.5, highestValue*1.5]);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar
hold on;
line([0, 0], [0 size(p_PBS_dF,2)], 'Color', 'red','LineStyle','--');
saveas(gcf,fullfile(pathName,'U_Process_PBS'));
% % % SOMA TMEV % % %
% ax(3) = subplot(2,2,3); % bottom left plot will be soma tmev
figure()
imagesc([-max_soma_time max_soma_time], [1 size(soma_TMEV_dF,1)], soma_TMEV_dF);
yticks([1:size(soma_TMEV_videos,1)]);
yticklabels(soma_TMEV_labels);
% Make new tick marks.
xlim([-round(max_soma_time) round(max_soma_time)]);
xl = xlim(); % Get existing range.
xTickLocations = linspace(xl(1), xl(2), numberOfXTickMarks);
set(gca,'XTick', xTickLocations);
% Make new labels for the new tick marks
set(gca,'XTickLabel', soma_time_label);
set(gca,'YDir','reverse');
title('TMEV Soma'); xlabel('Time (sec)'); ylabel('Video, DPI, ROI #');
set(gca,'TickLabelInterpreter','none');
myColorMap = jet(256);
caxis(gca,[lowestValue*1.5, highestValue*1.5]);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar
hold on;
line([0, 0], [0 size(soma_TMEV_dF,2)], 'Color', 'red','LineStyle','--');
saveas(gcf,fullfile(pathName,'U_Soma_TMEV'));
% % % PROCESS TMEV % % %
% ax(4) = subplot(2,2,4); % bottom right plot will be process tmev
figure()
h = imagesc([-max_process_time max_process_time], [1 size(p_TMEV_dF,1)], p_TMEV_dF);
h.AlphaData = ones(size(h.CData));
h.AlphaData(isnan(h.CData)) = 0;
yticks([1:size(p_TMEV_videos,1)]);
yticklabels(p_TMEV_labels);
% Make new tick marks.
xlim([-round(max_process_time) round(max_process_time)]);
xl = xlim(); % Get existing range.
xTickLocations = linspace(xl(1), xl(2), numberOfXTickMarks);
set(gca,'XTick', xTickLocations);
% Make new labels for the new tick marks
set(gca,'XTickLabel', process_time_label);
set(gca,'YDir','reverse');
title('TMEV Process'); xlabel('Time (sec)'); ylabel('Video, DPI, ROI #');
set(gca,'TickLabelInterpreter','none');
myColorMap = jet(256);
caxis(gca,[lowestValue*1.5, highestValue*1.5]);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar
hold on;
line([0, 0], [0 size(p_TMEV_dF,2)], 'Color', 'red','LineStyle','--');
saveas(gcf,fullfile(pathName, 'U_Process_TMEV'));
end
    