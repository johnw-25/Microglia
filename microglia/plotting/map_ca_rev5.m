function map_ca_rev5(soma_data, process_data)
% this function takes the output from the "excel_matching" function and
% creates a heatmap proportional to the number of cell arrays within
% excel_matching. Nested sorting built into function. Written by John Wagner. I do not claim that this is an
% optimal solution for our data, but I do ensure that it is interesting and
% a bit convoluted.
%% SOMA shifting: start by lining up Ca traces based on soma_timeHitRadius1
soma_frames = max(cellfun('size',soma_data(:,3),2));
fs = 0.97; % frames per second (sampling rate)
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
include = logical(cell2mat(soma_data(:,8)));
sorted_F = sorted_F(include,:);
% % create dF/F using first 50 seconds
dF_sorted = zeros(size(sorted_F));
soma_timeHitRadius1 = cell2mat(soma_data(:,6)); % soma_time in seconds
soma_timeHitRadius1 = soma_timeHitRadius1(include,:);
soma_time_shifts = round(soma_timeHitRadius1.*fs);
n = 0; % how much to move signal back in soma_time
% creating soma_time vectors because we'll use these for labeling horiz. axis
soma_time = [-soma_frames:soma_frames]./fs;
soma_time = round(soma_time-n-1,1);
soma_time = soma_time(rem(soma_time, 100) == 0);
soma_time_label = string([soma_time]);
dF_shifted_int = zeros(size(sorted_F,1),soma_frames);
before_soma_timeHitRadius1 = zeros(size(sorted_F,1),soma_frames);
before_soma_timeHitRadius1(:,:) = NaN;
for jj = 1:size(sorted_F,1)
% % calculate dF/F first
    F0 = mean(sorted_F(jj,1:50));
    dF_sorted(jj,:) = (sorted_F(jj,:)-F0)./F0;
    % % Do an initial shift to line up soma_timeHitRadius1's
    % % absolute value makes sense. think about circles
    dF_first_shift(jj,:) = circshift(dF_sorted(jj,:),-soma_time_shifts(jj)+n+1);
    % % extract position where the original beginning of trace got shifted
    % to
    test_index = abs(size(dF_sorted(jj,:),2)-soma_time_shifts(jj)+2):size(dF_sorted(jj,:),2);
    % % extract the stuff before soma_timeHitRadius1 for each trace
    before_soma_timeHitRadius1(jj,end-size(dF_first_shift(jj,test_index),2)+1:end) = dF_first_shift(jj,test_index);
    % % replace the shifted stuff with NaN so we don't repeat data when
    % visualizing
    dF_first_shift(jj,test_index) = NaN;
    % store the final shifted product(includes soma_timeHitRadius1 until end of
    % video)in temp variable (this might be a redundancy)
    temp_dF = dF_first_shift(jj,:);
    dF_shifted_int(jj,:) = [temp_dF];
end
% concactenate the before_soma_timeHitRadius1 and shifted trace to obtain the entire
% matrix of Ca traces, all lined up by soma_timeHitRadius1 (this took way too long)
soma_dF_shifted = [before_soma_timeHitRadius1, dF_shifted_int];
for kk = 1:size(soma_dF_shifted,1)
    soma_data{kk, 7} = soma_dF_shifted(kk,:);
end

%% PROCESS shifting
process_frames = max(cellfun('size',process_data(:,3),2));
fs = 0.97; % frames per second (sampling rate)
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
pinclude = logical(cell2mat(process_data(:,8)));
sorted_F = sorted_F(pinclude,:);
dF_sorted = zeros(size(sorted_F));
process_timeHitRadius1 = cell2mat(process_data(:,6)); % process_time in seconds
process_timeHitRadius1 = process_timeHitRadius1(pinclude,:);
process_time_shifts = round(process_timeHitRadius1.*fs);
n = 0; % how much to move signal back in process_time
% creating process_time vectors because we'll use these for labeling horiz. axis
process_time = [-process_frames:process_frames]./fs;
process_time = round(process_time-n-1,1);
process_time = process_time(rem(process_time, 100) == 0);
process_time_label = string([process_time]);
dF_shifted_int = zeros(size(sorted_F,1),process_frames);
before_process_timeHitRadius1 = zeros(size(sorted_F,1),process_frames);
before_process_timeHitRadius1(:,:) = NaN;
for jj = 1:size(sorted_F,1)
% % calculate dF/F first
    F0 = mean(sorted_F(jj,1:50));
    dF_sorted(jj,:) = (sorted_F(jj,:)-F0)./F0;
    % % Do an initial shift to line up process_timeHitRadius1's
    % % absolute value makes sense. think about circles
    dF_first_shift(jj,:) = circshift(dF_sorted(jj,:),-process_time_shifts(jj)+n+1);
    % % extract position where the original beginning of trace got shifted
    % to
    test_index = abs(size(dF_sorted(jj,:),2)-process_time_shifts(jj)+2):size(dF_sorted(jj,:),2);
    % % extract the stuff before process_timeHitRadius1 for each trace
    before_process_timeHitRadius1(jj,end-size(dF_first_shift(jj,test_index),2)+1:end) = dF_first_shift(jj,test_index);
    % % replace the shifted stuff with NaN so we don't repeat data when
    % visualizing
    dF_first_shift(jj,test_index) = NaN;
    % store the final shifted product(includes process_timeHitRadius1 until end of
    % video)in temp variable (this might be a redundancy)
    temp_dF = dF_first_shift(jj,:);
    dF_shifted_int(jj,:) = [temp_dF];
end
% concactenate the before_process_timeHitRadius1 and shifted trace to obtain the entire
% matrix of Ca traces, all lined up by process_timeHitRadius1 (this took way too long)
process_dF_shifted = [before_process_timeHitRadius1, dF_shifted_int];
for kk = 1:size(process_dF_shifted,1)
    process_data{kk, 7} = process_dF_shifted(kk,:);
end
%% this section just extracts the sorted traces into separate arrays based on soma/process and pathology
% % % soma PBS/TMEV
soma_TMEV_logicals = contains(soma_data(:,4),'TMEV');
somaTmevInclude = soma_TMEV_logicals(include);
% rows where TMEV occurs
soma_TMEV_dF = soma_dF_shifted(somaTmevInclude,:);
% extract main table data
soma_TMEV_extracted = soma_data(soma_TMEV_logicals,:);
for k = 1:size(soma_TMEV_dF,1)
soma_TMEV_extracted{k,3} = soma_TMEV_dF(k,:);
end
% repeat for PBS group
soma_PBS_logicals = contains(soma_data(:,4),'PBS');
somaPbsInclude = soma_PBS_logicals(include);
soma_PBS_dF = soma_dF_shifted(somaPbsInclude,:);
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
pTmevInclude = p_TMEV_logicals(pinclude);
p_TMEV_dF = process_dF_shifted(pTmevInclude,:);
% extract main table data
p_TMEV_extracted = process_data(p_TMEV_logicals,:);
for k = 1:size(p_TMEV_dF,1)
p_TMEV_extracted{k,3} = p_TMEV_dF(k,:);
end
% repeat for PBS group
p_PBS_logicals = contains(process_data(:,4),'PBS'); 
pPbsInclude = p_PBS_logicals(pinclude);
p_PBS_dF = process_dF_shifted(pPbsInclude,:);
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

% % Find max and mins for heatmap color axis % %
%outliers
% somaTMEVmax = max(soma_TMEV_dF);
% somaTMEVmin = min(soma_TMEV_dF);
% somaPBSmax = max(soma_PBS_dF);
% somaPBSmin = min(soma_PBS_dF);
% pTMEVmax = max(p_TMEV_dF);
% pTMEVmin = min(p_TMEV_dF);
% pPBSmax = max(p_PBS_dF);
% pPBSmin = min(p_PBS_dF);
%+3SDs above mean
somaTMEVmax = mean(soma_TMEV_dF) + 3*std(soma_TMEV_dF);
somaTMEVmin = mean(soma_TMEV_dF) - 1.8*std(soma_TMEV_dF);
somaPBSmax = mean(soma_PBS_dF) + 3*std(soma_PBS_dF);
somaPBSmin = mean(soma_PBS_dF) - 1.8*std(soma_PBS_dF);
pTMEVmax = mean(p_TMEV_dF) + 3*std(p_TMEV_dF);
pTMEVmin = mean(p_TMEV_dF) - 1.8*std(p_TMEV_dF);
pPBSmax = mean(p_PBS_dF) + 3*std(p_PBS_dF);
pPBSmin = mean(p_PBS_dF) - 1.8*std(p_PBS_dF);

highestValue = max([somaTMEVmax, somaPBSmax, pTMEVmax, pPBSmax]);
lowestValue = min([somaTMEVmin, somaPBSmin, pTMEVmin, pPBSmin]);

highestValue = 6;
lowestValue = -2;

highestYAxis = max([size(soma_PBS_dF,1),size(soma_TMEV_dF,1), size(p_PBS_dF,1), size(p_TMEV_dF,1)]); 
%% let the heatmapping begin!
% % % SOMA PBS % % %
somaPbs = figure();
h1 = imagesc([-max_soma_time max_soma_time], [1 size(soma_PBS_dF,1)],soma_PBS_dF);
axisLines = gca;
set(axisLines,'Parent',somaPbs,'YDir','normal', ...
    'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255], ...
    'FontSize',10,'FontName','Arial','Color','none');
% yticks([1:size(soma_PBS_videos,1)]);
% yticklabels(soma_PBS_labels);
% Make new tick marks.
xlim([-round(max_soma_time) round(max_soma_time)]);
xl = xlim(); % Get existing range.
ylim([0 highestYAxis]);
xTickLocations = [round(xl(1)) 0 round(xl(2))];
yTickLocations = [round(size(soma_PBS_dF,1)/4), round(size(soma_PBS_dF,1)/1.5)];
set(axisLines,'XTick', xTickLocations, ...
    'YTick',yTickLocations);
% Make new labels for the new tick marks
set(axisLines,'XTickLabel', {num2str(round(xl(1)/60)), num2str(0), num2str(round(xl(2)/60))});
title('PBS Soma','FontWeight','Normal'); %xlabel(labels,'Time (min)'); ylabel('ROI');
set(axisLines,'TickLabelInterpreter','none');
set(h1,'AlphaData',~isnan(soma_PBS_dF));
set(gcf, 'Units','inches','position',[4 4 1.5 1.5]);
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 1.5 1.5]);
myColorMap = parula(32);
caxis(gca,[lowestValue, highestValue]);
colormap(myColorMap);
colorbar
hold on
line(axisLines,[0, 0], [0 size(soma_PBS_dF,2)], 'Color', 'red','LineStyle','--','LineWidth',1);
savefig('Heatmap_Soma_PBS.fig');

% % % RROCESS PBS % % %
% ax(2) = subplot(2,2,2); % top right plot will be process pbs
figure()
h2 = imagesc([-max_process_time max_process_time], [1 size(p_PBS_dF,1)], p_PBS_dF);
% yticks([1:size(p_PBS_videos,1)]);
% yticklabels(p_PBS_labels);
% Make new tick marks.
xlim([-round(max_process_time) round(max_process_time)]);
xl = xlim(); % Get existing range.
xTickLocations = [round(xl(1)) 0 round(xl(2))];
set(gca,'XTick', xTickLocations);
% Make new labels for the new tick marks
set(gca,'XTickLabel', {num2str(round(xl(1)/60)), num2str(0), num2str(round(xl(2)/60))});
set(gca,'YDir','reverse');
title('PBS Process','FontWeight','Normal'); xlabel('Time (sec)'); ylabel('Video, DPI, ROI #');
set(gca,'TickLabelInterpreter','none');
set(h2,'AlphaData',~isnan(p_PBS_dF));
set(gcf, 'Units','inches','position',[4 4 1.5 1.5]);
set(gca,'FontSize',10,'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 1.5 1.5]);
caxis(gca,[lowestValue, highestValue]);
colormap(myColorMap);
colorbar
hold on;
line([0, 0], [0 size(p_PBS_dF,2)], 'Color', 'red','LineStyle','--','LineWidth',1);
savefig('Heatmap_Processes_PBS.fig');
% % % SOMA TMEV % % %
% ax(3) = subplot(2,2,3); % bottom left plot will be soma tmev
figure()
h3 = imagesc([-max_soma_time max_soma_time], [1 size(soma_TMEV_dF,1)], soma_TMEV_dF);
% yticks([1:size(soma_TMEV_videos,1)]);
% yticklabels(soma_TMEV_labels);
% Make new tick marks.
xlim([-round(max_soma_time) round(max_soma_time)]);
xl = xlim(); % Get existing range.
xTickLocations = [round(xl(1)) 0 round(xl(2))];
set(gca,'XTick', xTickLocations);
% Make new labels for the new tick marks
set(gca,'XTickLabel', {num2str(round(xl(1)/60)), num2str(0), num2str(round(xl(2)/60))});
set(gca,'YDir','reverse');
title('TMEV Soma','FontWeight','Normal'); xlabel('Time (sec)'); ylabel('Video, DPI, ROI #');
set(gca,'TickLabelInterpreter','none');
set(h3,'AlphaData',~isnan(soma_TMEV_dF));
set(gcf, 'Units','inches','position',[4 4 1.5 1.5]);
set(gca,'FontSize',10,'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 1.5 1.5]);
caxis(gca,[lowestValue, highestValue]);
colormap(myColorMap);
colorbar
hold on;
line([0, 0], [0 size(soma_TMEV_dF,2)], 'Color', 'red','LineStyle','--','LineWidth',1);
savefig('Heatmap_Soma_TMEV.fig');
% % % PROCESS TMEV % % %
% ax(4) = subplot(2,2,4); % bottom right plot will be process tmev
figure()
h4 = imagesc([-max_process_time max_process_time], [1 size(p_TMEV_dF,1)], p_TMEV_dF);
% yticks([1:size(p_TMEV_videos,1)]);
% yticklabels(p_TMEV_labels);
% Make new tick marks.
xlim([-round(max_process_time) round(max_process_time)]);
xl = xlim(); % Get existing range.
xTickLocations = [round(xl(1)) 0 round(xl(2))];
set(gca,'XTick', xTickLocations);
% Make new labels for the new tick marks
set(gca,'XTickLabel', {num2str(round(xl(1)/60)), num2str(0), num2str(round(xl(2)/60))});
set(gca,'YDir','reverse');
title('TMEV Process','FontWeight','Normal'); xlabel('Time (sec)'); ylabel('Video, DPI, ROI #');
set(gca,'TickLabelInterpreter','none');
set(h4,'AlphaData',~isnan(p_TMEV_dF));
set(gcf, 'Units','inches','position',[4 4 1.5 1.5]);
set(gca,'FontSize',10,'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 1.5 1.5]);
caxis(gca,[lowestValue, highestValue]);
colormap(myColorMap);
colorbar
hold on;
line([0, 0], [0 size(p_TMEV_dF,2)], 'Color', 'red','LineStyle','--','LineWidth',1);
savefig('Heatmap_Processes_TMEV.fig');
end
    