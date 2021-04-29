function map_ca_rev6(soma_data, process_data)
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
somaF = sorted_F;
% % create dF/F using first 50 seconds
soma_timeHitRadius1 = cell2mat(soma_data(:,6)); % soma_time in seconds
soma_timeHitRadius1 = soma_timeHitRadius1(include,:);

% % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % %
% % uncomment this to do just phase 2-3 % % % % %
sorted_F = (sorted_F - mean(sorted_F(:,1:50), 2))./mean(sorted_F(:,1:50), 2);
soma_timeHitRadius1 = soma_timeHitRadius1 - 90; %
sorted_F = sorted_F(:,91:end);  % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % %
% matrixSortBy = mean(sorted_F(:,500:end), 2);
% [~, sortIdx] = sort(matrixSortBy, 'descend');
% sorted_F = sorted_F(sortIdx,:);

dF_sorted = zeros(size(sorted_F));
soma_time_shifts = round(soma_timeHitRadius1.*fs);
n = 0; % how much to move signal back in soma_time
% creating soma_time vectors because we'll use these for labeling horiz. axis
soma_time = [-soma_frames:soma_frames]./fs;
soma_time = round(soma_time-n-1,1);
soma_time = soma_time(rem(soma_time, 100) == 0);
soma_time_label = string([soma_time]);
dF_shifted_int = zeros(size(sorted_F,1),soma_frames - 90);
before_soma_timeHitRadius1 = zeros(size(sorted_F,1),soma_frames);
before_soma_timeHitRadius1(:,:) = NaN;
for jj = 1:size(sorted_F,1)
% % calculate dF/F first
    F0 = mean(sorted_F(jj,1:50));
    dF_sorted(jj,:) = ((sorted_F(jj,:)));%-F0)./F0;
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
processF = sorted_F;

process_timeHitRadius1 = cell2mat(process_data(:,6)); % process_time in seconds
process_timeHitRadius1 = process_timeHitRadius1(pinclude,:);

% % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % %
% % uncomment this to do just phase 2-3 % % % % %
F0 = mean(sorted_F(:,1:50), 2);
sorted_F = (sorted_F - F0)./F0;
process_timeHitRadius1 = process_timeHitRadius1 - 90; %
sorted_F = sorted_F(:,91:end);  % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % %
% matrixSortBy = mean(sorted_F(:,500:1000), 2);
% [~, sortIdx] = sort(matrixSortBy, 'descend');
% sorted_F = sorted_F(sortIdx,:);

dF_sorted = zeros(size(sorted_F));
process_time_shifts = round(process_timeHitRadius1.*fs);
n = 0; % how much to move signal back in process_time
% creating process_time vectors because we'll use these for labeling horiz. axis
process_time = [-process_frames:process_frames]./fs;
process_time = round(process_time-n-1,1);
process_time = process_time(rem(process_time, 100) == 0);
process_time_label = string([process_time]);
dF_shifted_int = zeros(size(sorted_F,1),process_frames - 90);
before_process_timeHitRadius1 = zeros(size(sorted_F,1),process_frames);
before_process_timeHitRadius1(:,:) = NaN;
for jj = 1:size(sorted_F,1)
% % calculate dF/F first
    F0 = mean(sorted_F(jj,1:50));
    dF_sorted(jj,:) = ((sorted_F(jj,:)));%-F0)./F0;
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
allDPI = cell2mat(soma_TMEV_extracted(:,5));
soma_TMEV2_dF = soma_TMEV_dF(1:19,:);
soma_TMEV5_dF = soma_TMEV_dF(20:40,:);
soma_TMEV15_dF = soma_TMEV_dF(41:end,:);

% repeat for PBS group
soma_PBS_logicals = contains(soma_data(:,4),'PBS');
somaPbsInclude = soma_PBS_logicals(include);
soma_PBS_dF = soma_dF_shifted(somaPbsInclude,:);
soma_PBS_extracted = soma_data(soma_PBS_logicals,:);
for k = 1:size(soma_PBS_dF,1)
soma_PBS_extracted{k,3} = soma_PBS_dF(k,:);
end

% uncomment for df/f in phase0-1
for x = 1:size(somaF,1)
    somaF(x,:) = (somaF(x,:) - mean(somaF(x,1:50)))./mean(somaF(x,1:50));
end

% % Unshifted stuff % %
somaTmev = somaF(somaTmevInclude,:);
somaTmev2 = somaTmev(1:19,:);
somaTmev5 = somaTmev(20:40,:);
somaTmev15 = somaTmev(41:end,:);
somaPbs = somaF(somaPbsInclude,:);

% % sort shifted soma df/f arrays % %
% soma pbs
matrixSortBy = mean(somaPbs(:,500:1000), 2);
[~, sortIdx] = sort(matrixSortBy, 'descend');
soma_PBS_dF = soma_PBS_dF(sortIdx,:);
% soma tmev 2
matrixSortBy = mean(somaTmev2(:,500:1000), 2);
[~, sortIdx] = sort(matrixSortBy, 'descend');
soma_TMEV2_dF = soma_TMEV2_dF(sortIdx,:);
% soma tmev 5
matrixSortBy = mean(somaTmev5(:,500:1000), 2);
[~, sortIdx] = sort(matrixSortBy, 'descend');
soma_TMEV5_dF = soma_TMEV5_dF(sortIdx,:);
% soma tmev 15
matrixSortBy = mean(somaTmev15(:,500:1000), 2);
[~, sortIdx] = sort(matrixSortBy, 'descend');
soma_TMEV15_dF = soma_TMEV15_dF(sortIdx,:);

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
p_TMEV2_dF = p_TMEV_dF(1:64,:);
p_TMEV5_dF = p_TMEV_dF(65:80,:);
p_TMEV15_dF = p_TMEV_dF(81:end,:);
% repeat for PBS group
p_PBS_logicals = contains(process_data(:,4),'PBS'); 
pPbsInclude = p_PBS_logicals(pinclude);
p_PBS_dF = process_dF_shifted(pPbsInclude,:);
p_PBS_extracted = process_data(p_PBS_logicals,:);
for k = 1:size(p_PBS_dF,1)
p_PBS_extracted{k,3} = p_PBS_dF(k,:);
end

% uncomment for df/f in phase0-1
for x = 1:size(processF,1)
    processF(x,:) = (processF(x,:) - mean(processF(x,1:50)))./mean(processF(x,1:50));
end
% % Unshifted stuff % %
processTmev = processF(pTmevInclude,:);
processTmev2 = processTmev(1:64,:);
processTmev5 = processTmev(65:80,:);
processTmev15 = processTmev(81:end,:);
processPbs = processF(pPbsInclude,:);


% % sort shifted p df/f arrays % %
% p pbs
matrixSortBy = mean(processPbs(:,500:1000), 2);
[~, sortIdx] = sort(matrixSortBy, 'descend');
p_PBS_dF = p_PBS_dF(sortIdx,:);
% p tmev 2
matrixSortBy = mean(processTmev2(:,500:1000), 2);
[~, sortIdx] = sort(matrixSortBy, 'descend');
p_TMEV2_dF = p_TMEV2_dF(sortIdx,:);
% p tmev 5
matrixSortBy = mean(processTmev5(:,500:1000), 2);
[~, sortIdx] = sort(matrixSortBy, 'descend');
p_TMEV5_dF = p_TMEV5_dF(sortIdx,:);
% p tmev 15
matrixSortBy = mean(processTmev15(:,500:1000), 2);
[~, sortIdx] = sort(matrixSortBy, 'descend');
p_TMEV15_dF = p_TMEV15_dF(sortIdx,:);

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
lowestValue = -1;

% uncomment for raw f
% highestValue = 1600;
% lowestValue = 0;

highestYAxis = max([size(soma_PBS_dF,1),size(soma_TMEV_dF,1), size(p_PBS_dF,1), size(p_TMEV_dF,1)]); 
%% let the heatmapping begin!
time = [round(-1740/0.97), round(1740/0.97)]; %seconds
% myColorMap = parula(30);
myColorMap = parula(7);
myColorMap(3,:) = [159/256, 159/256, 95/256];
myColorMap(4,:) = [1, 1, 0];
myColorMap(5,:) = [255/256, 127/256, 0];
myColorMap(6,:) = [1, 0, 0];
myColorMap(7,:) = [139/256, 0, 0];
% very non saturated yellow: [159/256, 159/256, 95/256]
% % % SOMA PBS % % %
HMMatrix(time, size(soma_PBS_dF,1), soma_PBS_dF, highestYAxis,lowestValue, highestValue,'Soma PBS Phase 2-3 Sorted',myColorMap)
% % % RROCESS PBS % % %
HMMatrix(time, size(p_PBS_dF,1), p_PBS_dF, highestYAxis,lowestValue, highestValue,'Processes PBS 2-3 Sorted',myColorMap)
% % % SOMA TMEV % % %
HMMatrix(time, size(soma_TMEV2_dF,1), soma_TMEV2_dF, highestYAxis,lowestValue, highestValue,'Soma TMEV 2 DPI 2-3 Sorted',myColorMap)
HMMatrix(time, size(soma_TMEV5_dF,1), soma_TMEV5_dF, highestYAxis,lowestValue, highestValue,'Soma TMEV 5 DPI 2-3 Sorted',myColorMap)
HMMatrix(time, size(soma_TMEV15_dF,1), soma_TMEV15_dF, highestYAxis,lowestValue, highestValue,'Soma TMEV 15 DPI 2-3 Sorted',myColorMap)
% % % PROCESS TMEV % % %
HMMatrix(time, size(p_TMEV2_dF,1), p_TMEV2_dF, highestYAxis,lowestValue, highestValue,'Processes TMEV 2 DPI 2-3 Sorted',myColorMap)
HMMatrix(time, size(p_TMEV5_dF,1), p_TMEV5_dF, highestYAxis,lowestValue, highestValue,'Processes TMEV 5 DPI 2-3 Sorted',myColorMap)
HMMatrix(time, size(p_TMEV15_dF,1), p_TMEV15_dF, highestYAxis,lowestValue, highestValue,'Processes TMEV 15 DPI 2-3 Sorted',myColorMap)

%% Break up into phases
% time = [1, 90];
% %soma
% phase1SomaPbs = somaPbs(:,1:90);
% phase1SomaTmev2 = somaTmev2(:,1:90);
% phase1SomaTmev5 = somaTmev5(:,1:90);
% phase1SomaTmev15 = somaTmev15(:,1:90);
% 
% Phase01HeatMap(time, size(phase1SomaPbs,1), phase1SomaPbs, highestYAxis,lowestValue, highestValue,'Soma PBS Phase 0-1',myColorMap)
% Phase01HeatMap(time, size(phase1SomaTmev2,1), phase1SomaTmev2, highestYAxis,lowestValue, highestValue,'Soma TMEV2 Phase 0-1',myColorMap)
% Phase01HeatMap(time, size(phase1SomaTmev5,1), phase1SomaTmev5, highestYAxis,lowestValue, highestValue,'Soma TMEV5 Phase 0-1',myColorMap)
% Phase01HeatMap(time, size(phase1SomaTmev15,1), phase1SomaTmev15, highestYAxis,lowestValue, highestValue,'Soma TMEV15 Phase 0-1',myColorMap)
% 
% %process
% phase1ProcessPbs = processPbs(:,1:90);
% phase1ProcessTmev2 = processTmev2(:,1:90);
% phase1ProcessTmev5 = processTmev5(:,1:90);
% phase1ProcessTmev15 = processTmev15(:,1:90);
% 
% Phase01HeatMap(time, size(phase1ProcessPbs,1), phase1ProcessPbs, highestYAxis,lowestValue, highestValue,'Process PBS Phase 0-1',myColorMap)
% Phase01HeatMap(time, size(phase1ProcessTmev2,1), phase1ProcessTmev2, highestYAxis,lowestValue, highestValue,'Process TMEV2 Phase 0-1',myColorMap)
% Phase01HeatMap(time, size(phase1ProcessTmev5,1), phase1ProcessTmev5, highestYAxis,lowestValue, highestValue,'Process TMEV5 Phase 0-1',myColorMap)
% Phase01HeatMap(time, size(phase1ProcessTmev15,1), phase1ProcessTmev15, highestYAxis,lowestValue, highestValue,'Process TMEV15 Phase 0-1',myColorMap)
% 
% 
% % phase 2 - 3
% %soma
% time = [91, 1740];
% phase23SomaPbs = somaPbs(:,91:end);
% phase23SomaTmev2 = somaTmev2(:,91:end);
% phase23SomaTmev5 = somaTmev5(:,91:end);
% phase23SomaTmev15 = somaTmev15(:,91:end);
% 
% Phase23HeatMap(time, size(phase23SomaPbs,1), phase23SomaPbs, highestYAxis,lowestValue, highestValue,'Soma PBS Phase 2-3',myColorMap)
% Phase23HeatMap(time, size(phase23SomaTmev2,1), phase23SomaTmev2, highestYAxis,lowestValue, highestValue,'Soma TMEV2 Phase 2-3',myColorMap)
% Phase23HeatMap(time, size(phase23SomaTmev5,1), phase23SomaTmev5, highestYAxis,lowestValue, highestValue,'Soma TMEV5 Phase 2-3',myColorMap)
% Phase23HeatMap(time, size(phase23SomaTmev15,1), phase23SomaTmev15, highestYAxis,lowestValue, highestValue,'Soma TMEV15 Phase 2-3',myColorMap)
% 
% %process
% phase23ProcessPbs = processPbs(:,91:end);
% phase23ProcessTmev2 = processTmev2(:,91:end);
% phase23ProcessTmev5 = processTmev5(:,91:end);
% phase23ProcessTmev15 = processTmev15(:,91:end);
% 
% Phase23HeatMap(time, size(phase23ProcessPbs,1), phase23ProcessPbs, highestYAxis,lowestValue, highestValue,'Process PBS Phase 2-3',myColorMap)
% Phase23HeatMap(time, size(phase23ProcessTmev2,1), phase23ProcessTmev2, highestYAxis,lowestValue, highestValue,'Process TMEV2 Phase 2-3',myColorMap)
% Phase23HeatMap(time, size(phase23ProcessTmev5,1), phase23ProcessTmev5, highestYAxis,lowestValue, highestValue,'Process TMEV5 Phase 2-3',myColorMap)
% Phase23HeatMap(time, size(phase23ProcessTmev15,1), phase23ProcessTmev15, highestYAxis,lowestValue, highestValue,'Process TMEV15 Phase 2-3',myColorMap)
end
    