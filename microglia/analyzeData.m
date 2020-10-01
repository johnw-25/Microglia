function analyzeData(sortedData,process_data)
%%% This is a function to take in a sortedData structure and analyze the
%%% data within the groups.

% % % Remember we want to look at each subsequent metric in relative time

% dont have timehitR in struct yet but this is how we will do it:
%sortedData.Group.DPI_X.Video.RelativeLocs =
%sortedData.Group.DPI_X.Video.PeakLocs -
%sortedData.Group.DPI_X.Video.timehitR;

% bring groups up a layer to keep things less confusing
PBS = sortedData.PBS;
TMEV = sortedData.TMEV;

PBSfields = fieldnames(PBS);
allPBSLocs = []; %allocate space
allPBSPeaks = [];
allPeakDurs = [];
PBS_Alengths = [];
PBS_Blengths = [];
all_RAWauc = [];
all_auc = [];
allIEI = [];
allLocs = [];
before_0_auc = [];
after_0_auc = [];
before_0_peaks = [];
after_0_peaks = [];
before_0_durs = [];
after_0_durs = [];
frames = 1740;
allPBSTraces = {};
for k = 1:numel(PBSfields)
    % extract
    tempIEI = PBS.(PBSfields{k}).IEI;
    allIEI = [allIEI;tempIEI];
    
    temp_timehitR = PBS.(PBSfields{k}).timehitR;
    tempLocs = PBS.(PBSfields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    tempLocs = tempLocs(tempLocs > 70); % exclude burn
%     if ~isempty(tempLocs)
%     cutout = length(tempLocs);
%     tempLocs(cutout) = [];
%     tempLocs = tempLocs - temp_timehitR;
%     end
    allLocs = [allLocs;tempLocs];
    tempDurs = PBS.(PBSfields{k}).PeakDurs ./ 0.97;
    tempDurs = tempDurs(~isnan(tempDurs));
    tempA = PBS.(PBSfields{k}).LENGTH_A;
    tempB = PBS.(PBSfields{k}).LENGTH_B;
    tempRAWAUC = PBS.(PBSfields{k}).RAW_AUC;
    tempRAWAUC = abs(tempRAWAUC(~isnan(tempRAWAUC)));
    tempAUC = PBS.(PBSfields{k}).AUC;
    tempAUC = abs(tempAUC(~isnan(tempAUC)));
    tempPeaks = PBS.(PBSfields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    
    % update all trace names
    for i = 1:length(tempRAWAUC)
        allPBSTraces = [allPBSTraces; PBSfields(k)];
    end
    
    allPBSPeaks = [allPBSPeaks; tempPeaks];
    all_RAWauc = [all_RAWauc; tempRAWAUC];
    all_auc = [all_auc; tempAUC];
    shiftedLocs = tempLocs - temp_timehitR;
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));
    tempBefore0 = tempRAWAUC(shiftedLocs < 0 );
    tempAfter0 = tempRAWAUC(shiftedLocs >= 0 );
    B0Peaks = tempPeaks(shiftedLocs < 0);
    A0Peaks = tempPeaks(shiftedLocs >= 0);
    B0Durs = tempDurs(shiftedLocs < 0);
    A0Durs = tempDurs(shiftedLocs >= 0);
    before_0_durs = [before_0_durs; B0Durs];
    after_0_durs = [after_0_durs; A0Durs];
    before_0_auc = [before_0_auc; tempBefore0];
    after_0_auc = [after_0_auc; tempAfter0];
    before_0_peaks = [before_0_peaks; B0Peaks];
    after_0_peaks = [after_0_peaks; A0Peaks];
    PBS.(PBSfields{k}).shiftedLocs = shiftedLocs;
    
    allPBSLocs = [allPBSLocs; shiftedLocs];
    allPeakDurs = [allPeakDurs; tempDurs];
    PBS_Alengths = [PBS_Alengths; tempA];
    PBS_Blengths = [PBS_Blengths; tempB];
end

% % Stacked PBS Traces % %
figure()
xlabel('Relative time(s)'); ylabel('Trace Names');
for k = 1:numel(PBSfields)
    temp_timehitR = PBS.(PBSfields{k}).timehitR;
    R = round(temp_timehitR*0.97, 1);
    tempLocs = PBS.(PBSfields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    if temp_timehitR > frames 
        shiftedLocs = (tempLocs - R)./0.97;
    else
        shiftedLocs = (tempLocs - R)./0.97;
    end
    tempPeaks = PBS.(PBSfields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));   
    tempF = PBS.(PBSfields{k}).trace;
    F0 = mean(tempF(1:50));
    tempF = (tempF - F0) / F0;
    tempTime = (1:1740) - R;
    tempTime = tempTime./0.97;
    hold on
    plot(tempTime, tempF + k);
    hold on
    plot(shiftedLocs, tempF(tempLocs) + k, 'db', 'MarkerFaceColor','r','MarkerSize',6);
end
yticks([1:size(PBSfields)]);
yticklabels(PBSfields);
binEdges = [-1600:100:1600];
numBins = numel(binEdges);

% % Peak Magnitude vs Peak Duration % %
figure()
plot(allPeakDurs, allPBSPeaks,'bd','MarkerFaceColor','r');
xlabel('Peak durations (s)'); ylabel('Peak magnitudes (dF/F)'); title('PBS: Peak vs duration');

% % Peak Location Histogram % %
figure()
histogram(allPBSLocs,numBins, 'BinEdges', binEdges)
title('All Shifted PBS Peak Locations');
xlabel('Relative time (s)'); ylabel('Frequency of events');

% % Peak Duration Scatter % %
figure()
plot(allPBSLocs, allPeakDurs, 'dc','MarkerFaceColor','b')
title('PBS: Peak Durations in relative time'); xlabel('Relative time (s)'); ylabel('Peak Duration (s)');

% % Peak Magnitude Scatter % %
figure()
plot(allPBSLocs, allPBSPeaks, 'dc','MarkerFaceColor','b')
title('PBS: Peak Magnitudes in relative time'); xlabel('Relative time (s)'); ylabel('Peak Height (dF/F)');

figure()
% try out different data visulization of signal arc lengths
groupA = ones(size(PBS_Alengths));
groupB = ones(size(PBS_Blengths))+1;
dotA = [groupA,PBS_Alengths];
dotB = [groupB,PBS_Blengths];
plot(groupA, PBS_Alengths, 'MarkerSize', 12);
hold on
plot(groupB, PBS_Blengths);
legend('Length A', 'Length B');
title('PBS: Arc lengths of traces before and after timehitR'); xlabel('Group A / Group B'); ylabel('Arc length (dF/F)');
xticks([0, 1, 2, 3])
xlim([0 3])

% % filtered event auc scatter % %
figure()
plot(allPBSLocs, all_auc, 'dr' ,'MarkerFaceColor','r')
title('PBS: filtered signal auc'); xlabel('relative time(s)');

% % raw event auc scatter % %
fig = figure('DeleteFcn','doc datacursormode');
hold on
plot(allPBSLocs, all_RAWauc, 'dr' ,'MarkerFaceColor','b')
dcm_obj = datacursormode(fig);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,allPBSTraces})
title('PBS: raw signal auc'); xlabel('relative time(s)');

% % Plotting Inter-event Intervals % %
figure()
xlabel('Relative Time (s)'); ylabel('IEI (s) plotted @ first peak'); title('PBS: IEI');
for k = 1:length(PBSfields)
    tempIEI = PBS.(PBSfields{k}).IEI;
    tempLocs = PBS.(PBSfields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    tempLocs = tempLocs(1:end-1); % exclude last peak - plotting IEI on first peak in calculation
    tempLocs = tempLocs./0.97;
    R = PBS.(PBSfields{k}).timehitR;
    tempLocs = tempLocs - R;
    
    hold on
    plot(tempLocs, tempIEI, 'db','MarkerFaceColor','r');
end
    
    
    

% % boxplots of event auc before time zero and after % %
ALL = [before_0_auc; after_0_auc];
groups = [ones(size(before_0_auc)); ones(size(after_0_auc)) * 2];
figure()
boxplot(ALL,groups,'Whisker',4)
title('PBS Raw Event AUC before relative time zero and after time zero');
ylabel('AUC (dF/F)');

% % boxplots of event Peaks before time zero and after % %
ALL = [before_0_peaks; after_0_peaks];
figure()
boxplot(ALL,groups,'Whisker',4)
title('PBS Raw Event Peaks before relative time zero and after time zero');
ylabel('AUC (dF/F)');

% % boxplots of event durations before time zero and after % %
ALL = [before_0_durs; after_0_durs];
figure()
boxplot(ALL,groups,'Whisker',4)
title('PBS Raw Event Durations before relative time zero and after time zero');
ylabel('duration (s)');

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% TMEV  2 DPI %%%%%
%%%%%%%%%%%%%%%%%%%%%%%
TMEV_2DPI = TMEV.DPI_2;
TMEV2fields = fieldnames(TMEV_2DPI);
allTMEV2Locs = []; %allocate space
allTMEV2Peaks = [];
allPeakDurs = [];
TMEV2_Alengths = [];
TMEV2_Blengths = [];
all_RAWauc = [];
all_auc = [];
before_0_auc = [];
after_0_auc = [];
before_0_peaks = [];
after_0_peaks = [];
before_0_durs = [];
after_0_durs = [];
for k  = 1:numel(TMEV2fields)
    % extract
    temp_timehitR = TMEV_2DPI.(TMEV2fields{k}).timehitR;
    tempLocs = TMEV_2DPI.(TMEV2fields{k}).PeakLocs;
    tempDurs = TMEV_2DPI.(TMEV2fields{k}).PeakDurs .* 1.04;
    tempDurs = tempDurs(~isnan(tempDurs));
    tempA = TMEV_2DPI.(TMEV2fields{k}).LENGTH_A;
    tempB = TMEV_2DPI.(TMEV2fields{k}).LENGTH_B;
    tempRAWAUC = TMEV_2DPI.(TMEV2fields{k}).RAW_AUC;
    tempRAWAUC = abs(tempRAWAUC(~isnan(tempRAWAUC)));
    tempAUC = TMEV_2DPI.(TMEV2fields{k}).AUC;
    tempAUC = abs(tempAUC(~isnan(tempAUC)));
    tempPeaks = TMEV_2DPI.(TMEV2fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    
    allTMEV2Peaks = [allTMEV2Peaks; tempPeaks];
    all_RAWauc = [all_RAWauc; tempRAWAUC];
    all_auc = [all_auc; tempAUC];
    shiftedLocs = tempLocs - temp_timehitR;
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));
    tempBefore0 = tempRAWAUC(shiftedLocs < 0 );
    tempAfter0 = tempRAWAUC(shiftedLocs >= 0 );
    B0Peaks = tempPeaks(shiftedLocs < 0);
    A0Peaks = tempPeaks(shiftedLocs >= 0);
    before_0_auc = [before_0_auc; tempBefore0];
    after_0_auc = [after_0_auc; tempAfter0];
    before_0_peaks = [before_0_peaks; B0Peaks];
    after_0_peaks = [after_0_peaks; A0Peaks];
    TMEV_2DPI.(TMEV2fields{k}).shiftedLocs = shiftedLocs;
    B0Peaks = tempPeaks(shiftedLocs < 0);
    A0Peaks = tempPeaks(shiftedLocs >= 0);
    B0Durs = tempDurs(shiftedLocs < 0);
    A0Durs = tempDurs(shiftedLocs >= 0);
    before_0_durs = [before_0_durs; B0Durs];
    after_0_durs = [after_0_durs; A0Durs];
    
    allTMEV2Locs = [allTMEV2Locs; shiftedLocs];
    allPeakDurs = [allPeakDurs; tempDurs];
    TMEV2_Alengths = [TMEV2_Alengths; tempA];
    TMEV2_Blengths = [TMEV2_Blengths; tempB];
end

% % Stacked TMEV2 Traces % %
shiftedTime = [-1739:1740]./0.97;
figure()
xlabel('Relative time(s)'); ylabel('Trace Names');
for k = 1:numel(TMEV2fields)
    temp_timehitR = TMEV.(TMEV2fields{k}).timehitR;
    R = round(temp_timehitR*0.97, 1);
    tempLocs = TMEV.(TMEV2fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    if temp_timehitR > frames 
        shiftedLocs = (tempLocs - R)./0.97;
    else
        shiftedLocs = (tempLocs - R)./0.97;
    end
    tempPeaks = TMEV2.(TMEV2fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));   
    tempF = TMEV.(TMEV2fields{k}).trace;
    F0 = mean(tempF(1:50));
    tempF = (tempF - F0) / F0;
    tempTime = (1:1740) - R;
    tempTime = tempTime./0.97;
    hold on
    plot(tempTime, tempF + k);
    hold on
    plot(shiftedLocs, tempF(tempLocs) + k, 'db', 'MarkerFaceColor','r','MarkerSize',6);
end
yticks([1:size(TMEV2fields)]);
yticklabels(TMEV2fields);

binEdges = [-1600:100:1600];
numBins = numel(binEdges);
figure()
histogram(allTMEV2Locs,numBins, 'BinEdges', binEdges)
title('All Shifted TMEV 2 DPI Peak Locations');
xlabel('Relative time (s)'); ylabel('Frequency of events');

figure()
plot(allTMEV2Locs, allPeakDurs, 'dc','MarkerFaceColor','b')
title('TMEV 2 DPI: Peak Durations in relative time'); xlabel('Relative time (s)'); ylabel('Peak Duration (s)');

figure()
plot(allTMEV2Locs,all_auc,'d');
title('TMEV 2 DPI: AUC vs. Relative time: Filtered Signal');
xlabel('relative time (s)'); ylabel('AUC');

figure()
plot(allTMEV2Locs,all_RAWauc,'d');
title('TMEV 2 DPI: AUC vs. Relative time: Raw Signal');
xlabel('relative time (s)'); ylabel('AUC');

figure()
plot(allTMEV2Locs, allTMEV2Peaks, 'dc','MarkerFaceColor','b')
title('TMEV 2 DPI: Peak Magnitudes in relative time'); xlabel('Relative time (s)'); ylabel('Peak Height (dF/F)');

% % Plotting Inter-event Intervals % %
figure()
xlabel('Relative Time (s)'); ylabel('IEI (s) plotted @ first peak'); title('TMEV2: IEI');
for k = 1:length(TMEV2fields)
    tempIEI = TMEV.(TMEV2fields{k}).IEI;
    tempLocs = TMEV.(TMEV2fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    tempLocs = tempLocs(1:end-1); % exclude last peak - plotting IEI on first peak in calculation
    tempLocs = tempLocs./0.97;
    R = TMEV.(TMEV2fields{k}).timehitR;
    tempLocs = tempLocs - R;
    
    hold on
    plot(tempLocs, tempIEI, 'db','MarkerFaceColor','r');
end

% % boxplots of event auc before time zero and after % %
ALL = [before_0_auc; after_0_auc];
groups = [ones(size(before_0_auc)); ones(size(after_0_auc)) * 2];
figure()
boxplot(ALL,groups)
title('TMEV 2 DPI Raw Event AUC before relative time zero and after time zero');
ylabel('AUC (dF/F)');

% % boxplots of event Peaks before time zero and after % %
ALL = [before_0_peaks; after_0_peaks];
figure()
boxplot(ALL,groups,'Whisker',4)
title('TMEV 2 DPI: Raw Event Peaks before relative time zero and after time zero');
ylabel('AUC (dF/F)');

% % boxplots of event durations before time zero and after % %
ALL = [before_0_durs; after_0_durs];
figure()
boxplot(ALL,groups,'Whisker',4)
title('TMEV 2 DPI: Raw Event Durations before relative time zero and after time zero');
ylabel('duration (s)');

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% TMEV  5 DPI %%%%%
%%%%%%%%%%%%%%%%%%%%%%%
TMEV_5DPI = TMEV.DPI_5;
TMEV5fields = fieldnames(TMEV_5DPI);
allTMEV5Locs = []; %allocate space
allTMEV5Peaks = [];
allPeakDurs = [];
TMEV5_Alengths = [];
TMEV5_Blengths = [];
all_RAWauc = [];
all_auc = [];
before_0_auc = [];
after_0_auc = [];
before_0_peaks = [];
after_0_peaks = [];
before_0_durs = [];
after_0_durs = [];
for k  = 1:numel(TMEV5fields)
    % extract
    temp_timehitR = TMEV_5DPI.(TMEV5fields{k}).timehitR;
    tempLocs = TMEV_5DPI.(TMEV5fields{k}).PeakLocs;
    tempDurs = TMEV_5DPI.(TMEV5fields{k}).PeakDurs .* 1.04;
    tempDurs = tempDurs(~isnan(tempDurs));
    tempA = TMEV_5DPI.(TMEV5fields{k}).LENGTH_A;
    tempB = TMEV_5DPI.(TMEV5fields{k}).LENGTH_B;
    tempRAWAUC = TMEV_5DPI.(TMEV5fields{k}).RAW_AUC;
    tempRAWAUC = abs(tempRAWAUC(~isnan(tempRAWAUC)));
    tempAUC = TMEV_5DPI.(TMEV5fields{k}).AUC;
    tempAUC = abs(tempAUC(~isnan(tempAUC)));
    tempPeaks = TMEV_5DPI.(TMEV5fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    
    allTMEV5Peaks = [allTMEV5Peaks; tempPeaks];
    all_RAWauc = [all_RAWauc; tempRAWAUC];
    all_auc = [all_auc; tempAUC];
    shiftedLocs = tempLocs - temp_timehitR;
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));
    tempBefore0 = tempRAWAUC(shiftedLocs < 0 );
    tempAfter0 = tempRAWAUC(shiftedLocs >= 0 );
    B0Peaks = tempPeaks(shiftedLocs < 0);
    A0Peaks = tempPeaks(shiftedLocs >= 0);
    before_0_auc = [before_0_auc; tempBefore0];
    after_0_auc = [after_0_auc; tempAfter0];
    before_0_peaks = [before_0_peaks; B0Peaks];
    after_0_peaks = [after_0_peaks; A0Peaks];
    TMEV_5DPI.(TMEV5fields{k}).shiftedLocs = shiftedLocs;
    B0Peaks = tempPeaks(shiftedLocs < 0);
    A0Peaks = tempPeaks(shiftedLocs >= 0);
    B0Durs = tempDurs(shiftedLocs < 0);
    A0Durs = tempDurs(shiftedLocs >= 0);
    before_0_durs = [before_0_durs; B0Durs];
    after_0_durs = [after_0_durs; A0Durs];
    
    allTMEV5Locs = [allTMEV5Locs; shiftedLocs];
    allPeakDurs = [allPeakDurs; tempDurs];
    TMEV5_Alengths = [TMEV5_Alengths; tempA];
    TMEV5_Blengths = [TMEV5_Blengths; tempB];
end

% % Stacked TMEV5 Traces % %
shiftedTime = [-1739:1740]./0.97;
figure()
xlabel('Relative time(s)'); ylabel('Trace Names');
for k = 1:numel(TMEV5fields)
    temp_timehitR = TMEV.(TMEV5fields{k}).timehitR;
    R = round(temp_timehitR*0.97, 1);
    tempLocs = TMEV.(TMEV5fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    if temp_timehitR > frames 
        shiftedLocs = (tempLocs - R)./0.97;
    else
        shiftedLocs = (tempLocs - R)./0.97;
    end
    tempPeaks = TMEV5.(TMEV5fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));   
    tempF = TMEV.(TMEV5fields{k}).trace;
    F0 = mean(tempF(1:50));
    tempF = (tempF - F0) / F0;
    tempTime = (1:1740) - R;
    tempTime = tempTime./0.97;
    hold on
    plot(tempTime, tempF + k);
    hold on
    plot(shiftedLocs, tempF(tempLocs) + k, 'db', 'MarkerFaceColor','r','MarkerSize',6);
end
yticks([1:size(TMEV5fields)]);
yticklabels(TMEV5fields);

binEdges = [-1600:100:1600];
numBins = numel(binEdges);
figure()
histogram(allTMEV5Locs,numBins, 'BinEdges', binEdges)
title('All Shifted TMEV 5 DPI Peak Locations');
xlabel('Relative time (s)'); ylabel('Frequency of events');

figure()
plot(allTMEV5Locs, allPeakDurs, 'dc','MarkerFaceColor','b')
title('tmev 5 dpi: Peak Durations in relative time'); xlabel('Relative time (s)'); ylabel('Peak Duration (s)');

figure()
plot(allTMEV5Locs,all_auc,'d');
title('tmev 5 dpi: AUC vs. Relative time: Filtered Signal');
xlabel('relative time (s)'); ylabel('AUC');

figure()
plot(allTMEV5Locs,all_RAWauc,'d');
title('tmev 5 dpi: AUC vs. Relative time: Raw Signal');
xlabel('relative time (s)'); ylabel('AUC');

figure()
plot(allTMEV5Locs, allTMEV5Peaks, 'dc','MarkerFaceColor','b')
title('TMEV 5 DPI: Peak Magnitudes in relative time'); xlabel('Relative time (s)'); ylabel('Peak Height (dF/F)');

% % Plotting Inter-event Intervals % %
figure()
xlabel('Relative Time (s)'); ylabel('IEI (s) plotted @ first peak'); title('TMEV5: IEI');
for k = 1:length(TMEV5fields)
    tempIEI = TMEV.(TMEV5fields{k}).IEI;
    tempLocs = TMEV.(TMEV5fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    tempLocs = tempLocs(1:end-1); % exclude last peak - plotting IEI on first peak in calculation
    tempLocs = tempLocs./0.97;
    R = TMEV.(TMEV5fields{k}).timehitR;
    tempLocs = tempLocs - R;
    
    hold on
    plot(tempLocs, tempIEI, 'db','MarkerFaceColor','r');
end

% % boxplots of event auc before time zero and after % %
ALL = [before_0_auc; after_0_auc];
groups = [ones(size(before_0_auc)); ones(size(after_0_auc)) * 2];
figure()
boxplot(ALL,groups)
title('TMEV 5 DPI Raw Event AUC before relative time zero and after time zero');
ylabel('AUC (dF/F)');

% % boxplots of event Peaks before time zero and after % %
ALL = [before_0_peaks; after_0_peaks];
figure()
boxplot(ALL,groups,'Whisker',4)
title('TMEV 5 DPI: Raw Event Peaks before relative time zero and after time zero');
ylabel('AUC (dF/F)');

% % boxplots of event durations before time zero and after % %
ALL = [before_0_durs; after_0_durs];
figure()
boxplot(ALL,groups,'Whisker',4)
title('TMEV 5 DPI: Raw Event Durations before relative time zero and after time zero');
ylabel('duration (s)');

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% TMEV  14 DPI %%%%
%%%%%%%%%%%%%%%%%%%%%%%

TMEV_14DPI = TMEV.DPI_14;
TMEV14fields = fieldnames(TMEV_14DPI);
allTMEV14Locs = []; %allocate space
allTMEV14Peaks = [];
allPeakDurs = [];
TMEV14_Alengths = [];
TMEV14_Blengths = [];
all_RAWauc = [];
all_auc = [];
before_0_auc = [];
after_0_auc = [];
before_0_peaks = [];
after_0_peaks = [];
before_0_durs = [];
after_0_durs = [];
for k  = 1:numel(TMEV14fields)
    % extract
    temp_timehitR = TMEV_14DPI.(TMEV14fields{k}).timehitR;
    tempLocs = TMEV_14DPI.(TMEV14fields{k}).PeakLocs;
    tempDurs = TMEV_14DPI.(TMEV14fields{k}).PeakDurs .* 1.04;
    tempDurs = tempDurs(~isnan(tempDurs));
    tempA = TMEV_14DPI.(TMEV14fields{k}).LENGTH_A;
    tempB = TMEV_14DPI.(TMEV14fields{k}).LENGTH_B;
    tempRAWAUC = TMEV_14DPI.(TMEV14fields{k}).RAW_AUC;
    tempRAWAUC = abs(tempRAWAUC(~isnan(tempRAWAUC)));
    tempAUC = TMEV_14DPI.(TMEV14fields{k}).AUC;
    tempAUC = abs(tempAUC(~isnan(tempAUC)));
    tempPeaks = TMEV_14DPI.(TMEV14fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    
    allTMEV14Peaks = [allTMEV14Peaks; tempPeaks];
    all_RAWauc = [all_RAWauc; tempRAWAUC];
    all_auc = [all_auc; tempAUC];
    shiftedLocs = tempLocs - temp_timehitR;
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));
    tempBefore0 = tempRAWAUC(shiftedLocs < 0 );
    tempAfter0 = tempRAWAUC(shiftedLocs >= 0 );
    B0Peaks = tempPeaks(shiftedLocs < 0);
    A0Peaks = tempPeaks(shiftedLocs >= 0);
    before_0_auc = [before_0_auc; tempBefore0];
    after_0_auc = [after_0_auc; tempAfter0];
    before_0_peaks = [before_0_peaks; B0Peaks];
    after_0_peaks = [after_0_peaks; A0Peaks];
    TMEV_14DPI.(TMEV14fields{k}).shiftedLocs = shiftedLocs;
    B0Peaks = tempPeaks(shiftedLocs < 0);
    A0Peaks = tempPeaks(shiftedLocs >= 0);
    B0Durs = tempDurs(shiftedLocs < 0);
    A0Durs = tempDurs(shiftedLocs >= 0);
    before_0_durs = [before_0_durs; B0Durs];
    after_0_durs = [after_0_durs; A0Durs];
    
    allTMEV14Locs = [allTMEV14Locs; shiftedLocs];
    allPeakDurs = [allPeakDurs; tempDurs];
    TMEV14_Alengths = [TMEV14_Alengths; tempA];
    TMEV14_Blengths = [TMEV14_Blengths; tempB];
end

% % Stacked TMEV14 Traces % %
shiftedTime = [-1739:1740]./0.97;
figure()
xlabel('Relative time(s)'); ylabel('Trace Names');
for k = 1:numel(TMEV14fields)
    temp_timehitR = TMEV.(TMEV14fields{k}).timehitR;
    R = round(temp_timehitR*0.97, 1);
    tempLocs = TMEV.(TMEV14fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    if temp_timehitR > frames 
        shiftedLocs = (tempLocs - R)./0.97;
    else
        shiftedLocs = (tempLocs - R)./0.97;
    end
    tempPeaks = TMEV14.(TMEV14fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));   
    tempF = TMEV.(TMEV14fields{k}).trace;
    F0 = mean(tempF(1:50));
    tempF = (tempF - F0) / F0;
    tempTime = (1:1740) - R;
    tempTime = tempTime./0.97;
    hold on
    plot(tempTime, tempF + k);
    hold on
    plot(shiftedLocs, tempF(tempLocs) + k, 'db', 'MarkerFaceColor','r','MarkerSize',6);
end
yticks([1:size(TMEV14fields)]);
yticklabels(TMEV14fields);

binEdges = [-1600:100:1600];
numBins = numel(binEdges);
figure()
histogram(allTMEV14Locs,numBins, 'BinEdges', binEdges)
% plot(allPBSLocs,'d')
title('All Shifted TMEV 14 DPI Peak Locations');
xlabel('Relative time (s)'); ylabel('Frequency of events');

figure()
plot(allTMEV14Locs, allPeakDurs, 'dc','MarkerFaceColor','b')
title('tmev 14 dpi: Peak Durations in relative time'); xlabel('Relative time (s)'); ylabel('Peak Duration (s)');

figure()
plot(allTMEV14Locs,all_auc,'d');
title('tmev 14 dpi: AUC vs. Relative time: Filtered Signal');
xlabel('relative time (s)'); ylabel('AUC');

figure()
plot(allTMEV14Locs,all_RAWauc,'d');
title('tmev 14 dpi: AUC vs. Relative time: Raw Signal');
xlabel('relative time (s)'); ylabel('AUC');

figure()
plot(allTMEV14Locs, allTMEV14Peaks, 'dc','MarkerFaceColor','b')
title('TMEV 14 DPI: Peak Magnitudes in relative time'); xlabel('Relative time (s)'); ylabel('Peak Height (dF/F)');

% % Plotting Inter-event Intervals % %
figure()
xlabel('Relative Time (s)'); ylabel('IEI (s) plotted @ first peak'); title('TMEV14: IEI');
for k = 1:length(TMEV14fields)
    tempIEI = TMEV.(TMEV14fields{k}).IEI;
    tempLocs = TMEV.(TMEV14fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    tempLocs = tempLocs(1:end-1); % exclude last peak - plotting IEI on first peak in calculation
    tempLocs = tempLocs./0.97;
    R = TMEV.(TMEV14fields{k}).timehitR;
    tempLocs = tempLocs - R;
    
    hold on
    plot(tempLocs, tempIEI, 'db','MarkerFaceColor','r');
end

% % boxplots of event auc before time zero and after % %
ALL = [before_0_auc; after_0_auc];
groups = [ones(size(before_0_auc)); ones(size(after_0_auc)) * 2];
figure()
boxplot(ALL,groups)
title('TMEV 14 DPI Raw Event AUC before relative time zero and after time zero');
ylabel('AUC (dF/F)');

% % boxplots of event Peaks before time zero and after % %
ALL = [before_0_peaks; after_0_peaks];
figure()
boxplot(ALL,groups,'Whisker',4)
title('TMEV 14 DPI: Raw Event Peaks before relative time zero and after time zero');
ylabel('Peak (dF/F)');

% % boxplots of event durations before time zero and after % %
ALL = [before_0_durs; after_0_durs];
figure()
boxplot(ALL,groups,'Whisker',4)
title('TMEV 14 DPI: Raw Event Durations before relative time zero and after time zero');
ylabel('duration (s)');
%%
AL = [];
A_RMS = [];
B_RMS = [];
BL = [];
PBS_R = [];
%     TMEV_R = [];
allRMS{1,1} = 'Trace + ROI';
allRMS{1,2} = 'TMEV/PBS';
allRMS{1,3} = 'DPI';
allRMS{1,4} = 'A RMS^2';
allRMS{1,5} = 'B RMS^2';
filtered_traces = cell(size(process_data,1),3);
for ii = 1:size(process_data,1)
    ROInum = process_data{ii,2};
    ROInum = num2str(ROInum);
    ROI = strcat('  ROI: ', ROInum);
    traceName = process_data{ii,1};
    traceName = strcat(traceName,ROI);
    DPI = process_data{ii,5};
    if strcmp(process_data{ii,4},'TMEV')
        if process_data{ii,5} == 2 %%% TMEV DPI condition
            shifted_trace = double(process_data{ii,7});
            shifted_trace(isnan(shifted_trace)) = 0;
            
            timehitR = double(process_data{ii,6});
            timehitR_marker = round(timehitR);
            filtered_traces{ii,3} = timehitR_marker;
            PBS_R = [PBS_R; timehitR];
            shifted_trace(isnan(shifted_trace)) = 0;
            frames = length(shifted_trace)/2;
            traceA = shifted_trace(1:frames);
            traceA(traceA==0) = [];
            traceB = shifted_trace(frames+1:end);
            traceB(traceB==0) = [];
            
            %                 [LPF, ~] = noiseRemover(100, 0.055, 1.04);
            %                 if length(traceA) < 300
            %                     paddingL = 300 - length(traceA)+1;
            %                     padding = zeros(1,paddingL);
            %                     traceA = [padding, traceA];
            %                     traceA = filtfilt(LPF, traceA);
            %                     traceA(1:paddingL) = [];
            %
            %                 else
            %                     traceA = filtfilt(LPF, traceA);
            %                 end
            %
            %                 % % filtered trace A step % %
            %
            %
            %                 if length(traceB) < 300
            %                     paddingL = 300 - length(traceB)+1;
            %                     padding = zeros(1,paddingL);
            %                     traceB = [padding, traceB];
            %                     traceB = filtfilt(LPF, traceB);
            %                     traceB(1:paddingL) = [];
            %                 else
            %                     traceB = filtfilt(LPF, traceB);
            %                 end
            traceA(1:61) = [];
            % % envelopes % %
            %                 [Aupper, ~] = envelope(traceA, 35,'peak');
            %                 [Bupper, ~] = envelope(traceB, 35,'peak');
            
            
            %%%% CALCULATE RMS/Power and add to each trace a/b vectors
            filtered_traces{ii,1} = [traceA, traceB];
            filtered_traces{ii,2} = traceName;
            % % Find arc length of two halves of trace % %
            lengthA = signalLength(traceA,1.04)/length(traceA);
            lengthB = abs(signalLength(traceB,1.04)/length(traceB));
            AL = [AL;lengthA];
            BL = [BL;lengthB];
            
            % % Find RMS of A and B portions % %
            rmsA = rms(traceA);
            rmsB = rms(traceB);
            A_RMS = [A_RMS; rmsA];
            B_RMS = [B_RMS; rmsB];
            
            allRMS{ii+1,1} = traceName;
            allRMS{ii+1,2} = 'TMEV';
            allRMS{ii+1,3} = DPI;
            allRMS{ii+1,4} = rmsA;
            allRMS{ii+1,5} = rmsB;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% SUBPLOT for split trace and complete trace %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % filtered trace A step % %
            figure()
            subplot(2,3,1);
            plot(traceA)
            title(strcat('Trace A: ',traceName)); xlabel('frames'); ylabel('df/f');
            
            % % filtered trace B step % %
            subplot(2,3,2);
            plot(traceB)
            title(strcat('Trace B: ',traceName)); xlabel('frames'); ylabel('df/f');
            
            % % Original trace with marker for timehitR % %
            original_trace = shifted_trace;
            original_trace(original_trace == 0) = [];
            subplot(2,3,3);
            plot(original_trace)
            if timehitR_marker < length(original_trace)
                hold on
                plot(timehitR_marker,original_trace(timehitR_marker),'d','MarkerFaceColor','r');
                title(traceName)
            else
                %%% nothing %%%
            end
            % % Plot Power for trace A/trace B
            subplot(2,3,4)
            plot(1,rmsA,'d','MarkerFaceColor','b');
            hold on
            plot(2, rmsB,'d','MarkerFaceColor','b');
            xticks([0, 1, 2, 3])
            xlim([0 3])
            title('Trace A vs Trace B Power Analysis'); xlabel('Trace A / Trace B'); ylabel('RMS^2');
            
            subplot(2,3,5)
            plot(1,lengthA,'d','MarkerFaceColor','b');
            hold on
            plot(2,lengthB,'d','MarkerFaceColor','b');
            xticks([0, 1, 2, 3])
            xlim([0 3])
            title('Trace A vs Trace B Lengths'); xlabel('Trace A / Trace B'); ylabel('df/f');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Subplot for all traces(split/orig)COMPLETE %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %             end %%%% uncomment for TMEV DPI condition
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Plotting Results %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Paired RMS % %
figure()
xticks([0, 1, 2, 3])
xlim([0 3])
xlabel('Before time 0 / After time 0'); ylabel('RMS'); title('Soma PBS RMS');
for k = 1:length(A_RMS)
    hold on
    if A_RMS(k) < B_RMS(k)
        plot([1,2], [log10(A_RMS(k)), log10(B_RMS(k))], 'r-d');
    else
        plot([1,2], [A_RMS(k), B_RMS(k)], 'b-d');
        
    end
end

groupA = ones(size(AL));
groupB = ones(size(BL))+1;
% % PLOTTING RMS^2 (POWER)
figure()
plot(groupA, A_RMS, 'dk');
hold on
plot(groupB, B_RMS, 'dk');
legend('Length A', 'Length B');
title('PBS: RMS LEVELS traces before and after timehitR'); xlabel('Group A / Group B'); ylabel('RMS^2');
xticks([0, 1, 2, 3])
xlim([0 3])

% % PLOTTING TMEV POWER
figure()
for k = 1:length(A_RMS)
hold on
plot([1,log10(A_RMS(k))], [2, log10(B_RMS(k))], 'dk');
% hold on
% plot(groupB, B_RMS, 'dk');
% legend('Length A', 'Length B');
title('TMEV 2 DPI: RMS levels of filtered traces before and after timehitR'); xlabel('Group A / Group B'); ylabel('RMS');
xticks([0, 1, 2, 3])
xlim([0 3])
% A_AVERAGE = mean(A_RMS);
end

% % PLOTTING PBS VS TMEV TIME HIT R
TMEV_GROUP = ones(size(TMEV_R));
PBS_GROUP = ones(size(PBS_R))+1;
figure()
plot(TMEV_GROUP,TMEV_R,'dk');
hold on
plot(PBS_GROUP,PBS_R,'dk');
title('TMEV vs. PBS timehitsR''s'); xlabel('TMEV / PBS'); ylabel('timehitR (s)');
xticks([0, 1, 2, 3])
xlim([0 3])

% % PLOTTING TRACE A + TRACE B W/ timehitR MARKER
figure()
n = 0; % initialize counter
for k = 1:size(filtered_traces,1)
    if ~isempty(filtered_traces{k,1})
    trace = filtered_traces{k,1};
    split = filtered_traces{k,3};
    hold on
    plot(trace + n);
    if split < length(trace)
    hold on
    plot(split, trace(split) + n, 'd','MarkerFaceColor','r');
    end
    n = n + 1; % step counter
    end
end
trace_labels = filtered_traces(:,2);
trace_labels = trace_labels(~cellfun('isempty',trace_labels));
yticks([1:size(trace_labels,1)]);
yticklabels(trace_labels);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Labels = filtered_traces(:,2);
figure()
for x = 1:size(filtered_traces,1)
    if ~isempty(filtered_traces{x,1})
    hold on
    plot(filtered_traces{x,1} + x)
    end
end
% yticklabels(Labels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% HYPOTHESIS TESTING %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % parametric 2 sample t-test
format long e

[h2,p2,ci2,stats2] = ttest2(A_RMS,B_RMS,'Alpha',0.05,'Vartype','unequal')

[H, P, CI, STATS] = vartest2(A_RMS,B_RMS)

% % MANN-WHITNEY TEST FOR RMS^2 (POWER) % %
MWW_STATS = mwwtest(A_RMS', B_RMS')
% % END MWW TEST % %

%% SOMA %%
AL = [];
A_RMS = [];
B_RMS = [];
BL = [];
PBS_R = [];
%     TMEV_R = [];
allRMS{1,1} = 'Trace + ROI';
allRMS{1,2} = 'TMEV/PBS';
allRMS{1,3} = 'DPI';
allRMS{1,4} = 'A RMS^2';
allRMS{1,5} = 'B RMS^2';
filtered_traces = cell(size(soma_data,1),3);
for ii = 1:size(soma_data,1)
    ROInum = soma_data{ii,2};
    ROInum = num2str(ROInum);
    ROI = strcat('  ROI: ', ROInum);
    traceName = soma_data{ii,1};
    traceName = strcat(traceName,ROI);
    DPI = soma_data{ii,5};
    if strcmp(soma_data{ii,4},'PBS')
        %         if soma_data{ii,5} == 2 %%% TMEV DPI condition
        shifted_trace = double(soma_data{ii,7});
        shifted_trace(isnan(shifted_trace)) = 0;
        
        timehitR = double(soma_data{ii,6});
        timehitR_marker = round(timehitR);
        filtered_traces{ii,3} = timehitR_marker;
        PBS_R = [PBS_R; timehitR];
        shifted_trace(isnan(shifted_trace)) = 0;
        frames = length(shifted_trace)/2;
        traceA = shifted_trace(1:frames);
        traceA(traceA==0) = [];
        traceB = shifted_trace(frames+1:end);
        traceB(traceB==0) = [];
        traceA(1:61) = [];
        %%%% CALCULATE RMS/Power and add to each trace a/b vectors
        filtered_traces{ii,1} = [traceA, traceB];
        filtered_traces{ii,2} = traceName;
        % % Find arc length of two halves of trace % %
        lengthA = signalLength(traceA,1.04)/length(traceA);
        lengthB = abs(signalLength(traceB,1.04)/length(traceB));
        AL = [AL;lengthA];
        BL = [BL;lengthB];
        
        % % Find RMS of A and B portions % %
        rmsA = rms(traceA);
        rmsB = rms(traceB);
        A_RMS = [A_RMS; rmsA];
        B_RMS = [B_RMS; rmsB];
        
        allRMS{ii+1,1} = traceName;
        allRMS{ii+1,2} = 'TMEV';
        allRMS{ii+1,3} = DPI;
        allRMS{ii+1,4} = rmsA;
        allRMS{ii+1,5} = rmsB;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SUBPLOT for split trace and complete trace %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % filtered trace A step % %
        figure()
        subplot(2,3,1);
        plot(traceA)
        title(strcat('Trace A: ',traceName)); xlabel('frames'); ylabel('df/f');
        
        % % filtered trace B step % %
        subplot(2,3,2);
        plot(traceB)
        title(strcat('Trace B: ',traceName)); xlabel('frames'); ylabel('df/f');
        
        % % Original trace with marker for timehitR % %
        original_trace = shifted_trace;
        original_trace(original_trace == 0) = [];
        subplot(2,3,3);
        plot(original_trace)
        if timehitR_marker < length(original_trace)
            hold on
            plot(timehitR_marker,original_trace(timehitR_marker),'d','MarkerFaceColor','r');
            title(traceName)
        else
            %%% nothing %%%
        end
        % % Plot Power for trace A/trace B
        subplot(2,3,4)
        plot(1,rmsA,'d','MarkerFaceColor','b');
        hold on
        plot(2, rmsB,'d','MarkerFaceColor','b');
        xticks([0, 1, 2, 3])
        xlim([0 3])
        title('Trace A vs Trace B Power Analysis'); xlabel('Trace A / Trace B'); ylabel('RMS^2');
        
        subplot(2,3,5)
        plot(1,lengthA,'d','MarkerFaceColor','b');
        hold on
        plot(2,lengthB,'d','MarkerFaceColor','b');
        xticks([0, 1, 2, 3])
        xlim([0 3])
        title('Trace A vs Trace B Lengths'); xlabel('Trace A / Trace B'); ylabel('df/f');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Subplot for all traces(split/orig)COMPLETE %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %             end %%%% uncomment for TMEV DPI condition
    end
    
end

% % Paired RMS % %
figure()
xticks([0, 1, 2, 3])
xlim([0 3])
xlabel('Before time 0 / After time 0'); ylabel('RMS'); title('Soma PBS RMS');
for k = 1:length(A_RMS)
    hold on
    if A_RMS(k) < B_RMS(k)
        plot([1,2], [A_RMS(k), B_RMS(k)], 'r-d');
    else
        plot([1,2], [A_RMS(k), B_RMS(k)], 'b-d');
        
    end
end


end
