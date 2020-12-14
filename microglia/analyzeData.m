%%% This is a function to take in a sortedData structure and analyze the
%%% data within the groups.

% % % Remember we want to look at each subsequent metric in relative time

% dont have timehitR in struct yet but this is how we will do it:
%sortedData.Group.DPI_X.Video.RelativeLocs =
%sortedData.Group.DPI_X.Video.PeakLocs -
%sortedData.Group.DPI_X.Video.timehitR;
burnFrame = 90;
DATA = process_data;
% DATA = soma_data;
isburn = @(x) logical(x>90);
% bring groups up a layer to keep things less confusing
PBS = sortedData.PBS;
TMEV = sortedData.TMEV;
burnFlag = false;
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
allPBSAUC = {};
allPBSDurs = {};
allPBSMags = {};
allPBSROI = {};
allPBSpath = {};
allPBSDPI = {};
allPBSNormalLocs = {};
allPBSShiftedLocs = {};
allPBSEventGroup = {};
allPBSMice = {};
allPBSIEI = {};

Group1_AUC = [];
Group2_AUC = [];
Group3_AUC = [];

Group1_PeakMags = [];
Group2_PeakMags = [];
Group3_PeakMags = [];

Group1_PeakDurs = [];
Group2_PeakDurs = [];
Group3_PeakDurs = [];

Group1_Locs = [];
Group2_Locs = [];
Group3_Locs = [];

% % count number of mice for PBS % %
PBS_logicals = contains(DATA(:,4),'PBS'); 
PBS_extracted = DATA(PBS_logicals,:);
mice = PBS_extracted(:,7);
uniquePBS = unique(cellfun(@num2str,mice,'uni',0));
PBS_mice_n = numel(uniquePBS);
PBS_mice_str = ['(N = ' num2str(PBS_mice_n) 'mice)'];
for k = 1:numel(PBSfields)

    % extract
    tempIEI = PBS.(PBSfields{k}).IEI;
    allIEI = [allIEI;tempIEI];
     
    temp_timehitR = PBS.(PBSfields{k}).timehitR;
    tempLocs = PBS.(PBSfields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    
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
        allPBSAUC = [allPBSAUC; tempRAWAUC(i)];
        allPBSDurs = [allPBSDurs; tempDurs(i)];
        allPBSMags = [allPBSMags; tempPeaks(i)];
        allPBSROI = [allPBSROI; PBS.(PBSfields{k}).ROI];
        allPBSpath = [allPBSpath; PBS.(PBSfields{k}).Pathology];
        allPBSDPI = [allPBSDPI; PBS.(PBSfields{k}).DPI];
        allPBSNormalLocs = [allPBSNormalLocs; tempLocs(i)];
        allPBSShiftedLocs = [allPBSShiftedLocs; tempLocs(i)-round(temp_timehitR*0.97)];
        if i <= length(tempIEI)
            allPBSIEI = [allPBSIEI; tempIEI(i)];
        else
            allPBSIEI = [allPBSIEI; 'Skipped b/c last peak'];
        end
        if tempLocs(i) < burnFrame
            allPBSEventGroup = [allPBSEventGroup; 1];
        elseif tempLocs(i) >= round(temp_timehitR*0.97)
            allPBSEventGroup = [allPBSEventGroup; 3];
        else
            allPBSEventGroup = [allPBSEventGroup; 2];
        end
%         for j = 1:length(mice)
%             if strcmp(PBS_extracted{j},PBSfields{k})
%                 allPBSMice = [allPBSMice; mice(j)];
%             end
%         end

    end
 
        % % % Conditional if excluding burns % % %
    if burnFlag    
        noBurnLocs = tempLocs(tempLocs > burnFrame); % exclude burn
        tempPeaks = tempPeaks(tempLocs > burnFrame);
        tempRAWAUC = tempRAWAUC(tempLocs > burnFrame);
        tempDurs = tempDurs(tempLocs > burnFrame);
    else
        noBurnLocs = tempLocs;
        tempPeaks = tempPeaks;
        tempRAWAUC = tempRAWAUC;
        tempDurs = tempDurs;
        %%%%% Break up into 3 groups (burn, b/w burn and timehitR, post
        %%%%% timehitR)
        temp_Group1_AUC = tempRAWAUC(tempLocs <= burnFrame);
        temp_Group2_AUC = tempRAWAUC((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_AUC = tempRAWAUC(tempLocs >=  round(temp_timehitR*0.97));
        Group1_AUC = [Group1_AUC; temp_Group1_AUC];
        Group2_AUC = [Group2_AUC; temp_Group2_AUC];
        Group3_AUC = [Group3_AUC; temp_Group3_AUC];
        
        temp_Group1_PeakMags = tempPeaks(tempLocs <= burnFrame);
        temp_Group2_PeakMags = tempPeaks((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_PeakMags = tempPeaks(tempLocs >=  round(temp_timehitR*0.97));
        Group1_PeakMags = [Group1_PeakMags; temp_Group1_PeakMags];
        Group2_PeakMags = [Group2_PeakMags; temp_Group2_PeakMags];
        Group3_PeakMags = [Group3_PeakMags; temp_Group3_PeakMags];
        
        temp_Group1_PeakDurs = tempDurs(tempLocs <= burnFrame);
        temp_Group2_PeakDurs = tempDurs((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_PeakDurs = tempDurs(tempLocs >=  round(temp_timehitR*0.97));
        Group1_PeakDurs = [Group1_PeakDurs; temp_Group1_PeakDurs];
        Group2_PeakDurs = [Group2_PeakDurs; temp_Group2_PeakDurs];
        Group3_PeakDurs = [Group3_PeakDurs; temp_Group3_PeakDurs];
        
        temp_Group1_Locs = tempLocs(tempLocs <= burnFrame);
        temp_Group2_Locs = tempLocs((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_Locs = tempLocs(tempLocs >=  round(temp_timehitR*0.97));
        Group1_Locs = [Group1_Locs; temp_Group1_Locs];
        Group2_Locs = [Group2_Locs; temp_Group2_Locs];
        Group3_Locs = [Group3_Locs; temp_Group3_Locs];
    end
    
    allPBSPeaks = [allPBSPeaks; tempPeaks];
    all_RAWauc = [all_RAWauc; tempRAWAUC];
    all_auc = [all_auc; tempAUC];
    shiftedLocs = noBurnLocs - temp_timehitR;
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
for j = 1:length(mice)
    for k = 1:length(allPBSTraces)
        if strcmp(PBS_extracted{j,1},allPBSTraces{k})
            allPBSMice = [allPBSMice; mice(j)];
        end
    end
end
        
PBS_Data = cell(length(allPBSTraces)+1, 12);
PBS_Data(1,:) = {'Trace Name','ROI #','DPI','TMEV/PBS','Duration','Peak Amp','AUC','Peak Locations','Shifted Locations','Event Group','MouseID','IEI'};
PBS_Data(2:end,:) = [allPBSTraces,allPBSROI, allPBSDPI, allPBSpath, allPBSDurs,allPBSMags,allPBSAUC,allPBSNormalLocs,allPBSShiftedLocs,allPBSEventGroup,allPBSMice,allPBSIEI];
save('Process_PBS_Data.mat','PBS_Data');
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
savefig('PBS_stacked_traces_NoBurns.fig');

% % Shifted scatter of JUST peaks % %
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
    K = ones(size(shiftedLocs))*k;
    hold on
    plot(shiftedLocs, K, '.k', 'MarkerFaceColor','k','MarkerSize',6);
end
yticks([1:size(PBSfields)]);
yticklabels(PBSfields);

% % Normal traces/unshifted time % %
figure()
xlabel('time(s)'); ylabel('Trace Names');
for k = 1:numel(PBSfields)
    temp_timehitR = PBS.(PBSfields{k}).timehitR;
    R = round(temp_timehitR*0.97, 1);
    tempLocs = PBS.(PBSfields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    tempPeaks = PBS.(PBSfields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks)); 
    tempF = PBS.(PBSfields{k}).trace;
    F0 = mean(tempF(1:50));
    tempF = (tempF - F0) / F0;
    hold on
    plot(tempF + k);
    hold on
    plot(tempLocs, tempF(tempLocs) + k, 'db', 'MarkerFaceColor','r','MarkerSize',6);
end
yticks([1:size(PBSfields)]);
yticklabels(PBSfields);

%%%% PBS Plotting Zones %%%%
% % Zones - EVENT AUC % %
figure()
plot(Group1_Locs, Group1_AUC,'kd','MarkerFaceColor','g');
hold on
plot(Group2_Locs, Group2_AUC,'kd','MarkerFaceColor','b');
hold on
plot(Group3_Locs, Group3_AUC,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('AUC (dF/F)');
title('PBS: Zones 1-3 AUC');
legend('Zone 1','Zone 2','Zone 3');
savefig('PBS_Zones1-3_AUC.fig');

% % Zones - EVENT Duration % %
figure()
plot(Group1_Locs, Group1_PeakDurs,'kd','MarkerFaceColor','g');
hold on
plot(Group2_Locs, Group2_PeakDurs,'kd','MarkerFaceColor','b');
hold on
plot(Group3_Locs, Group3_PeakDurs,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('Peak durations (s)');
title('PBS: Zones 1-3 Peak Duration');
legend('Zone 1','Zone 2','Zone 3');
savefig('PBS_Zones1-3_EventDurations.fig');

% % Zones - EVENT Peak Magnitude % %
figure()
plot(Group1_Locs, Group1_PeakMags,'kd','MarkerFaceColor','g');
hold on
plot(Group2_Locs, Group2_PeakMags,'kd','MarkerFaceColor','b');
hold on
plot(Group3_Locs, Group3_PeakMags,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('Peak Amplitude (dF/F)');
title('PBS: Zones 1-3 Peak Amplitude');
legend('Zone 1','Zone 2','Zone 3');
savefig('PBS_Zones1-3_PeakAmp.fig');

% % Zones - Duration vs. Amplitude % %
figure()
plot(Group1_PeakMags, Group1_PeakDurs,'kd','MarkerFaceColor','g');
hold on
plot(Group2_PeakMags, Group2_PeakDurs,'kd','MarkerFaceColor','b');
hold on
plot(Group3_PeakMags, Group3_PeakDurs,'kd','MarkerFaceColor','r');
xlabel('Peak Amplitude (dF/F)'); ylabel('Peak Duration (s)');
title('PBS: Zones 1-3 Peak Duration vs Amplitude');
legend('Zone 1','Zone 2','Zone 3');
savefig('PBS_Zones1-3_DurationVsAmp.fig');

% % Zones - sqrt(AUC) vs. Amplitude % %
figure()
plot(Group1_PeakMags, sqrt(Group1_AUC),'kd','MarkerFaceColor','g');
hold on
plot(Group2_PeakMags, sqrt(Group2_AUC),'kd','MarkerFaceColor','b');
hold on
plot(Group3_PeakMags, sqrt(Group3_AUC),'kd','MarkerFaceColor','r');
xlabel('Peak Amplitude (dF/F)'); ylabel('sqrt(AUC) (sqrt(dF/F))');
title('PBS: Zones 1-3 Event sqrt(AUC) vs Amplitude');
legend('Zone 1','Zone 2','Zone 3');
savefig('PBS_Zones1-3_sqrtAUCVsAmp.fig');
%%%% End Zone Plotting %%%%

% % Peak Location Histogram % %
figure()
binEdges = [-1600:100:1600];
numBins = numel(binEdges);
h = histogram(allPBSLocs,numBins, 'BinEdges', binEdges);
title('All Shifted PBS Peak Locations');
xlabel('Relative time (s)'); ylabel('Frequency of events');
savefig('PBS_PeakLocation_Histogram_NoBurns.fig');
histHeights = h.Values./PBS_mice_n;
figure()
binEdges = [-1500:100:1600];
bg = bar(binEdges, histHeights);
title(strcat('PBS: Event frequency normalized to # of animals', PBS_mice_str)); ylabel('Events / Animal'); xlabel('Relative time binned (s)');
ylim([0 5]);
set(bg ,'FaceColor','b');

% % Peak Duration Scatter % %
figure()
plot(allPBSLocs, allPeakDurs, 'dc','MarkerFaceColor','b')
title('PBS: Peak Durations in relative time'); xlabel('Relative time (s)'); ylabel('Peak Duration (s)');
ylim([0,300]);
savefig('PBS_PeakDurations_NoBurns.fig');

% % Peak Magnitude Scatter % %
figure()
plot(allPBSLocs, allPBSPeaks, 'dc','MarkerFaceColor','b')
title('PBS: Peak Magnitudes in relative time'); xlabel('Relative time (s)'); ylabel('Peak Height (dF/F)');
ylim([-0.2, 4]);
savefig('PBS_PeakMags_scatter_NoBurns.fig');

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

% % raw event auc scatter % %
fig = figure('DeleteFcn','doc datacursormode');
hold on
plot(allPBSLocs, all_RAWauc, 'dr' ,'MarkerFaceColor','b')
title('PBS: raw signal auc'); xlabel('relative time(s)');
ylim([0, 400]);
savefig('PBS_PeakAUC_scatter_NoBurns.fig');

% % Peak Magnitude vs Peak Duration % %
figure()
plot(allPeakDurs, allPBSPeaks,'bd','MarkerFaceColor','r');
xlabel('Peak durations (s)'); ylabel('Peak magnitudes (dF/F)'); title('PBS: Peak vs duration');
ylim([-0.2, 4]);
[c,p1] = polyfit(allPeakDurs,allPBSPeaks,1);
y_est = polyval(c,allPeakDurs);
hold on
plot(allPeakDurs, y_est, 'r--','LineWidth',2);
savefig('PBS_PeakMags_scatter_NoBurns.fig');

% % Amplitude vs. sqrt(AUC) % %
figure()
plot(sqrt(all_RAWauc), allPBSPeaks, 'd');
xlabel('sqrt(AUC)'); ylabel('Peak Amplitude'); title('PBS: amp vs auc');
[c2,p2] = polyfit(sqrt(all_RAWauc),allPBSPeaks,1);
y_est2 = polyval(c2,sqrt(all_RAWauc));
hold on
plot(sqrt(all_RAWauc), y_est2,'--r');

% % Plotting both linear regressions % %
figure()
plot(sqrt(all_RAWauc), y_est2,'--r','LineWidth',2);
hold on
plot(allPeakDurs, y_est, 'b--','LineWidth',2);
grid on;
xlabel('duration and auc'); ylabel('Peak amp'); title('PBS: both linear regressions');
legend('Amp vs AUC','Amp vs Duration');

% % Plotting Inter-event Intervals % %
figure()
xlabel('Relative Time (s)'); ylabel('IEI (s) plotted @ first peak'); title('PBS: IEI');
for k = 1:length(PBSfields)

    tempLocs = PBS.(PBSfields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    IEILocs = tempLocs(1:end);
    tempIEI = PBS.(PBSfields{k}).IEI;
    if burnFlag
        tempLocs = tempLocs(tempLocs > burnFrame);
        tempIEI = tempIEI(IEILocs(1:end-1) > burnFrame);
    else
        tempLocs = tempLocs;
        tempIEI = tempIEI;
    end
    tempLocs = tempLocs(1:end-1); % exclude last peak - plotting IEI on first peak in calculation
    tempLocs = tempLocs./0.97;
    R = PBS.(PBSfields{k}).timehitR;
    tempLocs = tempLocs - R;  
    hold on
    plot(tempLocs, tempIEI, 'db','MarkerFaceColor','r');
end
ylim([0,1500]);
savefig('PBS_IEI_scatter_NoBurns.fig');
    
    

% % boxplots of event auc before time zero and after % %
ALL = [Group1_AUC; Group2_AUC; Group3_AUC];
groups = [ones(size(Group1_AUC)); ones(size(Group2_AUC)) * 2; ones(size(Group3_AUC)) * 3];
figure()
boxplot(ALL,groups,'Whisker',4)
title('PBS Event AUC (Burns//Between burn and timehitR//Post timehitR)');
ylabel('AUC (dF/F)');
ylim([0, 400]);
savefig('PBS_boxplot_AUC_3groups.fig');

% % boxplots of event Peaks before time zero and after % %
ALL = [Group1_PeakMags; Group2_PeakMags; Group3_PeakMags];
figure()
boxplot(ALL,groups,'Whisker',4)
title('PBS Event Peaks (Burns//Between burn and timehitR//Post timehitR)');
ylabel('Peak (dF/F)');
ylim([-0.2, 2]);
savefig('PBS_boxplot_PeakMags_3groups.fig');

% % boxplots of event durations before time zero and after % %
ALL = [Group1_PeakDurs; Group2_PeakDurs; Group3_PeakDurs];
figure()
boxplot(ALL,groups,'Whisker',4)
title('PBS Event Durations (Burns//Between burn and timehitR//Post timehitR)');
ylabel('duration (s)');
ylim([0,300]);
savefig('PBS_boxplot_PeakDurations_3groups.fig');

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
allIEI = [];

allTMEV2Traces = {};
allTMEV2AUC = {};
allTMEV2Durs = {};
allTMEV2Mags = {};
allTMEV2ROI = {};
allTMEV2path = {};
allTMEV2DPI = {};
allTMEV2NormalLocs = {};
allTMEV2ShiftedLocs = {};
allTMEV2EventGroup = {};
allTMEV2Mice = {};
allTMEV2IEI = {};

Group1_AUC = [];
Group2_AUC = [];
Group3_AUC = [];

Group1_PeakMags = [];
Group2_PeakMags = [];
Group3_PeakMags = [];

Group1_PeakDurs = [];
Group2_PeakDurs = [];
Group3_PeakDurs = [];

Group1_Locs = [];
Group2_Locs = [];
Group3_Locs = [];

% % count number of mice for TMEV2 % %
TMEV_logicals = contains(DATA(:,4),'TMEV'); 
TMEV_extracted = DATA(TMEV_logicals,:);
% TMEV2_extracted = TMEV_extracted(1:19,:); % soma
TMEV2_extracted = TMEV_extracted(1:64,:); % process
mice = TMEV2_extracted(:,7);
uniqueTMEV2 = unique(cellfun(@num2str,mice,'uni',0));
TMEV2_mice_n = numel(uniqueTMEV2);
TMEV2_mice_str = ['(N = ' num2str(TMEV2_mice_n) 'mice)'];
for k  = 1:numel(TMEV2fields)
    % extract
    tempIEI = TMEV_2DPI.(TMEV2fields{k}).IEI;
    allIEI = [allIEI;tempIEI];
    temp_timehitR = TMEV_2DPI.(TMEV2fields{k}).timehitR;
    tempLocs = TMEV_2DPI.(TMEV2fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
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
    % % % Conditional if excluding burns % % %
    if burnFlag
        
        noBurnLocs = tempLocs(tempLocs > burnFrame); % exclude burn
        tempPeaks = tempPeaks(tempLocs > burnFrame);
        tempRAWAUC = tempRAWAUC(tempLocs > burnFrame);
        tempDurs = tempDurs(tempLocs > burnFrame);
    else
        noBurnLocs = tempLocs;
        tempPeaks = tempPeaks;
        tempRAWAUC = tempRAWAUC;
        tempDurs = tempDurs;
        
        %%%%% Break up into 3 groups (burn, b/w burn and timehitR, post
        %%%%% timehitR)
        temp_Group1_AUC = tempRAWAUC(tempLocs <= burnFrame);
        temp_Group2_AUC = tempRAWAUC((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_AUC = tempRAWAUC(tempLocs >=  round(temp_timehitR*0.97));
        Group1_AUC = [Group1_AUC; temp_Group1_AUC];
        Group2_AUC = [Group2_AUC; temp_Group2_AUC];
        Group3_AUC = [Group3_AUC; temp_Group3_AUC];
        
        temp_Group1_PeakMags = tempPeaks(tempLocs <= burnFrame);
        temp_Group2_PeakMags = tempPeaks((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_PeakMags = tempPeaks(tempLocs >=  round(temp_timehitR*0.97));
        Group1_PeakMags = [Group1_PeakMags; temp_Group1_PeakMags];
        Group2_PeakMags = [Group2_PeakMags; temp_Group2_PeakMags];
        Group3_PeakMags = [Group3_PeakMags; temp_Group3_PeakMags];
        
        temp_Group1_PeakDurs = tempDurs(tempLocs <= burnFrame);
        temp_Group2_PeakDurs = tempDurs((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_PeakDurs = tempDurs(tempLocs >=  round(temp_timehitR*0.97));
        Group1_PeakDurs = [Group1_PeakDurs; temp_Group1_PeakDurs];
        Group2_PeakDurs = [Group2_PeakDurs; temp_Group2_PeakDurs];
        Group3_PeakDurs = [Group3_PeakDurs; temp_Group3_PeakDurs];
        
        temp_Group1_Locs = tempLocs(tempLocs <= burnFrame);
        temp_Group2_Locs = tempLocs((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_Locs = tempLocs(tempLocs >=  round(temp_timehitR*0.97));
        Group1_Locs = [Group1_Locs; temp_Group1_Locs];
        Group2_Locs = [Group2_Locs; temp_Group2_Locs];
        Group3_Locs = [Group3_Locs; temp_Group3_Locs];
    end
    
    % update all trace names
    for i = 1:length(tempRAWAUC)
        allTMEV2Traces = [allTMEV2Traces; TMEV2fields(k)];
        allTMEV2AUC = [allTMEV2AUC; tempRAWAUC(i)];
        allTMEV2Durs = [allTMEV2Durs; tempDurs(i)];
        allTMEV2Mags = [allTMEV2Mags; tempPeaks(i)];
        allTMEV2ROI = [allTMEV2ROI; TMEV_2DPI.(TMEV2fields{k}).ROI];
        allTMEV2path = [allTMEV2path; TMEV_2DPI.(TMEV2fields{k}).Pathology];
        allTMEV2DPI = [allTMEV2DPI; TMEV_2DPI.(TMEV2fields{k}).DPI];
        allTMEV2NormalLocs = [allTMEV2NormalLocs; tempLocs(i)];
        allTMEV2ShiftedLocs = [allTMEV2ShiftedLocs; tempLocs(i)-round(temp_timehitR*0.97)];
        if i <= length(tempIEI)
            allTMEV2IEI = [allTMEV2IEI; tempIEI(i)];
        else
            allTMEV2IEI = [allTMEV2IEI; 'Skipped b/c last peak'];
        end
        if tempLocs(i) < burnFrame
            allTMEV2EventGroup = [allTMEV2EventGroup; 1];
        elseif tempLocs(i) >= round(temp_timehitR*0.97)
            allTMEV2EventGroup = [allTMEV2EventGroup; 3];
        else
            allTMEV2EventGroup = [allTMEV2EventGroup; 2];
        end

    end
    
    allTMEV2Peaks = [allTMEV2Peaks; tempPeaks];
    all_RAWauc = [all_RAWauc; tempRAWAUC];
    all_auc = [all_auc; tempAUC];
    shiftedLocs = noBurnLocs - temp_timehitR;
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
for j = 1:length(mice)
    for k = 1:length(allTMEV2Traces)
        if strcmp(TMEV2_extracted{j,1},allTMEV2Traces{k})
            allTMEV2Mice = [allTMEV2Mice; mice(j)];
        end
    end
end
TMEV2_Data = cell(length(allTMEV2Traces)+1, 12);
TMEV2_Data(1,:) = {'Trace Name','ROI #','DPI','TMEV/PBS','Duration','Peak Amp','AUC','Peak Locations','Shifted Locations','Event Group','MouseID','IEI'};
TMEV2_Data(2:end,:) = [allTMEV2Traces,allTMEV2ROI, allTMEV2DPI, allTMEV2path, allTMEV2Durs,allTMEV2Mags,allTMEV2AUC, allTMEV2NormalLocs, allTMEV2ShiftedLocs,allTMEV2EventGroup,allTMEV2Mice,allTMEV2IEI];
save('Process_TMEV_2DPI_Data.mat','TMEV2_Data');
% % Stacked TMEV2 Traces % %
shiftedTime = [-1739:1740]./0.97;
figure()
xlabel('Relative time(s)'); ylabel('Trace Names');
for k = 1:numel(TMEV2fields)
    temp_timehitR = TMEV_2DPI.(TMEV2fields{k}).timehitR;
    R = round(temp_timehitR*0.97, 1);
    tempLocs = TMEV_2DPI.(TMEV2fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    if temp_timehitR > frames 
        shiftedLocs = (tempLocs - R);
    else
        shiftedLocs = (tempLocs - R);
    end
    tempPeaks = TMEV_2DPI.(TMEV2fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));
    shiftedLocs = shiftedLocs *0.97;
    tempF = TMEV_2DPI.(TMEV2fields{k}).trace;
    F0 = mean(tempF(1:50));
    tempF = (tempF - F0) / F0;
    tempTime = [1:numel(tempF)];
    tempTime = tempTime./0.97 - R;
    hold on
    plot(tempTime, tempF + k);
    hold on
    plot(shiftedLocs, tempF(tempLocs) + k, 'db', 'MarkerFaceColor','r','MarkerSize',6);
end
yticks([1:size(TMEV2fields)]);
yticklabels(TMEV2fields);
savefig('TMEV_2DPI_stacked_traces_NoBurns.fig');

% % Shifted scatter of JUST peaks % %
figure()
xlabel('Relative time(s)'); ylabel('Trace Names');
for k = 1:numel(TMEV2fields)
    temp_timehitR = TMEV_2DPI.(TMEV2fields{k}).timehitR;
    R = round(temp_timehitR*0.97, 1);
    tempLocs = TMEV_2DPI.(TMEV2fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    if temp_timehitR > frames 
        shiftedLocs = (tempLocs - R)./0.97;
    else
        shiftedLocs = (tempLocs - R)./0.97;
    end
    tempPeaks = TMEV_2DPI.(TMEV2fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));   
    tempF = TMEV_2DPI.(TMEV2fields{k}).trace;
    F0 = mean(tempF(1:50));
    tempF = (tempF - F0) / F0;
    tempTime = (1:1740) - R;
    tempTime = tempTime./0.97;
    K = ones(size(shiftedLocs))*k;
    hold on
    plot(shiftedLocs, K, '.k', 'MarkerFaceColor','k','MarkerSize',6);
end
yticks([1:size(TMEV2fields)]);
yticklabels(TMEV2fields);

% % normal traces % %
figure()
xlabel('Relative time(s)'); ylabel('Trace Names');
for k = 1:numel(TMEV2fields)
    temp_timehitR = TMEV_2DPI.(TMEV2fields{k}).timehitR;
    R = round(temp_timehitR*0.97, 1);
    tempLocs = TMEV_2DPI.(TMEV2fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    tempPeaks = TMEV_2DPI.(TMEV2fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    tempF = TMEV_2DPI.(TMEV2fields{k}).trace;
    F0 = mean(tempF(1:50));
    tempF = (tempF - F0) / F0;
    hold on
    plot(tempF + k);
    hold on
    plot(tempLocs, tempF(tempLocs) + k, 'db', 'MarkerFaceColor','r','MarkerSize',6);
end
yticks([1:size(TMEV2fields)]);
yticklabels(TMEV2fields);

%%%% TMEV 2 DPI Plotting Zones %%%%
% % Zones - EVENT AUC % %
figure()
plot(Group1_Locs, Group1_AUC,'kd','MarkerFaceColor','g');
hold on
plot(Group2_Locs, Group2_AUC,'kd','MarkerFaceColor','b');
hold on
plot(Group3_Locs, Group3_AUC,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('AUC (dF/F)');
title('TMEV 2 DPI: Zones 1-3 AUC');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV2DPI_Zones1-3_AUC.fig');

% % Zones - EVENT Duration % %
figure()
plot(Group1_Locs, Group1_PeakDurs,'kd','MarkerFaceColor','g');
hold on
plot(Group2_Locs, Group2_PeakDurs,'kd','MarkerFaceColor','b');
hold on
plot(Group3_Locs, Group3_PeakDurs,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('Peak durations (s)');
title('TMEV 2 DPI: Zones 1-3 Peak Duration');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV2DPI_Zones1-3_PeakDuration.fig');
% % Zones - EVENT Peak Magnitude % %
figure()
plot(Group1_Locs, Group1_PeakMags,'kd','MarkerFaceColor','g');
hold on
plot(Group2_Locs, Group2_PeakMags,'kd','MarkerFaceColor','b');
hold on
plot(Group3_Locs, Group3_PeakMags,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('Peak Amplitude (dF/F)');
title('TMEV 2 DPI: Zones 1-3 Peak Amplitude');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV2DPI_Zones1-3_PeakAmp.fig');

% % Zones - Duration vs. Amplitude % %
figure()
plot(Group1_PeakMags, Group1_PeakDurs,'kd','MarkerFaceColor','g');
hold on
plot(Group2_PeakMags, Group2_PeakDurs,'kd','MarkerFaceColor','b');
hold on
plot(Group3_PeakMags, Group3_PeakDurs,'kd','MarkerFaceColor','r');
xlabel('Peak Amplitude (dF/F)'); ylabel('Peak Duration (s)');
title('TMEV 2 DPI: Zones 1-3 Peak Duration vs Amplitude');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV2DPI_Zones1-3_DurationvsAmp.fig');

% % Zones - sqrt(AUC) vs. Amplitude % %
figure()
plot(Group1_PeakMags, sqrt(Group1_AUC),'kd','MarkerFaceColor','g');
hold on
plot(Group2_PeakMags, sqrt(Group2_AUC),'kd','MarkerFaceColor','b');
hold on
plot(Group3_PeakMags, sqrt(Group3_AUC),'kd','MarkerFaceColor','r');
xlabel('Peak Amplitude (dF/F)'); ylabel('sqrt(AUC) (sqrt(dF/F))');
title('TMEV 2 DPI: Zones 1-3 Event sqrt(AUC) vs Amplitude');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV2DPI_Zones1-3_sqrtAUCVsAmp.fig');

%%%% End Zone Plotting %%%%

binEdges = [-1600:100:1600];
numBins = numel(binEdges);
figure()
h = histogram(allTMEV2Locs,numBins, 'BinEdges', binEdges);
title('All Shifted TMEV 2 DPI Peak Locations');
xlabel('Relative time (s)'); ylabel('Frequency of events');
savefig('TMEV_2DPI_LocsHistogram_NoBurns.fig');
histHeights = h.Values./PBS_mice_n;
figure()
binEdges = [-1500:100:1600];
bg = bar(binEdges, histHeights);
title(strcat('TMEV 2DPI: Event frequency normalized to # of animals',TMEV2_mice_str)); ylabel('Events / Animal'); xlabel('Relative time binned (s)');
ylim([0 5]);
set(bg, 'FaceColor','r');


figure()
plot(allTMEV2Locs, allPeakDurs, 'dc','MarkerFaceColor','b')
title('TMEV 2 DPI: Peak Durations in relative time'); xlabel('Peak locs(s)'); ylabel('Peak Duration (s)');
ylim([0 300])
savefig('TMEV_2DPI_PeakDurations_NoBurns.fig');

figure()
plot(allTMEV2Locs,all_RAWauc,'d');
title('TMEV 2 DPI: AUC vs. Relative time: Raw Signal');
xlabel('relative time (s)'); ylabel('AUC');
ylim([0, 400]);
savefig('TMEV_2DPI_AUC_scatter_NoBurns.fig');

% % amp vs duration % %
figure()
plot(allPeakDurs, allTMEV2Peaks,'bd','MarkerFaceColor','r');
xlabel('Peak durations (s)'); ylabel('Peak magnitudes (dF/F)'); title('TMEV 2 DPI: Peak vs duration');
ylim([-0.2, 4]);
[c,p1] = polyfit(allPeakDurs,allTMEV2Peaks,1);
y_est = polyval(c,allPeakDurs);
hold on
plot(allPeakDurs, y_est, 'r--','LineWidth',2);

% % Amplitude vs. sqrt(AUC) % %
figure()
plot(sqrt(all_RAWauc), allTMEV2Peaks, 'd');
xlabel('sqrt(AUC)'); ylabel('Peak Amplitude'); title('TMEV 2 DPI: amp vs auc');
[c2,p2] = polyfit(sqrt(all_RAWauc),allTMEV2Peaks,1);
y_est2 = polyval(c2,sqrt(all_RAWauc));
hold on
plot(sqrt(all_RAWauc), y_est2,'--r');

% % Plotting both linear regressions % %
figure()
plot(sqrt(all_RAWauc), y_est2,'--r','LineWidth',2);
hold on
plot(allPeakDurs, y_est, 'b--','LineWidth',2);
xlabel('duration and auc'); ylabel('Peak amp'); title('TMEV 2 DPI: two linear regressions');
legend('Amp vs AUC','Amp vs Duration'); grid on;


figure()
plot(allTMEV2Locs, allTMEV2Peaks, 'dc','MarkerFaceColor','b')
title('TMEV 2 DPI: Peak Magnitudes in relative time'); xlabel('Relative time (s)'); ylabel('Peak Height (dF/F)');
savefig('TMEV_2DPI_PeakMags_scatter_NoBurns.fig');
ylim([-0.2, 4]);
% % Plotting Inter-event Intervals % %
figure()
xlabel('Relative Time (s)'); ylabel('IEI (s) plotted @ first peak'); title('TMEV2: IEI');
for k = 1:length(TMEV2fields)
    tempLocs = TMEV_2DPI.(TMEV2fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    IEILocs = tempLocs(1:end);
    tempIEI = TMEV_2DPI.(TMEV2fields{k}).IEI;
    if burnFlag
        tempLocs = tempLocs(tempLocs > 85);
        tempIEI = tempIEI(IEILocs(1:end-1) > 85);
    else
        tempLocs = tempLocs;
        tempIEI = tempIEI;
    end
    tempLocs = tempLocs(1:end-1); % exclude last peak - plotting IEI on first peak in calculation
    tempLocs = tempLocs./0.97;
    R = TMEV_2DPI.(TMEV2fields{k}).timehitR;
    tempLocs = tempLocs - R;
    
    hold on
    plot(tempLocs, tempIEI, 'db','MarkerFaceColor','r');
end
ylim([0,1500]);
savefig('TMEV_2DPI_IEI_scatter_NoBurns.fig');
% % boxplots of event auc before time zero and after % %
ALL = [Group1_AUC; Group2_AUC; Group3_AUC];
groups = [ones(size(Group1_AUC)); ones(size(Group2_AUC)) * 2; ones(size(Group3_AUC)) * 3];
figure()
boxplot(ALL,groups)
title('TMEV 2 DPI: Event AUC (Burns//Between burn and timehitR//Post timehitR)');
ylabel('AUC (dF/F)');
ylim([0, 400]);
savefig('TMEV_2DPI_AUC_boxplot_3groups.fig');
% % boxplots of event Peaks before time zero and after % %
ALL = [Group1_PeakMags; Group2_PeakMags; Group3_PeakMags];
figure()
boxplot(ALL,groups,'Whisker',4)
title('TMEV 2 DPI: Event Peaks (Burns//Between burn and timehitR//Post timehitR)');
ylabel('Peak (dF/F)');
ylim([-0.2, 2.0]);
savefig('TMEV_2DPI_PeakMags_boxplot_3groups.fig');
% % boxplots of event durations before time zero and after % %
ALL = [Group1_PeakDurs; Group2_PeakDurs; Group3_PeakDurs];
figure()
boxplot(ALL,groups,'Whisker',4)
title('TMEV 2 DPI: Event Durations (Burns//Between burn and timehitR//Post timehitR)');
ylabel('duration (s)');
ylim([0, 300]);
savefig('TMEV_2DPI_EventDurations_boxplot_3groups.fig');
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

allTMEV5Traces = {};
allTMEV5AUC = {};
allTMEV5Durs = {};
allTMEV5Mags = {};
allTMEV5ROI = {};
allTMEV5path = {};
allTMEV5DPI = {};
allTMEV5NormalLocs = {};
allTMEV5ShiftedLocs = {};
allTMEV5EventGroup = {};
allTMEV5Mice = {};
allTMEV5IEI = {};

Group1_AUC = [];
Group2_AUC = [];
Group3_AUC = [];

Group1_PeakMags = [];
Group2_PeakMags = [];
Group3_PeakMags = [];

Group1_PeakDurs = [];
Group2_PeakDurs = [];
Group3_PeakDurs = [];

Group1_Locs = [];
Group2_Locs = [];
Group3_Locs = [];

% % count number of mice for TMEV2 % %
% TMEV5_extracted = TMEV_extracted(20:40,:); % soma
TMEV5_extracted = TMEV_extracted(65:80,:); % process
mice = TMEV5_extracted(:,7);
uniqueTMEV5 = unique(cellfun(@num2str,mice,'uni',0));
TMEV5_mice_n = numel(uniqueTMEV5);
TMEV5_mice_str = ['(N = ' num2str(TMEV5_mice_n) 'mice)'];
for k  = 1:numel(TMEV5fields)
    % extract
    temp_timehitR = TMEV_5DPI.(TMEV5fields{k}).timehitR;
    tempLocs = TMEV_5DPI.(TMEV5fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
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
    % % % Conditional if excluding burns % % %
    if burnFlag
        
        noBurnLocs = tempLocs(tempLocs > burnFrame); % exclude burn
        tempPeaks = tempPeaks(tempLocs > burnFrame);
        tempRAWAUC = tempRAWAUC(tempLocs > burnFrame);
        tempDurs = tempDurs(tempLocs > burnFrame);
    else
        noBurnLocs = tempLocs;
        tempPeaks = tempPeaks;
        tempRAWAUC = tempRAWAUC;
        tempDurs = tempDurs;
        %%%%% Break up into 3 groups (burn, b/w burn and timehitR, post
        %%%%% timehitR)
        temp_Group1_AUC = tempRAWAUC(tempLocs <= burnFrame);
        temp_Group2_AUC = tempRAWAUC((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_AUC = tempRAWAUC(tempLocs >=  round(temp_timehitR*0.97));
        Group1_AUC = [Group1_AUC; temp_Group1_AUC];
        Group2_AUC = [Group2_AUC; temp_Group2_AUC];
        Group3_AUC = [Group3_AUC; temp_Group3_AUC];
        
        temp_Group1_PeakMags = tempPeaks(tempLocs <= burnFrame);
        temp_Group2_PeakMags = tempPeaks((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_PeakMags = tempPeaks(tempLocs >=  round(temp_timehitR*0.97));
        Group1_PeakMags = [Group1_PeakMags; temp_Group1_PeakMags];
        Group2_PeakMags = [Group2_PeakMags; temp_Group2_PeakMags];
        Group3_PeakMags = [Group3_PeakMags; temp_Group3_PeakMags];
        
        temp_Group1_PeakDurs = tempDurs(tempLocs <= burnFrame);
        temp_Group2_PeakDurs = tempDurs((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_PeakDurs = tempDurs(tempLocs >=  round(temp_timehitR*0.97));
        Group1_PeakDurs = [Group1_PeakDurs; temp_Group1_PeakDurs];
        Group2_PeakDurs = [Group2_PeakDurs; temp_Group2_PeakDurs];
        Group3_PeakDurs = [Group3_PeakDurs; temp_Group3_PeakDurs];
        
        temp_Group1_Locs = tempLocs(tempLocs <= burnFrame);
        temp_Group2_Locs = tempLocs((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_Locs = tempLocs(tempLocs >=  round(temp_timehitR*0.97));
        Group1_Locs = [Group1_Locs; temp_Group1_Locs];
        Group2_Locs = [Group2_Locs; temp_Group2_Locs];
        Group3_Locs = [Group3_Locs; temp_Group3_Locs];
    end
    
    % update all trace names
    for i = 1:length(tempRAWAUC)
        allTMEV5Traces = [allTMEV5Traces; TMEV5fields(k)];
        allTMEV5AUC = [allTMEV5AUC; tempRAWAUC(i)];
        allTMEV5Durs = [allTMEV5Durs; tempDurs(i)];
        allTMEV5Mags = [allTMEV5Mags; tempPeaks(i)];
        allTMEV5ROI = [allTMEV5ROI; TMEV_5DPI.(TMEV5fields{k}).ROI];
        allTMEV5path = [allTMEV5path; TMEV_5DPI.(TMEV5fields{k}).Pathology];
        allTMEV5DPI = [allTMEV5DPI; TMEV_5DPI.(TMEV5fields{k}).DPI];
        allTMEV5NormalLocs = [allTMEV5NormalLocs; tempLocs(i)];
        allTMEV5ShiftedLocs = [allTMEV5ShiftedLocs; tempLocs(i)-round(temp_timehitR*0.97)];
        if i <= length(tempIEI)
            allTMEV5IEI = [allTMEV5IEI; tempIEI(i)];
        else
            allTMEV5IEI = [allTMEV5IEI; 'Skipped b/c last peak'];
        end
        if tempLocs(i) < burnFrame
            allTMEV5EventGroup = [allTMEV5EventGroup; 1];
        elseif tempLocs(i) >= round(temp_timehitR*0.97)
            allTMEV5EventGroup = [allTMEV5EventGroup; 3];
        else
            allTMEV5EventGroup = [allTMEV5EventGroup; 2];
        end

    end
    
    allTMEV5Peaks = [allTMEV5Peaks; tempPeaks];
    all_RAWauc = [all_RAWauc; tempRAWAUC];
    all_auc = [all_auc; tempAUC];
    shiftedLocs = noBurnLocs - temp_timehitR;
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
for j = 1:length(mice)
    for k = 1:length(allTMEV5Traces)
        if strcmp(TMEV5_extracted{j,1},allTMEV5Traces{k})
            allTMEV5Mice = [allTMEV5Mice; mice(j)];
        end
    end
end
TMEV5_Data = cell(length(allTMEV5Traces)+1, 12);
TMEV5_Data(1,:) = {'Trace Name','ROI #','DPI','TMEV/PBS','Duration','Peak Amp','AUC','Peak Locations','Shifted Locations','Event Group','MouseID','IEI'};
TMEV5_Data(2:end,:) = [allTMEV5Traces,allTMEV5ROI, allTMEV5DPI, allTMEV5path, allTMEV5Durs,allTMEV5Mags,allTMEV5AUC, allTMEV5NormalLocs, allTMEV5ShiftedLocs,allTMEV5EventGroup, allTMEV5Mice,allTMEV5IEI];
save('Process_TMEV_5DPI_Data.mat','TMEV5_Data');
% % Stacked TMEV5 Traces % %
shiftedTime = [-1739:1740]./0.97;
figure()
xlabel('Relative time(s)'); ylabel('Trace Names');
for k = 1:numel(TMEV5fields)
    temp_timehitR = TMEV_5DPI.(TMEV5fields{k}).timehitR;
    R = round(temp_timehitR*0.97, 1);
    tempLocs = TMEV_5DPI.(TMEV5fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    if temp_timehitR > frames 
        shiftedLocs = (tempLocs - R)./0.97;
    else
        shiftedLocs = (tempLocs - R)./0.97;
    end
    tempPeaks = TMEV_5DPI.(TMEV5fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));   
    tempF = TMEV_5DPI.(TMEV5fields{k}).trace;
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
savefig('TMEV_5DPI_stacked_traces_NoBurns.fig');

% % Shifted scatter of JUST peaks % %
figure()
xlabel('Relative time(s)'); ylabel('Trace Names');
for k = 1:numel(TMEV5fields)
    temp_timehitR = TMEV_5DPI.(TMEV5fields{k}).timehitR;
    R = round(temp_timehitR*0.97, 1);
    tempLocs = TMEV_5DPI.(TMEV5fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    if temp_timehitR > frames 
        shiftedLocs = (tempLocs - R)./0.97;
    else
        shiftedLocs = (tempLocs - R)./0.97;
    end
    tempPeaks = TMEV_5DPI.(TMEV5fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));   
    tempF = TMEV_5DPI.(TMEV5fields{k}).trace;
    F0 = mean(tempF(1:50));
    tempF = (tempF - F0) / F0;
    tempTime = (1:1740) - R;
    tempTime = tempTime./0.97;
    K = ones(size(shiftedLocs))*k;
    hold on
    plot(shiftedLocs, K, '.k', 'MarkerFaceColor','k','MarkerSize',6);
end
yticks([1:size(TMEV5fields)]);
yticklabels(TMEV5fields);

% % normal traces % %
figure()
xlabel('Relative time(s)'); ylabel('Trace Names');
for k = 1:numel(TMEV5fields)
    temp_timehitR = TMEV_5DPI.(TMEV5fields{k}).timehitR;
    R = round(temp_timehitR*0.97, 1);
    tempLocs = TMEV_5DPI.(TMEV5fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    tempPeaks = TMEV_5DPI.(TMEV5fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));   
    tempF = TMEV_5DPI.(TMEV5fields{k}).trace;
    F0 = mean(tempF(1:50));
    tempF = (tempF - F0) / F0;
    hold on
    plot(tempF + k);
    hold on
    plot(tempLocs, tempF(tempLocs) + k, 'db', 'MarkerFaceColor','r','MarkerSize',6);
end
yticks([1:size(TMEV5fields)]);
yticklabels(TMEV5fields);

%%%% TMEV 5 DPI Plotting Zones %%%%
% % Zones - EVENT AUC % %
figure()
plot(Group1_Locs, Group1_AUC,'kd','MarkerFaceColor','g');
hold on
plot(Group2_Locs, Group2_AUC,'kd','MarkerFaceColor','b');
hold on
plot(Group3_Locs, Group3_AUC,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('AUC (dF/F)');
title('TMEV 5 DPI: Zones 1-3 AUC');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV5DPI_Zones1-3_AUC.fig');

% % Zones - EVENT Duration % %
figure()
plot(Group1_Locs, Group1_PeakDurs,'kd','MarkerFaceColor','g');
hold on
plot(Group2_Locs, Group2_PeakDurs,'kd','MarkerFaceColor','b');
hold on
plot(Group3_Locs, Group3_PeakDurs,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('Peak durations (s)');
title('TMEV 5 DPI: Zones 1-3 Peak Duration');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV5DPI_Zones1-3_PeakDuration.fig');

% % Zones - EVENT Peak Magnitude % %
figure()
plot(Group1_Locs, Group1_PeakMags,'kd','MarkerFaceColor','g');
hold on
plot(Group2_Locs, Group2_PeakMags,'kd','MarkerFaceColor','b');
hold on
plot(Group3_Locs, Group3_PeakMags,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('Peak Amplitude (dF/F)');
title('TMEV 5 DPI: Zones 1-3 Peak Amplitude');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV5DPI_Zones1-3_PeakAmp.fig');

% % Zones - Duration vs. Amplitude % %
figure()
plot(Group1_PeakMags, Group1_PeakDurs,'kd','MarkerFaceColor','g');
hold on
plot(Group2_PeakMags, Group2_PeakDurs,'kd','MarkerFaceColor','b');
hold on
plot(Group3_PeakMags, Group3_PeakDurs,'kd','MarkerFaceColor','r');
xlabel('Peak Amplitude (dF/F)'); ylabel('Peak Duration (s)');
title('TMEV 5 DPI: Zones 1-3 Peak Duration vs Amplitude');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV5DPI_Zones1-3_DurationVsAmp.fig');

% % Zones - sqrt(AUC) vs. Amplitude % %
figure()
plot(Group1_PeakMags, sqrt(Group1_AUC),'kd','MarkerFaceColor','g');
hold on
plot(Group2_PeakMags, sqrt(Group2_AUC),'kd','MarkerFaceColor','b');
hold on
plot(Group3_PeakMags, sqrt(Group3_AUC),'kd','MarkerFaceColor','r');
xlabel('Peak Amplitude (dF/F)'); ylabel('sqrt(AUC) (sqrt(dF/F))');
title('TMEV 5 DPI: Zones 1-3 Event sqrt(AUC) vs Amplitude');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV5DPI_Zones1-3_sqrtAUCVsAmp.fig');
%%%% End Zone Plotting %%%%

binEdges = [-1600:100:1600];
numBins = numel(binEdges);
figure()
h = histogram(allTMEV5Locs,numBins, 'BinEdges', binEdges);
title('All Shifted TMEV 5 DPI Peak Locations');
xlabel('Relative time (s)'); ylabel('Frequency of events');
savefig('TMEV_5DPI_LocationHistogram_NoBurns.fig');
histHeights = h.Values./TMEV5_mice_n;
figure()
binEdges = [-1600:100:1500];
bg = bar(binEdges, histHeights);
title(strcat('TMEV 5DPI: Event frequency normalized to # of animals',TMEV5_mice_str)); ylabel('Events / Animal'); xlabel('Relative time binned (s)');
ylim([0 8.5]);
set(bg,'FaceColor','r');

figure()
plot(allTMEV5Locs, allPeakDurs, 'dc','MarkerFaceColor','b')
title('tmev 5 dpi: Peak Durations in relative time'); xlabel('Relative time (s)'); ylabel('Peak Duration (s)');
ylim([0 300]);
savefig('TMEV_5DPI_PeakDurations_NoBurns.fig');

figure()
plot(allTMEV5Locs,all_RAWauc,'d');
title('tmev 5 dpi: AUC vs. Relative time: Raw Signal');
xlabel('relative time (s)'); ylabel('AUC');
savefig('TMEV_5DPI_AUC_scatter_NoBurns.fig');
ylim([0, 400]);

% % Amp vs duration % %
figure()
plot(allPeakDurs, allTMEV5Peaks,'bd','MarkerFaceColor','r');
xlabel('Peak durations (s)'); ylabel('Peak magnitudes (dF/F)'); title('TMEV 5 DPI: Peak vs duration');
ylim([-0.2, 4]);
[c,p1] = polyfit(allPeakDurs,allTMEV5Peaks,1);
y_est = polyval(c,allPeakDurs);
hold on
plot(allPeakDurs, y_est, 'r--','LineWidth',2);

% % Amplitude vs. sqrt(AUC) % %
figure()
plot(sqrt(all_RAWauc), allTMEV5Peaks, 'd');
xlabel('sqrt(AUC)'); ylabel('Peak Amplitude'); title('TMEV 5 DPI: amp vs auc');
[c2,p2] = polyfit(sqrt(all_RAWauc),allTMEV5Peaks,1);
y_est2 = polyval(c2,sqrt(all_RAWauc));
hold on
plot(sqrt(all_RAWauc), y_est2,'--r');

% % Plotting both linear regressions % %
figure()
plot(sqrt(all_RAWauc), y_est2,'--r','LineWidth',2);
hold on
plot(allPeakDurs, y_est, 'b--','LineWidth',2);
xlabel('duration and auc'); ylabel('Peak amp'); title('TMEV 5 DPI: two linear regressions');
legend('Amp vs AUC','Amp vs Duration');
grid on

figure()
plot(allTMEV5Locs, allTMEV5Peaks, 'dc','MarkerFaceColor','b')
title('TMEV 5 DPI: Peak Magnitudes in relative time'); xlabel('Relative time (s)'); ylabel('Peak Height (dF/F)');
savefig('TMEV_5DPI_AUC_scatter_NoBurns.fig');
ylim([-0.2, 4]);
% % Plotting Inter-event Intervals % %
figure()
xlabel('Relative Time (s)'); ylabel('IEI (s) plotted @ first peak'); title('TMEV5: IEI');
for k = 1:length(TMEV5fields)
    tempLocs = TMEV_5DPI.(TMEV5fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    IEILocs = tempLocs(1:end);
    tempIEI = TMEV_5DPI.(TMEV5fields{k}).IEI;
    if burnFlag
        tempLocs = tempLocs(tempLocs > burnFrame);
        tempIEI = tempIEI(IEILocs(1:end-1) > burnFrame);
    else
        tempLocs = tempLocs;
        tempIEI = tempIEI;
    end
    tempLocs = tempLocs(1:end-1); % exclude last peak - plotting IEI on first peak in calculation
    tempLocs = tempLocs./0.97;
    R = TMEV_5DPI.(TMEV5fields{k}).timehitR;
    tempLocs = tempLocs - R;    
    hold on
    plot(tempLocs, tempIEI, 'db','MarkerFaceColor','r');
end
ylim([0,1500]);
savefig('TMEV_5DPI_IEI_scatter_NoBurns.fig');
% % boxplots of event auc before time zero and after % %
ALL = [Group1_AUC; Group2_AUC; Group3_AUC];
groups = [ones(size(Group1_AUC)); ones(size(Group2_AUC)) * 2; ones(size(Group3_AUC)) * 3];
figure()
boxplot(ALL,groups)
title('TMEV 5 DPI: Event AUC (Burns//Between burn and timehitR//Post timehitR)');
ylabel('AUC (dF/F)');
ylim([0, 400]);
savefig('TMEV_5DPI_AUC_boxplot_3groups.fig');
% % boxplots of event Peaks before time zero and after % %
ALL = [Group1_PeakMags; Group2_PeakMags; Group3_PeakMags];
figure()
boxplot(ALL,groups,'Whisker',4)
title('TMEV 5 DPI: Event Peaks (Burns//Between burn and timehitR//Post timehitR)');
ylabel('Peak (dF/F)');
ylim([-0.2, 2]);
savefig('TMEV_5DPI_PeakMags_boxplot_3groups.fig');
% % boxplots of event durations before time zero and after % %
ALL = [Group1_PeakDurs; Group2_PeakDurs; Group3_PeakDurs];
figure()
boxplot(ALL,groups,'Whisker',4)
title('TMEV 5 DPI: Event Durations (Burns//Between burn and timehitR//Post timehitR)');
ylabel('duration (s)');
ylim([0, 300]);
savefig('TMEV_5DPI_PeakDurations_boxplot_3groups.fig');
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

allTMEV14Traces = {};
allTMEV14AUC = {};
allTMEV14Durs = {};
allTMEV14Mags = {};
allTMEV14ROI = {};
allTMEV14path = {};
allTMEV14DPI = {};
allTMEV14NormalLocs = {};
allTMEV14ShiftedLocs = {};
allTMEV14EventGroup = {};
allTMEV14Mice = {};
allTMEV14IEI = {};

Group1_AUC = [];
Group2_AUC = [];
Group3_AUC = [];

Group1_PeakMags = [];
Group2_PeakMags = [];
Group3_PeakMags = [];

Group1_PeakDurs = [];
Group2_PeakDurs = [];
Group3_PeakDurs = [];

Group1_Locs = [];
Group2_Locs = [];
Group3_Locs = [];

% TMEV14_extracted = TMEV_extracted(41:end,:); % soma
TMEV14_extracted = TMEV_extracted(81:end,:); % process
mice = TMEV14_extracted(:,7);
uniqueTMEV14 = unique(cellfun(@num2str,mice,'uni',0));
TMEV14_mice_n = numel(uniqueTMEV14);
TMEV14_mice_str = ['(N = ' num2str(TMEV14_mice_n) 'mice)'];
for k  = 1:numel(TMEV14fields)
    % extract
    temp_timehitR = TMEV_14DPI.(TMEV14fields{k}).timehitR;
    tempLocs = TMEV_14DPI.(TMEV14fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
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
    % % % Conditional if excluding burns % % %
    if burnFlag
        
        noBurnLocs = tempLocs(tempLocs > burnFrame); % exclude burn
        tempPeaks = tempPeaks(tempLocs > burnFrame);
        tempRAWAUC = tempRAWAUC(tempLocs > burnFrame);
        tempDurs = tempDurs(tempLocs > burnFrame);
    else
        noBurnLocs = tempLocs;
        tempPeaks = tempPeaks;
        tempRAWAUC = tempRAWAUC;
        tempDurs = tempDurs;
        %%%%% Break up into 3 groups (burn, b/w burn and timehitR, post
        %%%%% timehitR)
        temp_Group1_AUC = tempRAWAUC(tempLocs <= burnFrame);
        temp_Group2_AUC = tempRAWAUC((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_AUC = tempRAWAUC(tempLocs >=  round(temp_timehitR*0.97));
        Group1_AUC = [Group1_AUC; temp_Group1_AUC];
        Group2_AUC = [Group2_AUC; temp_Group2_AUC];
        Group3_AUC = [Group3_AUC; temp_Group3_AUC];
        
        temp_Group1_PeakMags = tempPeaks(tempLocs <= burnFrame);
        temp_Group2_PeakMags = tempPeaks((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_PeakMags = tempPeaks(tempLocs >=  round(temp_timehitR*0.97));
        Group1_PeakMags = [Group1_PeakMags; temp_Group1_PeakMags];
        Group2_PeakMags = [Group2_PeakMags; temp_Group2_PeakMags];
        Group3_PeakMags = [Group3_PeakMags; temp_Group3_PeakMags];
        
        temp_Group1_PeakDurs = tempDurs(tempLocs <= burnFrame);
        temp_Group2_PeakDurs = tempDurs((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_PeakDurs = tempDurs(tempLocs >=  round(temp_timehitR*0.97));
        Group1_PeakDurs = [Group1_PeakDurs; temp_Group1_PeakDurs];
        Group2_PeakDurs = [Group2_PeakDurs; temp_Group2_PeakDurs];
        Group3_PeakDurs = [Group3_PeakDurs; temp_Group3_PeakDurs];
        
        temp_Group1_Locs = tempLocs(tempLocs <= burnFrame);
        temp_Group2_Locs = tempLocs((isburn(tempLocs)) & (tempLocs < round(temp_timehitR*0.97)));
        temp_Group3_Locs = tempLocs(tempLocs >=  round(temp_timehitR*0.97));
        Group1_Locs = [Group1_Locs; temp_Group1_Locs];
        Group2_Locs = [Group2_Locs; temp_Group2_Locs];
        Group3_Locs = [Group3_Locs; temp_Group3_Locs];
    end
    
    % update all trace names
    for i = 1:length(tempRAWAUC)
        allTMEV14Traces = [allTMEV14Traces; TMEV14fields(k)];
        allTMEV14AUC = [allTMEV14AUC; tempRAWAUC(i)];
        allTMEV14Durs = [allTMEV14Durs; tempDurs(i)];
        allTMEV14Mags = [allTMEV14Mags; tempPeaks(i)];
        allTMEV14ROI = [allTMEV14ROI; TMEV_14DPI.(TMEV14fields{k}).ROI];
        allTMEV14path = [allTMEV14path; TMEV_14DPI.(TMEV14fields{k}).Pathology];
        allTMEV14DPI = [allTMEV14DPI; TMEV_14DPI.(TMEV14fields{k}).DPI];
        allTMEV14NormalLocs = [allTMEV14NormalLocs; tempLocs(i)];
        allTMEV14ShiftedLocs = [allTMEV14ShiftedLocs; tempLocs(i)-round(temp_timehitR*0.97)];
        if i <= length(tempIEI)
            allTMEV14IEI = [allTMEV14IEI; tempIEI(i)];
        else
            allTMEV14IEI = [allTMEV14IEI; 'Skipped b/c last peak'];
        end
        if tempLocs(i) < burnFrame
            allTMEV14EventGroup = [allTMEV14EventGroup; 1];
        elseif tempLocs(i) >= round(temp_timehitR*0.97)
            allTMEV14EventGroup = [allTMEV14EventGroup; 3];
        else
            allTMEV14EventGroup = [allTMEV14EventGroup; 2];
        end

    end
    
    allTMEV14Peaks = [allTMEV14Peaks; tempPeaks];
    all_RAWauc = [all_RAWauc; tempRAWAUC];
    all_auc = [all_auc; tempAUC];
    shiftedLocs = noBurnLocs - temp_timehitR;
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
for j = 1:length(mice)
    for k = 1:length(allTMEV14Traces)
        if strcmp(TMEV14_extracted{j,1},allTMEV14Traces{k})
            allTMEV14Mice = [allTMEV14Mice; mice(j)];
        end
    end
end
TMEV14_Data = cell(length(allTMEV14Traces)+1, 12);
TMEV14_Data(1,:) = {'Trace Name','ROI #','DPI','TMEV/PBS','Duration','Peak Amp','AUC','Peak Locations','Shifted Locations','Event Group','MouseID','IEI'};
TMEV14_Data(2:end,:) = [allTMEV14Traces,allTMEV14ROI, allTMEV14DPI, allTMEV14path, allTMEV14Durs,allTMEV14Mags,allTMEV14AUC, allTMEV14NormalLocs, allTMEV14ShiftedLocs,allTMEV14EventGroup, allTMEV14Mice,allTMEV14IEI];
save('Process_TMEV_14DPI_Data.mat','TMEV14_Data');

% % Stacked TMEV14 Traces % %
shiftedTime = [-1739:1740]./0.97;
figure()
xlabel('Relative time(s)'); ylabel('Trace Names');
for k = 1:numel(TMEV14fields)
    temp_timehitR = TMEV_14DPI.(TMEV14fields{k}).timehitR;
    R = round(temp_timehitR*0.97, 1);
    tempLocs = TMEV_14DPI.(TMEV14fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    if temp_timehitR > frames 
        shiftedLocs = (tempLocs - R)./0.97;
    else
        shiftedLocs = (tempLocs - R)./0.97;
    end
    tempPeaks = TMEV_14DPI.(TMEV14fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));   
    tempF = TMEV_14DPI.(TMEV14fields{k}).trace;
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
savefig('TMEV_14DPI_stacked_traces_NoBurns.fig');

% % normal traces % %
figure()
xlabel('Relative time(s)'); ylabel('Trace Names');
for k = 1:numel(TMEV14fields)
    temp_timehitR = TMEV_14DPI.(TMEV14fields{k}).timehitR;
    R = round(temp_timehitR*0.97, 1);
    tempLocs = TMEV_14DPI.(TMEV14fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    tempPeaks = TMEV_14DPI.(TMEV14fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));   
    tempF = TMEV_14DPI.(TMEV14fields{k}).trace;
    F0 = mean(tempF(1:50));
    tempF = (tempF - F0) / F0;
    tempTime = tempTime./0.97;
    hold on
    plot(tempF + k);
    hold on
    plot(tempLocs, tempF(tempLocs) + k, 'db', 'MarkerFaceColor','r','MarkerSize',6);
end
yticks([1:size(TMEV14fields)]);
yticklabels(TMEV14fields);

%%%% TMEV 14 DPI Plotting Zones %%%%
% % Zones - EVENT AUC % %
figure()
plot(Group1_Locs, Group1_AUC,'kd','MarkerFaceColor','g');
hold on
plot(Group2_Locs, Group2_AUC,'kd','MarkerFaceColor','b');
hold on
plot(Group3_Locs, Group3_AUC,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('AUC (dF/F)');
title('TMEV 14 DPI: Zones 1-3 AUC');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV14DPI_Zones1-3_AUC.fig');

% % Shifted scatter of JUST peaks % %
figure()
xlabel('Relative time(s)'); ylabel('Trace Names');
for k = 1:numel(TMEV14fields)
    temp_timehitR = TMEV_14DPI.(TMEV14fields{k}).timehitR;
    R = round(temp_timehitR*0.97, 1);
    tempLocs = TMEV_14DPI.(TMEV14fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    if temp_timehitR > frames 
        shiftedLocs = (tempLocs - R)./0.97;
    else
        shiftedLocs = (tempLocs - R)./0.97;
    end
    tempPeaks = TMEV_14DPI.(TMEV14fields{k}).PeakMags;
    tempPeaks = tempPeaks(~isnan(tempPeaks));
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));   
    tempF = TMEV_14DPI.(TMEV14fields{k}).trace;
    F0 = mean(tempF(1:50));
    tempF = (tempF - F0) / F0;
    tempTime = (1:1740) - R;
    tempTime = tempTime./0.97;
    K = ones(size(shiftedLocs))*k;
    hold on
    plot(shiftedLocs, K, '.k', 'MarkerFaceColor','k','MarkerSize',6);
end
yticks([1:size(TMEV14fields)]);
yticklabels(TMEV14fields);

% % Zones - EVENT Duration % %
figure()
plot(Group1_Locs, Group1_PeakDurs,'kd','MarkerFaceColor','g');
hold on
plot(Group2_Locs, Group2_PeakDurs,'kd','MarkerFaceColor','b');
hold on
plot(Group3_Locs, Group3_PeakDurs,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('Peak durations (s)');
title('TMEV 14 DPI: Zones 1-3 Peak Duration');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV14DPI_Zones1-3_PeakDuration.fig');

% % Zones - EVENT Peak Magnitude % %
figure()
plot(Group1_Locs, Group1_PeakMags,'kd','MarkerFaceColor','g');
hold on
plot(Group2_Locs, Group2_PeakMags,'kd','MarkerFaceColor','b');
hold on
plot(Group3_Locs, Group3_PeakMags,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('Peak Amplitude (dF/F)');
title('TMEV 14 DPI: Zones 1-3 Peak Amplitude');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV14DPI_Zones1-3_PeakAmp.fig');

% % Zones - Duration vs. Amplitude % %
figure()
plot(Group1_PeakMags, Group1_PeakDurs,'kd','MarkerFaceColor','g');
hold on
plot(Group2_PeakMags, Group2_PeakDurs,'kd','MarkerFaceColor','b');
hold on
plot(Group3_PeakMags, Group3_PeakDurs,'kd','MarkerFaceColor','r');
xlabel('Peak Amplitude (dF/F)'); ylabel('Peak Duration (s)');
title('TMEV 14 DPI: Zones 1-3 Peak Duration vs Amplitude');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV14DPI_Zones1-3_DurationVsAmp.fig');

% % Zones - sqrt(AUC) vs. Amplitude % %
figure()
plot(Group1_PeakMags, sqrt(Group1_AUC),'kd','MarkerFaceColor','g');
hold on
plot(Group2_PeakMags, sqrt(Group2_AUC),'kd','MarkerFaceColor','b');
hold on
plot(Group3_PeakMags, sqrt(Group3_AUC),'kd','MarkerFaceColor','r');
xlabel('Peak Amplitude (dF/F)'); ylabel('sqrt(AUC) (sqrt(dF/F))');
title('TMEV 14 DPI: Zones 1-3 Event sqrt(AUC) vs Amplitude');
legend('Zone 1','Zone 2','Zone 3');
savefig('TMEV14DPI_Zones1-3_sqrtAUCVsAmp.fig');

%%%% End Zone Plotting %%%%

binEdges = [-1600:100:1600];
numBins = numel(binEdges);
figure()
h = histogram(allTMEV14Locs,numBins, 'BinEdges', binEdges);
title('All Shifted TMEV 14 DPI Peak Locations');
xlabel('Relative time (s)'); ylabel('Frequency of events');
savefig('TMEV_14DPI_LocationHistogram_NoBurns.fig');
histHeights = h.Values./TMEV5_mice_n;
figure()
binEdges = [-1600:100:1500];
bg = bar(binEdges, histHeights);
title(strcat('TMEV 14DPI: Event frequency normalized to # of animals',TMEV14_mice_str)); ylabel('Events / Animal'); xlabel('Relative time binned (s)');
ylim([0 8.5]);
set(bg,'FaceColor','r');


figure()
plot(allTMEV14Locs, allPeakDurs, 'dc','MarkerFaceColor','b')
title('tmev 14 dpi: Peak Durations in relative time'); xlabel('Relative time (s)'); ylabel('Peak Duration (s)');
ylim([0 300]);



figure()
plot(allTMEV14Locs,all_RAWauc,'d');
title('tmev 14 dpi: AUC vs. Relative time: Raw Signal');
xlabel('relative time (s)'); ylabel('AUC');
savefig('TMEV_14DPI_AUC_scatter_NoBurns.fig');
ylim([0, 400]);

% % Amp vs duration % %
figure()
plot(allPeakDurs, allTMEV14Peaks,'bd','MarkerFaceColor','r');
xlabel('Peak durations (s)'); ylabel('Peak magnitudes (dF/F)'); title('TMEV 14 DPI: Peak vs duration');
ylim([-0.2, 4]);
[c,p1] = polyfit(allPeakDurs,allTMEV14Peaks,1);
y_est = polyval(c,allPeakDurs);
hold on
plot(allPeakDurs, y_est, 'r--','LineWidth',2);

% % Amplitude vs. sqrt(AUC) % %
figure()
plot(sqrt(all_RAWauc), allTMEV14Peaks, 'd');
xlabel('sqrt(AUC)'); ylabel('Peak Amplitude'); title('TMEV 14 DPI: amp vs auc');
[c2,p2] = polyfit(sqrt(all_RAWauc),allTMEV14Peaks,1);
y_est2 = polyval(c2,sqrt(all_RAWauc));
hold on
plot(sqrt(all_RAWauc), y_est2,'--r');

% % Plotting both linear regressions % %
figure()
plot(sqrt(all_RAWauc), y_est2,'--r','LineWidth',2);
hold on
plot(allPeakDurs, y_est, 'b--','LineWidth',2);
xlabel('duration and auc'); ylabel('Peak amp'); title('TMEV 14 DPI: two linear regressions');
legend('Amp vs AUC','Amp vs Duration');
grid on

figure()
plot(allTMEV14Locs, allTMEV14Peaks, 'dc','MarkerFaceColor','b')
title('TMEV 14 DPI: Peak Magnitudes in relative time'); xlabel('Relative time (s)'); ylabel('Peak Height (dF/F)');
savefig('TMEV_14DPI_PeakMags_scatter_NoBurns.fig');
ylim([-0.2, 2.0]);

% % Plotting Inter-event Intervals % %
figure()
xlabel('Relative Time (s)'); ylabel('IEI (s) plotted @ first peak'); title('TMEV14: IEI');
savefig('TMEV_14DPI_IEI_scatter_NoBurns.fig');

for k = 1:length(TMEV14fields)
    tempLocs = TMEV_14DPI.(TMEV14fields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    IEILocs = tempLocs(1:end);
    tempIEI = TMEV_14DPI.(TMEV14fields{k}).IEI;
    if burnFlag
        tempLocs = tempLocs(tempLocs > burnFrame);
        tempIEI = tempIEI(IEILocs(1:end-1) > burnFrame);
    else
        tempLocs = tempLocs;
        tempIEI = tempIEI;
    end
    tempLocs = tempLocs(1:end-1); % exclude last peak - plotting IEI on first peak in calculation
    tempLocs = tempLocs./0.97;
    R = TMEV_14DPI.(TMEV14fields{k}).timehitR;
    tempLocs = tempLocs - R;
    
    hold on
    plot(tempLocs, tempIEI, 'db','MarkerFaceColor','r');
end
ylim([0,1500]);
savefig('TMEV_14DPI_IEI_scatter_NoBurns.fig');
% % boxplots of event auc before time zero and after % %
ALL = [Group1_AUC; Group2_AUC; Group3_AUC];
groups = [ones(size(Group1_AUC)); ones(size(Group2_AUC)) * 2; ones(size(Group3_AUC)) * 3];
figure()
boxplot(ALL,groups)
title('TMEV 14 DPI: Event AUC (Burns//Between burn and timehitR//Post timehitR)');
ylabel('AUC (dF/F)');
ylim([0, 400]);
savefig('TMEV_14DPI_AUC_boxplot_3groups.fig');

% % boxplots of event Peaks before time zero and after % %
ALL = [Group1_PeakMags; Group2_PeakMags; Group3_PeakMags];
figure()
boxplot(ALL,groups,'Whisker',4)
title('TMEV 14 DPI: Event Peaks (Burns//Between burn and timehitR//Post timehitR)');
ylabel('Peak (dF/F)');
ylim([-0.2, 2]);
savefig('TMEV_14DPI_PeakMags_boxplot_3groups.fig');

% % boxplots of event durations before time zero and after % %
ALL = [Group1_PeakDurs; Group2_PeakDurs; Group3_PeakDurs];
figure()
boxplot(ALL,groups,'Whisker',4)
title('TMEV 14 DPI: Event Durations (Burns//Between burn and timehitR//Post timehitR)');
ylabel('duration (s)');
ylim([0, 300]);
savefig('TMEV_14DPI_PeakDurationsws_boxplot_3groups.fig');
%%
AL = [];
PBS_A_RMS = [];
PBS_B_RMS = [];
PBS_C_RMS = [];

TMEV2_A_RMS = [];
TMEV2_B_RMS = [];
TMEV2_C_RMS = [];

TMEV5_A_RMS = [];
TMEV5_B_RMS = [];
TMEV5_C_RMS = [];

TMEV14_A_RMS = [];
TMEV14_B_RMS = [];
TMEV14_C_RMS = [];

BL = [];
PBS_R = [];
%     TMEV_R = [];
allRMS{1,1} = 'Trace + ROI';
allRMS{1,2} = 'TMEV/PBS';
allRMS{1,3} = 'DPI';
allRMS{1,4} = 'A(1 to 200) RMS^2';
allRMS{1,5} = 'B(201 to timehitR) RMS^2';
allRMS{1,6} = 'C(timehitR to end) RMS^2';
allRMS{1,7} = 'timehitR Marker/Index';
allRMS{1,8} = 'PBS/TMEV';
filtered_traces = cell(size(process_data,1),3);
for ii = 1:size(process_data,1)
    ROInum = process_data{ii,2};
    ROInum = num2str(ROInum);
    ROI = strcat('  ROI: ', ROInum);
    traceName = process_data{ii,1};
    traceName = strcat(traceName,ROI);
    DPI = process_data{ii,5};
    timehitR = double(process_data{ii,6});
    timehitR_marker = round(timehitR);
    currentTrace = process_data{ii,3};
    F0 = mean(currentTrace(1:50));
    currentTrace = (currentTrace-F0)./F0;
    if timehitR_marker > 200 && timehitR_marker < length(currentTrace)
        
        rmsA = rms(currentTrace(1:200));
        rmsB = rms(currentTrace(201:timehitR_marker-1));
        rmsC = rms(currentTrace(timehitR_marker:end));
        allRMS{ii+1,1} = traceName;
        allRMS{ii+1,2} = process_data{ii,4};
        allRMS{ii+1,3} = DPI;
        allRMS{ii+1,4} = rmsA;
        allRMS{ii+1,5} = rmsB;
        allRMS{ii+1,6} = rmsC;
        allRMS{ii+1,7} = timehitR_marker;
        if strcmp(process_data{ii,4},'PBS')
            % % Find RMS of A and B and C portions % %
            PBS_A_RMS = [PBS_A_RMS; rmsA];
            PBS_B_RMS = [PBS_B_RMS; rmsB];
            PBS_C_RMS = [PBS_C_RMS; rmsC];
            allRMS{ii+1,8} = 'PBS';
        else
            allRMS{ii+1,8} = 'TMEV';
            if DPI == 2
                TMEV2_A_RMS = [TMEV2_A_RMS; rmsA];
                TMEV2_B_RMS = [TMEV2_B_RMS; rmsB];
                TMEV2_C_RMS = [TMEV2_C_RMS; rmsC];
            elseif DPI == 4 || DPI == 5 || DPI == 6
                TMEV5_A_RMS = [TMEV5_A_RMS; rmsA];
                TMEV5_B_RMS = [TMEV5_B_RMS; rmsB];
                TMEV5_C_RMS = [TMEV5_C_RMS; rmsC];
            else
                TMEV14_A_RMS = [TMEV14_A_RMS; rmsA];
                TMEV14_B_RMS = [TMEV14_B_RMS; rmsB];
                TMEV14_C_RMS = [TMEV14_C_RMS; rmsC];
            end            
        end
    else
        rmsA = 'Skipped';
        rmsB = 'Skipped';
        rmsC = 'Skipped';
    end 
end
    
% end
PBS_A_MEAN = mean(PBS_A_RMS);
PBS_B_MEAN = mean(PBS_B_RMS);
PBS_C_MEAN = mean(PBS_C_RMS);
PBSy = [PBS_A_MEAN, PBS_B_MEAN, PBS_C_MEAN];
PBS_sd_vct = [std(PBS_A_RMS)/sqrt(length(PBS_A_RMS)), std(PBS_B_RMS)/sqrt(length(PBS_B_RMS)), std(PBS_C_RMS)/sqrt(length(PBS_C_RMS))];

TMEV2_A_MEAN = mean(TMEV2_A_RMS);
TMEV2_B_MEAN = mean(TMEV2_B_RMS);
TMEV2_C_MEAN = mean(TMEV2_C_RMS);
TMEV2y = [TMEV2_A_MEAN, TMEV2_B_MEAN, TMEV2_C_MEAN];
TMEV2_sd_vct = [std(TMEV2_A_RMS)/sqrt(length(TMEV2_A_RMS)), std(TMEV2_B_RMS)/sqrt(length(TMEV2_B_RMS)), std(TMEV2_C_RMS)/sqrt(length(TMEV2_C_RMS))];

TMEV5_A_MEAN = mean(TMEV5_A_RMS);
TMEV5_B_MEAN = mean(TMEV5_B_RMS);
TMEV5_C_MEAN = mean(TMEV5_C_RMS);
TMEV5y = [TMEV5_A_MEAN, TMEV5_B_MEAN, TMEV5_C_MEAN];
TMEV5_sd_vct = [std(TMEV5_A_RMS)/sqrt(length(TMEV5_A_RMS)), std(TMEV5_B_RMS)/sqrt(length(TMEV5_B_RMS)), std(TMEV5_C_RMS)/sqrt(length(TMEV5_C_RMS))];

TMEV14_A_MEAN = mean(TMEV14_A_RMS);
TMEV14_B_MEAN = mean(TMEV14_B_RMS);
TMEV14_C_MEAN = mean(TMEV14_C_RMS);
TMEV14y = [TMEV14_A_MEAN, TMEV14_B_MEAN, TMEV14_C_MEAN];
TMEV14_sd_vct = [std(TMEV14_A_RMS)/sqrt(length(TMEV14_A_RMS)), std(TMEV14_B_RMS)/sqrt(length(TMEV14_B_RMS)), std(TMEV14_C_RMS)/sqrt(length(TMEV14_C_RMS))];

x = [0, 1, 2];
figure()
h(1) = errorbar(x,PBSy,PBS_sd_vct,'k');
hold on
h(2) = plot(x,PBSy,'kd','MarkerFaceColor','k');
hold on
h(3) = errorbar(x,TMEV2y,TMEV2_sd_vct,'r');
hold on
h(4) = plot(x,TMEV2y,'rs','MarkerFaceColor','r');
hold on
h(5) = errorbar(x,TMEV5y,TMEV5_sd_vct,'g');
hold on
h(6) = plot(x,TMEV5y,'rd','MarkerFaceColor','r');
hold on
h(7) = errorbar(x,TMEV14y,TMEV14_sd_vct,'b');
hold on
h(8) = plot(x,TMEV14y,'r*','MarkerFaceColor','r');

title('RMS: 1-200 frames, 201-timehitR, timehitR-end');
xlabel('Groups'); ylabel('RMS'); grid on;
xlim([-1 4]);
legend(h([2,4,6,8]),'PBS','TMEV 2 DPI','TMEV 5 DPI','TMEV 14 DPI');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Plotting Results %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Paired RMS % %
figure()
xticks([0, 1, 2, 3])
xlim([0 3])
xlabel('Before time 0 / After time 0'); ylabel('RMS'); title('Process PBS RMS');
for k = 1:length(A_RMS)
    hold on
    if A_RMS(k) < B_RMS(k)
        plot([1,2], [A_RMS(k), B_RMS(k)], 'r-d');
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

%% SOMA %%
%%
AL = [];
PBS_A_RMS = [];
PBS_B_RMS = [];
PBS_C_RMS = [];

TMEV2_A_RMS = [];
TMEV2_B_RMS = [];
TMEV2_C_RMS = [];

TMEV5_A_RMS = [];
TMEV5_B_RMS = [];
TMEV5_C_RMS = [];

TMEV14_A_RMS = [];
TMEV14_B_RMS = [];
TMEV14_C_RMS = [];

BL = [];
PBS_R = [];
%     TMEV_R = [];
allRMS{1,1} = 'Trace + ROI';
allRMS{1,2} = 'TMEV/PBS';
allRMS{1,3} = 'DPI';
allRMS{1,4} = 'A(1 to 200) RMS^2';
allRMS{1,5} = 'B(201 to timehitR) RMS^2';
allRMS{1,6} = 'C(timehitR to end) RMS^2';
allRMS{1,7} = 'timehitR Marker/Index';
allRMS{1,8} = 'PBS/TMEV';
filtered_traces = cell(size(soma_data,1),3);
for ii = 1:size(soma_data,1)
    ROInum = soma_data{ii,2};
    ROInum = num2str(ROInum);
    ROI = strcat('  ROI: ', ROInum);
    traceName = soma_data{ii,1};
    traceName = strcat(traceName,ROI);
    DPI = soma_data{ii,5};
    timehitR = double(soma_data{ii,6});
    timehitR_marker = round(timehitR);
    currentTrace = soma_data{ii,3};
    F0 = mean(currentTrace(1:50));
    currentTrace = (currentTrace-F0)./F0;
    if timehitR_marker > 200 && timehitR_marker < length(currentTrace)
        
        rmsA = rms(currentTrace(1:200));
        rmsB = rms(currentTrace(201:timehitR_marker-1));
        rmsC = rms(currentTrace(timehitR_marker:end));
        allRMS{ii+1,1} = traceName;
        allRMS{ii+1,2} = soma_data{ii,4};
        allRMS{ii+1,3} = DPI;
        allRMS{ii+1,4} = rmsA;
        allRMS{ii+1,5} = rmsB;
        allRMS{ii+1,6} = rmsC;
        allRMS{ii+1,7} = timehitR_marker;
        if strcmp(soma_data{ii,4},'PBS')
            % % Find RMS of A and B and C portions % %
            PBS_A_RMS = [PBS_A_RMS; rmsA];
            PBS_B_RMS = [PBS_B_RMS; rmsB];
            PBS_C_RMS = [PBS_C_RMS; rmsC];
            allRMS{ii+1,8} = 'PBS';
        else
            allRMS{ii+1,8} = 'TMEV';
            if DPI == 2
                TMEV2_A_RMS = [TMEV2_A_RMS; rmsA];
                TMEV2_B_RMS = [TMEV2_B_RMS; rmsB];
                TMEV2_C_RMS = [TMEV2_C_RMS; rmsC];
            elseif DPI == 4 || DPI == 5 || DPI == 6
                TMEV5_A_RMS = [TMEV5_A_RMS; rmsA];
                TMEV5_B_RMS = [TMEV5_B_RMS; rmsB];
                TMEV5_C_RMS = [TMEV5_C_RMS; rmsC];
            else
                TMEV14_A_RMS = [TMEV14_A_RMS; rmsA];
                TMEV14_B_RMS = [TMEV14_B_RMS; rmsB];
                TMEV14_C_RMS = [TMEV14_C_RMS; rmsC];
            end            
        end
    else
        rmsA = 'Skipped';
        rmsB = 'Skipped';
        rmsC = 'Skipped';
    end 
end
    
% end
PBS_A_MEAN = mean(PBS_A_RMS);
PBS_B_MEAN = mean(PBS_B_RMS);
PBS_C_MEAN = mean(PBS_C_RMS);
PBSy = [PBS_A_MEAN, PBS_B_MEAN, PBS_C_MEAN];
PBS_sd_vct = [std(PBS_A_RMS)/sqrt(length(PBS_A_RMS)), std(PBS_B_RMS)/sqrt(length(PBS_B_RMS)), std(PBS_C_RMS)/sqrt(length(PBS_C_RMS))];

TMEV2_A_MEAN = mean(TMEV2_A_RMS);
TMEV2_B_MEAN = mean(TMEV2_B_RMS);
TMEV2_C_MEAN = mean(TMEV2_C_RMS);
TMEV2y = [TMEV2_A_MEAN, TMEV2_B_MEAN, TMEV2_C_MEAN];
TMEV2_sd_vct = [std(TMEV2_A_RMS)/sqrt(length(TMEV2_A_RMS)), std(TMEV2_B_RMS)/sqrt(length(TMEV2_B_RMS)), std(TMEV2_C_RMS)/sqrt(length(TMEV2_C_RMS))];

TMEV5_A_MEAN = mean(TMEV5_A_RMS);
TMEV5_B_MEAN = mean(TMEV5_B_RMS);
TMEV5_C_MEAN = mean(TMEV5_C_RMS);
TMEV5y = [TMEV5_A_MEAN, TMEV5_B_MEAN, TMEV5_C_MEAN];
TMEV5_sd_vct = [std(TMEV5_A_RMS)/sqrt(length(TMEV5_A_RMS)), std(TMEV5_B_RMS)/sqrt(length(TMEV5_B_RMS)), std(TMEV5_C_RMS)/sqrt(length(TMEV5_C_RMS))];

TMEV14_A_MEAN = mean(TMEV14_A_RMS);
TMEV14_B_MEAN = mean(TMEV14_B_RMS);
TMEV14_C_MEAN = mean(TMEV14_C_RMS);
TMEV14y = [TMEV14_A_MEAN, TMEV14_B_MEAN, TMEV14_C_MEAN];
TMEV14_sd_vct = [std(TMEV14_A_RMS)/sqrt(length(TMEV14_A_RMS)), std(TMEV14_B_RMS)/sqrt(length(TMEV14_B_RMS)), std(TMEV14_C_RMS)/sqrt(length(TMEV14_C_RMS))];

x = [0, 1, 2];
figure()
h(1) = errorbar(x,PBSy,PBS_sd_vct,'k');
hold on
h(2) = plot(x,PBSy,'kd','MarkerFaceColor','k');
hold on
h(3) = errorbar(x,TMEV2y,TMEV2_sd_vct,'r');
hold on
h(4) = plot(x,TMEV2y,'rs','MarkerFaceColor','r');
hold on
h(5) = errorbar(x,TMEV5y,TMEV5_sd_vct,'g');
hold on
h(6) = plot(x,TMEV5y,'rd','MarkerFaceColor','r');
hold on
h(7) = errorbar(x,TMEV14y,TMEV14_sd_vct,'b');
hold on
h(8) = plot(x,TMEV14y,'r*','MarkerFaceColor','r');

title('RMS: 1-200 frames, 201-timehitR, timehitR-end');
xlabel('Groups'); ylabel('RMS'); grid on;
xlim([-1 4]);
legend(h([2,4,6,8]),'PBS','TMEV 2 DPI','TMEV 5 DPI','TMEV 14 DPI');
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