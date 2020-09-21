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
allPeakDurs = [];
PBS_Alengths = [];
PBS_Blengths = [];
all_RAWauc = [];
all_auc = [];
allIEI = [];
allLocs = [];
before_0_auc = [];
after_0_auc = [];
for k = 1:numel(PBSfields)
    % extract
    tempIEI = PBS.(PBSfields{k}).IEI;
    allIEI = [allIEI;tempIEI];
    
    temp_timehitR = PBS.(PBSfields{k}).timehitR;
    tempLocs = PBS.(PBSfields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
%     if ~isempty(tempLocs)
%     cutout = length(tempLocs);
%     tempLocs(cutout) = [];
%     tempLocs = tempLocs - temp_timehitR;
%     end
    allLocs = [allLocs;tempLocs];
    tempDurs = PBS.(PBSfields{k}).PeakDurs .* 1.04;
    tempDurs = tempDurs(~isnan(tempDurs));
    tempA = PBS.(PBSfields{k}).LENGTH_A;
    tempB = PBS.(PBSfields{k}).LENGTH_B;
    tempRAWAUC = PBS.(PBSfields{k}).RAW_AUC;
    tempRAWAUC = tempRAWAUC(~isnan(tempRAWAUC));
    tempAUC = PBS.(PBSfields{k}).AUC;
    tempAUC = tempAUC(~isnan(tempAUC));
    
    
    all_RAWauc = [all_RAWauc; tempRAWAUC];
    all_auc = [all_auc; tempAUC];
    shiftedLocs = tempLocs - temp_timehitR;
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));
    tempBefore0 = tempRAWAUC(shiftedLocs < 0 );
    tempAfter0 = tempRAWAUC(shiftedLocs >= 0 );
    before_0_auc = [before_0_auc; tempBefore0];
    after_0_auc = [after_0_auc; tempAfter0];
    PBS.(PBSfields{k}).shiftedLocs = shiftedLocs;
    
    allPBSLocs = [allPBSLocs; shiftedLocs];
    allPeakDurs = [allPeakDurs; tempDurs];
    PBS_Alengths = [PBS_Alengths; tempA];
    PBS_Blengths = [PBS_Blengths; tempB];
end
binEdges = [-1600:100:1600];
numBins = numel(binEdges);
figure()
histogram(allPBSLocs,numBins, 'BinEdges', binEdges)
% plot(allPBSLocs,'d')
title('All Shifted PBS Peak Locations');
xlabel('Relative time (s)'); ylabel('Frequency of events');

figure()
plot(allPBSLocs, allPeakDurs, 'dc','MarkerFaceColor','b')
title('Peak Durations in relative time'); xlabel('Relative time (s)'); ylabel('Peak Duration (s)');

figure()
% try out different data visulization of signal arc lengths
groupA = ones(size(PBS_Alengths));
groupB = ones(size(PBS_Blengths))+1;
dotA = [groupA,PBS_Alengths];
dotB = [groupB,PBS_Blengths];
plot(groupA, PBS_Alengths, '.k', 'MarkerSize', 12);
hold on
plot(groupB, PBS_Blengths, 'dk');
legend('Length A', 'Length B');
title('Arc lengths of traces before and after timehitR'); xlabel('Group A / Group B'); ylabel('Arc length (dF/F)');
xticks([0, 1, 2, 3])
xlim([0 3])

figure()
plot(allPBSLocs, all_auc, 'dr' ,'MarkerFaceColor','r')
title('filtered signal auc'); xlabel('relative time(s)');

figure()
plot(allPBSLocs, all_RAWauc, 'dr' ,'MarkerFaceColor','b')
title('raw signal auc'); xlabel('relative time(s)');

% figure()
% plot(allLocs(2:2:end),allIEI, 'dr','MarkerFaceColor','b')
% xlabel('relative time (s)'); ylabel('IEI (s)');
% title('plotting interval at First Peak');

% % boxplots of event auc before time zero and after % %
ALL = [before_0_auc; after_0_auc];
groups = [ones(size(before_0_auc)); ones(size(after_0_auc)) * 2];
figure()
boxplot(ALL,groups,'Whisker',4)
title('PBS Raw Event AUC before relative time zero and after time zero');
ylabel('AUC (dF/F)');

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% TMEV  2 DPI %%%%%
%%%%%%%%%%%%%%%%%%%%%%%
TMEV_2DPI = TMEV.DPI_2;
TMEV2fields = fieldnames(TMEV_2DPI);
allTMEV2Locs = []; %allocate space
allPeakDurs = [];
TMEV2_Alengths = [];
TMEV2_Blengths = [];
all_RAWauc = [];
all_auc = [];
before_0_auc = [];
after_0_auc = [];
for k  = 1:numel(TMEV2fields)
    % extract
    temp_timehitR = TMEV_2DPI.(TMEV2fields{k}).timehitR;
    tempLocs = TMEV_2DPI.(TMEV2fields{k}).PeakLocs;
    tempDurs = TMEV_2DPI.(TMEV2fields{k}).PeakDurs .* 1.04;
    tempDurs = tempDurs(~isnan(tempDurs));
    tempA = TMEV_2DPI.(TMEV2fields{k}).LENGTH_A;
    tempB = TMEV_2DPI.(TMEV2fields{k}).LENGTH_B;
    tempRAWAUC = TMEV_2DPI.(TMEV2fields{k}).RAW_AUC;
    tempRAWAUC = tempRAWAUC(~isnan(tempRAWAUC));
    tempAUC = TMEV_2DPI.(TMEV2fields{k}).AUC;
    tempAUC = tempAUC(~isnan(tempAUC));
    
    all_RAWauc = [all_RAWauc; tempRAWAUC];
    all_auc = [all_auc; tempAUC];
    shiftedLocs = tempLocs - temp_timehitR;
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));
    tempBefore0 = tempRAWAUC(shiftedLocs < 0 );
    tempAfter0 = tempRAWAUC(shiftedLocs >= 0 );
    before_0_auc = [before_0_auc; tempBefore0];
    after_0_auc = [after_0_auc; tempAfter0];
    TMEV_2DPI.(TMEV2fields{k}).shiftedLocs = shiftedLocs;
    
    allTMEV2Locs = [allTMEV2Locs; shiftedLocs];
    allPeakDurs = [allPeakDurs; tempDurs];
    TMEV2_Alengths = [TMEV2_Alengths; tempA];
    TMEV2_Blengths = [TMEV2_Blengths; tempB];
end

binEdges = [-1600:100:1600];
numBins = numel(binEdges);
figure()
histogram(allTMEV2Locs,numBins, 'BinEdges', binEdges)
% plot(allPBSLocs,'d')
title('All Shifted TMEV 2 DPI Peak Locations');
xlabel('Relative time (s)'); ylabel('Frequency of events');

figure()
plot(allTMEV2Locs, allPeakDurs, 'dc','MarkerFaceColor','b')
title('Peak Durations in relative time'); xlabel('Relative time (s)'); ylabel('Peak Duration (s)');

figure()
plot(allTMEV2Locs,all_auc,'d');
title('AUC vs. Relative time: Raw Signal');
xlabel('relative time (s)'); ylabel('AUC');

figure()
plot(allTMEV2Locs,all_RAWauc,'d');
title('AUC vs. Relative time: Raw Signal');
xlabel('relative time (s)'); ylabel('AUC');

% % boxplots of event auc before time zero and after % %
ALL = [before_0_auc; after_0_auc];
groups = [ones(size(before_0_auc)); ones(size(after_0_auc)) * 2];
figure()
boxplot(ALL,groups)
title('TMEV 2 DPI Raw Event AUC before relative time zero and after time zero');
ylabel('AUC (dF/F)');

    
%%%%%%%%%%%%%%%%%%%%%%%
%%%%% TMEV  5 DPI %%%%%
%%%%%%%%%%%%%%%%%%%%%%%
TMEV_5DPI = TMEV.DPI_5;
TMEV5fields = fieldnames(TMEV_5DPI);
allTMEV5Locs = []; %allocate space
allPeakDurs = [];
TMEV5_Alengths = [];
TMEV5_Blengths = [];
all_RAWauc = [];
all_auc = [];
before_0_auc = [];
after_0_auc = [];
for k  = 1:numel(TMEV5fields)
    % extract
    temp_timehitR = TMEV_5DPI.(TMEV5fields{k}).timehitR;
    tempLocs = TMEV_5DPI.(TMEV5fields{k}).PeakLocs;
    tempDurs = TMEV_5DPI.(TMEV5fields{k}).PeakDurs .* 1.04;
    tempDurs = tempDurs(~isnan(tempDurs));
    tempA = TMEV_5DPI.(TMEV5fields{k}).LENGTH_A;
    tempB = TMEV_5DPI.(TMEV5fields{k}).LENGTH_B;
    tempRAWAUC = TMEV_5DPI.(TMEV5fields{k}).RAW_AUC;
    tempRAWAUC = tempRAWAUC(~isnan(tempRAWAUC));
    tempAUC = TMEV_5DPI.(TMEV5fields{k}).AUC;
    tempAUC = tempAUC(~isnan(tempAUC));
    
    all_RAWauc = [all_RAWauc; tempRAWAUC];
    all_auc = [all_auc; tempAUC];
    shiftedLocs = tempLocs - temp_timehitR;
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));
    tempBefore0 = tempRAWAUC(shiftedLocs < 0 );
    tempAfter0 = tempRAWAUC(shiftedLocs >= 0 );
    before_0_auc = [before_0_auc; tempBefore0];
    after_0_auc = [after_0_auc; tempAfter0];
    TMEV_5DPI.(TMEV5fields{k}).shiftedLocs = shiftedLocs;
    
    allTMEV5Locs = [allTMEV5Locs; shiftedLocs];
    allPeakDurs = [allPeakDurs; tempDurs];
    TMEV5_Alengths = [TMEV5_Alengths; tempA];
    TMEV5_Blengths = [TMEV5_Blengths; tempB];
end

binEdges = [-1600:100:1600];
numBins = numel(binEdges);
figure()
histogram(allTMEV5Locs,numBins, 'BinEdges', binEdges)
% plot(allPBSLocs,'d')
title('All Shifted TMEV 5 DPI Peak Locations');
xlabel('Relative time (s)'); ylabel('Frequency of events');

figure()
plot(allTMEV5Locs, allPeakDurs, 'dc','MarkerFaceColor','b')
title('tmev 5 dpi: Peak Durations in relative time'); xlabel('Relative time (s)'); ylabel('Peak Duration (s)');

figure()
plot(allTMEV5Locs,all_auc,'d');
title('tmev 5 dpi: AUC vs. Relative time: Raw Signal');
xlabel('relative time (s)'); ylabel('AUC');

figure()
plot(allTMEV5Locs,all_RAWauc,'d');
title('tmev 5 dpi: AUC vs. Relative time: Raw Signal');
xlabel('relative time (s)'); ylabel('AUC');

% % boxplots of event auc before time zero and after % %
ALL = [before_0_auc; after_0_auc];
groups = [ones(size(before_0_auc)); ones(size(after_0_auc)) * 2];
figure()
boxplot(ALL,groups)
title('TMEV 5 DPI Raw Event AUC before relative time zero and after time zero');
ylabel('AUC (dF/F)');
    %%
    AL = [];
    A_RMS = [];
    B_RMS = [];
    BL = [];
    PBS_R = [];
%     TMEV_R = [];
    filtered_traces = cell(size(process_data,1),3);
    for ii = 1:size(process_data,1)
        ROInum = process_data{ii,2};
        ROInum = num2str(ROInum);
        ROI = strcat('  ROI: ', ROInum);
        traceName = process_data{ii,1};
        traceName = strcat(traceName,ROI);
        if strcmp(process_data{ii,4},'TMEV')
%             if process_data{ii,5} == 2 %%% TMEV DPI condition
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
                rmsA = rms(traceA)^2;
                rmsB = rms(traceB)^2;
                A_RMS = [A_RMS; rmsA];
                B_RMS = [B_RMS; rmsB];

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
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Plotting Results %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
plot(groupA, A_RMS, 'dk');
hold on
plot(groupB, B_RMS, 'dk');
legend('Length A', 'Length B');
title('TMEV 2 DPI: RMS levels of filtered traces before and after timehitR'); xlabel('Group A / Group B'); ylabel('RMS');
xticks([0, 1, 2, 3])
xlim([0 3])
A_AVERAGE = mean(A_RMS);

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


end



