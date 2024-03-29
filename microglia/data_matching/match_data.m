%% use this script to call excel_matching %%
% obtain file and path for excel workbook to pull data from
[fileName, pathName] = uigetfile('*.xlsx');
mainS2p = uigetdir();
velocityCell = ExtractVelocityData(pathName,fileName);
sortedCell = sortrows(velocityCell,[4,5]);
[NEW_SOMA, SOMA_TABLE] = soma_matching(fileName, pathName, mainS2p);
[NEW_PROCESS, PROCESS_TABLE] = pro_matching(fileName, pathName, mainS2p);
NEW_SOMA = sortrows(NEW_SOMA,9);
NEW_PROCESS = sortrows(NEW_PROCESS,9);
somaVideos = cell2mat(NEW_SOMA(:,9));
processVideos = cell2mat(NEW_PROCESS(:,9));
somaCells = cell2mat(NEW_SOMA(:,8));
processCells = cell2mat(NEW_PROCESS(:,8));
[uniqueSomaVideos, somaIdx] = unique(somaVideos);
[uniqueProcessVideos, processIdx] = unique(processVideos);
somaRoiCellFreq = [];
for i = 1:length(somaIdx)
    currentVideo = uniqueSomaVideos(i);
    if i < length(somaIdx)
        cellSlice = sort(somaCells(somaIdx(i):somaIdx(i+1)-1));
    else
        cellSlice = sort(somaCells(somaIdx(i):length(somaCells)));
    end
    cellFreq = [];
    uniqueSomaCells = unique(cellSlice);
    for index = 1:length(uniqueSomaCells)
        cellFreq(index) = sum(cellSlice == uniqueSomaCells(index));
    end
    cellFreq = cellFreq(cellFreq > 0);
    for k = 1:length(cellFreq)
        somaRoiCellFreq = [somaRoiCellFreq; currentVideo, uniqueSomaCells(k), cellFreq(k)];
    end
end
processRoiCellFreq = [];
for i = 1:length(processIdx)
    currentVideo = uniqueProcessVideos(i);
    if i < length(processIdx)
        cellSlice = sort(processCells(processIdx(i):processIdx(i+1)-1));
    else
        cellSlice = sort(processCells(processIdx(i):length(processCells)));
    end
    cellFreq = [];
    uniqueProcessCells = unique(cellSlice);
    for index = 1:length(uniqueProcessCells)
        cellFreq(index) = sum(cellSlice == uniqueProcessCells(index));
    end
    cellFreq = cellFreq(cellFreq > 0);
    for k = 1:length(cellFreq)
        processRoiCellFreq = [processRoiCellFreq; currentVideo, uniqueProcessCells(k), cellFreq(k)];
    end
end
% organize data for two way anova
groups = sortedCell(:,4);
DPI = sortedCell(:,5);
DPI(1:310) = {0};
velocityDistance = [1,2];
velocity = cell2mat(sortedCell(:,3));
distToR = cell2mat(sortedCell(:,7));
velGrouping = ones(size(velocity));
distGrouping = ones(size(distToR));

anovaData = [velocity; distToR];
velDistGrouping = [velGrouping; distGrouping];
groupGrouping = [groups; groups];
dpiGrouping = [DPI; DPI];
[p, tbl, stats, terms] = anovan(anovaData,{velDistGrouping, groupGrouping, cell2mat(dpiGrouping)});
% first query is an integer for excel sheet number
% second is the directory for Fall.mat files
% [sorted_data,excel_table1] = excel_matching(fileName, pathName);
% 
% [sorted_data2,excel_table2] = second_best_matching(fileName, pathName);
% 
% [sorted_data3,excel_table3] = third_best_matching(fileName, pathName);

%% matching paired ROIs



%% shifted
load('soma_data.mat');
load('process_data.mat');
map_ca_rev6(soma_data, process_data);
% unshifted_map_rev3(sorted_data);
%% unshifted

%% find duplicate roi's %%
allImage = SOMA_TABLE.Image;
[val,idx] = unique(allImage);
markForReview = zeros(size(SOMA_TABLE,1),1);
for i = 1:length(idx)-1
    sub_table = SOMA_TABLE(idx(i):idx(i+1)-1,:);
    subReview = markForReview(idx(i):idx(i+1)-1);
    tempROIs = sub_table.Final_calcium_soma;
    [dupeROIs,roiIdx] = unique(tempROIs);
    dupeInd = setdiff(1:size(tempROIs,1),roiIdx);
    subReview(dupeInd) = 1;
    markForReview(idx(i):idx(i+1)-1) = subReview;
end

Tracks = SOMA_TABLE.TrackID;
CellID = SOMA_TABLE.Nr;
timehitR1 = SOMA_TABLE.TimeHitRadius1;
Suite2p_ROIs = SOMA_TABLE.Final_calcium_soma;

Soma_Marked = table(Tracks,CellID,Suite2p_ROIs,markForReview,timehitR1);

%% Match paired ROIs and Event Data
% load event data from matlab folder and paired data
SomaProcessedEvents = load('ProcessedEventsSOMA.mat');
SomaProcessedEvents = SomaProcessedEvents.ProcessedEvents;
SomaProcessedEvents.T2096_15.ROI = 23+1;
ProcessProcessedEvents = load('ProcessedEventsPROCESSES.mat');
ProcessProcessedEvents = ProcessProcessedEvents.ProcessedEvents;
set(0,'defaultfigurecolor',[1 1 1])
[fileName, pathName] = uigetfile('*.xlsx');
mainS2p = uigetdir();
[TRACES, notPaired] = PairedROI(mainS2p,fileName, pathName);

% outer loop is number of fields in paired data
pairedFields = fieldnames(TRACES);
somaFields = fieldnames(SomaProcessedEvents);
processFields = fieldnames(ProcessProcessedEvents);

% count unpaired events
unpairedCount = [];
for k = 1:length(notPaired)
    for x = 1:length(somaFields)
        if strcmp(notPaired{k},somaFields{x})
            sLocs = SomaProcessedEvents.(somaFields{x}).PeakLocs;
            sLocs = sLocs(~isnan(sLocs));
            unpairedCount = unpairedCount + length(sLocs);
        end
    end
    for y = 1:length(processFields)
        if strcmp(notPaired{k},processFields{y})
            pLocs = ProcessProcessedEvents.(processFields{y}).PeakLocs;
            pLocs = pLocs(~isnan(pLocs));
            unpairedCount = unpairedCount + length(pLocs);
        end
    end
end
% Match data to insert events location data
for k = 1:numel(pairedFields)
    subPair = TRACES.(pairedFields{k});
    subFields = fieldnames(subPair);
    for i = 1:numel(subFields)
        currentVideo = subFields{i};
        currentVideo = currentVideo(1:5);
        for j = 1:numel(processFields)
            processVid = processFields{j};
            processVid = processVid(1:5);
            for x = 1:numel(somaFields)
                somaVid = somaFields{x};
                somaVid = somaVid(1:5);
                if ~strcmp(subFields{i},'Pathology') && ~strcmp(subFields{i},'DPI')
                    if isfield(TRACES.(pairedFields{k}).(subFields{i}),'Soma')
                        if strcmp(currentVideo,somaVid)%strcmp(subFields{i},somaFields{x}) %|| strcmp(subFields{i},processFields{j})
                            if TRACES.(pairedFields{k}).(subFields{i}).Soma.ROI == SomaProcessedEvents.(somaFields{x}).ROI-1
                                TRACES.(pairedFields{k}).(subFields{i}).Soma.PeakLocs = SomaProcessedEvents.(somaFields{x}).PeakLocs;
                                TRACES.(pairedFields{k}).(subFields{i}).Soma.PeakMags = SomaProcessedEvents.(somaFields{x}).PeakMags;
                                TRACES.(pairedFields{k}).(subFields{i}).Soma.PeakDurs = SomaProcessedEvents.(somaFields{x}).PeakDurs;
                                TRACES.(pairedFields{k}).(subFields{i}).Soma.RAW_AUC = SomaProcessedEvents.(somaFields{x}).RAW_AUC;
                                TRACES.(pairedFields{k}).Pathology = SomaProcessedEvents.(somaFields{x}).Pathology;
                                TRACES.(pairedFields{k}).DPI = SomaProcessedEvents.(somaFields{x}).DPI;
                            end
                        end
                    end
                end
            end
            
            if ~strcmp(subFields{i},'Pathology') && ~strcmp(subFields{i},'DPI')
                if isfield(TRACES.(pairedFields{k}).(subFields{i}),'Process')
                    if strcmp(currentVideo,processVid)%if strcmp(subFields{i},processFields{j}) %|| strcmp(subFields{i},somaFields{x})
                        if TRACES.(pairedFields{k}).(subFields{i}).Process.ROI == ProcessProcessedEvents.(processFields{j}).ROI-1
                            TRACES.(pairedFields{k}).(subFields{i}).Process.PeakLocs = ProcessProcessedEvents.(processFields{j}).PeakLocs;
                            TRACES.(pairedFields{k}).(subFields{i}).Process.PeakMags = ProcessProcessedEvents.(processFields{j}).PeakMags;
                            TRACES.(pairedFields{k}).(subFields{i}).Process.PeakDurs = ProcessProcessedEvents.(processFields{j}).PeakDurs;
                            TRACES.(pairedFields{k}).(subFields{i}).Process.RAW_AUC = ProcessProcessedEvents.(processFields{j}).RAW_AUC;
                        end
                    end
                end
            end
        end
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % TIME WINDOW % % % % % % % % % % % % % % % % % % % % % %
timeCheck = round(15*0.97); %threshold for matched events % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Pair events based on how close they are together in time
allMins = [];
allPLocs = [];

allPBSMins = [];
allPBSPLocs = [];
allPbsunpaired = [];

allTmev2Mins = [];
allTmev2PLocs = [];
allTmev2unpaired = [];

allTmev5Mins = [];
allTmev5PLocs = [];
allTmev5unpaired = [];

allTmev15Mins = [];
allTmev15PLocs = [];
allTmev15unpaired = [];

somaPbsPeakLocs = [];
somatmev2PeakLocs = [];
somatmev5PeakLocs = [];
somatmev15PeakLocs = [];
somaPbsMags = [];
somatmev2Mags = [];
somatmev5Mags = [];
somatmev15Mags = [];
somaPbsAuc = [];
somatmev2Auc = [];
somatmev5Auc = [];
somatmev15Auc = [];
somaPbsDurs = [];
somatmev2Durs = [];
somatmev5Durs = [];
somatmev15Durs = [];

processPbsPeakLocs = [];
processtmev2PeakLocs = [];
processtmev5PeakLocs = [];
processtmev15PeakLocs = [];
processPbsMags = [];
processtmev2Mags = [];
processtmev5Mags = [];
processtmev15Mags = [];
processPbsAuc = [];
processtmev2Auc = [];
processtmev5Auc = [];
processtmev15Auc = [];
processPbsDurs = [];
processtmev2Durs = [];
processtmev5Durs = [];
processtmev15Durs = [];

PairedPbsAuc = [];
PairedTmev2Auc = [];
PairedTmev5Auc = [];
PairedTmev15Auc = [];

PairedPbsDurs = [];
PairedTmev2Durs = [];
PairedTmev5Durs = [];
PairedTmev15Durs = [];

PairedPbsPhase3 = {};
somaPairedF = [];
processPairedF = [];
somaTraceNames = {};
processTraceNames = {};
for k = 1:numel(pairedFields)
    subPair = TRACES.(pairedFields{k});
    subFields = fieldnames(subPair);
    for i = 1:numel(subFields)
        PairedEvents = [];
        tempPaired = [];
        if isfield(TRACES.(pairedFields{k}).(subFields{i}), 'Soma') && isfield(TRACES.(pairedFields{k}).(subFields{i}), 'Process')
            if isfield(TRACES.(pairedFields{k}).(subFields{i}).Soma,'PeakLocs') && isfield(TRACES.(pairedFields{k}).(subFields{i}).Process,'PeakLocs')
                tempSomaLocs = TRACES.(pairedFields{k}).(subFields{i}).Soma.PeakLocs;
                tempProcessLocs = TRACES.(pairedFields{k}).(subFields{i}).Process.PeakLocs;
                tempSomaDurs = TRACES.(pairedFields{k}).(subFields{i}).Soma.PeakDurs;
                tempProcessDurs = TRACES.(pairedFields{k}).(subFields{i}).Process.PeakDurs;
                tempSomaDurs = tempSomaDurs(~isnan(tempSomaDurs));
                tempProcessDurs = tempProcessDurs(~isnan(tempProcessDurs));
                tempSomaMags = TRACES.(pairedFields{k}).(subFields{i}).Soma.PeakMags;
                tempProcessMags = TRACES.(pairedFields{k}).(subFields{i}).Process.PeakMags;
                tempSomaAuc = TRACES.(pairedFields{k}).(subFields{i}).Soma.RAW_AUC;
                tempProcessAuc = TRACES.(pairedFields{k}).(subFields{i}).Process.RAW_AUC;
                tempSomaLocs = tempSomaLocs(~isnan(tempSomaLocs));
                tempProcessLocs = tempProcessLocs(~isnan(tempProcessLocs));
                tempSomaMags = tempSomaMags(~isnan(tempSomaMags));
                tempProcessMags = tempProcessMags(~isnan(tempProcessMags));
                tempSomaAuc = tempSomaAuc(~isnan(tempSomaAuc));
                tempProcessAuc = tempProcessAuc(~isnan(tempProcessAuc));
                somaRoi = TRACES.(pairedFields{k}).(subFields{i}).Soma.ROI;
                processRoi = TRACES.(pairedFields{k}).(subFields{i}).Process.ROI;
                tempSomaF = TRACES.(pairedFields{k}).(subFields{i}).Soma.traceF;
                tempProcessF = TRACES.(pairedFields{k}).(subFields{i}).Process.traceF;
                for x = 1:length(tempProcessLocs)
                    for y = 1:length(tempSomaLocs)
                        % if tempPaired is < -15, process event happened
                        % outisde of 15s? if > 15 soma event happened
                        % outisde of 15s?
                        tempPaired = [tempPaired; tempProcessLocs(x)-tempSomaLocs(y)];
                        if abs(tempSomaLocs(y) - tempProcessLocs(x)) <= timeCheck && abs(tempSomaLocs(y) - tempProcessLocs(x)) >= 0
                            % query tree for event trace
                            PairedEvents = [PairedEvents; tempSomaLocs(y), tempProcessLocs(x)];
                            if tempSomaLocs(y)+30 <= length(tempSomaF)
                                somaPairedF = [somaPairedF; tempSomaF(tempSomaLocs(y)-30:tempSomaLocs(y)+30)];
                                somaTraceNames = [somaTraceNames; strcat(subFields{i},': ',num2str(somaRoi))];
                            else
%                                 somaPairedF = [somaPairedF; tempSomaF(tempSomaLocs(y)-30:end)];
                            end
                            
                            if tempProcessLocs(x)+30 <= length(tempProcessF)
                                processPairedF = [processPairedF; tempProcessF(tempProcessLocs(x)-30:tempProcessLocs(x)+30)];
                                processTraceNames = [processTraceNames; strcat(subFields{i},': ',num2str(processRoi))];
                            else
%                                 processPairedF = [processPairedF; tempProcessF(tempProcessLocs(x)-30:end)];
                            end 
                        else
                            % nothing
                        end
                    end
                    if ~isempty(tempPaired)
                        closestToZero = abs(tempPaired);
                        [~,I] = min(closestToZero);
                        allMins = [allMins; tempPaired(I), tempProcessLocs(x)];
                        if strcmp(TRACES.(pairedFields{k}).Pathology, 'PBS')
                            allPBSMins = [allPBSMins; tempPaired(I), tempProcessLocs(x)];
                            allPbsunpaired = allPbsunpaired + length(tempPaired(~I));
                            somaPbsPeakLocs = [somaPbsPeakLocs; tempSomaLocs];
                            processPbsPeakLocs = [processPbsPeakLocs; tempProcessLocs];
                            somaPbsMags = [somaPbsMags; tempSomaMags];
                            processPbsMags = [processPbsMags; tempProcessMags];
                            somaPbsAuc = [somaPbsAuc; tempSomaAuc];
                            processPbsAuc = [processPbsAuc; tempProcessAuc];
                            somaPbsDurs = [somaPbsDurs; tempSomaDurs];
                            processPbsDurs = [processPbsDurs; tempProcessDurs];
                            
                                for y = 1:length(tempSomaLocs)
                                    tempPaired = [tempPaired; tempProcessLocs(x)-tempSomaLocs(y)];
                                    if abs(tempSomaLocs(y) - tempProcessLocs(x)) <= timeCheck && tempSomaLocs(y) - tempProcessLocs(x) > 0
                                        %%% sliceIdx = max(somaPbsDurs(y),processPbsDurs(x));
                                        %%% 
                                        PairedPbsAuc = [PairedPbsAuc; tempSomaAuc(y), tempProcessAuc(x)];
                                        PairedPbsDurs = [PairedPbsDurs; tempSomaDurs(y), tempProcessDurs(x)];
                                    end
                                end
                        else %TMEV
                            if TRACES.(pairedFields{k}).DPI == 2
                                allTmev2Mins = [allTmev2Mins; tempPaired(I), tempProcessLocs(x)];
                                allTmev2unpaired = allTmev2unpaired + length(tempPaired(~I));
                                somatmev2PeakLocs = [somatmev2PeakLocs; tempSomaLocs];
                                processtmev2PeakLocs = [processtmev2PeakLocs; tempProcessLocs];
                                somatmev2Mags = [somatmev2Mags; tempSomaMags];
                                processtmev2Mags = [processtmev2Mags; tempProcessMags];
                                somatmev2Auc = [somatmev2Auc; tempSomaAuc];
                                processtmev2Auc = [processtmev2Auc; tempProcessAuc];
                                for y = 1:length(tempSomaLocs)
                                    tempPaired = [tempPaired; tempProcessLocs(x)-tempSomaLocs(y)];
                                    if abs(tempSomaLocs(y) - tempProcessLocs(x)) <= timeCheck && tempSomaLocs(y) - tempProcessLocs(x) > 0
                                        PairedTmev2Auc = [PairedTmev2Auc; tempSomaAuc(y), tempProcessAuc(x)];
                                        PairedTmev2Durs = [PairedTmev2Durs; tempSomaDurs(y), tempProcessDurs(x)];
                                    end
                                end
                            elseif TRACES.(pairedFields{k}).DPI >= 14
                                allTmev15Mins = [allTmev15Mins; tempPaired(I), tempProcessLocs(x)];
                                allTmev15unpaired = allTmev15unpaired + length(tempPaired(~I));
                                somatmev15PeakLocs = [somatmev15PeakLocs; tempSomaLocs];
                                processtmev15PeakLocs = [processtmev15PeakLocs; tempProcessLocs];
                                somatmev15Mags = [somatmev15Mags; tempSomaMags];
                                processtmev15Mags = [processtmev15Mags; tempProcessMags];
                                somatmev15Auc = [somatmev15Auc; tempSomaAuc];
                                processtmev15Auc = [processtmev15Auc; tempProcessAuc];
                                for y = 1:length(tempSomaLocs)
                                    tempPaired = [tempPaired; tempProcessLocs(x)-tempSomaLocs(y)];
                                    if abs(tempSomaLocs(y) - tempProcessLocs(x)) <= timeCheck && tempSomaLocs(y) - tempProcessLocs(x) >= 0
                                        PairedTmev15Auc = [PairedTmev15Auc; tempSomaAuc(y), tempProcessAuc(x)];
                                        PairedTmev15Durs = [PairedTmev15Durs; tempSomaDurs(y), tempProcessDurs(x)];
                                    end
                                end
                            else
                                allTmev5Mins = [allTmev5Mins; tempPaired(I), tempProcessLocs(x)];
                                allTmev5unpaired = allTmev5unpaired + length(tempPaired(~I));
                                somatmev5PeakLocs = [somatmev5PeakLocs; tempSomaLocs];
                                processtmev5PeakLocs = [processtmev5PeakLocs; tempProcessLocs];
                                somatmev5Mags = [somatmev5Mags; tempSomaMags];
                                processtmev5Mags = [processtmev5Mags; tempProcessMags];
                                somatmev5Auc = [somatmev5Auc; tempSomaAuc];
                                processtmev5Auc = [processtmev5Auc; tempProcessAuc];
                                for y = 1:length(tempSomaLocs)
                                    tempPaired = [tempPaired; tempProcessLocs(x)-tempSomaLocs(y)];
                                    if abs(tempSomaLocs(y) - tempProcessLocs(x)) <= timeCheck && tempSomaLocs(y) - tempProcessLocs(x) > 0
                                        PairedTmev5Auc = [PairedTmev5Auc; tempSomaAuc(y), tempProcessAuc(x)];
                                        PairedTmev5Durs = [PairedTmev5Durs; tempSomaDurs(y), tempProcessDurs(x)];
                                    end
                                end
                            end
                        end
                    end
                end
                TRACES.(pairedFields{k}).(subFields{i}).PairedEvents = PairedEvents;
                allPLocs = [allPLocs; tempProcessLocs];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%all paired events plotted%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracesFields = fieldnames(TRACES);
for i = 1:numel(tracesFields)
    trackFields = fieldnames(TRACES.(tracesFields{i}));
    for j = 1:numel(trackFields)
        if isfield(TRACES.(tracesFields{i}).(trackFields{j}),'Soma') && isfield(TRACES.(tracesFields{i}).(trackFields{j}),'Process')
            figTitle = trackFields{j};
            somaF = TRACES.(tracesFields{i}).(trackFields{j}).Soma.traceF;
            somaRoi = strcat('ROI: ',num2str(TRACES.(tracesFields{i}).(trackFields{j}).Soma.ROI));
            processRoi = strcat('ROI: ',num2str(TRACES.(tracesFields{i}).(trackFields{j}).Process.ROI));
            processF = TRACES.(tracesFields{i}).(trackFields{j}).Process.traceF;
            somaLocs = removeNaN(TRACES.(tracesFields{i}).(trackFields{j}).Soma.PeakLocs);
            processLocs = removeNaN(TRACES.(tracesFields{i}).(trackFields{j}).Process.PeakLocs);
            paired = TRACES.(tracesFields{i}).(trackFields{j}).PairedEvents;
            figure();
            hold on
            s = plot(somaF);
            plot(somaLocs,somaF(somaLocs),'o')
            p = plot(processF);
            plot(processLocs,processF(processLocs),'o')
            title(figTitle,'Interpreter','none')
            if ~isempty(paired)
                plot(paired,somaF(paired(:,1)),'o','MarkerFaceColor',[0,0,0])
                plot(paired,processF(paired(:,2)),'o','MarkerFaceColor',[0,0,0])
            end
            legend([s,p],somaRoi,processRoi);
        end
    end
    
end


% % % compare paired auc ratios to non paired
pbsSomaPairedAuc = PairedPbsAuc(:,1);
pbsProcessPairedAuc = PairedPbsAuc(:,2);
pairedPbsRatio = pbsSomaPairedAuc./pbsProcessPairedAuc;
pbsSomaRandomAuc = randsample(somaPbsAuc, 45);
pbsProcessRandomAuc = randsample(processPbsAuc, 45);
pbsSomaRandomDurs = randsample(somaPbsDurs, 45);
pbsProcessRandomDurs = randsample(processPbsDurs, 45);
randomPbsRatio = pbsSomaRandomAuc./pbsProcessRandomAuc;

pbsLocs = allPBSMins(:,1);
windowPairedPbs = [(pbsLocs <= 14 & pbsLocs >= -14), (pbsLocs <= 14 & pbsLocs >= -14), (pbsLocs <= 14 & pbsLocs >= -14)];
unpairedPBS = (pbsLocs < -15 | pbsLocs > 15);
pbsBefore = (pbsLocs < -1 & pbsLocs >= -15);
pbsDuring = (pbsLocs >= -1 & pbsLocs <= 1);
pbsAfter = (pbsLocs <= 15 & pbsLocs > 1);

tmev2Locs = allTmev2Mins(:,1);
windowPairedtmev2 = [(tmev2Locs <= 14 & tmev2Locs >= -14), (tmev2Locs <= 14 & tmev2Locs >= -14), (tmev2Locs <= 14 & tmev2Locs >= -14)];
unpairedTmev2 = (tmev2Locs < -15 | tmev2Locs > 15);
tmev2Before = (tmev2Locs < -1 & tmev2Locs >= -15);
tmev2During = (tmev2Locs >= -1 & tmev2Locs <= 1);
tmev2After = (tmev2Locs <= 15 & tmev2Locs > 0);


tmev5Locs = allTmev5Mins(:,1);
windowPairedtmev5 = [(tmev5Locs <= 14 & tmev5Locs >= -14),(tmev5Locs <= 14 & tmev5Locs >= -14),(tmev5Locs <= 14 & tmev5Locs >= -14)];
unpairedTmev5 = (tmev5Locs < -15 | tmev5Locs > 15);
tmev5Before = (tmev5Locs < -1 & tmev5Locs >= -15);
tmev5During = (tmev5Locs >= -1 & tmev5Locs <= 1);
tmev5After = (tmev5Locs <= 15 & tmev5Locs > 1);

tmev15Locs = allTmev15Mins(:,1);
windowPairedtmev15 = [(tmev15Locs <= 14 & tmev15Locs >= -14),(tmev15Locs <= 14 & tmev15Locs >= -14),(tmev15Locs <= 14 & tmev15Locs >= -14)];
unpairedTmev15 = (tmev15Locs < -15 | tmev15Locs > 15);
tmev15Before = (tmev15Locs < -1 & tmev15Locs >= -15);
tmev15During = (tmev15Locs >= -1 & tmev15Locs <= 1);
tmev15After = (tmev15Locs <= 15 & tmev15Locs > 1);

allPLocs = allPLocs(~isnan(allPLocs));
% allMins = allMins(~isnan(allMins));
%% all events
figure()
plot(allMins(:,1)./0.97,allMins(:,2)./0.97,'ok','MarkerSize',4,'Color',[84/256, 109/256, 238/256]);
ylimit = ylim;
ylim([0, ylimit(2)]);
yticks(ylimit(2));
yTickL = {ChangeTextColor(num2str(ylimit(2)), [0 0 0])};
xlimit = xlim;
xlim([-xlimit(2), xlimit(2)]);
xticks([-xlimit(2), 0, xlimit(2)]);
xTickL = {ChangeTextColor(num2str(-xlimit(2)), [0 0 0]), ...
    ChangeTextColor(num2str(0), [0 0 0]), ...
    ChangeTextColor(num2str(xlimit(2)), [0 0 0])};
hold on
line([-15, -15], [0, ylimit(2)],'Color','r')
line([15, 15], [0, ylimit(2)],'Color','r')
hold off
box off
set(gca,'YTickLabel',yTickL','XTickLabel',xTickL');
set(gca,'FontSize',10,'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 2 2.5]);
set(gcf, 'Units','inches','position',[4 4 2 2.5]);
print('all_events.png', '-r900','-dpng')
diffs = allMins(:,1)./0.97;
locs = allMins(:,2)./0.97;
noExtremes = [diffs(diffs < 15 & diffs > -15), locs(diffs < 15 & diffs > -15)];
h = figure();
hold on
% plot(noExtremes(:,1),noExtremes(:,2),'ok','MarkerSize',4,'Color',[84/256, 109/256, 238/256]);
plot(allMins(:,1)./0.97,allMins(:,2)./0.97,'ok','MarkerSize',4,'Color',[84/256, 109/256, 238/256]);
ylimit = ylim;
ylim([0, ylimit(2)]);
yticks(ylimit(2));
xlimit = xlim;
xlim([-30, 30]);
xticks([-30, 0, 30]);
xTickL = {ChangeTextColor(num2str(-30), [0 0 0]), ...
    ChangeTextColor(num2str(0), [0 0 0]), ...
    ChangeTextColor(num2str(30), [0 0 0])};
line([-15, -15], [0, ylimit(2)],'Color','r')
line([15, 15], [0, ylimit(2)],'Color','r')
box off
set(gca,'YTickLabel',yTickL','XTickLabel',xTickL');
set(gca, 'TickLabelInterpreter', 'tex');
set(gca,'FontSize',10,'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 2 4]);
set(gcf, 'Units','inches','position',[4 4 2 4]);
print('zoomed in all events.png', '-r900','-dpng')

figure()
hold on
plot(TRACES.T2096.T2096_10.Soma.traceF)
plot(TRACES.T2096.T2096_10.Process.traceF)
legend('soma','process')
%% PBS Plot

pbsFit = fitdist(allPBSMins(:,1)./0.97,'kernel');
index = linspace(min(allPBSMins(:,1)./0.97), max(allPBSMins(:,1)./0.97), 1000);
figure()
plot(index, pdf(pbsFit, index))
figure()
plot(allPBSMins(:,1)./0.97,allPBSMins(:,2)./0.97,'ok','MarkerSize',4,'Color',[84/256, 109/256, 238/256]);
ylimit = ylim;
ylim([0, ylimit(2)]);
yticks(ylimit(2));
line([-15, -15], [0, ylimit(2)],'Color','r')
line([15, 15], [0, ylimit(2)],'Color','r')
yTickL = {ChangeTextColor(num2str(ylimit(2)), [0 0 0])};
xlimit = xlim;
xlim([-xlimit(2), xlimit(2)]);
% xticks([-xlimit(2):5:xlimit(2)]);
xTickL = {ChangeTextColor(num2str(-xlimit(2)), [0 0 0]), ...
    ChangeTextColor(num2str(-10), [0 0 0]), ...
    ChangeTextColor(num2str(-5), [0 0 0]), ...
    ChangeTextColor(num2str(0), [0 0 0]), ...
    ChangeTextColor(num2str(5), [0 0 0]), ...
    ChangeTextColor(num2str(10), [0 0 0]), ...
    ChangeTextColor(num2str(xlimit(2)), [0 0 0])};
box off
% set(gca,'YTickLabel',yTickL','XTickLabel',xTickL');
set(gca,'FontSize',10,'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 2 4]);
set(gcf, 'Units','inches','position',[4 4 2 4]);
print('pbs_events.png', '-r900','-dpng')
diffs = allPBSMins(:,1)./0.97;
locs = allPBSMins(:,2)./0.97;
noExtremes = [diffs(diffs < 15 & diffs > -15), locs(diffs < 15 & diffs > -15)];
h = figure();
hold on
plot(noExtremes(:,1),noExtremes(:,2),'ok','MarkerSize',4,'Color',[84/256, 109/256, 238/256]);
ylimit = ylim;
ylim([0, ylimit(2)]);
yticks(ylimit(2));
xlimit = xlim;
xlim([-xlimit(2), xlimit(2)]);
xticks([-xlimit(2):5:xlimit(2)]);
xTickL = {ChangeTextColor(num2str(-xlimit(2)), [0 0 0]), ...
    ChangeTextColor(num2str(-10), [0 0 0]), ...
    ChangeTextColor(num2str(-5), [0 0 0]), ...
    ChangeTextColor(num2str(0), [0 0 0]), ...
    ChangeTextColor(num2str(5), [0 0 0]), ...
    ChangeTextColor(num2str(10), [0 0 0]), ...
    ChangeTextColor(num2str(xlimit(2)), [0 0 0])};
line([-5, -5], [0, ylimit(2)],'Color','r')
line([5, 5], [0, ylimit(2)],'Color','r')
box off
% set(gca,'YTickLabel',yTickL','XTickLabel',xTickL');
set(gca, 'TickLabelInterpreter', 'tex');
set(gca,'FontSize',10,'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 2 2.5]);
set(gcf, 'Units','inches','position',[4 4 2 2.5]);
print('cutoffPBS_events.png', '-r900','-dpng')

%% TMEV 2 dpi plot
figure()
plot(allTmev2Mins(:,1)./0.97,allTmev2Mins(:,2)./0.97,'ok','MarkerSize',4,'Color',[84/256, 109/256, 238/256]);
ylimit = ylim;
ylim([0, ylimit(2)]);
yticks(ylimit(2));
yTickL = {ChangeTextColor(num2str(ylimit(2)), [0 0 0])};
xlimit = xlim;
xlim([-xlimit(2), xlimit(2)]);
xticks([-xlimit(2), 0, xlimit(2)]);
xTickL = {ChangeTextColor(num2str(-xlimit(2)), [0 0 0]), ...
    ChangeTextColor(num2str(0), [0 0 0]), ...
    ChangeTextColor(num2str(xlimit(2)), [0 0 0])};
box off
set(gca,'YTickLabel',yTickL','XTickLabel',xTickL');
set(gca,'FontSize',10,'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 2 4]);
set(gcf, 'Units','inches','position',[4 4 2 4]);
print('tmev2_events.png', '-r900','-dpng')
diffs = allTmev2Mins(:,1)./0.97;
locs = allTmev2Mins(:,2)./0.97;
noExtremes = [diffs(diffs < 15 & diffs > -15), locs(diffs < 15 & diffs > -15)];
h = figure();
hold on
plot(noExtremes(:,1),noExtremes(:,2),'ok','MarkerSize',4,'Color',[84/256, 109/256, 238/256]);
ylimit = ylim;
ylim([0, ylimit(2)]);
yticks(ylimit(2));
xlimit = xlim;
xlim([-xlimit(2), xlimit(2)]);
xlim([-15, 15]);
xlimit = xlim;
xticks([-xlimit(2):5:xlimit(2)]);
xTickL = {ChangeTextColor(num2str(-xlimit(2)), [0 0 0]), ...
    ChangeTextColor(num2str(-10), [0 0 0]), ...
    ChangeTextColor(num2str(-5), [0 0 0]), ...
    ChangeTextColor(num2str(0), [0 0 0]), ...
    ChangeTextColor(num2str(5), [0 0 0]), ...
    ChangeTextColor(num2str(10), [0 0 0]), ...
    ChangeTextColor(num2str(xlimit(2)), [0 0 0])};
line([-5, -5], [0, ylimit(2)],'Color','r')
line([5, 5], [0, ylimit(2)],'Color','r')
box off
set(gca,'YTickLabel',yTickL','XTickLabel',xTickL');
set(gca, 'TickLabelInterpreter', 'tex');
set(gca,'FontSize',10,'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 2 2.5]);
set(gcf, 'Units','inches','position',[4 4 2 2.5]);
print('cutoffTmev2_events.png', '-r900','-dpng')

%% TMEV 5 dpi
figure()
plot(allTmev5Mins(:,1)./0.97,allTmev5Mins(:,2)./0.97,'ok','MarkerSize',4,'Color',[84/256, 109/256, 238/256]);
ylimit = ylim;
ylim([0, ylimit(2)]);
yticks(ylimit(2));
yTickL = {ChangeTextColor(num2str(ylimit(2)), [0 0 0])};
xlimit = xlim;
xlim([-xlimit(2), xlimit(2)]);
xticks([-xlimit(2), 0, xlimit(2)]);
xTickL = {ChangeTextColor(num2str(-xlimit(2)), [0 0 0]), ...
    ChangeTextColor(num2str(0), [0 0 0]), ...
    ChangeTextColor(num2str(xlimit(2)), [0 0 0])};
box off
set(gca,'YTickLabel',yTickL','XTickLabel',xTickL');
set(gca,'FontSize',10,'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 2 4]);
set(gcf, 'Units','inches','position',[4 4 2 4]);
print('tmev5_events.png', '-r900','-dpng')
diffs = allTmev5Mins(:,1)./0.97;
locs = allTmev5Mins(:,2)./0.97;
noExtremes = [diffs(diffs < 15 & diffs > -15), locs(diffs < 15 & diffs > -15)];
h = figure();
hold on
plot(noExtremes(:,1),noExtremes(:,2),'ok','MarkerSize',4,'Color',[84/256, 109/256, 238/256]);
ylimit = ylim;
ylim([0, ylimit(2)]);
yticks(ylimit(2));
xlimit = xlim;
xlim([-15, 15]);
xlimit = xlim;
xticks([-xlimit(2):5:xlimit(2)]);
xTickL = {ChangeTextColor(num2str(-xlimit(2)), [0 0 0]), ...
    ChangeTextColor(num2str(-10), [0 0 0]), ...
    ChangeTextColor(num2str(-5), [0 0 0]), ...
    ChangeTextColor(num2str(0), [0 0 0]), ...
    ChangeTextColor(num2str(5), [0 0 0]), ...
    ChangeTextColor(num2str(10), [0 0 0]), ...
    ChangeTextColor(num2str(xlimit(2)), [0 0 0])};
line([-5, -5], [0, ylimit(2)],'Color','r')
line([5, 5], [0, ylimit(2)],'Color','r')
box off
set(gca,'YTickLabel',yTickL','XTickLabel',xTickL');
set(gca, 'TickLabelInterpreter', 'tex');
set(gca,'FontSize',10,'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 2 2.5]);
set(gcf, 'Units','inches','position',[4 4 2 2.5]);
print('cutoffTmev5_events.png', '-r900','-dpng')

%% TMEV 15 dpi
figure()
plot(allTmev15Mins(:,1)./0.97,allTmev15Mins(:,2)./0.97,'ok','MarkerSize',4,'Color',[84/256, 109/256, 238/256]);
ylimit = ylim;
ylim([0, ylimit(2)]);
yticks(ylimit(2));
yTickL = {ChangeTextColor(num2str(ylimit(2)), [0 0 0])};
xlimit = xlim;
xlim([-xlimit(2), xlimit(2)]);
xticks([-xlimit(2), 0, xlimit(2)]);
xTickL = {ChangeTextColor(num2str(-xlimit(2)), [0 0 0]), ...
    ChangeTextColor(num2str(0), [0 0 0]), ...
    ChangeTextColor(num2str(xlimit(2)), [0 0 0])};
box off
set(gca,'YTickLabel',yTickL','XTickLabel',xTickL');
set(gca,'FontSize',10,'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 2 4]);
set(gcf, 'Units','inches','position',[4 4 2 4]);
print('tmev15_events.png', '-r900','-dpng')
diffs = allTmev15Mins(:,1)./0.97;
locs = allTmev15Mins(:,2)./0.97;
noExtremes = [diffs(diffs < 15 & diffs > -15), locs(diffs < 15 & diffs > -15)];
h = figure();
hold on
plot(noExtremes(:,1),noExtremes(:,2),'ok','MarkerSize',4,'Color',[84/256, 109/256, 238/256]);
ylimit = ylim;
ylim([0, ylimit(2)]);
yticks(ylimit(2));
xlimit = xlim;
xlim([-15, 15]);
xlimit = xlim;
xticks([-xlimit(2):5:xlimit(2)]);
xTickL = {ChangeTextColor(num2str(-xlimit(2)), [0 0 0]), ...
    ChangeTextColor(num2str(-10), [0 0 0]), ...
    ChangeTextColor(num2str(-5), [0 0 0]), ...
    ChangeTextColor(num2str(0), [0 0 0]), ...
    ChangeTextColor(num2str(5), [0 0 0]), ...
    ChangeTextColor(num2str(10), [0 0 0]), ...
    ChangeTextColor(num2str(xlimit(2)), [0 0 0])};
line([-5, -5], [0, ylimit(2)],'Color','r')
line([5, 5], [0, ylimit(2)],'Color','r')
box off
set(gca,'YTickLabel',yTickL','XTickLabel',xTickL');
set(gca, 'TickLabelInterpreter', 'tex');
set(gca,'FontSize',10,'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 2 2.5]);
set(gcf, 'Units','inches','position',[4 4 2 2.5]);
print('cutoffTmev15_events.png', '-r900','-dpng')
%%
% Plot traces; first column is soma second column is process
OrderEvents = {'Trace','Group','DPI','SomaROI','ProcesROI','SameTime','SomaFirst','ProcessFirst'};
for k = 1:numel(pairedFields)
    subPair = TRACES.(pairedFields{k});
    subFields = fieldnames(subPair);
    for i = 1:numel(subFields)
        if isfield(TRACES.(pairedFields{k}).(subFields{i}),'PairedEvents')            
            somaTrace = TRACES.(pairedFields{k}).(subFields{i}).Soma.traceF;
            somaTrace = (somaTrace-mean(somaTrace(1:50)))./mean(somaTrace(1:50));
            processTrace = TRACES.(pairedFields{k}).(subFields{i}).Process.traceF;
            processTrace = (processTrace-mean(processTrace(1:50)))./mean(processTrace(1:50));
            pairedLocs = TRACES.(pairedFields{k}).(subFields{i}).PairedEvents;
            if ~isempty(pairedLocs)                               
                somaLocs = pairedLocs(:,1);
                somaLocs = somaLocs(~isnan(somaLocs));
                processLocs = pairedLocs(:,2);
                processLocs = processLocs(~isnan(processLocs));
                eventOrder = somaLocs - processLocs;
                tempOrderEvents = cell(length(eventOrder),7);
                for x = 1:length(eventOrder)
                    tempOrderEvents{x,1} = subFields{i};
                    tempOrderEvents{x,2} = TRACES.(pairedFields{k}).Pathology;
                    tempOrderEvents{x,3} = TRACES.(pairedFields{k}).DPI;
                    tempOrderEvents{x,4} = TRACES.(pairedFields{k}).(subFields{i}).Soma.ROI;
                    tempOrderEvents{x,5} = TRACES.(pairedFields{k}).(subFields{i}).Process.ROI;
                    if eventOrder(x) > 0
                        tempOrderEvents{x,6} = 0;
                        tempOrderEvents{x,7} = 1;
                        tempOrderEvents{x,8} = 0;
                    elseif eventOrder(x) < 0
                        tempOrderEvents{x,6} = 0;
                        tempOrderEvents{x,7} = 0;
                        tempOrderEvents{x,8} = 1;
                    else
                        tempOrderEvents{x,6} = 1;
                        tempOrderEvents{x,7} = 0;
                        tempOrderEvents{x,8} = 0;
                    end
                end
                group = strcat(TRACES.(subFields{i}).Pathology, num2str(TRACES.(subFields{i}).DPI));
                OrderEvents = [OrderEvents; tempOrderEvents];
                time = (1:length(somaTrace))./0.97;
                figure()
                hold on
                p1 = plot(time,somaTrace);
                p2 = plot(somaLocs./0.97,somaTrace(somaLocs),'v','MarkerFaceColor','b');
                p3 = plot(time,processTrace);
                p4 = plot(processLocs./0.97,processTrace(processLocs),'v','MarkerFaceColor','r','MarkerSize',12);
                title(strcat(subFields{i},group),'interpreter','none');
                legend([p1(1), p3(1)],'Soma','Process');
            end
        end
    end
end

pbsSoma = TRACES.T2096.T2096_9.Soma.traceF;
pbsSoma = (pbsSoma - mean(pbsSoma(1:50)))./mean(pbsSoma(1:50));
pbsLocs = TRACES.T2096.T2096_9.PairedEvents;
pbsProcess = TRACES.T2096.T2096_9.Process.traceF;
pbsProcess = (pbsProcess - mean(pbsProcess(1:50)))./mean(pbsProcess(1:50));
pbsSomaLocs = pbsLocs(:,1);
pbsProcessLocs = pbsLocs(:,2);
time = (1:1740)./0.97;
figure()
hold on
plot(time,pbsProcess,'Color',[255/256, 200/256, 0],'LineWidth',1.5)
plot(time, pbsSoma,'Color',[67/256, 94/256, 237/256],'LineWidth',1.5)
plot(pbsSomaLocs./0.97, pbsSoma(pbsSomaLocs)*1.05,'vk','MarkerFaceColor',[0 0 0])
plot(pbsProcessLocs./0.97, pbsProcess(pbsProcessLocs).*1.05,'vk','MarkerFaceColor',[0 0 0])
set(gcf, 'Units','inches','position',[4 4 5 3]);
set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [1 1 1],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 5 3]);
print('example_pairs.png', '-r900','-dpng')