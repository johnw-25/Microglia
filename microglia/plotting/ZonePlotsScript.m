%% Script summary
% This script will be employ custom data organizing and plotting functions
% for the sole purpose of visualizing data for the Wilcox SMBB lab - John W

%% Data organization of Events Data
% Load in data from MATLAB folder
SOMA_sortedData = load('SortedTracesSOMA_1012.mat');
SOMA_sortedData = SOMA_sortedData.sortedData;
PROCESSES_sortedData = load('SortedTracesPROCESSES_1012.mat');
PROCESSES_sortedData = PROCESSES_sortedData.sortedData;
SOMA_DATA = load('soma_data_04-12-21.mat');
SOMA_DATA = SOMA_DATA.NEW_SOMA;
PROCESSES_DATA = load('process_data_04-12-21.mat');
PROCESSES_DATA = PROCESSES_DATA.NEW_PROCESS;

PBS_logicals = contains(SOMA_DATA(:,4),'PBS'); 
PBS_SOMA_DATA = SOMA_DATA(PBS_logicals,:);
PBS_logicals = contains(PROCESSES_DATA(:,4),'PBS'); 
PBS_PROCESSES_DATA = PROCESSES_DATA(PBS_logicals,:);


TMEV_logicals = contains(SOMA_DATA(:,4),'TMEV'); 
TMEV_extracted = SOMA_DATA(TMEV_logicals,:);
SOMA_TMEV2_DATA = TMEV_extracted(1:19,:); % soma
SOMA_TMEV5_DATA = TMEV_extracted(20:40,:); % soma
SOMA_TMEV14_DATA = TMEV_extracted(41:end,:); % soma
TMEV_logicals = contains(PROCESSES_DATA(:,4),'TMEV');
TMEV_extracted = PROCESSES_DATA(TMEV_logicals,:);
PROCESSES_TMEV2_DATA = TMEV_extracted(1:64,:); % process
PROCESSES_TMEV5_DATA = TMEV_extracted(65:80,:); % process
PROCESSES_TMEV14_DATA = TMEV_extracted(81:end,:); % process
% Organize data for plotting
% Soma first
SOMA_PBS = SOMA_sortedData.PBS;
SOMA_TMEV = SOMA_sortedData.TMEV;
SOMA_TMEV2 = SOMA_TMEV.DPI_2;
SOMA_TMEV5 = SOMA_TMEV.DPI_5;
SOMA_TMEV14 = SOMA_TMEV.DPI_14;

% Soma PBS:
%   - T2042_12, ROI: 11, Loc 79, Phase 1
%   - T2313_5, ROI: 22, LOC 68, PHASE 1
T2042_12F = SOMA_PBS.T2042_12.trace; % valleys: 46, 104
T2042_12F = (T2042_12F-mean(T2042_12F(1:50)))/mean(T2042_12F(1:50));
durs = SOMA_PBS.T2042_12.PeakDurs;
durs(1) = 104-46;
SOMA_PBS.T2042_12.PeakDurs = durs;
auc = SOMA_PBS.T2042_12.RAW_AUC;
auc(1) = sum(abs(T2042_12F(46:104)));
SOMA_PBS.T2042_12.RAW_AUC = auc;

% Recalculate event features that are 0 when they shouldn't be
T2313_5F = SOMA_PBS.T2313_5.trace; % valleys: 58,141
T2313_5F = (T2313_5F-mean(T2313_5F(1:50)))/mean(T2313_5F(1:50));
durs = SOMA_PBS.T2313_5.PeakDurs;
durs(1) = 141-58;
SOMA_PBS.T2313_5.PeakDurs = durs;
auc = SOMA_PBS.T2313_5.RAW_AUC;
auc(1) = sum(abs(T2313_5F(58:141)));
SOMA_PBS.T2313_5.RAW_AUC = auc;

% Processes next
PROCESSES_PBS = PROCESSES_sortedData.PBS;
PROCESSES_TMEV = PROCESSES_sortedData.TMEV;
PROCESSES_TMEV2 = PROCESSES_TMEV.DPI_2;
PROCESSES_TMEV5 = PROCESSES_TMEV.DPI_5;
PROCESSES_TMEV14 = PROCESSES_TMEV.DPI_14;


% Process PBS:
%   - T2112_11, ROI: 9, Loc 68, Phase 1
T2112_11F = PROCESSES_PBS.T2112_11.trace; % valleys: 45,95
T2112_11F = (T2112_11F-mean(T2112_11F(1:50)))/mean(T2112_11F(1:50));
durs = PROCESSES_PBS.T2112_11.PeakDurs;
durs(1) = 95-45;
PROCESSES_PBS.T2112_11.PeakDurs = durs;
auc = PROCESSES_PBS.T2112_11.RAW_AUC;
auc(1) = sum(abs(T2112_11F(45:95)));
PROCESSES_PBS.T2112_11.RAW_AUC = auc;
% Process Tmev 2:
%   - T2052_1, ROI: 1, Loc 1501, Phase 3
T2052_1F = PROCESSES_TMEV2.T2052_1.trace; % valleys: 1476,1558
T2052_1F = (T2052_1F-mean(T2052_1F(1:50)))/mean(T2052_1F(1:50));
durs = PROCESSES_TMEV2.T2052_1.PeakDurs;
durs(end) = 1558-1476;
PROCESSES_TMEV2.T2052_1.PeakDurs = durs;
auc = PROCESSES_TMEV2.T2052_1.RAW_AUC;
auc(end) = sum(abs(T2052_1F(1476:1558)));
PROCESSES_TMEV2.T2052_1.RAW_AUC = auc;

% Processes Tmev 14:
%   - T2129_7, ROi: 21, Loc 74, Phase 1
T2129_7F = PROCESSES_TMEV14.T2129_7.trace; % valleys: 56,125
T2129_7F = (T2129_7F-mean(T2129_7F(1:50)))/mean(T2129_7F(1:50));
durs = PROCESSES_TMEV14.T2129_7.PeakDurs;
durs(1) = 125-56;
PROCESSES_TMEV14.T2129_7.PeakDurs = durs;
auc = PROCESSES_TMEV14.T2129_7.RAW_AUC;
auc(1) = sum(abs(T2129_7F(56:125)));
PROCESSES_TMEV14.T2129_7.RAW_AUC = auc;

% Organize Soma and Processes
burnFrame = 90;
[SOMA_PBS_Data,SOMA_PBS_TimeZones,SOMA_PBS_F] = DataOrganizer(SOMA_PBS, PBS_SOMA_DATA, burnFrame);
[SOMA_TMEV2_Data,SOMA_TMEV2_TimeZones,SOMA_TMEV2_F] = DataOrganizer(SOMA_TMEV2, SOMA_TMEV2_DATA, burnFrame);
[SOMA_TMEV5_Data,SOMA_TMEV5_TimeZones,SOMA_TMEV5_F] = DataOrganizer(SOMA_TMEV5, SOMA_TMEV5_DATA, burnFrame);
[SOMA_TMEV14_Data,SOMA_TMEV14_TimeZones,SOMA_TMEV15_F] = DataOrganizer(SOMA_TMEV14, SOMA_TMEV14_DATA, burnFrame);
[PROCESSES_PBS_Data,PROCESSES_PBS_TimeZones,PROCESSES_PBS_F] = DataOrganizer(PROCESSES_PBS, PBS_PROCESSES_DATA, burnFrame);
[PROCESSES_TMEV2_Data,PROCESSES_TMEV2_TimeZones,PROCESSES_TMEV2_F] = DataOrganizer(PROCESSES_TMEV2, PROCESSES_TMEV2_DATA, burnFrame);
[PROCESSES_TMEV5_Data,PROCESSES_TMEV5_TimeZones,PROCESSES_TMEV5_F] = DataOrganizer(PROCESSES_TMEV5, PROCESSES_TMEV5_DATA, burnFrame);
[PROCESSES_TMEV14_Data,PROCESSES_TMEV14_TimeZones,PROCESSES_TMEV15_F] = DataOrganizer(PROCESSES_TMEV14, PROCESSES_TMEV14_DATA, burnFrame);



% Process data to obtain RMS
SOMA_PBS_allRMS = ProcessRMS(PBS_SOMA_DATA,'LPF');
SOMA_TMEV2_allRMS = ProcessRMS(SOMA_TMEV2_DATA,'LPF');
SOMA_TMEV5_allRMS = ProcessRMS(SOMA_TMEV5_DATA,'LPF');
SOMA_TMEV14_allRMS = ProcessRMS(SOMA_TMEV14_DATA,'LPF');

PROCESSES_PBS_allRMS = ProcessRMS(PBS_PROCESSES_DATA,'LPF');
PROCESSES_TMEV2_allRMS = ProcessRMS(PROCESSES_TMEV2_DATA,'LPF');
PROCESSES_TMEV5_allRMS = ProcessRMS(PROCESSES_TMEV5_DATA,'LPF');
PROCESSES_TMEV14_allRMS = ProcessRMS(PROCESSES_TMEV14_DATA,'LPF');

% yMax
yMax = numel(fieldnames(PROCESSES_PBS));

%% SOMA EVENTS/image
sPbsZones  = cell2mat(SOMA_PBS_Data(2:end,10));
somaTraces = SOMA_PBS_Data(2:end,1);
% somaTraces = PBS_SOMA_DATA(1:end,1);
sPbsMice = length(unique(cell2mat(PBS_SOMA_DATA(2:end,7))));
for k = 1:length(somaTraces)
    trace = somaTraces{k};
    somaTraces(k) = {trace(2:5)};
end
somaPbsTraces = str2double(cellstr(somaTraces));
somaPbsTraces = somaPbsTraces(sPbsZones == 3);
[uniqSPBS, idx] = unique(somaPbsTraces);
sPbsFreqs = diff([sort(idx); length(somaPbsTraces)+1]);


sTmev2Zones  = cell2mat(SOMA_TMEV2_Data(2:end,10));
somaTmev2Traces = SOMA_TMEV2_Data(2:end,1);
% somaTmev2Traces = SOMA_TMEV2_DATA(1:end,1);
sTmev2Mice = length(unique(cell2mat(SOMA_TMEV2_DATA(2:end,7))));
for k = 1:length(somaTmev2Traces)
    trace = somaTmev2Traces{k};
    somaTmev2Traces(k) = {trace(2:5)};
end
somaTmev2Traces = str2double(cellstr(somaTmev2Traces));
somaTmev2Traces = somaTmev2Traces(sTmev2Zones == 3);
[uniqSTmev2, idx] = unique(somaTmev2Traces);
sTmev2Freqs = diff([sort(idx); length(somaTmev2Traces)+1]);

sTmev5Zones  = cell2mat(SOMA_TMEV5_Data(2:end,10));
somaTmev5Traces = SOMA_TMEV5_Data(2:end,1);
% somaTmev5Traces = SOMA_TMEV5_DATA(1:end,1);
sTmev5Mice = length(unique(cell2mat(SOMA_TMEV5_DATA(2:end,7))));
for k = 1:length(somaTmev5Traces)
    trace = somaTmev5Traces{k};
    somaTmev5Traces(k) = {trace(2:5)};
end
somaTmev5Traces = str2double(cellstr(somaTmev5Traces));
somaTmev5Traces = somaTmev5Traces(sTmev5Zones == 3);
[uniqSTmev5, idx] = unique(somaTmev5Traces);
sTmev5Freqs = diff([sort(idx); length(somaTmev5Traces)+1]);
% sTmev5Zone1Mice = unique(sTmev5Mice(sTmev5Zones == 1));
% sTmev5Zone2Mice = unique(sTmev5Mice(sTmev5Zones == 2));
% sTmev5Zone3Mice = unique(sTmev5Mice(sTmev5Zones == 3));

sTmev14Zones  = cell2mat(SOMA_TMEV14_Data(2:end,10));
somaTmev14Traces = SOMA_TMEV14_Data(2:end,1);
% somaTmev14Traces = SOMA_TMEV14_Data(1:end,1);
sTmev14Mice = length(unique(cell2mat(SOMA_TMEV14_DATA(2:end,7))));
for k = 1:length(somaTmev14Traces)
    trace = somaTmev14Traces{k};
    somaTmev14Traces(k) = {trace(2:5)};
end
somaTmev14Traces = str2double(cellstr(somaTmev14Traces));
somaTmev14Traces = somaTmev14Traces(sTmev14Zones == 3);
[uniqSTmev14, idx] = unique(somaTmev14Traces);
sTmev14Freqs = diff([sort(idx); length(somaTmev14Traces)+1]);
% sTmev14Zone1Mice = unique(sTmev14Mice(sTmev14Zones == 1));
% sTmev14Zone2Mice = unique(sTmev14Mice(sTmev14Zones == 2));
% sTmev14Zone3Mice = unique(sTmev14Mice(sTmev14Zones == 3));

somaPbsEvents = [uniqSPBS, sPbsFreqs];
somaTmev2Events = [uniqSTmev2, sTmev2Freqs];
somaTmev5Events = [uniqSTmev5, sTmev5Freqs];
somaTmev14Events = [uniqSTmev14, sTmev14Freqs];

%% Processes - Did track reach burn/radius?

pbs = cell2mat(PBS_PROCESSES_DATA(:,11)) ./ cell2mat(PBS_PROCESSES_DATA(:,10));
tmev2 = cell2mat(PROCESSES_TMEV2_DATA(:,11)) ./ cell2mat(PROCESSES_TMEV2_DATA(:,10));
tmev5 = cell2mat(PROCESSES_TMEV5_DATA(:,11)) ./ cell2mat(PROCESSES_TMEV5_DATA(:,10));
tmev15 = cell2mat(PROCESSES_TMEV14_DATA(:,11)) ./ cell2mat(PROCESSES_TMEV14_DATA(:,10));

%% Process events/image
pPbsZones  = cell2mat(PROCESSES_PBS_Data(2:end,10));
processTraces = PROCESSES_PBS_Data(2:end,1);
% processTraces = PBS_PROCESSES_DATA(1:end,1);
pPbsMice = length(unique(cell2mat(PBS_PROCESSES_DATA(2:end,7))));
for k = 1:length(processTraces)
    trace = processTraces{k};
    processTraces(k) = {trace(2:5)};
end
processPbsTraces = str2double(cellstr(processTraces));
processPbsTraces = processPbsTraces(pPbsZones == 3);
[uniqSPBS, idx] = unique(processPbsTraces);
sPbsFreqs = diff([sort(idx); length(processPbsTraces)+1]);
% pPbsZone1Mice = unique(pPbsMice(pPbsZones == 1));
% pPbsZone2Mice = unique(pPbsMice(pPbsZones == 2));
% pPbsZone3Mice = unique(pPbsMice(pPbsZones == 3));

pTmev2Zones  = cell2mat(PROCESSES_TMEV2_Data(2:end,10));
processTmev2Traces = PROCESSES_TMEV2_Data(2:end,1);
% processTmev2Traces = PROCESSES_TMEV2_DATA(1:end,1);
pTmev2Mice = length(unique(cell2mat(PROCESSES_TMEV2_DATA(2:end,7))));
for k = 1:length(processTmev2Traces)
    trace = processTmev2Traces{k};
    processTmev2Traces(k) = {trace(2:5)};
end
processTmev2Traces = str2double(cellstr(processTmev2Traces));
processTmev2Traces = processTmev2Traces(pTmev2Zones == 3);
[uniqSTmev2, idx] = unique(processTmev2Traces);
sTmev2Freqs = diff([sort(idx); length(processTmev2Traces)+1]);
% pTmev2Zone1Mice = unique(pTmev2Mice(pTmev2Zones == 1));
% pTmev2Zone2Mice = unique(pTmev2Mice(pTmev2Zones == 2));
% pTmev2Zone3Mice = unique(pTmev2Mice(pTmev2Zones == 3));

pTmev5Zones  = cell2mat(PROCESSES_TMEV5_Data(2:end,10));
processTmev5Traces = PROCESSES_TMEV5_Data(2:end,1);
% processTmev5Traces = PROCESSES_TMEV5_DATA(1:end,1);
pTmev5Mice = cell2mat(PROCESSES_TMEV5_DATA(2:end,7));
for k = 1:length(processTmev5Traces)
    trace = processTmev5Traces{k};
    processTmev5Traces(k) = {trace(2:5)};
end
processTmev5Traces = str2double(cellstr(processTmev5Traces));
processTmev5Traces = processTmev5Traces(pTmev5Zones == 3);
[uniqSTmev5, idx] = unique(processTmev5Traces);
sTmev5Freqs = diff([sort(idx); length(processTmev5Traces)+1]);
% pTmev5Zone1Mice = unique(pTmev5Mice(pTmev5Zones == 1));
% pTmev5Zone2Mice = unique(pTmev5Mice(pTmev5Zones == 2));
% pTmev5Zone3Mice = unique(pTmev5Mice(pTmev5Zones == 3));

pTmev14Zones  = cell2mat(PROCESSES_TMEV14_Data(2:end,10));
processTmev14Traces = PROCESSES_TMEV14_Data(2:end,1);
% processTmev14Traces = PROCESSES_TMEV14_DATA(1:end,1);
pTmev14Mice = length(unique(cell2mat(PROCESSES_TMEV14_DATA(2:end,7))));
for k = 1:length(processTmev14Traces)
    trace = processTmev14Traces{k};
    processTmev14Traces(k) = {trace(2:5)};
end
processTmev14Traces = str2double(cellstr(processTmev14Traces));
processTmev14Traces = processTmev14Traces(pTmev14Zones == 3);
[uniqSTmev14, idx] = unique(processTmev14Traces);
sTmev14Freqs = diff([sort(idx); length(processTmev14Traces)+1]);

processPbsEvents = [uniqSPBS, sPbsFreqs];
processTmev2Events = [uniqSTmev2, sTmev2Freqs];
processTmev5Events = [uniqSTmev5, sTmev5Freqs];
processTmev14Events = [uniqSTmev14, sTmev14Freqs];

%% SOMA AUC Zones Plotting
limits = [0 250];
PlotZonesAUC(SOMA_PBS_TimeZones, 'Soma PBS',limits);
PlotZonesAUC(SOMA_TMEV2_TimeZones, 'Soma TMEV 2 DPI',limits);
PlotZonesAUC(SOMA_TMEV5_TimeZones, 'Soma TMEV 5 DPI',limits);
PlotZonesAUC(SOMA_TMEV14_TimeZones, 'Soma TMEV 15 DPI',limits);

% % Statistics % %

%% SOMA Peak Durations Zones Plotting
limits = [0 250];
PlotZonesPeakDurs(SOMA_PBS_TimeZones, 'Soma PBS',limits);
PlotZonesPeakDurs(SOMA_TMEV2_TimeZones, 'Soma TMEV 2 DPI',limits);
PlotZonesPeakDurs(SOMA_TMEV5_TimeZones, 'Soma TMEV 5 DPI',limits);
PlotZonesPeakDurs(SOMA_TMEV14_TimeZones, 'Soma TMEV 15 DPI',limits);

%% SOMA Peak Magnitudes Zones Plotting
limits = [0 1.5];
PlotZonesPeakMags(SOMA_PBS_TimeZones, 'Soma PBS',limits);
PlotZonesPeakMags(SOMA_TMEV2_TimeZones, 'Soma TMEV 2 DPI',limits);
PlotZonesPeakMags(SOMA_TMEV5_TimeZones, 'Soma TMEV 5 DPI',limits);
PlotZonesPeakMags(SOMA_TMEV14_TimeZones, 'Soma TMEV 15 DPI',limits);

%% SOMA Peak AUC vs Amp Zones Plotting
limits = [0 16];
PlotZonesSqrtAUCVsAmp(SOMA_PBS_TimeZones, 'Soma PBS', limits);
PlotZonesSqrtAUCVsAmp(SOMA_TMEV2_TimeZones, 'Soma TMEV 2 DPI', limits);
PlotZonesSqrtAUCVsAmp(SOMA_TMEV5_TimeZones, 'Soma TMEV 5 DPI', limits);
PlotZonesSqrtAUCVsAmp(SOMA_TMEV14_TimeZones, 'Soma TMEV 15 DPI', limits);

%% SOMA Peak Duration vs Amp Zones Plotting
limits = [0 250];
PlotZonesDurationVsAmp(SOMA_PBS_TimeZones, 'Soma PBS',limits);
PlotZonesDurationVsAmp(SOMA_TMEV2_TimeZones, 'Soma TMEV 2 DPI',limits);
PlotZonesDurationVsAmp(SOMA_TMEV5_TimeZones, 'Soma TMEV 5 DPI',limits);
PlotZonesDurationVsAmp(SOMA_TMEV14_TimeZones, 'Soma TMEV 15 DPI',limits);

%% SOMA Linear Regression AUC vs Amplitude
varNames = {'Peak Amplitude','AUC','AUC pred'};
somaPbsMagsZone1 = [SOMA_PBS_TimeZones.Zone1_PeakMags]; %SOMA_PBS_TimeZones.Zone2_PeakMags; SOMA_PBS_TimeZones.Zone3_PeakMags];
somaPbsAucZone1 = [SOMA_PBS_TimeZones.Zone1_AUC]; %SOMA_PBS_TimeZones.Zone2_AUC; SOMA_PBS_TimeZones.Zone3_AUC];
[Rsq, PBS1fitLine,PBS1Model] = LinearAnalysis(somaPbsMagsZone1, somaPbsAucZone1,'PBS Soma Auc vs amp Zone1',varNames,180,1.4);

somaTmev2MagsZone1 = [SOMA_TMEV2_TimeZones.Zone1_PeakMags]; %SOMA_TMEV2_TimeZones.Zone2_PeakMags; SOMA_TMEV2_TimeZones.Zone3_PeakMags];
somaTmev2AucZone1 = [SOMA_TMEV2_TimeZones.Zone1_AUC]; %SOMA_TMEV2_TimeZones.Zone2_AUC; SOMA_TMEV2_TimeZones.Zone3_AUC];
[Rsq, TMEV2P1fitLine] = LinearAnalysis(somaTmev2MagsZone1, somaTmev2AucZone1,'TMEV 2 Soma Zone 1 ',varNames);

somaTmev5MagsZone1 = [SOMA_TMEV5_TimeZones.Zone1_PeakMags]; %SOMA_TMEV5_TimeZones.Zone2_PeakMags; SOMA_TMEV5_TimeZones.Zone3_PeakMags];
somaTmev5AucZone1 = [SOMA_TMEV5_TimeZones.Zone1_AUC]; %SOMA_TMEV5_TimeZones.Zone2_AUC; SOMA_TMEV5_TimeZones.Zone3_AUC];
[Rsq, TME52P1fitLine] = LinearAnalysis(somaTmev5MagsZone1, somaTmev5AucZone1,'TMEV 5 Soma Zone 1 ',varNames);

somaTmev15MagsZone1 = [SOMA_TMEV14_TimeZones.Zone1_PeakMags]; %SOMA_TMEV14_TimeZones.Zone2_PeakMags; SOMA_TMEV14_TimeZones.Zone3_PeakMags];
somaTmev15AucZone1 = [SOMA_TMEV14_TimeZones.Zone1_AUC]; %SOMA_TMEV14_TimeZones.Zone2_AUC; SOMA_TMEV14_TimeZones.Zone3_AUC];
[Rsq, TMEV15P1fitLine] = LinearAnalysis(somaTmev15MagsZone1, somaTmev15AucZone1,'TMEV 15 Soma Auc vs amp Zone 1 ',varNames,180,1.4);

somaPbsMagsZone2 = [SOMA_PBS_TimeZones.Zone2_PeakMags]; 
somaPbsAucZone2 = [SOMA_PBS_TimeZones.Zone2_AUC]; 
[Rsq, PBS2fitLine] = LinearAnalysis(somaPbsMagsZone2, somaPbsAucZone2,'PBS Soma Auc vs amp Zone 2 ',varNames,180,1.4);

somaTmev2MagsZone2 = [SOMA_TMEV2_TimeZones.Zone2_PeakMags]; 
somaTmev2AucZone2 = [SOMA_TMEV2_TimeZones.Zone2_AUC]; 
[Rsq, TMEV2P2fitLine] = LinearAnalysis(somaTmev2MagsZone2, somaTmev2AucZone2,'TMEV 2 Soma Zone 2 ',varNames);

% somaTmev5MagsZone2 = [SOMA_TMEV5_TimeZones.Zone2_PeakMags]; 
% somaTmev5AucZone2 = [SOMA_TMEV5_TimeZones.Zone2_AUC];
% [Rsq, TMEV5P2fitLine] = LinearAnalysis(somaTmev5MagsZone2, somaTmev5AucZone2,'TMEV 5 Soma');

somaTmev15MagsZone2 = [SOMA_TMEV14_TimeZones.Zone2_PeakMags]; 
somaTmev15AucZone2 = [SOMA_TMEV14_TimeZones.Zone2_AUC]; 
[Rsq, TMEV15P2fitLine] = LinearAnalysis(somaTmev15MagsZone2, somaTmev15AucZone2,'TMEV 15 Soma Auc vs amp Zone 2',varNames,180,1.4);

somaPbsMagsZone3 = [SOMA_PBS_TimeZones.Zone3_PeakMags]; 
somaPbsAucZone3 = [SOMA_PBS_TimeZones.Zone3_AUC]; 
[Rsq, PBS3fitLine] = LinearAnalysis(somaPbsMagsZone3, somaPbsAucZone3,'PBS Soma Auc vs amp Zone 3',varNames,180,1.4);

somaTmev2MagsZone3 = [SOMA_TMEV2_TimeZones.Zone3_PeakMags]; 
somaTmev2AucZone3 = [SOMA_TMEV2_TimeZones.Zone3_AUC]; 
[Rsq, TMEV2P3fitLine] = LinearAnalysis(somaTmev2MagsZone3, somaTmev2AucZone3,'TMEV 2 Soma Zone 3',varNames);

somaTmev5MagsZone3 = [SOMA_TMEV5_TimeZones.Zone3_PeakMags]; 
somaTmev5AucZone3 = [SOMA_TMEV5_TimeZones.Zone3_AUC]; 
[Rsq, TMEV5P3fitLine] = LinearAnalysis(somaTmev5MagsZone3, somaTmev5AucZone3,'TMEV 5 Soma Zone 3',varNames);

somaTmev15MagsZone3 = [SOMA_TMEV14_TimeZones.Zone3_PeakMags]; 
somaTmev15AucZone3 = [SOMA_TMEV14_TimeZones.Zone3_AUC]; 
[Rsq, TMEV15P3fitLine] = LinearAnalysis(somaTmev15MagsZone3, somaTmev15AucZone3,'TMEV 15 Soma Auc vs amp Zone 3',varNames,180,1.4);

figure()
title('PBS Soma Fit lines'); xlabel('Duration'); ylabel('AUC');


% figure()
h = CustomFigure();
hold on
plot(somaPbsMagsZone2, PBS2fitLine,'-b')
plot(somaPbsMagsZone3, PBS3fitLine,'-r')
plot(somaPbsMagsZone2, somaPbsAucZone2,'ob');
plot(somaPbsMagsZone3, somaPbsAucZone3,'or');
xlabel('Duration (s)'); ylabel('AUC (dF/F)'); title('Soma PBS Phase 2 vs Phase 3'); %legend('Pbs Phase 2','Pbs Phase 3','phase2','phase3'); 
grid on

h = CustomFigure();
hold on
plot(somaTmev15MagsZone2, TMEV15P2fitLine,'-b')
plot(somaTmev15MagsZone3, TMEV15P3fitLine,'-r')
plot(somaTmev15MagsZone2, somaTmev15AucZone2,'ob');
plot(somaTmev15MagsZone3, somaTmev15AucZone3,'or');
xlabel('Duration (s)'); ylabel('AUC (dF/F)'); title('Soma TMEV15 Phase 2 vs Phase 3'); %legend('Pbs Phase 2','Pbs Phase 3','phase2','phase3'); 
grid on

figure()
hold on
plot(somaTmev2PbsMagsZone2, TMEV2P2fitLine,'-b')
plot(somaTmev2PbsMagsZone3, TMEV2P3fitLine,'-r')
legend('Pbs Phase 2','Pbs Phase 3'); xlabel('Duration (s)'); ylabel('AUC (dF/F)'); title('Soma PBS Phase 2 vs Phase 3');
grid on


figure()
hold on
scatter(somaPbsMagsZone3,somaPbsAucZone3,'or','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaPbsMagsZone2,somaPbsAucZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaPbsMagsZone1,somaPbsAucZone1,'ok','MarkerFaceColor',[0 0.75 0.25],'MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',0.2)
plot(somaPbsMagsZone2,PBS2fitLine,'-b');
plot(somaPbsMagsZone3,PBS3fitLine,'-r');
plot(somaPbsMagsZone1,PBS1fitLine,'-g');
legend('phase 3','phase 2','phase 1'); xlabel('duration'); ylabel('AUC'); title('pbs SOMA');

figure()
hold on
scatter(somaTmev15MagsZone3,somaTmev15AucZone3,'or','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaTmev15MagsZone2,somaTmev15AucZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaTmev15MagsZone1,somaTmev15AucZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.2)
plot(somaTmev15MagsZone2,TMEV15P2fitLine,'-b');
plot(somaTmev15MagsZone3,TMEV15P3fitLine,'-r');
plot(somaTmev15MagsZone1,TMEV15P1fitLine,'-','Color',[0.9 0.9 0]);
legend('phase 3','phase 2','phase 1'); xlabel('duration'); ylabel('AUC'); title('Tmev15 SOMA');

figure()
hold on
scatter(somaPbsMagsZone3,somaPbsAucZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaPbsMagsZone2,somaPbsAucZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaPbsMagsZone1,somaPbsAucZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.2)
plot(somaPbsMagsZone2,PBS2fitLine,'-b');
plot(somaPbsMagsZone3,PBS3fitLine,'-r');
plot(somaPbsMagsZone1,PBS1fitLine,'-','Color',[0.9 0.9 0]);
legend('phase 3','phase 2','phase 1'); xlabel('amplitude'); ylabel('AUC'); title('pbs soma');
xlim([0 1.4]);
figure()
hold on
scatter(somaTmev15MagsZone3,somaTmev15AucZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaTmev15MagsZone2,somaTmev15AucZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaTmev15MagsZone1,somaTmev15AucZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.2)
plot(somaTmev15MagsZone2,TMEV15P2fitLine,'-b');
plot(somaTmev15MagsZone3,TMEV15P3fitLine,'-r');
plot(somaTmev15MagsZone1,TMEV15P1fitLine,'-','Color',[0.9 0.9 0]);
legend('phase 3','phase 2','phase 1'); xlabel('amplitude'); ylabel('AUC'); title('Tmev15 soma');
xlim([0 1.4]);

%% SOMA Linear Regression Duration vs Amplitude
varNames = {'Peak Amplitude','Peak Duration','duration pred'};
SOMAPbsAmpZone1 = [SOMA_PBS_TimeZones.Zone1_PeakMags]; 
SOMAPbsPeakDursZone1 = [SOMA_PBS_TimeZones.Zone1_PeakDurs];
[Rsq, PBS1fitLine] = LinearAnalysis(SOMAPbsAmpZone1, SOMAPbsPeakDursZone1,'PBS SOMA Durs vs amp Zone 1',varNames,250,2);

SOMAPbsAmpZone2 = [SOMA_PBS_TimeZones.Zone2_PeakMags]; 
SOMAPbsPeakDursZone2 = [SOMA_PBS_TimeZones.Zone2_PeakDurs];
[Rsq, PBS2fitLine] = LinearAnalysis(SOMAPbsAmpZone2, SOMAPbsPeakDursZone2,'PBS SOMA Durs vs amp Zone 2',varNames,250,2);

SOMAPbsAmpZone3 = [SOMA_PBS_TimeZones.Zone3_PeakMags]; 
SOMAPbsPeakDursZone3 = [SOMA_PBS_TimeZones.Zone3_PeakDurs];
[Rsq, PBS3fitLine] = LinearAnalysis(SOMAPbsAmpZone3, SOMAPbsPeakDursZone3,'PBS SOMA Durs vs amp Zone 3',varNames,250,2);

SOMATmev2AmpZone1 = [SOMA_TMEV2_TimeZones.Zone1_PeakMags]; 
SOMATmev2PeakDursZone1 = [SOMA_TMEV2_TimeZones.Zone1_PeakDurs];
[Rsq, fitLine] = LinearAnalysis(SOMATmev2AmpZone1, SOMATmev2PeakDursZone1,'Tmev2 SOMA Zone 1',varNames);

SOMATmev2AmpZone2 = [SOMA_PBS_TimeZones.Zone2_PeakMags]; 
SOMATmev2PeakDursZone2 = [SOMA_PBS_TimeZones.Zone2_PeakDurs];
[Rsq, fitLine] = LinearAnalysis(SOMAPbsAmpZone2, SOMAPbsPeakDursZone2,'Tmev 2 SOMA Zone 2',varNames);

SOMATmev2AmpZone3 = [SOMA_PBS_TimeZones.Zone3_PeakMags]; 
SOMATmev2PeakDursZone3 = [SOMA_PBS_TimeZones.Zone3_PeakDurs];
[Rsq, fitLine] = LinearAnalysis(SOMAPbsAmpZone3, SOMAPbsPeakDursZone3,'Tmev2 SOMA Zone 3',varNames);

SOMATmev5AmpZone1 = [SOMA_TMEV5_TimeZones.Zone1_PeakMags];
SOMATmev5PeakDursZone1 = [SOMA_TMEV5_TimeZones.Zone1_PeakDurs];
[Rsq, fitLine] = LinearAnalysis(SOMATmev5AmpZone1, SOMATmev5PeakDursZone1,'TMEV 5 SOMA Zone 1',varNames);

SOMATmev5AmpZone2 = [SOMA_TMEV5_TimeZones.Zone2_PeakMags];
SOMATmev5PeakDursZone2 = [SOMA_TMEV5_TimeZones.Zone2_PeakDurs];
[Rsq, fitLine] = LinearAnalysis(SOMATmev5AmpZone2, SOMATmev5PeakDursZone2,'TMEV 5 SOMA Zone 2',varNames);

SOMATmev5AmpZone3 = [SOMA_TMEV5_TimeZones.Zone3_PeakMags];
SOMATmev5PeakDursZone3 = [SOMA_TMEV5_TimeZones.Zone3_PeakDurs];
[Rsq, fitLine] = LinearAnalysis(SOMATmev5AmpZone3, SOMATmev5PeakDursZone3,'TMEV 5 SOMA Zone 3',varNames);

SOMATmev15AmpZone1 = [SOMA_TMEV14_TimeZones.Zone1_PeakMags];
SOMATmev15PeakDursZone1 = [SOMA_TMEV14_TimeZones.Zone1_PeakDurs];
[Rsq, TMEV15P1fitLine] = LinearAnalysis(SOMATmev15AmpZone1, SOMATmev15PeakDursZone1,'TMEV 15 SOMA Durs vs amp Zone 1',varNames,250,2);

SOMATmev15AmpZone2 = [SOMA_TMEV14_TimeZones.Zone2_PeakMags];
SOMATmev15PeakDursZone2 = [SOMA_TMEV14_TimeZones.Zone2_PeakDurs];
[Rsq, TMEV15P2fitLine] = LinearAnalysis(SOMATmev15AmpZone2, SOMATmev15PeakDursZone2,'TMEV 15 SOMA Durs vs amp Zone 2',varNames,250,2);

SOMATmev15AmpZone3 = [SOMA_TMEV14_TimeZones.Zone3_PeakMags];
SOMATmev15PeakDursZone3 = [SOMA_TMEV14_TimeZones.Zone3_PeakDurs];
[Rsq, TMEV15P3fitLine] = LinearAnalysis(SOMATmev15AmpZone3, SOMATmev15PeakDursZone3,'TMEV 15 SOMA Durs vs amp Zone 3',varNames,250,2);

figure()
hold on
scatter(SOMAPbsAmpZone3,SOMAPbsPeakDursZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(SOMAPbsAmpZone2,SOMAPbsPeakDursZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(SOMAPbsAmpZone1,SOMAPbsPeakDursZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.2)
plot(SOMAPbsAmpZone2,PBS2fitLine,'-b');
plot(SOMAPbsAmpZone3,PBS3fitLine,'-r');
plot(SOMAPbsAmpZone1,PBS1fitLine,'-','Color',[0.9 0.9 0]);
legend('phase 3','phase 2','phase 1'); xlabel('amplitude'); ylabel('PeakDurs'); title('pbs SOMA');
xlim([0 1.4]); ylim([0 250]);
figure()
hold on
scatter(SOMATmev15AmpZone3,SOMATmev15PeakDursZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(SOMATmev15AmpZone2,SOMATmev15PeakDursZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(SOMATmev15AmpZone1,SOMATmev15PeakDursZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.2)
plot(SOMATmev15AmpZone2,TMEV15P2fitLine,'-b');
plot(SOMATmev15AmpZone3,TMEV15P3fitLine,'-r');
plot(SOMATmev15AmpZone1,TMEV15P1fitLine,'-','Color',[0.9 0.9 0]);
legend('phase 3','phase 2','phase 1'); xlabel('amplitude'); ylabel('PeakDurs'); title('Tmev15 SOMA');
xlim([0 1.4]); ylim([0 250]);

%% SOMA Linear Regression AUC vs Duration
varNames = {'Peak Duration','AUC','AUC pred'};
somaPbsDursZone1 = [SOMA_PBS_TimeZones.Zone1_PeakDurs]; %SOMA_PBS_TimeZones.Zone2_PeakDurs; SOMA_PBS_TimeZones.Zone3_PeakDurs];
somaPbsAucZone1 = [SOMA_PBS_TimeZones.Zone1_AUC]; %SOMA_PBS_TimeZones.Zone2_AUC; SOMA_PBS_TimeZones.Zone3_AUC];
[Rsq, PBS1fitLine,PBS1Model] = LinearAnalysis(somaPbsDursZone1, somaPbsAucZone1,'PBS Soma Zone1 AUC vs Dur',varNames,180,250);

somaTmev2DursZone1 = [SOMA_TMEV2_TimeZones.Zone1_PeakDurs]; %SOMA_TMEV2_TimeZones.Zone2_PeakDurs; SOMA_TMEV2_TimeZones.Zone3_PeakDurs];
somaTmev2AucZone1 = [SOMA_TMEV2_TimeZones.Zone1_AUC]; %SOMA_TMEV2_TimeZones.Zone2_AUC; SOMA_TMEV2_TimeZones.Zone3_AUC];
[Rsq, TMEV2P1fitLine] = LinearAnalysis(somaTmev2DursZone1, somaTmev2AucZone1,'TMEV 2 Soma Zone 1 AUC vs Dur',varNames);

somaTmev5DursZone1 = [SOMA_TMEV5_TimeZones.Zone1_PeakDurs]; %SOMA_TMEV5_TimeZones.Zone2_PeakDurs; SOMA_TMEV5_TimeZones.Zone3_PeakDurs];
somaTmev5AucZone1 = [SOMA_TMEV5_TimeZones.Zone1_AUC]; %SOMA_TMEV5_TimeZones.Zone2_AUC; SOMA_TMEV5_TimeZones.Zone3_AUC];
[Rsq, TME52P1fitLine] = LinearAnalysis(somaTmev5DursZone1, somaTmev5AucZone1,'TMEV 5 Soma Zone 1 AUC vs Dur',varNames);

somaTmev15DursZone1 = [SOMA_TMEV14_TimeZones.Zone1_PeakDurs]; %SOMA_TMEV14_TimeZones.Zone2_PeakDurs; SOMA_TMEV14_TimeZones.Zone3_PeakDurs];
somaTmev15AucZone1 = [SOMA_TMEV14_TimeZones.Zone1_AUC]; %SOMA_TMEV14_TimeZones.Zone2_AUC; SOMA_TMEV14_TimeZones.Zone3_AUC];
[Rsq, TMEV15P1fitLine] = LinearAnalysis(somaTmev15DursZone1, somaTmev15AucZone1,'TMEV 15 Soma Zone 1 AUC vs Dur',varNames,180,250);

somaPbsDursZone2 = [SOMA_PBS_TimeZones.Zone2_PeakDurs]; 
somaPbsAucZone2 = [SOMA_PBS_TimeZones.Zone2_AUC]; 
[Rsq, PBS2fitLine] = LinearAnalysis(somaPbsDursZone2, somaPbsAucZone2,'PBS Soma Zone 2 AUC vs Dur',varNames,180,250);

somaTmev2DursZone2 = [SOMA_TMEV2_TimeZones.Zone2_PeakDurs]; 
somaTmev2AucZone2 = [SOMA_TMEV2_TimeZones.Zone2_AUC]; 
[Rsq, TMEV2P2fitLine] = LinearAnalysis(somaTmev2DursZone2, somaTmev2AucZone2,'TMEV 2 Soma Zone 2 AUC vs Dur',varNames);
% 
% somaTmev5DursZone2 = [SOMA_TMEV5_TimeZones.Zone2_PeakDurs]; 
% somaTmev5AucZone2 = [SOMA_TMEV5_TimeZones.Zone2_AUC];
% [Rsq, TMEV5P2fitLine] = LinearAnalysis(somaTmev5DursZone2, somaTmev5AucZone2,'TMEV 5 Soma');

somaTmev15DursZone2 = [SOMA_TMEV14_TimeZones.Zone2_PeakDurs]; 
somaTmev15AucZone2 = [SOMA_TMEV14_TimeZones.Zone2_AUC]; 
[Rsq, TMEV15P2fitLine] = LinearAnalysis(somaTmev15DursZone2, somaTmev15AucZone2,'TMEV 15 Soma Zone 2 AUC vs Dur',varNames,180,250);

somaPbsDursZone3 = [SOMA_PBS_TimeZones.Zone3_PeakDurs]; 
somaPbsAucZone3 = [SOMA_PBS_TimeZones.Zone3_AUC]; 
[Rsq, PBS3fitLine] = LinearAnalysis(somaPbsDursZone3, somaPbsAucZone3,'PBS Soma Zone 3 AUC vs Dur',varNames,180,250);

somaTmev2DursZone3 = [SOMA_TMEV2_TimeZones.Zone3_PeakDurs]; 
somaTmev2AucZone3 = [SOMA_TMEV2_TimeZones.Zone3_AUC]; 
[Rsq, TMEV2P3fitLine] = LinearAnalysis(somaTmev2DursZone3, somaTmev2AucZone3,'TMEV 2 Soma Zone 3 AUC vs Dur',varNames);

somaTmev5DursZone3 = [SOMA_TMEV5_TimeZones.Zone3_PeakDurs]; 
somaTmev5AucZone3 = [SOMA_TMEV5_TimeZones.Zone3_AUC]; 
[Rsq, TMEV5P3fitLine] = LinearAnalysis(somaTmev5DursZone3, somaTmev5AucZone3,'TMEV 5 Soma Zone 3 AUC vs Dur',varNames);

somaTmev15DursZone3 = [SOMA_TMEV14_TimeZones.Zone3_PeakDurs]; 
somaTmev15AucZone3 = [SOMA_TMEV14_TimeZones.Zone3_AUC]; 
[Rsq, TMEV15P3fitLine] = LinearAnalysis(somaTmev15DursZone3, somaTmev15AucZone3,'TMEV 15 Soma Zone 3 AUC vs Dur',varNames,180,250);

figure()
title('PBS Soma Fit lines'); xlabel('Duration'); ylabel('AUC');


% figure()
h = CustomFigure();
hold on
plot(somaPbsDursZone2, PBS2fitLine,'-b')
plot(somaPbsDursZone3, PBS3fitLine,'-r')
plot(somaPbsDursZone2, somaPbsAucZone2,'ob');
plot(somaPbsDursZone3, somaPbsAucZone3,'or');
xlabel('Duration (s)'); ylabel('AUC (dF/F)'); title('Soma PBS Phase 2 vs Phase 3'); %legend('Pbs Phase 2','Pbs Phase 3','phase2','phase3'); 
grid on


figure()
hold on
scatter(somaPbsDursZone3,somaPbsAucZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaPbsDursZone2,somaPbsAucZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaPbsDursZone1,somaPbsAucZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.2)
plot(somaPbsDursZone2,PBS2fitLine,'-b');
plot(somaPbsDursZone3,PBS3fitLine,'-r');
plot(somaPbsDursZone1,PBS1fitLine,'-','Color',[0.9 0.9 0]);
legend('phase 3','phase 2','phase 1'); xlabel('duration'); ylabel('AUC'); title('pbs SOMA');
xlim([0 250]);
figure()
hold on
scatter(somaTmev15DursZone3,somaTmev15AucZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaTmev15DursZone2,somaTmev15AucZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaTmev15DursZone1,somaTmev15AucZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.2)
plot(somaTmev15DursZone2,TMEV15P2fitLine,'-b');
plot(somaTmev15DursZone3,TMEV15P3fitLine,'-r');
plot(somaTmev15DursZone1,TMEV15P1fitLine,'-','Color',[0.9 0.9 0]);
legend('phase 3','phase 2','phase 1'); xlabel('duration'); ylabel('AUC'); title('Tmev15 SOMA');
xlim([0 250]);

figure()
hold on
scatter(somaPbsDursZone3,somaPbsAucZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaPbsDursZone2,somaPbsAucZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaPbsDursZone1,somaPbsAucZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.2)
plot(somaPbsDursZone2,PBS2fitLine,'-b');
plot(somaPbsDursZone3,PBS3fitLine,'-r');
plot(somaPbsDursZone1,PBS1fitLine,'-g');
legend('phase 3','phase 2','phase 1'); xlabel('duration'); ylabel('AUC'); title('pbs soma');

figure()
hold on
scatter(somaTmev15DursZone3,somaTmev15AucZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaTmev15DursZone2,somaTmev15AucZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(somaTmev15DursZone1,somaTmev15AucZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.2)
plot(somaTmev15DursZone2,TMEV15P2fitLine,'-b');
plot(somaTmev15DursZone3,TMEV15P3fitLine,'-r');
plot(somaTmev15DursZone1,TMEV15P1fitLine,'-','Color',[0.9 0.9 0]);
legend('phase 3','phase 2','phase 1'); xlabel('duration'); ylabel('AUC'); title('Tmev15 soma');

%% SOMA Events Dot Plot
time = round((1740/0.97)/60);
yMax = numel(fieldnames(PROCESSES_PBS));
xlimit = [-time time];
EventDotPlot(SOMA_PBS, 'Soma PBS Phase 2-3', xlimit, yMax);
EventDotPlot(SOMA_TMEV2, 'Soma TMEV 2 DPI Phase 2-3', xlimit, yMax);
EventDotPlot(SOMA_TMEV5, 'Soma TMEV 5 DPI Phase 2-3', xlimit, yMax);
EventDotPlot(SOMA_TMEV14, 'Soma TMEV 15 DPI Phase 2-3', xlimit, yMax);

%% SOMA Phase0-1 Events Dot Plot
time = round((1740/0.97)/60);
yMax = numel(fieldnames(PROCESSES_PBS));
xlimit = [-1 0.5];
Phase01DotPlot(SOMA_PBS, 'Soma PBS Phase 0-1', xlimit, yMax);
Phase01DotPlot(SOMA_TMEV2, 'Soma TMEV 2 DPI Phase 0-1', xlimit, yMax);
Phase01DotPlot(SOMA_TMEV5, 'Soma TMEV 5 DPI Phase 0-1', xlimit, yMax);
Phase01DotPlot(SOMA_TMEV14, 'Soma TMEV 15 DPI Phase 0-1', xlimit, yMax);

%% SOMA Grouped Bar Charts
barTitle = 'Soma';
xDim = 2;
yDim = 1.5;
PairedBarGraph(barTitle, xDim, yDim, SOMA_PBS_TimeZones, SOMA_TMEV2_TimeZones, SOMA_TMEV5_TimeZones, SOMA_TMEV14_TimeZones);

%% SOMA Zone derivs
barTitle = 'Soma';
xDim = 2;
yDim = 1.5;
PairedZonesDeriv(barTitle, xDim, yDim, PBS_SOMA_DATA, SOMA_TMEV2_DATA, SOMA_TMEV5_DATA, SOMA_TMEV14_DATA);

%% SOMA UnGrouped Bar Charts
barTitle = 'Soma';
xDim = 1.0;
yDim = 1.5;
SingleZoneBar('Soma PBS', xDim, yDim,SOMA_PBS_TimeZones, [0, 0, 0])
SingleZoneBar('Soma TMEV 2 DPI', xDim, yDim,SOMA_TMEV2_TimeZones, [1 0 0])
SingleZoneBar('Soma TMEV 5 DPI', xDim, yDim,SOMA_TMEV5_TimeZones, [1 0 0])
SingleZoneBar('Soma TMEV 15 DPI', xDim, yDim,SOMA_TMEV14_TimeZones, [1 0 0])

%% SOMA Plot RMS Data
ScatterRMS(SOMA_PBS_allRMS,'PBS'); % use T2114_1 as example
ScatterRMS(SOMA_TMEV2_allRMS,'TMEV 2 DPI'); % T2190_10
ScatterRMS(SOMA_TMEV5_allRMS,'TMEV 5 DPI'); % T2102_12
ScatterRMS(SOMA_TMEV14_allRMS,'TMEV 15 DPI'); % T2326_1

%% SOMA RMS Examples
LPF = noiseRemover(100, 0.01, 0.97);
% PBS Example
pbsMarker = 250;
pbsTrace = PBS_SOMA_DATA{51,3};
pbsTrace = (pbsTrace - mean(pbsTrace(1:50)))/mean(pbsTrace(1:50));
% pbsTrace = filtfilt(LPF, double(pbsTrace));
figure()
subplot(1,3,1)
plot(pbsTrace(1:200))
title(strcat('RMS: ',num2str(rms(pbsTrace(1:200)))));
subplot(1,3,2)
plot(pbsTrace(201:pbsMarker-1))
title(strcat('RMS: ',num2str(rms(pbsTrace(201:pbsMarker-1)))));
subplot(1,3,3)
plot(pbsTrace(pbsMarker:end))
title(strcat('RMS: ',num2str(rms(pbsTrace(pbsMarker:end)))));
suptitle('Soma PBS: T2114_1');
% TMEV 2 Example
tmev2Marker = 1508;
tmev2Trace = SOMA_TMEV2_DATA{16,3};
tmev2Trace = (tmev2Trace - mean(tmev2Trace(1:50)))/mean(tmev2Trace(1:50));
figure()
subplot(1,3,1)
plot(tmev2Trace(1:200))
title(strcat('RMS: ',num2str(rms(tmev2Trace(1:200)))));
subplot(1,3,2)
plot(tmev2Trace(201:tmev2Marker-1))
title(strcat('RMS: ',num2str(rms(tmev2Trace(201:tmev2Marker-1)))));
subplot(1,3,3)
plot(tmev2Trace(tmev2Marker:end))
title(strcat('RMS: ',num2str(rms(tmev2Trace(tmev2Marker:end)))));
suptitle('Soma TMEV 2 DPI: T2190_10');
% TMEV 5 Example
tmev5Marker = 858;
tmev5Trace = SOMA_TMEV5_DATA{5,3};
tmev5Trace = (tmev5Trace - mean(tmev5Trace(1:50)))/mean(tmev5Trace(1:50));
figure()
subplot(1,3,1)
plot(tmev5Trace(1:200))
title(strcat('RMS: ',num2str(rms(tmev5Trace(1:200)))));
subplot(1,3,2)
plot(tmev5Trace(201:tmev5Marker-1))
title(strcat('RMS: ',num2str(rms(tmev5Trace(201:tmev5Marker-1)))));
subplot(1,3,3)
plot(tmev5Trace(tmev5Marker:end))
title(strcat('RMS: ',num2str(rms(tmev5Trace(tmev5Marker:end)))));
suptitle('Soma TMEV 5 DPI: T2102_12');
% TMEV 15 Example
tmev15Marker = 250;
tmev15Trace = SOMA_TMEV14_DATA{7,3};
tmev15Trace = (tmev15Trace - mean(tmev15Trace(1:50)))/mean(tmev15Trace(1:50));
figure()
subplot(1,3,1)
plot(tmev15Trace(1:200))
title(strcat('RMS: ',num2str(rms(tmev15Trace(1:200)))));
subplot(1,3,2)
plot(tmev15Trace(201:tmev15Marker-1))
title(strcat('RMS: ',num2str(rms(tmev15Trace(201:tmev15Marker-1)))));
subplot(1,3,3)
plot(tmev15Trace(tmev15Marker:end))
title(strcat('RMS: ',num2str(rms(tmev15Trace(tmev15Marker:end)))));
suptitle('Soma TMEV 15 DPI: T2326_1');

%% SOMA LENGTH Examples
pbsTraceA = sum(pbsTrace(1:200))/length(pbsTrace(1:200));
pbsTraceB = sum(pbsTrace(201:pbsMarker-1))/length(pbsTrace(201:pbsMarker-1));
pbsTraceC = sum(pbsTrace(pbsMarker:end))/length(pbsTrace(pbsMarker:end));

tmev2TraceA = sum(tmev2Trace(1:200))/length(tmev2Trace(1:200));
tmev2TraceB = sum(tmev2Trace(201:tmev2Marker-1))/length(tmev2Trace(201:tmev2Marker-1));
tmev2TraceC = sum(tmev2Trace(tmev2Marker:end))/length(tmev2Trace(tmev2Marker:end));

tmev5TraceA = sum(tmev5Trace(1:200))/length(tmev5Trace(1:200));
tmev5TraceB = sum(tmev5Trace(201:tmev5Marker-1))/length(tmev5Trace(201:tmev5Marker-1));
tmev5TraceC = sum(tmev5Trace(tmev5Marker:end))/length(tmev5Trace(tmev5Marker:end));

tmev15TraceA = sum(tmev15Trace(1:200))/length(tmev15Trace(1:200));
tmev15TraceB = sum(tmev15Trace(201:tmev15Marker-1))/length(tmev15Trace(201:tmev15Marker-1));
tmev15TraceC = sum(tmev15Trace(tmev15Marker:end))/length(tmev15Trace(tmev15Marker:end));

figure()
subplot(1,3,1)
plot(pbsTrace(1:200))
title(strcat('Length: ',num2str(pbsTraceA)));
subplot(1,3,2)
plot(pbsTrace(201:pbsMarker-1))
title(strcat('Length: ',num2str(pbsTraceB)));
subplot(1,3,3)
plot(pbsTrace(pbsMarker:end))
title(strcat('Length: ',num2str(pbsTraceC)));
suptitle('Soma PBS: T2114_1');
% TMEV 2 Example
figure()
subplot(1,3,1)
plot(tmev2Trace(1:200))
title(strcat('Length: ',num2str(tmev2TraceA)));
subplot(1,3,2)
plot(tmev2Trace(201:tmev2Marker-1))
title(strcat('Length: ',num2str(tmev2TraceB)));
subplot(1,3,3)
plot(tmev2Trace(tmev2Marker:end))
title(strcat('Length: ',num2str(tmev2TraceC)));
suptitle('Soma TMEV 2 DPI: T2190_10');
% TMEV 5 Example
figure()
subplot(1,3,1)
plot(tmev5Trace(1:200))
title(strcat('Length: ',num2str(tmev5TraceA)));
subplot(1,3,2)
plot(tmev5Trace(201:tmev5Marker-1))
title(strcat('Length: ',num2str(tmev5TraceB)));
subplot(1,3,3)
plot(tmev5Trace(tmev5Marker:end))
title(strcat('Length: ',num2str(tmev5TraceC)));
suptitle('Soma TMEV 5 DPI: T2102_12');
% TMEV 15 Example
figure()
subplot(1,3,1)
plot(tmev15Trace(1:200))
title(strcat('Length: ',num2str(tmev15TraceA)));
subplot(1,3,2)
plot(tmev15Trace(201:tmev15Marker-1))
title(strcat('Length: ',num2str(tmev15TraceB)));
subplot(1,3,3)
plot(tmev15Trace(tmev15Marker:end))
title(strcat('Length: ',num2str(tmev15TraceC)));
suptitle('Soma TMEV 15 DPI: T2326_1');

%% SOMA Approx. Derivative Examples
LPF = noiseRemover(100, 0.02, 0.97);
pbsTraceA = sum(abs(diff(pbsTrace(1:200))))/length(pbsTrace(1:200));
pbsTraceA = filtfilt(LPF,double(pbsTraceA));
pbsTraceB = sum(abs(diff(pbsTrace(201:pbsMarker-1))))/length(pbsTrace(201:pbsMarker-1));
pbsTraceB = filtfilt(LPF,double(pbsTraceB));
pbsTraceC = sum(abs(diff(pbsTrace(pbsMarker:end))))/length(pbsTrace(pbsMarker:end));
pbsTraceC = filtfilt(LPF,double(pbsTraceC));


tmev2TraceA = sum(abs(diff(tmev2Trace(1:200))))/length(tmev2Trace(1:200));
tmev2TraceA = filtfilt(LPF,double(tmev2TraceA));
tmev2TraceB = sum(abs(diff(tmev2Trace(201:tmev2Marker-1))))/length(tmev2Trace(201:tmev2Marker-1));
tmev2TraceB = filtfilt(LPF,double(tmev2TraceB));
tmev2TraceC = sum(abs(diff(tmev2Trace(tmev2Marker:end))))/length(tmev2Trace(tmev2Marker:end));
tmev2TraceC = filtfilt(LPF,double(tmev2TraceC));

tmev5TraceA = sum(abs(diff(tmev5Trace(1:200))))/length(tmev5Trace(1:200));
tmev5TraceA = filtfilt(LPF,double(tmev5TraceA));
tmev5TraceB = sum(abs(diff(tmev5Trace(201:tmev5Marker-1))))/length(tmev5Trace(201:tmev5Marker-1));
tmev5TraceB = filtfilt(LPF,double(tmev5TraceB));
tmev5TraceC = sum(abs(diff(tmev5Trace(tmev5Marker:end))))/length(tmev5Trace(tmev5Marker:end));
tmev5TraceC = filtfilt(LPF,double(tmev5TraceC));

tmev15TraceA = sum(abs(diff(tmev15Trace(1:200))))/length(tmev15Trace(1:200));
tmev15TraceA = filtfilt(LPF,double(tmev15TraceA));
tmev15TraceB = sum(abs(diff(tmev15Trace(201:tmev15Marker-1))))/length(tmev15Trace(201:tmev15Marker-1));
tmev15TraceB = filtfilt(LPF,double(tmev15TraceB));
tmev15TraceC = sum(abs(diff(tmev15Trace(tmev15Marker:end))))/length(tmev15Trace(tmev15Marker:end));
tmev15TraceC = filtfilt(LPF,double(tmev15TraceC));

figure()
subplot(1,3,1)
plot(pbsTrace(1:200))
title(strcat('Approx. Deriv: ',num2str(pbsTraceA)));
subplot(1,3,2)
plot(pbsTrace(201:pbsMarker-1))
title(strcat('Approx. Deriv: ',num2str(pbsTraceB)));
subplot(1,3,3)
plot(pbsTrace(pbsMarker:end))
title(strcat('Approx. Deriv: ',num2str(pbsTraceC)));
suptitle('Soma PBS: T2114_1');
% TMEV 2 Example
figure()
subplot(1,3,1)
plot(tmev2Trace(1:200))
title(strcat('Approx. Deriv: ',num2str(tmev2TraceA)));
subplot(1,3,2)
plot(tmev2Trace(201:tmev2Marker-1))
title(strcat('Approx. Deriv: ',num2str(tmev2TraceB)));
subplot(1,3,3)
plot(tmev2Trace(tmev2Marker:end))
title(strcat('Approx. Deriv: ',num2str(tmev2TraceC)));
suptitle('Soma TMEV 2 DPI: T2190_10');
% TMEV 5 Example
figure()
subplot(1,3,1)
plot(tmev5Trace(1:200))
title(strcat('Approx. Deriv: ',num2str(tmev5TraceA)));
subplot(1,3,2)
plot(tmev5Trace(201:tmev5Marker-1))
title(strcat('Approx. Deriv: ',num2str(tmev5TraceB)));
subplot(1,3,3)
plot(tmev5Trace(tmev5Marker:end))
title(strcat('Approx. Deriv: ',num2str(tmev5TraceC)));
suptitle('Soma TMEV 5 DPI: T2102_12');
% TMEV 15 Example
figure()
subplot(1,3,1)
plot(tmev15Trace(1:200))
title(strcat('Approx. Deriv: ',num2str(tmev15TraceA)));
subplot(1,3,2)
plot(tmev15Trace(201:tmev15Marker-1))
title(strcat('Approx. Deriv: ',num2str(tmev15TraceB)));
subplot(1,3,3)
plot(tmev15Trace(tmev15Marker:end))
title(strcat('Approx. Deriv: ',num2str(tmev15TraceC)));
suptitle('Soma TMEV 15 DPI: T2326_1');

%% SOMA Plot Signal Length Data
PairedSignalLength(PBS_SOMA_DATA,'PBS');
PairedSignalLength(SOMA_TMEV2_DATA,'TMEV 2 DPI');
PairedSignalLength(SOMA_TMEV5_DATA,'TMEV 5 DPI');
PairedSignalLength(SOMA_TMEV14_DATA,'TMEV 15 DPI');

%% SOMA Plot Sum of Approx Deriv Data
xDim = 1.5;
yDim = 2.5;
yMax = 0.03;
[somaPBSzones, SOMA_PBS_PVALUES, somaPbsVideos] = PairedDeriv(PBS_SOMA_DATA,'Soma PBS',xDim, yDim, [0, 0, 0], yMax);
[somaTMEV2zones, SOMA_TMEV2_PVALUES, somaTmev2Videos] = PairedDeriv(SOMA_TMEV2_DATA,'Soma TMEV 2 DPI',xDim, yDim, [1 0 0], yMax);
[somaTMEV5zones, SOMA_TMEV5_PVALUES, somaTmev5Videos] = PairedDeriv(SOMA_TMEV5_DATA,'Soma TMEV 5 DPI',xDim, yDim, [1 0 0], yMax);
[somaTMEV15zones, SOMA_TMEV15_PVALUES, somaTmev15Videos] = PairedDeriv(SOMA_TMEV14_DATA,'Soma TMEV 15 DPI',xDim, yDim, [1 0 0], yMax);

%% SOMA Plot All zones Deriv Data
xDim = 2.5;
yDim = 2.5;
yMax = 0.03;
SOMA_P = AllDerivZones(PBS_SOMA_DATA,SOMA_TMEV2_DATA,SOMA_TMEV5_DATA,SOMA_TMEV14_DATA,'Soma', yMax, xDim, yDim);

%% SOAM Plot Interevent Interval %%
somaPbsIEIZone1 = SOMA_PBS_TimeZones.Zone1_IEI;
somaPbsIEIZone2 = SOMA_PBS_TimeZones.Zone2_IEI;
somaPbsIEIZone3 = SOMA_PBS_TimeZones.Zone3_IEI;

somaTmev2IEIZone1 = SOMA_TMEV2_TimeZones.Zone1_IEI;
somaTmev2IEIZone2 = SOMA_TMEV2_TimeZones.Zone2_IEI;
somaTmev2IEIZone3 = SOMA_TMEV2_TimeZones.Zone3_IEI;

somaTmev5IEIZone1 = SOMA_TMEV5_TimeZones.Zone1_IEI;
somaTmev5IEIZone2 = SOMA_TMEV5_TimeZones.Zone2_IEI;
somaTmev5IEIZone3 = SOMA_TMEV5_TimeZones.Zone3_IEI;

somaTmev14IEIZone1 = SOMA_TMEV14_TimeZones.Zone1_IEI;
somaTmev14IEIZone2 = SOMA_TMEV14_TimeZones.Zone2_IEI;
somaTmev14IEIZone3 = SOMA_TMEV14_TimeZones.Zone3_IEI;

IEIBoxPlotScatter(somaPbsIEIZone1(:,2),somaTmev2IEIZone1(:,2),somaTmev5IEIZone1(:,2),somaTmev14IEIZone1(:,2));
IEIBoxPlotScatter(somaPbsIEIZone2(:,2),somaTmev2IEIZone2(:,2),somaTmev5IEIZone2(:,2),somaTmev14IEIZone2(:,2));
IEIBoxPlotScatter(somaPbsIEIZone3(:,2),somaTmev2IEIZone3(:,2),somaTmev5IEIZone3(:,2),somaTmev14IEIZone3(:,2));
figure()
plot(somaPbsIEIZone1(:,1),somaPbsIEIZone1(:,2),'o')
hold on
plot(somaPbsIEIZone2(:,1),somaPbsIEIZone2(:,2),'go')
plot(somaPbsIEIZone3(:,1),somaPbsIEIZone3(:,2),'ro')

figure()
plot(somaTmev2IEIZone1(:,1),somaTmev2IEIZone1(:,2),'o')
hold on
plot(somaTmev2IEIZone2(:,1),somaTmev2IEIZone2(:,2),'go')
plot(somaTmev2IEIZone3(:,1),somaTmev2IEIZone3(:,2),'ro')

%% PROCESSES AUC Zones Plotting
limits = [0 1000];
PlotZonesAUC(PROCESSES_PBS_TimeZones, 'PROCESSES PBS',limits);
PlotZonesAUC(PROCESSES_TMEV2_TimeZones, 'PROCESSES TMEV 2 DPI',limits);
PlotZonesAUC(PROCESSES_TMEV5_TimeZones, 'PROCESSES TMEV 5 DPI',limits);
PlotZonesAUC(PROCESSES_TMEV14_TimeZones, 'PROCESSES TMEV 15 DPI',limits);

%% PROCESSES Peak Durations Zones Plotting
limits = [0 300];
PlotZonesPeakDurs(PROCESSES_PBS_TimeZones, 'PROCESSES PBS',limits);
PlotZonesPeakDurs(PROCESSES_TMEV2_TimeZones, 'PROCESSES TMEV 2 DPI',limits);
PlotZonesPeakDurs(PROCESSES_TMEV5_TimeZones, 'PROCESSES TMEV 5 DPI',limits);
PlotZonesPeakDurs(PROCESSES_TMEV14_TimeZones, 'PROCESSES TMEV 15 DPI',limits);

%% PROCESSES Peak Magnitudes Zones Plotting
limits = [0 7];
PlotZonesPeakMags(PROCESSES_PBS_TimeZones, 'PROCESSES PBS',limits);
PlotZonesPeakMags(PROCESSES_TMEV2_TimeZones, 'PROCESSES TMEV 2 DPI',limits);
PlotZonesPeakMags(PROCESSES_TMEV5_TimeZones, 'PROCESSES TMEV 5 DPI',limits);
PlotZonesPeakMags(PROCESSES_TMEV14_TimeZones, 'PROCESSES TMEV 14 DPI',limits);

%% PROCESSES Peak AUC vs Amp Zones Plotting
limits = [0 45];
PlotZonesSqrtAUCVsAmp(PROCESSES_PBS_TimeZones, 'PROCESSES PBS',limits);
PlotZonesSqrtAUCVsAmp(PROCESSES_TMEV2_TimeZones, 'PROCESSES TMEV 2 DPI',limits);
PlotZonesSqrtAUCVsAmp(PROCESSES_TMEV5_TimeZones, 'PROCESSES TMEV 5 DPI',limits);
PlotZonesSqrtAUCVsAmp(PROCESSES_TMEV14_TimeZones, 'PROCESSES TMEV 15 DPI',limits);

%% PROCESSES Peak Duration vs Amp Zones Plotting
limits = [0 300];
PlotZonesDurationVsAmp(PROCESSES_PBS_TimeZones, 'Processes PBS',limits);
PlotZonesDurationVsAmp(PROCESSES_TMEV2_TimeZones, 'Processes TMEV 2 DPI',limits);
PlotZonesDurationVsAmp(PROCESSES_TMEV5_TimeZones, 'Processes TMEV 5 DPI',limits);
PlotZonesDurationVsAmp(PROCESSES_TMEV14_TimeZones, 'Processes TMEV 15 DPI',limits);

%% PROCESSES Peak AUC vs Duration
varNames = {'Duration','AUC','AUC pred'};
yMax = 400;
PROCESSESPbsDursZone1 = [PROCESSES_PBS_TimeZones.Zone1_PeakDurs]; %PROCESSES_PBS_TimeZones.Zone2_PeakDurs; PROCESSES_PBS_TimeZones.Zone3_PeakDurs];
PROCESSESPbsAucZone1 = [PROCESSES_PBS_TimeZones.Zone1_AUC]; %PROCESSES_PBS_TimeZones.Zone2_AUC; PROCESSES_PBS_TimeZones.Zone3_AUC];
PROCESSESPbsPeakMagsZone1 = [PROCESSES_PBS_TimeZones.Zone1_PeakMags]; 
[Rsq, PBS1fitLine] = LinearAnalysis(PROCESSESPbsDursZone1, log10(PROCESSESPbsAucZone1),'PBS PROCESSES Auc vs dur Zone 1', varNames,3,250);

PROCESSESTmev2DursZone1 = [PROCESSES_TMEV2_TimeZones.Zone1_PeakDurs]; %PROCESSES_TMEV2_TimeZones.Zone2_PeakDurs; PROCESSES_TMEV2_TimeZones.Zone3_PeakDurs];
PROCESSESTmev2AucZone1 = [PROCESSES_TMEV2_TimeZones.Zone1_AUC]; %PROCESSES_TMEV2_TimeZones.Zone2_AUC; PROCESSES_TMEV2_TimeZones.Zone3_AUC];
[Rsq, TMEV2P1fitLine] = LinearAnalysis(PROCESSESTmev2DursZone1, PROCESSESTmev2AucZone1,'TMEV 2 PROCESSES Zone 1', varNames,yMax);

PROCESSESTmev5DursZone1 = [PROCESSES_TMEV5_TimeZones.Zone1_PeakDurs]; %PROCESSES_TMEV5_TimeZones.Zone2_PeakDurs; PROCESSES_TMEV5_TimeZones.Zone3_PeakDurs];
PROCESSESTmev5AucZone1 = [PROCESSES_TMEV5_TimeZones.Zone1_AUC]; %PROCESSES_TMEV5_TimeZones.Zone2_AUC; PROCESSES_TMEV5_TimeZones.Zone3_AUC];
[Rsq, TMEV5P1fitLine] = LinearAnalysis(PROCESSESTmev5DursZone1, PROCESSESTmev5AucZone1,'TMEV 5 PROCESSES Zone 1', varNames,yMax);

PROCESSESTmev15DursZone1 = [PROCESSES_TMEV14_TimeZones.Zone1_PeakDurs]; %PROCESSES_TMEV14_TimeZones.Zone2_PeakDurs; PROCESSES_TMEV14_TimeZones.Zone3_PeakDurs];
PROCESSESTmev15AucZone1 = [PROCESSES_TMEV14_TimeZones.Zone1_AUC]; %PROCESSES_TMEV14_TimeZones.Zone2_AUC; PROCESSES_TMEV14_TimeZones.Zone3_AUC];
[Rsq, TMEV15P1fitLine] = LinearAnalysis(PROCESSESTmev15DursZone1, PROCESSESTmev15AucZone1,'TMEV 15 PROCESSES Auc vs dur Zone 1', varNames,1000,400);

PROCESSESPbsDursZone2 = [PROCESSES_PBS_TimeZones.Zone2_PeakDurs]; 
PROCESSESPbsAucZone2 = [PROCESSES_PBS_TimeZones.Zone2_AUC]; 
PROCESSESPbsPeakMagsZone2 = [PROCESSES_PBS_TimeZones.Zone2_PeakMags]; 
[Rsq, PBS2fitLine] = LinearAnalysis(PROCESSESPbsDursZone2, PROCESSESPbsAucZone2,'PBS PROCESSES Auc vs dur Zone 2', varNames,1000,400);

PROCESSESTmev2DursZone2 = [PROCESSES_TMEV2_TimeZones.Zone2_PeakDurs]; 
PROCESSESTmev2AucZone2 = [PROCESSES_TMEV2_TimeZones.Zone2_AUC];
PROCESSESTmev2PeakMagsZone2 = [PROCESSES_TMEV2_TimeZones.Zone2_PeakMags];
[Rsq, TMEV2P2fitLine] = LinearAnalysis(PROCESSESTmev2DursZone2, PROCESSESTmev2AucZone2,'TMEV 2 PROCESSES Zone 2', varNames,yMax);

PROCESSESTmev5DursZone2 = [PROCESSES_TMEV5_TimeZones.Zone2_PeakDurs]; 
PROCESSESTmev5AucZone2 = [PROCESSES_TMEV5_TimeZones.Zone2_AUC]; 
[Rsq, TMEV5P2fitLine] = LinearAnalysis(PROCESSESTmev5DursZone2, PROCESSESTmev5AucZone2,'TMEV 5 PROCESSES Zone 2', varNames,yMax);

PROCESSESTmev15DursZone2 = [PROCESSES_TMEV14_TimeZones.Zone2_PeakDurs]; 
PROCESSESTmev15AucZone2 = [PROCESSES_TMEV14_TimeZones.Zone2_AUC];
PROCESSESTmev15PeakMagsZone2 = [PROCESSES_TMEV14_TimeZones.Zone2_PeakMags];
[Rsq, TMEV15P2fitLine] = LinearAnalysis(PROCESSESTmev15DursZone2, PROCESSESTmev15AucZone2,'TMEV 15 PROCESSES Auc vs dur Zone 2', varNames,1000,400);

PROCESSESPbsDursZone3 = [PROCESSES_PBS_TimeZones.Zone3_PeakDurs]; 
PROCESSESPbsAucZone3 = [PROCESSES_PBS_TimeZones.Zone3_AUC];
PROCESSESPbsPeakMagsZone3 = [PROCESSES_PBS_TimeZones.Zone3_PeakMags];
[Rsq, PBS3fitLine] = LinearAnalysis(PROCESSESPbsDursZone3, PROCESSESPbsAucZone3,'PBS PROCESSES Auc vs dur Zone 3', varNames,1000,400);

PROCESSESTmev2DursZone3 = [PROCESSES_TMEV2_TimeZones.Zone3_PeakDurs]; 
PROCESSESTmev2AucZone3 = [PROCESSES_TMEV2_TimeZones.Zone3_AUC];
PROCESSESTmev2PeakMagsZone3 = [PROCESSES_TMEV2_TimeZones.Zone3_PeakMags];
[Rsq, TMEV2P3fitLine] = LinearAnalysis(PROCESSESTmev2DursZone3, PROCESSESTmev2AucZone3,'TMEV 2 PROCESSES Zone 3', varNames,yMax);

PROCESSESTmev5DursZone3 = [PROCESSES_TMEV5_TimeZones.Zone3_PeakDurs]; 
PROCESSESTmev5AucZone3 = [PROCESSES_TMEV5_TimeZones.Zone3_AUC]; 
[Rsq, TMEV5P3fitLine] = LinearAnalysis(PROCESSESTmev5DursZone3, PROCESSESTmev5AucZone3,'TMEV 5 PROCESSES Zone 3', varNames,yMax);

PROCESSESTmev15DursZone3 = [PROCESSES_TMEV14_TimeZones.Zone3_PeakDurs]; 
PROCESSESTmev15AucZone3 = [PROCESSES_TMEV14_TimeZones.Zone3_AUC]; 
PROCESSESTmev15PeakMagsZone3 = [PROCESSES_TMEV14_TimeZones.Zone3_PeakMags];
[Rsq, TMEV15P3fitLine] = LinearAnalysis(PROCESSESTmev15DursZone3, PROCESSESTmev15AucZone3,'TMEV 15 PROCESSES Auc vs dur Zone 3', varNames,1000,400);

h = CustomFigure();
hold on
plot(PROCESSESPbsDursZone2, PBS2fitLine,'-b')
plot(PROCESSESPbsDursZone3, PBS3fitLine,'-r')
plot(PROCESSESPbsDursZone2, PROCESSESPbsAucZone2,'ob')
plot(PROCESSESPbsDursZone3, PROCESSESPbsAucZone3,'or')
xlabel('Duration (s)'); ylabel('AUC (dF/F)'); title('Processes PBS Phase 2 vs Phase 3'); %legend('Pbs Phase 2','Pbs Phase 3');
grid on

h = CustomFigure();
hold on
plot(PROCESSESTmev15DursZone2, TMEV15P2fitLine,'-b')
plot(PROCESSESTmev15DursZone3, TMEV15P3fitLine,'-r')
plot(PROCESSESTmev15DursZone2, PROCESSESTmev15AucZone2,'ob')
plot(PROCESSESTmev15DursZone3, PROCESSESTmev15AucZone3,'or')
xlabel('Duration (s)'); ylabel('AUC (dF/F)'); title('Processes Tmev15 Phase 2 vs Phase 3'); %legend('Tmev15 Phase 2','Tmev15 Phase 3');
grid on

figure()
hold on
plot(PROCESSESTmev2DursZone2, TMEV2P2fitLine,'-b')
plot(PROCESSESTmev2DursZone3, TMEV2P3fitLine,'-r')
legend('Tmev2 Phase 2','Tmev2 Phase 3'); xlabel('Duration (s)'); ylabel('AUC (dF/F)'); title('Soma Tmev2 Phase 2 vs Phase 3');
grid on

% % 3D plots: AUC vs Amp vs Duration
figure()
hold on
plot3(PROCESSESPbsDursZone2,PROCESSESPbsAucZone2,PROCESSESPbsPeakMagsZone2,'o',PROCESSESPbsDursZone3,PROCESSESPbsAucZone3,PROCESSESPbsPeakMagsZone3,'or');
% plot3(PROCESSESPbsDursZone3,PROCESSESPbsAucZone3,PROCESSESPbsPeakMagsZone3,'or');
xlabel('Duration (s)'); ylabel('AUC'); zlabel('Amplitude'); title('Processes Pbs Phase 2');
legend('Phase 2','Phase 3')
grid on

SpherePlot(PROCESSESPbsDursZone2,PROCESSESPbsAucZone2,PROCESSESPbsPeakMagsZone2)
figure()
hold on
plot3(PROCESSESTmev15DursZone2,PROCESSESTmev15AucZone2,PROCESSESTmev15PeakMagsZone2,'o',PROCESSESTmev15DursZone3,PROCESSESTmev15AucZone3,PROCESSESTmev15PeakMagsZone3,'or');
% plot3(PROCESSESTmev15DursZone3,PROCESSESTmev15AucZone3,PROCESSESTmev15PeakMagsZone3,'or');
xlabel('Duration (s)'); ylabel('AUC'); zlabel('Amplitude'); title('Processes Tmev15 Phase 2');
legend('Phase 2','Phase 3')
grid on


[X, Y] = meshgrid(PROCESSESPbsDursZone3,PROCESSESPbsAucZone3);
tri = delaunay(PROCESSESPbsDursZone3,PROCESSESPbsAucZone3);
trisurf(tri,PROCESSESPbsDursZone3,PROCESSESPbsAucZone3,PROCESSESPbsPeakMagsZone3);
l = light('Position',[-50 -15 29]);
set(gca,'CameraPosition',[208 -50 7687])
lighting phong
shading interp
colorbar EastOutside

figure()
plot3(PROCESSESPbsDursZone3,PROCESSESPbsAucZone3,PROCESSESPbsPeakMagsZone3,'o');
surf([PROCESSESPbsDursZone3,PROCESSESPbsAucZone3,PROCESSESPbsPeakMagsZone3])
hold on
figure()
FitEllipsoid([PROCESSESPbsDursZone3,PROCESSESPbsAucZone3,PROCESSESPbsPeakMagsZone3], [1,0,0], [1 0 0],[1,1,1]);
hold on
FitEllipsoid([PROCESSESPbsDursZone2,PROCESSESPbsAucZone2,PROCESSESPbsPeakMagsZone2], [0,0,1], [0,0,1],[0,0,0]);
xlabel('Duration (s)'); ylabel('AUC'); zlabel('Amplitude'); title('Processes Pbs Phase 3');


figure()
hold on
scatter(PROCESSESPbsDursZone3,PROCESSESPbsAucZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESPbsDursZone2,PROCESSESPbsAucZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESPbsDursZone1,PROCESSESPbsAucZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',0.2)
plot(PROCESSESPbsDursZone2,PBS2fitLine,'-b');
plot(PROCESSESPbsDursZone3,PBS3fitLine,'-r');
plot(PROCESSESPbsDursZone1,PBS1fitLine,'-','Color',[0.9 0.9 0]);
xlim([0 250]); ylim([0 1000]);
legend('phase 3','phase 2','phase 1'); xlabel('duration'); ylabel('AUC'); title('pbs processes');

figure()
hold on
scatter(PROCESSESTmev15DursZone3,PROCESSESTmev15AucZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESTmev15DursZone2,PROCESSESTmev15AucZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESTmev15DursZone1,PROCESSESTmev15AucZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',0.2)
plot(PROCESSESTmev15DursZone2,TMEV15P2fitLine,'-b');
plot(PROCESSESTmev15DursZone3,TMEV15P3fitLine,'-r');
plot(PROCESSESTmev15DursZone1,TMEV15P1fitLine,'-','Color',[0.9 0.9 0]);
xlim([0 250]); ylim([0 1000]);
legend('phase 3','phase 2','phase 1'); xlabel('duration'); ylabel('AUC'); title('Tmev15 processes');

figure()
hold on
plot(PROCESSESPbsPeakMagsZone3,PROCESSESPbsAucZone3,'or')
plot(PROCESSESPbsPeakMagsZone2,PROCESSESPbsAucZone2,'ob','MarkerFaceColor','blue')
plot(PROCESSESPbsPeakMagsZone1,PROCESSESPbsAucZone1,'ok','MarkerFaceColor',[0.5 0.5 0])
% plot(PROCESSESPbsDursZone2,PBS2fitLine,'-b');
% plot(PROCESSESPbsDursZone3,PBS3fitLine,'-r');
% plot(PROCESSESPbsDursZone1,PBS1fitLine,'-k');
legend('phase 3','phase 2','phase 1'); xlabel('amplitude'); ylabel('AUC'); title('pbs processes duration vs amplitude');
Z = [PROCESSESPbsDursZone2,PROCESSESPbsAucZone2,PROCESSESPbsPeakMagsZone2];
[coeff,score,latent,~,explained] = pca(Z);

figure()
plot3(PROCESSESTmev2DursZone2,PROCESSESTmev2AucZone2,PROCESSESTmev2PeakMagsZone2,'o');
xlabel('Duration (s)'); ylabel('AUC'); zlabel('Amplitude'); title('Processes TMEV 2 Phase 2');
grid on

figure()
plot3(PROCESSESTmev2DursZone3,PROCESSESTmev2AucZone3,PROCESSESTmev2PeakMagsZone3,'o');
xlabel('Duration (s)'); ylabel('AUC'); zlabel('Amplitude'); title('Processes TMEV 2 Phase 3');
grid on

%% PROCESSES Plot RMS Data
ScatterRMS(PROCESSES_PBS_allRMS,'PBS');
ScatterRMS(PROCESSES_TMEV2_allRMS,'TMEV 2 DPI');
ScatterRMS(PROCESSES_TMEV5_allRMS,'TMEV 5 DPI');
ScatterRMS(PROCESSES_TMEV14_allRMS,'TMEV 15 DPI');

%% PROCESSES Plot Signal Length Data
PairedSignalLength(PBS_PROCESSES_DATA,'PBS');
PairedSignalLength(PROCESSES_TMEV2_DATA,'TMEV 2 DPI');
PairedSignalLength(PROCESSES_TMEV5_DATA,'TMEV 5 DPI');
PairedSignalLength(PROCESSES_TMEV14_DATA,'TMEV 15 DPI');

%% PROCESSES Events Dot Plot
yMax = numel(fieldnames(PROCESSES_PBS));
time = round((1740/0.97)/60);
xlimit = [-time time];
EventDotPlot(PROCESSES_PBS, 'Processes PBS Phase 2-3', xlimit, yMax);
EventDotPlot(PROCESSES_TMEV2, 'Processes TMEV 2 DPI Phase 2-3', xlimit, yMax);
EventDotPlot(PROCESSES_TMEV5, 'Processes TMEV 5 DPI Phase 2-3', xlimit, yMax);
EventDotPlot(PROCESSES_TMEV14, 'Processes TMEV 15 DPI Phase 2-3', xlimit, yMax);

%% PROCESSES Events Phase 0-1 Dot Plot
yMax = numel(fieldnames(PROCESSES_PBS));
time = round((1740/0.97)/60);
xlimit = [-1 0.5];
Phase01DotPlot(PROCESSES_PBS, 'Processes PBS Phase 0-1', xlimit, yMax);
Phase01DotPlot(PROCESSES_TMEV2, 'Processes TMEV 2 DPI Phase 0-1', xlimit, yMax);
Phase01DotPlot(PROCESSES_TMEV5, 'Processes TMEV 5 DPI Phase 0-1', xlimit, yMax);
Phase01DotPlot(PROCESSES_TMEV14, 'Processes TMEV 15 DPI Phase 0-1', xlimit, yMax);

%% PROCESSES Plot Sum of Approx Deriv Data
xDim = 1.5;
yDim = 2.5;
yMax = 0.03;
[processPBSzones, PROCESS_PBS_PVALUES, processPbsVideos] = PairedDeriv(PBS_PROCESSES_DATA,'Process PBS', xDim, yDim, [0 0 0], yMax);
[processTMEV2zones, PROCESS_TMEV2_PVALUES, processTmev2Videos] = PairedDeriv(PROCESSES_TMEV2_DATA,'Process TMEV 2 DPI', xDim, yDim, [1 0 0], yMax);
[processTMEV5zones, PROCESS_TMEV5_PVALUES, processTmev5Videos] = PairedDeriv(PROCESSES_TMEV5_DATA,'Process TMEV 5 DPI', xDim, yDim, [1 0 0], yMax);
[processTMEV15zones, PROCESS_TMEV15_PVALUES, processTmev15Videos] = PairedDeriv(PROCESSES_TMEV14_DATA,'Process TMEV 15 DPI', xDim, yDim, [1 0 0], yMax);

%% PROCESSES Plot All zones Deriv Data
xDim = 2.5;
yDim = 2.5;
yMax = 0.03;
PROCESS_P = AllDerivZones(PBS_PROCESSES_DATA,PROCESSES_TMEV2_DATA,PROCESSES_TMEV5_DATA,PROCESSES_TMEV14_DATA,'Processes', yMax, xDim, yDim);

%% PROCESSES Plot Interevent Interval %%
processPbsIEIZone1 = PROCESSES_PBS_TimeZones.Zone1_IEI;
processPbsIEIZone2 = PROCESSES_PBS_TimeZones.Zone2_IEI;
processPbsIEIZone3 = PROCESSES_PBS_TimeZones.Zone3_IEI;

processTmev2IEIZone1 = PROCESSES_TMEV2_TimeZones.Zone1_IEI;
processTmev2IEIZone2 = PROCESSES_TMEV2_TimeZones.Zone2_IEI;
processTmev2IEIZone3 = PROCESSES_TMEV2_TimeZones.Zone3_IEI;

processTmev5IEIZone1 = PROCESSES_TMEV5_TimeZones.Zone1_IEI;
processTmev5IEIZone2 = PROCESSES_TMEV5_TimeZones.Zone2_IEI;
processTmev5IEIZone3 = PROCESSES_TMEV5_TimeZones.Zone3_IEI;

processTmev14IEIZone1 = PROCESSES_TMEV14_TimeZones.Zone1_IEI;
processTmev14IEIZone2 = PROCESSES_TMEV14_TimeZones.Zone2_IEI;
processTmev14IEIZone3 = PROCESSES_TMEV14_TimeZones.Zone3_IEI;

IEIBoxPlotScatter(processPbsIEIZone1(:,2),processTmev2IEIZone1(:,2),processTmev5IEIZone1(:,2),processTmev14IEIZone1(:,2));
IEIBoxPlotScatter(processPbsIEIZone2(:,2),processTmev2IEIZone2(:,2),processTmev5IEIZone2(:,2),processTmev14IEIZone2(:,2));
IEIBoxPlotScatter(processPbsIEIZone3(:,2),processTmev2IEIZone3(:,2),processTmev5IEIZone3(:,2),processTmev14IEIZone3(:,2));

figure()
plot(processPbsIEIZone1(:,1),processPbsIEIZone1(:,2),'o')
hold on
plot(processPbsIEIZone2(:,1),processPbsIEIZone2(:,2),'go')
plot(processPbsIEIZone3(:,1),processPbsIEIZone3(:,2),'ro')

figure()
plot(processTmev2IEIZone1(:,1),processTmev2IEIZone1(:,2),'o')
hold on
plot(processTmev2IEIZone2(:,1),processTmev2IEIZone2(:,2),'go')
plot(processTmev2IEIZone3(:,1),processTmev2IEIZone3(:,2),'ro')

%% PROCESSES Grouped Bar Charts
barTitle = 'PROCESSES';
xDim = 1.5;
yDim = 1.5;
PairedBarGraph(barTitle, xDim, yDim, PROCESSES_PBS_TimeZones, PROCESSES_TMEV2_TimeZones, PROCESSES_TMEV5_TimeZones, PROCESSES_TMEV14_TimeZones);

%% PROCESSES UnGrouped Bar Charts
barTitle = 'PROCESSES';
xDim = 1.0;
yDim = 1.5;
SingleZoneBar('Processes', xDim, yDim,PROCESSES_PBS_TimeZones, [0.2, 0.2, 0.2])
SingleZoneBar('Processes 2 DPI', xDim, yDim,PROCESSES_TMEV2_TimeZones, [1 0 0])
SingleZoneBar('Processes 5 DPI', xDim, yDim,PROCESSES_TMEV5_TimeZones, [1 0 0])
SingleZoneBar('Processes 15 DPI', xDim, yDim,PROCESSES_TMEV14_TimeZones, [1 0 0])

%% PROCESSES Plot RMS Data
ScatterRMS(PROCESSES_PBS_allRMS,'PBS', [0.2 0.2 0.2]); % use T2114_1 as example
ScatterRMS(PROCESSES_TMEV2_allRMS,'TMEV 2 DPI', [1 0 0]); % T2190_10
ScatterRMS(PROCESSES_TMEV5_allRMS,'TMEV 5 DPI', [1 0 0]); % T2102_12
ScatterRMS(PROCESSES_TMEV14_allRMS,'TMEV 15 DPI', [1 0 0]); % T2326_1

%% PROCESSES Zone derivs
barTitle = 'Processes';
xDim = 2;
yDim = 1.5;
PairedZonesDeriv(barTitle, xDim, yDim, PBS_PROCESSES_DATA, PROCESSES_TMEV2_DATA, PROCESSES_TMEV5_DATA, PROCESSES_TMEV14_DATA);

%% PROCESSES Linear Regression AUC vs Amplitude
varNames = {'Peak Amplitude','AUC','auc pred'};
PROCESSESPbsAmpZone1 = [PROCESSES_PBS_TimeZones.Zone1_PeakMags]; %PROCESSES_PBS_TimeZones.Zone2_PeakMags; PROCESSES_PBS_TimeZones.Zone3_PeakMags];
PROCESSESPbsAucZone1 = [PROCESSES_PBS_TimeZones.Zone1_AUC]; %PROCESSES_PBS_TimeZones.Zone2_AUC; PROCESSES_PBS_TimeZones.Zone3_AUC];
[Rsq, PBS1fitLine] = LinearAnalysis(PROCESSESPbsAmpZone1, PROCESSESPbsAucZone1,'PBS PROCESSES Zone 1 AUC vs Amp',varNames,1000,4);

PROCESSESPbsAmpZone2 = [PROCESSES_PBS_TimeZones.Zone2_PeakMags];
PROCESSESPbsAucZone2 = [PROCESSES_PBS_TimeZones.Zone2_AUC];
[Rsq, PBS2fitLine] = LinearAnalysis(PROCESSESPbsAmpZone2, PROCESSESPbsAucZone2,'PBS PROCESSES Zone 2 AUC vs Amp',varNames,1000,4);

PROCESSESPbsAmpZone3 = [PROCESSES_PBS_TimeZones.Zone3_PeakMags];
PROCESSESPbsAucZone3 = [PROCESSES_PBS_TimeZones.Zone3_AUC]; 
[Rsq, PBS3fitLine] = LinearAnalysis(PROCESSESPbsAmpZone3, PROCESSESPbsAucZone3,'PBS PROCESSES Zone 3 AUC vs Amp',varNames,1000,4);

PROCESSESTmev2AmpZone1 = [PROCESSES_TMEV2_TimeZones.Zone1_PeakMags];
PROCESSESTmev2AucZone1 = [PROCESSES_TMEV2_TimeZones.Zone1_AUC];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev2AmpZone1, PROCESSESTmev2AucZone1,'TMEV 2 PROCESSES Zone 1 AUC vs Amp',varNames);

PROCESSESTmev2AmpZone2 = [PROCESSES_TMEV2_TimeZones.Zone2_PeakMags];
PROCESSESTmev2AucZone2 = [PROCESSES_TMEV2_TimeZones.Zone2_AUC];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev2AmpZone2, PROCESSESTmev2AucZone2,'TMEV 2 PROCESSES Zone 2',varNames);

PROCESSESTmev2AmpZone3 = [PROCESSES_TMEV2_TimeZones.Zone3_PeakMags];
PROCESSESTmev2AucZone3 = [PROCESSES_TMEV2_TimeZones.Zone3_AUC];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev2AmpZone3, PROCESSESTmev2AucZone3,'TMEV 2 PROCESSES Zone 3',varNames);

PROCESSESTmev5AmpZone1 = [PROCESSES_TMEV5_TimeZones.Zone1_PeakMags];
PROCESSESTmev5AucZone1 = [PROCESSES_TMEV5_TimeZones.Zone1_AUC];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev5AmpZone1, PROCESSESTmev5AucZone1,'TMEV 5 PROCESSES Zone 1',varNames);

PROCESSESTmev5AmpZone2 = [PROCESSES_TMEV5_TimeZones.Zone2_PeakMags];
PROCESSESTmev5AucZone2 = [PROCESSES_TMEV5_TimeZones.Zone2_AUC];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev5AmpZone2, PROCESSESTmev5AucZone2,'TMEV 5 PROCESSES Zone 2',varNames);

PROCESSESTmev5AmpZone3 = [PROCESSES_TMEV5_TimeZones.Zone3_PeakMags];
PROCESSESTmev5AucZone3 = [PROCESSES_TMEV5_TimeZones.Zone3_AUC];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev5AmpZone3, PROCESSESTmev5AucZone3,'TMEV 5 PROCESSES Zone 3',varNames);

PROCESSESTmev15AmpZone1 = [PROCESSES_TMEV14_TimeZones.Zone1_PeakMags];
PROCESSESTmev15AucZone1 = [PROCESSES_TMEV14_TimeZones.Zone1_AUC];
[Rsq, TMEV15P1fitLine] = LinearAnalysis(PROCESSESTmev15AmpZone1, PROCESSESTmev15AucZone1,'TMEV 15 PROCESSES Zone 1 AUC vs Amp',varNames,1000,4);

PROCESSESTmev15AmpZone2 = [PROCESSES_TMEV14_TimeZones.Zone2_PeakMags];
PROCESSESTmev15AucZone2 = [PROCESSES_TMEV14_TimeZones.Zone2_AUC];
[Rsq, TMEV15P2fitLine] = LinearAnalysis(PROCESSESTmev15AmpZone2, PROCESSESTmev15AucZone2,'TMEV 15 PROCESSES Zone 2 AUC vs Amp',varNames,1000,4);

PROCESSESTmev15AmpZone3 = [PROCESSES_TMEV14_TimeZones.Zone3_PeakMags];
PROCESSESTmev15AucZone3 = [PROCESSES_TMEV14_TimeZones.Zone3_AUC];
[Rsq, TMEV15P3fitLine] = LinearAnalysis(PROCESSESTmev15AmpZone3, PROCESSESTmev15AucZone3,'TMEV 15 PROCESSES Zone 3 AUC vs Amp',varNames,1000,4);

figure()
hold on
plot3(log10(PROCESSESPbsAmpZone1),log10(PROCESSESPbsAucZone1),log10(PROCESSESPbsDursZone1),'o')
plot3(log10(PROCESSESPbsAmpZone2),log10(PROCESSESPbsAucZone2),log10(PROCESSESPbsDursZone2),'o')
plot3(log10(PROCESSESPbsAmpZone3),log10(PROCESSESPbsAucZone3),log10(PROCESSESPbsDursZone3),'o')
xlabel('amplitude'); ylabel('auc'); zlabel('duration'); title('Log10 Pbs Processes');
legend('phase1','phase2','phase3');
% xlim([0 4]); ylim([0 1000]); zlim([0 250]);

figure()
hold on
plot3(log10(PROCESSESTmev15AmpZone1),log10(PROCESSESTmev15AucZone1),log10(PROCESSESTmev15DursZone1),'o')
plot3(log10(PROCESSESTmev15AmpZone2),log10(PROCESSESTmev15AucZone2),log10(PROCESSESTmev15DursZone2),'o')
plot3(log10(PROCESSESTmev15AmpZone3),log10(PROCESSESTmev15AucZone3),log10(PROCESSESTmev15DursZone3),'o')
xlabel('amplitude'); ylabel('auc'); zlabel('duration'); title('Log10 Tmev 15 Processes');
legend('phase1','phase2','phase3');
% xlim([0 4]); ylim([0 1000]); zlim([0 250]);

figure()
hold on
scatter(PROCESSESPbsAmpZone3,PROCESSESPbsAucZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESPbsAmpZone2,PROCESSESPbsAucZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESPbsAmpZone1,PROCESSESPbsAucZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.2)
plot(PROCESSESPbsAmpZone2,PBS2fitLine,'-b');
plot(PROCESSESPbsAmpZone3,PBS3fitLine,'-r');
plot(PROCESSESPbsAmpZone1,PBS1fitLine,'-','Color',[0.9 0.9 0]);
legend('phase 3','phase 2','phase 1'); xlabel('amplitude'); ylabel('AUC'); title('pbs PROCESSES');
xlim([0 4]); ylim([0 1000]);
figure()
hold on
scatter(PROCESSESTmev15AmpZone3,PROCESSESTmev15AucZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESTmev15AmpZone2,PROCESSESTmev15AucZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESTmev15AmpZone1,PROCESSESTmev15AucZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.2)
plot(PROCESSESTmev15AmpZone2,TMEV15P2fitLine,'-b');
plot(PROCESSESTmev15AmpZone3,TMEV15P3fitLine,'-r');
plot(PROCESSESTmev15AmpZone1,TMEV15P1fitLine,'-','Color',[0.9 0.9 0]);
legend('phase 3','phase 2','phase 1'); xlabel('amplitude'); ylabel('AUC'); title('Tmev15 PROCESSES');
xlim([0 4]); ylim([0 1000]);

% soma and processes
pbsPhase2Auc = [somaPbsAucZone2; PROCESSESPbsAucZone2];
tmev15Phase2Auc = [somaTmev15AucZone2; PROCESSESTmev15AucZone2];
pbsPhase2Amp = [somaPbsMagsZone2; PROCESSESPbsAmpZone2];
tmev15Phase2Amp = [somaTmev15MagsZone2; PROCESSESTmev15AmpZone2];
[Rsq, PBSP2fitLine] = LinearAnalysis(pbsPhase2Amp, pbsPhase2Auc,'PBS Zone 2 (processes + soma) AUC vs Amp',varNames,200);
[Rsq, PBSP2fitLine] = LinearAnalysis(tmev15Phase2Amp, tmev15Phase2Auc,'TMEV 15 Zone 2 (processes + soma) AUC vs Amp',varNames,200);

%% PROCESSES Linear Regression AUC vs Duration
varNames = {'Peak Durslitude','AUC','auc pred'};
PROCESSESPbsDursZone1 = [PROCESSES_PBS_TimeZones.Zone1_PeakDurs]; %PROCESSES_PBS_TimeZones.Zone2_PeakDurs; PROCESSES_PBS_TimeZones.Zone3_PeakDurs];
PROCESSESPbsAucZone1 = [PROCESSES_PBS_TimeZones.Zone1_AUC]; %PROCESSES_PBS_TimeZones.Zone2_AUC; PROCESSES_PBS_TimeZones.Zone3_AUC];
[Rsq, PBS1fitLine] = LinearAnalysis(PROCESSESPbsDursZone1, PROCESSESPbsAucZone1,'PBS PROCESSES Zone 1 AUC vs Durs',varNames,1000,250);

PROCESSESPbsDursZone2 = [PROCESSES_PBS_TimeZones.Zone2_PeakDurs];
PROCESSESPbsAucZone2 = [PROCESSES_PBS_TimeZones.Zone2_AUC];
[Rsq, PBS2fitLine] = LinearAnalysis(PROCESSESPbsDursZone2, PROCESSESPbsAucZone2,'PBS PROCESSES Zone 2 AUC vs Durs',varNames,1000,250);

PROCESSESPbsDursZone3 = [PROCESSES_PBS_TimeZones.Zone3_PeakDurs];
PROCESSESPbsAucZone3 = [PROCESSES_PBS_TimeZones.Zone3_AUC]; 
[Rsq, PBS3fitLine] = LinearAnalysis(PROCESSESPbsDursZone3, PROCESSESPbsAucZone3,'PBS PROCESSES Zone 3 AUC vs Durs',varNames,1000,250);

PROCESSESTmev2DursZone1 = [PROCESSES_TMEV2_TimeZones.Zone1_PeakDurs];
PROCESSESTmev2AucZone1 = [PROCESSES_TMEV2_TimeZones.Zone1_AUC];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev2DursZone1, PROCESSESTmev2AucZone1,'TMEV 2 PROCESSES Zone 1 AUC vs Durs',varNames);

PROCESSESTmev2DursZone2 = [PROCESSES_TMEV2_TimeZones.Zone2_PeakDurs];
PROCESSESTmev2AucZone2 = [PROCESSES_TMEV2_TimeZones.Zone2_AUC];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev2DursZone2, PROCESSESTmev2AucZone2,'TMEV 2 PROCESSES Zone 2',varNames);

PROCESSESTmev2DursZone3 = [PROCESSES_TMEV2_TimeZones.Zone3_PeakDurs];
PROCESSESTmev2AucZone3 = [PROCESSES_TMEV2_TimeZones.Zone3_AUC];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev2DursZone3, PROCESSESTmev2AucZone3,'TMEV 2 PROCESSES Zone 3',varNames);

PROCESSESTmev5DursZone1 = [PROCESSES_TMEV5_TimeZones.Zone1_PeakDurs];
PROCESSESTmev5AucZone1 = [PROCESSES_TMEV5_TimeZones.Zone1_AUC];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev5DursZone1, PROCESSESTmev5AucZone1,'TMEV 5 PROCESSES Zone 1',varNames);

PROCESSESTmev5DursZone2 = [PROCESSES_TMEV5_TimeZones.Zone2_PeakDurs];
PROCESSESTmev5AucZone2 = [PROCESSES_TMEV5_TimeZones.Zone2_AUC];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev5DursZone2, PROCESSESTmev5AucZone2,'TMEV 5 PROCESSES Zone 2',varNames);

PROCESSESTmev5DursZone3 = [PROCESSES_TMEV5_TimeZones.Zone3_PeakDurs];
PROCESSESTmev5AucZone3 = [PROCESSES_TMEV5_TimeZones.Zone3_AUC];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev5DursZone3, PROCESSESTmev5AucZone3,'TMEV 5 PROCESSES Zone 3',varNames);

PROCESSESTmev15DursZone1 = [PROCESSES_TMEV14_TimeZones.Zone1_PeakDurs];
PROCESSESTmev15AucZone1 = [PROCESSES_TMEV14_TimeZones.Zone1_AUC];
[Rsq, TMEV15P1fitLine] = LinearAnalysis(PROCESSESTmev15DursZone1, PROCESSESTmev15AucZone1,'TMEV 15 PROCESSES Zone 1 AUC vs Durs',varNames,1000,250);

PROCESSESTmev15DursZone2 = [PROCESSES_TMEV14_TimeZones.Zone2_PeakDurs];
PROCESSESTmev15AucZone2 = [PROCESSES_TMEV14_TimeZones.Zone2_AUC];
[Rsq, TMEV15P2fitLine] = LinearAnalysis(PROCESSESTmev15DursZone2, PROCESSESTmev15AucZone2,'TMEV 15 PROCESSES Zone 2 AUC vs Durs',varNames,1000,250);

PROCESSESTmev15DursZone3 = [PROCESSES_TMEV14_TimeZones.Zone3_PeakDurs];
PROCESSESTmev15AucZone3 = [PROCESSES_TMEV14_TimeZones.Zone3_AUC];
[Rsq, TMEV15P3fitLine] = LinearAnalysis(PROCESSESTmev15DursZone3, PROCESSESTmev15AucZone3,'TMEV 15 PROCESSES Zone 3 AUC vs Durs',varNames,1000,250);

figure()
hold on
scatter(PROCESSESPbsDursZone3,PROCESSESPbsAucZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESPbsDursZone2,PROCESSESPbsAucZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESPbsDursZone1,PROCESSESPbsAucZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.2)
plot(PROCESSESPbsDursZone2,PBS2fitLine,'-b');
plot(PROCESSESPbsDursZone3,PBS3fitLine,'-r');
plot(PROCESSESPbsDursZone1,PBS1fitLine,'-','Color',[0.9 0.9 0]);
legend('phase 3','phase 2','phase 1'); xlabel('duration'); ylabel('AUC'); title('pbs PROCESSES');
xlim([0 250]); ylim([0 1000]);
figure()
hold on
scatter(PROCESSESTmev15DursZone3,PROCESSESTmev15AucZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESTmev15DursZone2,PROCESSESTmev15AucZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESTmev15DursZone1,PROCESSESTmev15AucZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.2)
plot(PROCESSESTmev15DursZone2,TMEV15P2fitLine,'-b');
plot(PROCESSESTmev15DursZone3,TMEV15P3fitLine,'-r');
plot(PROCESSESTmev15DursZone1,TMEV15P1fitLine,'-','Color',[0.9 0.9 0]);
legend('phase 3','phase 2','phase 1'); xlabel('duration'); ylabel('AUC'); title('Tmev15 PROCESSES');
xlim([0 250]); ylim([0 1000]);

%% PROCESSES Linear Regression Duration vs Amplitude
varNames = {'Peak Amplitude','Peak Duration','duration pred'};
PROCESSESPbsAmpZone1 = [PROCESSES_PBS_TimeZones.Zone1_PeakMags]; 
PROCESSESPbsPeakDursZone1 = [PROCESSES_PBS_TimeZones.Zone1_PeakDurs];
[Rsq, PBS1fitLine] = LinearAnalysis(PROCESSESPbsAmpZone1, PROCESSESPbsPeakDursZone1,'PBS PROCESSES Durs vs amp Zone 1',varNames,250,4);

PROCESSESPbsAmpZone2 = [PROCESSES_PBS_TimeZones.Zone2_PeakMags]; 
PROCESSESPbsPeakDursZone2 = [PROCESSES_PBS_TimeZones.Zone2_PeakDurs];
[Rsq, PBS2fitLine] = LinearAnalysis(PROCESSESPbsAmpZone2, PROCESSESPbsPeakDursZone2,'PBS PROCESSES Durs vs amp Zone 2',varNames,250,4);

PROCESSESPbsAmpZone3 = [PROCESSES_PBS_TimeZones.Zone3_PeakMags]; 
PROCESSESPbsPeakDursZone3 = [PROCESSES_PBS_TimeZones.Zone3_PeakDurs];
[Rsq, PBS3fitLine] = LinearAnalysis(log10(PROCESSESPbsAmpZone3), log10(PROCESSESPbsPeakDursZone3),'PBS PROCESSES Durs vs amp Zone 3',varNames,3,1);

PROCESSESTmev2AmpZone1 = [PROCESSES_TMEV2_TimeZones.Zone1_PeakMags]; 
PROCESSESTmev2PeakDursZone1 = [PROCESSES_TMEV2_TimeZones.Zone1_PeakDurs];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev2AmpZone1, PROCESSESTmev2PeakDursZone1,'Tmev2 PROCESSES Zone 1',varNames);

PROCESSESTmev2AmpZone2 = [PROCESSES_PBS_TimeZones.Zone2_PeakMags]; 
PROCESSESTmev2PeakDursZone2 = [PROCESSES_PBS_TimeZones.Zone2_PeakDurs];
[Rsq, fitLine] = LinearAnalysis(PROCESSESPbsAmpZone2, PROCESSESPbsPeakDursZone2,'Tmev 2 PROCESSES Zone 2',varNames);

PROCESSESTmev2AmpZone3 = [PROCESSES_PBS_TimeZones.Zone3_PeakMags]; 
PROCESSESTmev2PeakDursZone3 = [PROCESSES_PBS_TimeZones.Zone3_PeakDurs];
[Rsq, fitLine] = LinearAnalysis(PROCESSESPbsAmpZone3, PROCESSESPbsPeakDursZone3,'Tmev2 PROCESSES Zone 3',varNames);

PROCESSESTmev5AmpZone1 = [PROCESSES_TMEV5_TimeZones.Zone1_PeakMags];
PROCESSESTmev5PeakDursZone1 = [PROCESSES_TMEV5_TimeZones.Zone1_PeakDurs];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev5AmpZone1, PROCESSESTmev5PeakDursZone1,'TMEV 5 PROCESSES Zone 1',varNames);

PROCESSESTmev5AmpZone2 = [PROCESSES_TMEV5_TimeZones.Zone2_PeakMags];
PROCESSESTmev5PeakDursZone2 = [PROCESSES_TMEV5_TimeZones.Zone2_PeakDurs];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev5AmpZone2, PROCESSESTmev5PeakDursZone2,'TMEV 5 PROCESSES Zone 2',varNames);

PROCESSESTmev5AmpZone3 = [PROCESSES_TMEV5_TimeZones.Zone3_PeakMags];
PROCESSESTmev5PeakDursZone3 = [PROCESSES_TMEV5_TimeZones.Zone3_PeakDurs];
[Rsq, fitLine] = LinearAnalysis(PROCESSESTmev5AmpZone3, PROCESSESTmev5PeakDursZone3,'TMEV 5 PROCESSES Zone 3',varNames);

PROCESSESTmev15AmpZone1 = [PROCESSES_TMEV14_TimeZones.Zone1_PeakMags];
PROCESSESTmev15PeakDursZone1 = [PROCESSES_TMEV14_TimeZones.Zone1_PeakDurs];
[Rsq, TMEV15P1fitLine] = LinearAnalysis(PROCESSESTmev15AmpZone1, PROCESSESTmev15PeakDursZone1,'TMEV 15 PROCESSES Durs vs amp Zone 1',varNames,250,4);

PROCESSESTmev15AmpZone2 = [PROCESSES_TMEV14_TimeZones.Zone2_PeakMags];
PROCESSESTmev15PeakDursZone2 = [PROCESSES_TMEV14_TimeZones.Zone2_PeakDurs];
[Rsq, TMEV15P2fitLine] = LinearAnalysis(PROCESSESTmev15AmpZone2, PROCESSESTmev15PeakDursZone2,'TMEV 15 PROCESSES Durs vs amp Zone 2',varNames,250,4);

PROCESSESTmev15AmpZone3 = [PROCESSES_TMEV14_TimeZones.Zone3_PeakMags];
PROCESSESTmev15PeakDursZone3 = [PROCESSES_TMEV14_TimeZones.Zone3_PeakDurs];
[Rsq, TMEV15P3fitLine] = LinearAnalysis(log10(PROCESSESTmev15AmpZone3), log10(PROCESSESTmev15PeakDursZone3),'TMEV 15 PROCESSES Durs vs amp Zone 3',varNames,3,1);

figure()
hold on
scatter(PROCESSESPbsAmpZone3,PROCESSESPbsPeakDursZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESPbsAmpZone2,PROCESSESPbsPeakDursZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESPbsAmpZone1,PROCESSESPbsPeakDursZone1,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.2)
plot(PROCESSESPbsAmpZone2,PBS2fitLine,'-b');
plot(PROCESSESPbsAmpZone3,PBS3fitLine,'-r');
plot(PROCESSESPbsAmpZone1,PBS1fitLine,'-','Color',[0.9 0.9 0]);
legend('phase 3','phase 2','phase 1'); xlabel('amplitude'); ylabel('duration'); title('pbs PROCESSES');
xlim([-0.5 4]); ylim([0 250]);
figure()
hold on
scatter(PROCESSESTmev15AmpZone3,PROCESSESTmev15PeakDursZone3,'or','MarkerFaceColor','red','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESTmev15AmpZone2,PROCESSESTmev15PeakDursZone2,'ob','MarkerFaceColor','blue','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
scatter(PROCESSESTmev15AmpZone2,PROCESSESTmev15PeakDursZone2,'ok','MarkerFaceColor',[0.9 0.9 0],'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.2)
plot(PROCESSESTmev15AmpZone2,TMEV15P2fitLine,'-b');
plot(PROCESSESTmev15AmpZone3,TMEV15P3fitLine,'-r');
plot(PROCESSESTmev15AmpZone1,TMEV15P1fitLine,'-','Color',[0.9 0.9 0]);
legend('phase 3','phase 2','phase 1'); xlabel('amplitude'); ylabel('duration'); title('Tmev15 PROCESSES');
xlim([-0.5 4]); ylim([0 250]);

%% trace example
exampleTrace = PROCESSES_PBS.T2042_3.trace;
exampleTrace = (exampleTrace - mean(exampleTrace(1:50)))./mean(exampleTrace(1:50));
yLimit = max(exampleTrace);
yMin = min(exampleTrace);
exampleLocs = PROCESSES_PBS.T2042_3.PeakLocs;
exampleLocs = exampleLocs(~isnan(exampleLocs));
time = (1:1740)./0.97;
figA = figure();
plot(time(1:200), exampleTrace(1:200), 'LineWidth', 2, 'Color',[67/256, 94/256, 237/256]);
set(gca,'XColor',[1 1 1],'YColor',[1 1 1])
hold on
plot([0, 0], [-0.25, 0.75], '-k',[0, 60], [-0.25 -0.25], '-k','LineWidth',2)
ylim([yMin yLimit])
hold off
set(gca,'Visible','off')
set(gcf, 'PaperPosition', [4 4 1.0 1.5]);
print(figA,'exampleA.png', '-r900','-dpng');

figB = figure();
plot(time(201:1016), exampleTrace(201:1016), 'LineWidth', 1, 'Color',[67/256, 94/256, 237/256])
ylim([yMin yLimit])
set(gca,'XColor',[1 1 1],'YColor',[1 1 1])
set(gca,'Visible','off')
set(gcf, 'PaperPosition', [4 4 1.0 1.5]);
print(figB,'exampleB.png', '-r900','-dpng');
figC = figure();
plot(time(1017:end), exampleTrace(1017:end), 'LineWidth', 1, 'Color',[67/256, 94/256, 237/256])
ylim([yMin yLimit])
set(gca,'XColor',[1 1 1],'YColor',[1 1 1])
set(gca,'Visible','off')
set(gcf, 'PaperPosition', [4 4 1.5 2.0]);
print(figC,'exampleC.png', '-r900','-dpng');

figure()
plot(exampleTrace)
hold on
plot(exampleLocs,exampleTrace(exampleLocs),'v')
hold off

%% events per image

