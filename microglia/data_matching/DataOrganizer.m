function [PBS_Data,TimeZones] = DataOrganizer(PBS, DATA, burnFrame) 
isburn = @(x) logical(x>90);
burnFlag = false;
PBSfields = fieldnames(PBS);
allPBSLocs = []; %allocate space
allPBSPeaks = [];
allPeakDurs = [];
all_RAWauc = [];
all_auc = [];
allIEI = [];
allLocs = [];
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

Group1_IEI = [];
Group2_IEI = [];
Group3_IEI = [];
% % count number of mice for PBS % %
mice = DATA(:,7);
for k = 1:numel(PBSfields)

    % extract
    tempIEI = PBS.(PBSfields{k}).IEI;
    temp_timehitR = PBS.(PBSfields{k}).timehitR;
    tempLocs = PBS.(PBSfields{k}).PeakLocs;
    tempLocs = tempLocs(~isnan(tempLocs));
    tempIEILocs = tempLocs(1:end-1);
    allIEI = [allIEI; tempIEI, tempIEILocs];
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
        
        temp_Group1_IEILocs = tempIEILocs(tempIEILocs <= burnFrame);
        temp_Group2_IEILocs = tempIEILocs((isburn(tempIEILocs)) & (tempIEILocs < round(temp_timehitR*0.97)));
        temp_Group3_IEILocs = tempIEILocs(tempIEILocs >=  round(temp_timehitR*0.97));
        temp_Group1_IEI = tempIEI(tempIEILocs <= burnFrame);
        temp_Group2_IEI = tempIEI((isburn(tempIEILocs)) & (tempIEILocs < round(temp_timehitR*0.97)));
        temp_Group3_IEI = tempIEI(tempIEILocs >=  round(temp_timehitR*0.97));
        Group1_IEI = [Group1_IEI; temp_Group1_IEILocs, temp_Group1_IEI];
        Group2_IEI = [Group2_IEI; temp_Group2_IEILocs, temp_Group2_IEI];
        Group3_IEI = [Group3_IEI; temp_Group3_IEILocs, temp_Group3_IEI];
    end
    
    allPBSPeaks = [allPBSPeaks; tempPeaks];
    all_RAWauc = [all_RAWauc; tempRAWAUC];
    all_auc = [all_auc; tempAUC];
    shiftedLocs = noBurnLocs - temp_timehitR;
    shiftedLocs = shiftedLocs(~isnan(shiftedLocs));
    PBS.(PBSfields{k}).shiftedLocs = shiftedLocs;
    
    allPBSLocs = [allPBSLocs; shiftedLocs];
    allPeakDurs = [allPeakDurs; tempDurs];
end
for j = 1:length(mice)
    for k = 1:length(allPBSTraces)
        if strcmp(DATA{j,1},allPBSTraces{k})
            allPBSMice = [allPBSMice; mice(j)];
        end
    end
end
TimeZones.Zone1_Locs = Group1_Locs;
TimeZones.Zone2_Locs = Group2_Locs;
TimeZones.Zone3_Locs = Group3_Locs;

TimeZones.Zone1_PeakDurs = Group1_PeakDurs;
TimeZones.Zone2_PeakDurs = Group2_PeakDurs;
TimeZones.Zone3_PeakDurs = Group3_PeakDurs;

TimeZones.Zone1_PeakMags = Group1_PeakMags;
TimeZones.Zone2_PeakMags = Group2_PeakMags;
TimeZones.Zone3_PeakMags = Group3_PeakMags;

TimeZones.Zone1_AUC = Group1_AUC;
TimeZones.Zone2_AUC = Group2_AUC;
TimeZones.Zone3_AUC = Group3_AUC;

TimeZones.Zone1_IEI = Group1_IEI;
TimeZones.Zone2_IEI = Group2_IEI;
TimeZones.Zone3_IEI = Group3_IEI;

PBS_Data = cell(length(allPBSTraces)+1, 12);
PBS_Data(1,:) = {'Trace Name','ROI #','DPI','TMEV/PBS','Duration','Peak Amp','AUC','Peak Locations','Shifted Locations','Event Group','MouseID','IEI'};
PBS_Data(2:end,:) = [allPBSTraces,allPBSROI, allPBSDPI, allPBSpath, allPBSDurs,allPBSMags,allPBSAUC,allPBSNormalLocs,allPBSShiftedLocs,allPBSEventGroup,allPBSMice,allPBSIEI];
end