function allRMS = ProcessRMS(DATA, filterSelect)
allRMS{1,1} = 'Trace + ROI';
allRMS{1,2} = 'TMEV/PBS';
allRMS{1,3} = 'DPI';
allRMS{1,4} = 'A(1 to 200) RMS^2';
allRMS{1,5} = 'B(201 to timehitR) RMS^2';
allRMS{1,6} = 'C(timehitR to end) RMS^2';
allRMS{1,7} = 'timehitR Marker/Index';
if strcmp(filterSelect, 'LPF')
%     LPF = noiseRemover(100,0.012,0.97);
    LPF = noiseRemover(100,0.1,0.97);
end
for ii = 1:size(DATA,1)
    % check if trace should be deleted/not included
%     if DATA{ii,8}
        ROInum = DATA{ii,2};
        ROInum = num2str(ROInum);
        ROI = strcat('  ROI: ', ROInum);
        traceName = DATA{ii,1};
        traceName = strcat(traceName,ROI);
        DPI = DATA{ii,5};
        timehitR = double(DATA{ii,6});
        timehitR_marker = round(timehitR);
        currentTrace = DATA{ii,3};
        F0 = mean(currentTrace(1:50));
        currentTrace = (currentTrace-F0)./F0;
%         figure()
%         subplot(1,2,1)
%         plot(currentTrace)
        if strcmp(filterSelect, 'LPF')
            currentTrace = filtfilt(LPF,double(currentTrace));
%             subplot(1,2,2)
%             plot(currentTrace)
        end
        if timehitR_marker > 200 && timehitR_marker < length(currentTrace)
            
            rmsA = rms(currentTrace(1:200));
            rmsB = rms(currentTrace(201:timehitR_marker-1));
            rmsC = rms(currentTrace(timehitR_marker:end));
            allRMS{ii+1,1} = traceName;
            allRMS{ii+1,2} = DATA{ii,4};
            allRMS{ii+1,3} = DPI;
            allRMS{ii+1,4} = rmsA;
            allRMS{ii+1,5} = rmsB;
            allRMS{ii+1,6} = rmsC;
            allRMS{ii+1,7} = timehitR_marker;
        else
            allRMS{ii+1,4} = NaN;
            allRMS{ii+1,5} = NaN;
            allRMS{ii+1,6} = NaN;
        end
%     else
%         continue
%     end
end
end
