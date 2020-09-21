function PEAK_METRICS = MainAnalysis(process_data, valleySD,peakSD,filterOrder,Fc,Fs,polyDeg,bp, varargin)

if nargin > 8
        % continue where last run instance left off
        lastTrace = length(fieldnames(varargin{1}));
        % truncate portion that has already been processed
        process_data = process_data(lastTrace+1:end,:);
end

for ii = 1:size(process_data,1)
    ROInum = process_data{ii,2};
    ROInum = num2str(ROInum);
    ROI = strcat('  ROI: ', ROInum);
    traceName = process_data{ii,1};
    traceName = strcat(traceName,ROI);
    
    shifted_trace = double(process_data{ii,7});
    timehitR = double(process_data{ii,6});
    shifted_trace(isnan(shifted_trace)) = 0;
    frames = length(shifted_trace)/2;
    traceA = shifted_trace(1:frames);
    traceB = shifted_trace(frames+1:end);
    % % Find arc length of two halves of trace % %
    lengthA = signalLength(traceA,Fs)/timehitR;
    lengthB = signalLength(traceB,Fs)/(frames*Fs-timehitR);
    test_trace = double(process_data{ii,3});
    F0 = mean(test_trace(1:50));
    test_trace = (test_trace - F0)./ F0;
    BP = 1:bp:frames;
    detrended_test = detrend(test_trace(:), polyDeg, BP, 'Continuous', false);
    % detrended_test = abs(detrended_test);
    % detrended_test(detrended_test < 0) = 0;
    % % Plot detrended and original to see difference
    TraceFigure = figure();
    plot(test_trace)
    hold on
    plot(detrended_test)
    legend('original trace','detrended trace');
    grid on
    
    [LPF, response] = noiseRemover(filterOrder, Fc, Fs);
    test_filter = filtfilt(LPF, detrended_test);
    
    [raw_auc, filtered_auc, peakDurations, peakWidth, interEventIntervals, pks, locs, filtFig, removed,removedValleys,addedValleys] = peakMeasure(test_trace,detrended_test,test_filter, valleySD, peakSD, traceName);
    PEAK_METRICS.(process_data{ii,1}).AUC = filtered_auc;
    PEAK_METRICS.(process_data{ii,1}).RAW_AUC = raw_auc;
    PEAK_METRICS.(process_data{ii,1}).PeakDurs = peakDurations;
    PEAK_METRICS.(process_data{ii,1}).PeakWidths = peakWidth;
    PEAK_METRICS.(process_data{ii,1}).IEI = interEventIntervals;
    PEAK_METRICS.(process_data{ii,1}).Figure = TraceFigure;
    PEAK_METRICS.(process_data{ii,1}).LENGTH_A = lengthA;
    PEAK_METRICS.(process_data{ii,1}).LENGTH_B = lengthB;
    PEAK_METRICS.(process_data{ii,1}).FreqResponse = response;
    PEAK_METRICS.(process_data{ii,1}).PeakFig = filtFig;
    PEAK_METRICS.(process_data{ii,1}).PeakLocs = locs;
    PEAK_METRICS.(process_data{ii,1}).PeakMags = pks;
    PEAK_METRICS.(process_data{ii,1}).RemovedPks = removed;
    PEAK_METRICS.(process_data{ii,1}).RemovedValleys = removedValleys;
    PEAK_METRICS.(process_data{ii,1}).AddedValleys = addedValleys;
    PEAK_METRICS.(process_data{ii,1}).Pathology = process_data{ii,4};
    PEAK_METRICS.(process_data{ii,1}).DPI = process_data{ii,5};
    PEAK_METRICS.(process_data{ii,1}).ROI = process_data{ii,2};
    PEAK_METRICS.(process_data{ii,1}).timehitR = timehitR;
    
    set(gcf, 'NumberTitle','off');
    set(gcf, 'Name', process_data{ii,1});
    % pause;
    close all
    nexTrace = questdlg('Continue to the next trace or exit the program? Current progress will be saved.', 'Program Query', 'Continue', 'Exit', 'Continue');
    
    switch nexTrace
        case 'Continue'
            continue % proceed to next trace
        case 'Exit'
            return; % return to program that called MainAnalysis
    end
end

end