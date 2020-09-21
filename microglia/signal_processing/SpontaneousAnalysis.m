function PEAK_METRICS = SpontaneousAnalysis(spontaneousData, valleySD,peakSD, varargin)

if nargin > 3
    % continue where last run instance left off
    lastTrace = length(varargin);
    % truncate portion that has already been processed
    spontaneousData = spontaneousData(lastTrace+1:end,:);
end

for ii = 1:size(spontaneousData,1)
    F_matrix = double(spontaneousData{ii,3});
    videoROI = spontaneousData{ii,4};
    for jj = 1:size(F_matrix,1)
        currentROI = videoROI(jj);
        test_trace = F_matrix(jj,:);
        frames = length(test_trace)/2;
        % % Find arc length of two halves of trace % %
        test_trace = double(spontaneousData{ii,3});
        F0 = mean(test_trace(1:50));
        test_trace = (test_trace - F0)./ F0;
        bp = 1:400:frames;
        detrended_test = detrend(test_trace(:), 2, bp, 'Continuous', false);
        % detrended_test = abs(detrended_test);
        % detrended_test(detrended_test < 0) = 0;
        % % Plot detrended and original to see difference
        TraceFigure = figure();
        plot(test_trace)
        hold on
        plot(detrended_test)
        legend('original trace','detrended trace');
        grid on
        
        [LPF, response] = noiseRemover(100, 0.06, 1.04);
        test_filter = filtfilt(LPF, detrended_test);
        Length = signalLength(test_filter, 1.04);
        [raw_auc, filtered_auc, peakDurations, peakWidth, interEventIntervals, pks, locs, filtFig, removed,removedValleys,addedValleys] = peakMeasure(test_trace,detrended_test,test_filter, valleySD, peakSD);
        PEAK_METRICS.(spontaneousData{ii,1}).(currentROI).AUC = filtered_auc;
        PEAK_METRICS.(spontaneousData{ii,1}).(currentROI).RAW_AUC = raw_auc;
        PEAK_METRICS.(spontaneousData{ii,1}).(currentROI).PeakDurs = peakDurations;
        PEAK_METRICS.(spontaneousData{ii,1}).(currentROI).PeakWidths = peakWidth;
        PEAK_METRICS.(spontaneousData{ii,1}).(currentROI).IEI = interEventIntervals;
        PEAK_METRICS.(spontaneousData{ii,1}).(currentROI).Figure = TraceFigure;
        PEAK_METRICS.(spontaneousData{ii,1}).(currentROI).LENGTH = Length;
        PEAK_METRICS.(spontaneousData{ii,1}).(currentROI).FreqResponse = response;
        PEAK_METRICS.(spontaneousData{ii,1}).(currentROI).PeakFig = filtFig;
        PEAK_METRICS.(spontaneousData{ii,1}).(currentROI).PeakLocs = locs;
        PEAK_METRICS.(spontaneousData{ii,1}).(currentROI).PeakMags = pks;
        PEAK_METRICS.(spontaneousData{ii,1}).(currentROI).RemovedPks = removed;
        PEAK_METRICS.(spontaneousData{ii,1}).(currentROI).RemovedValleys = removedValleys;
        PEAK_METRICS.(spontaneousData{ii,1}).(currentROI).AddedValleys = addedValleys;
        set(gcf, 'NumberTitle','off');
        set(gcf, 'Name', spontaneousData{ii,1});
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
end