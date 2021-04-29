function peak_data = AutoSignalPeak(signal,Peak_SDs,filterOrder,Fc,bp)


Fs = 1.04;
LPF = noiseRemover(filterOrder, Fc, Fs);
% measure length of signal
frames = length(signal);
% set breakpoints in signal for detrending
BP = 1:bp:frames;
% detrend signal
dtTrace = detrend(signal, 2, BP, 'Continuous', false);
[~, ~, ~, proms] = findpeaks(dtTrace);
MeanProms = mean(proms);%mean prominence of all peaks for individual trace
StdProms = std(proms);%standard deviation of prominence of all peaks for individual trace
UseProms = MeanProms + StdProms*Peak_SDs;%determine N standard deviations

filteredTempF = filtfilt(LPF, double(dtTrace));
[~, locs, w, p] = findpeaks(filteredTempF, 'MinPeakProminence', UseProms, 'MinPeakDistance', 10);

peak_data.locs = locs;
peak_data.proms = p;
peak_data.peaks = signal(locs);
peak_data.widths = w;
peak_data.filtered_signal = filteredTempF;
end
