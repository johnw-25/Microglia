function [process_data, all_F, peak_data, binary] = plot_burns_processes(process_data)

%% discriminate peaks
peak_data = cell(size(process_data,1),1);
peak_data(:,1) = process_data(:,1);
Fs = 1.04;
lag = 30; %window
threshold = 5.0; % threshold in standard deviations
influence = 0.99;
videos = size(process_data,1);
process_frames = max(cellfun('size',process_data(:,3),2));
fs = 1.04; % frames per second (sampling rate)

process_time = [-process_frames:process_frames]./fs;
process_time = round(process_time-1,1);
process_time = process_time(rem(process_time, 100) == 0);

all_F = zeros(videos, process_frames);
binary = zeros(videos, process_frames);
for j = 1:size(process_data,1)
    tempF = process_data{j,3};    
    % insert row numbers and convert to suite2p roi identifiers in 4th
    % columns
    F0 = mean(tempF(1,1:50));
%     F0 = mean(tempF(1,end-50:end)); % burn F0
    tempF(1,:) = (tempF(1,:) - F0)./F0;
    all_F(j,1:length(tempF)) = tempF;
%     [~, ~, ~, raw_p] = FindPeaks_Stim_v2(tempF(1,:) ,Fs, SDs, 'Frames');
    [signals,avgFilter,stdFilter] = peak_detect(tempF,lag,threshold,influence);
    [pks, locs, w, p] = findpeaks(signals);
    binary(j,1:length(signals)) = signals;
    % use p1_logicals as a lower limit in case peaks are detected in
    % eventless videos
    p1_check = any(p);
    height_check = any(tempF(locs) > 0.40);
    if  p1_check == 1 %&& height_check == 1 %&& any(raw_p) == 1
        % % % DO NOTHING
        all_F(j,:) = all_F(j,:);
        temp_locs = locs;
        temp_proms = p;
        temp_peaks = tempF(locs);
        temp_widths = w;
        
        peak_data{j,2} = struct('locs',temp_locs);
        peak_data{j,2}.proms = temp_proms;
        peak_data{j,2}.peaks = temp_peaks;
        peak_data{j,2}.widths = temp_widths;
    else
%         all_F(j,:) = NaN;
%         ROI(j) = NaN;
        temp_locs = [];
        temp_proms = [];
        temp_peaks = [];
        temp_widths = [];
        
        peak_data{j,2} = struct('locs',temp_locs);
        peak_data{j,2}.proms = temp_proms;
        peak_data{j,2}.peaks = temp_peaks;
        peak_data{j,2}.widths = temp_widths;
    end


end
peak_data(:,3) = process_data(:,4);
peak_data(:,4) = process_data(:,5);
end
