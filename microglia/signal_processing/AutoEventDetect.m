function peak_data = AutoEventDetect(MyData,videos,Peak_SDs,filterOrder,Fc,bp)

% Loop on each folder
% MyData = cell(size(fallList,1),1);
% for i = 1:length(fallList) % 1 to number of Fall.mat files
%   filetoread = fullfile(fallList(i).folder,'Fall.mat'); % extract suite2p data from video subfolders
%   % store each Fall.mat file into MyData
%   fileID = fopen(filetoread);
%   MyData{i} = load(filetoread); % the format depends of your files
%   fclose(fileID);
% end

Fs = 1.04;
LPF = noiseRemover(filterOrder, Fc, Fs);
peak_data = cell(size(MyData,1),3);
for j = 1:length(MyData)
    fovName = str2double(videos(j).name);
    peak_data{j,1} = fovName;
    tempF = MyData{j}.F;
    % insert row numbers and convert to suite2p roi identifiers in 4th
    % columns
    ROI = [1:size(tempF,1)]-1;
    ROI = ROI';
    temp_locs = cell(size(tempF,1), 1);
    temp_proms = cell(size(tempF,1), 1);
    temp_peaks = cell(size(tempF,1), 1);
    temp_widths = cell(size(tempF,1), 1);
    for k = 1:size(tempF,1)
        % set our constant variables
        F0 = mean(tempF(k,1:50)); %F0
        tempF(k,:) = (tempF(k,:) - F0)./F0; % df/f
        trace = tempF(k,:);
        % measure length of signal
        frames = length(tempF(k,:));
        % set breakpoints in signal for detrending
        BP = 1:bp:frames;
        % detrend signal
        dtTrace = detrend(tempF(k,:), 2, BP, 'Continuous', false);
        [~, ~, ~, proms] = findpeaks(dtTrace);
        MeanProms = mean(proms);%mean prominence of all peaks for individual trace
        StdProms = std(proms);%standard deviation of prominence of all peaks for individual trace
        UseProms = MeanProms + StdProms*Peak_SDs;%determine N standard deviations
        
        filteredTempF = filtfilt(LPF, double(dtTrace));
        [~, locs, w, p] = findpeaks(filteredTempF, 'MinPeakProminence', UseProms, 'MinPeakDistance', 10);
        % use p1_logicals as a lower limit in case peaks are detected in
        % eventless videos
        p1_check = any(p);
        if  p1_check == 1
            % % % DO NOTHING
            tempF(k,:) = tempF(k,:);
            temp_locs{k,1} = locs;
            temp_proms{k,1} = p;
            temp_peaks{k,1} = trace(locs);
            temp_locs{k,2} = ROI(k);
            temp_widths{k,2} = w;
        else
            tempF(k,:) = NaN;
            ROI(k) = NaN;
        end
    end
    peak_data{j,2} = struct('locs',{temp_locs});
    peak_data{j,2}.proms = temp_proms;
    peak_data{j,2}.peaks = temp_peaks;
    peak_data{j,2}.ROI = ROI;
    peak_data{j,2}.widths = temp_widths;
    peak_data{j,3} = tempF;
end

end