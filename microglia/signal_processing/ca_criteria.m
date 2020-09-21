%% The big question: what is a signal?
% might be useful to look at F0 over many signals from spont. suite2p
% outputs to determine how far above the baseline a signal can be a signal
%
%% Load spont. suite2p data
[data_names, MyData, iscell_list] = loads2p();
%% discriminate peaks
TN = cell(size(data_names,1),1);
peak_data = cell(size(data_names,1),1);
peak_data(:,1) = data_names(:,1);
SDs = 2.5;
Fs = 1.04;
lag = 30; %window
threshold = 5.0; % threshold in standard deviations
influence = 0.99;
for j = 1:length(MyData)
    tempF = data_names{j,2}.F;
    % insert row numbers and convert to suite2p roi identifiers in 4th
    % columns
    ROI = [1:size(tempF,1)]-1;
    ROI = ROI';
    temp_locs = cell(size(tempF,1), 1);
    temp_proms = cell(size(tempF,1), 1);
    temp_peaks = cell(size(tempF,1), 1);
    temp_widths = cell(size(tempF,1), 1);
    for k = 1:size(tempF,1)
        F0 = mean(tempF(k,1:50)); % spont. F0
%         F0 = mean(tempF(k,end-50:end)); % burn F0
        tempF(k,:) = (tempF(k,:) - F0)./F0;
        [~, ~, ~, raw_p] = FindPeaks_Stim_v2(tempF(k,:) ,Fs, SDs, 'Frames');
        data = tempF(k,:);
        [signals,avgFilter,stdFilter] = peak_detect(data,lag,threshold,influence);
        [pks, locs, w, p] = findpeaks(signals);
        % use p1_logicals as a lower limit in case peaks are detected in
        % eventless videos
        p1_check = any(p);
        height_check = any(data(locs) > 0.40);
        if  p1_check == 1 && height_check == 1 %&& any(raw_p) == 1
            % % % DO NOTHING
            tempF(k,:) = tempF(k,:);
            temp_locs{k,1} = locs;
            temp_proms{k,1} = p;
            temp_peaks{k,1} = data(locs);
            temp_locs{k,2} = ROI(k);
            temp_widths{k,2} = w;
        else
%             tempF(k,:) = NaN;
%             ROI(k) = NaN;
        end
    end
    peak_data{j,2} = struct('locs',{temp_locs});
    peak_data{j,2}.proms = temp_proms;
    peak_data{j,2}.peaks = temp_peaks;
    peak_data{j,2}.ROI = ROI;
    peak_data{j,2}.widths = temp_widths;
    temp_iscell = iscell_list{1,j};
    temp_TN = sum(isnan(ROI));
    TN{j,1} = temp_TN;
    %     data_names{j,3} = tempF(logical(temp_iscell(:,1)),:);
    %     data_names{j,4} = ROI(logical(temp_iscell(:,1)),:);
    data_names{j,3} = tempF;
    data_names{j,4} = ROI;
end
%% plotting
all_skew = []; % init. step
all_prom = [];
for jj = 1:size(peak_data,1)
    ROIs = peak_data{jj,2}.ROI;
    ROIs = ROIs(~isnan(ROIs));
    plot_F = data_names{jj,3};
    figure()
    for x = 1:size(data_names{jj,4},1)
        signal = plot_F(x,:) + x;
        plot(signal);
        hold on
    end
    grid on
    title(num2str(data_names{jj,1}));
    legend(num2str(ROIs))
    for x = 1:size(data_names{jj,4},1)
        pks = peak_data{jj,2}.peaks{x,1} + x;
        locs = peak_data{jj,2}.locs{x,1};
        plot(locs,pks,'x');
        hold on
    end
end
%% PEAK DETECT
lag = 75; %window
threshold = 6; % threshold in standard deviations
influence = 0.99;
tempF = data_names{30,3};
data = tempF(64,:);
[signals,avgFilter,stdFilter] = peak_detect(data,lag,threshold,influence);
[~, ~, ~, raw_p] = FindPeaks_Stim_v2(data ,Fs, SDs, 'Frames');
[~,locs,~,~] = findpeaks(signals);
height_check = any(data(locs) > 0.40);
SDs = 2.5;
if height_check == 1 %&& any(raw_p) == 1
x = 1:length(data);
figure()
plot(x,data)
hold on
plot(x,signals)
legend('unfiltered', 'filtered')
% figure()
% plot(avgFilter)
% else
%     % nothing
end
%% testing individual traces
tempF = data_names{1,3};
data = tempF(40,:);
skew = skewness(tempF,0,2);
% skew(41) = 0;
test_sd = std(skew);
test_mean = mean(skew);
testthresh = test_mean + test_sd.*0.1;
figure()
histogram(skew)
figure()
histogram(pks)
[peaks, locs, w, p] = findpeaks(data,'MinPeakDistance', 25);
% 306 peaks before applying any additional findpeaks  arguments
% MinPeakDistance 20-25 seems to appropriate for 2195 ROI 15 (activity
% group 2 with quick, short peaks, high frequency, "talkative")
figure()
plot(data)
findpeaks(data,'MinPeakDistance', 25,'Annotate','extents');
title('2104 ROI 24: Talkative signal')
ylabel('\DeltaF/F');
figure()
histogram(w,'BinWidth',0.2)
title('2104 ROI 24: Width Frequencies')
figure()
histogram(p,'BinWidth',0.2)
title('2104 ROI 24: Prominence Frequencies')
MeanProms = mean(p);%mean prominence of all peaks for individual trace
StdProms = std(p);%standard deviation of prominence of all peaks for individual trace
UseProms = MeanProms + StdProms*2;%determine N standard deviations
MeanW = mean(w);%mean prominence of all peaks for individual trace
StdW = std(w);%standard deviation of prominence of all peaks for individual trace
UseW = MeanW + StdW*2;%determine N standard deviations
%% trying signal filtering in this section
SDs = 2.5;
tempF = data_names{1,3};
data = double(tempF(11,:));

window = 11; % arbitrary window width
kernel = ones(1,window)/window;

avg_filter = filter(kernel, 1, data);
cubic_filter = sgolayfilt(data, 5, 9);

[avg_pks, locs, ~, ~] = FindPeaks_Stim_v2(avg_filter ,Fs, SDs, 'Frames');
figure()
plot(data)
hold on
plot(avg_filter, 'LineWidth', 1.5)
hold on
plot(locs, avg_pks,'d','MarkerFaceColor', 'r')
legend('unfiltered', 'filtered');

% figure()
% plot(data)
% hold on
% plot(cubic_filter)
% legend('unfiltered', 'cubic');
decay_data = avg_filter(locs(1):locs(1)+30);
figure()
plot(decay_data)