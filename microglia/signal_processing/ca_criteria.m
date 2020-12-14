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
lag = 50; %window
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
        if  p1_check == 1 && height_check == 1 && any(raw_p)
            % % % DO NOTHING
            tempF(k,:) = tempF(k,:);
            temp_locs{k,1} = locs;
            temp_proms{k,1} = p;
            temp_peaks{k,1} = data(locs);
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
%% roi plotting
set(0,'defaultfigurecolor',[1 1 1])
T2314_ROIs = [3; 66];
T2314_useROI = 1:length(T2314_ROIs);
T2314plot = figure();
ROI_plot_general_purpose(T2314_ROIs+1, 'PBS')

locs2314 = peak_data{35,2}.locs;
T2314_F = data_names{35,3};
time2314 = (1:940)./1.04;
fig = figure();
hold on
for k = 1:length(T2314_ROIs)
    tempLocs = locs2314{T2314_ROIs(k)+1};
    tempF = T2314_F(T2314_ROIs(k)+1,:);
    plot(time2314, tempF+k*2,'LineWidth',1.0,'Color',[0.9882, 0.8, 0])
    plot(tempLocs./1.04, (tempF(tempLocs)+k*2).*1.15,'o','Color',[0, 0, 0],'MarkerFaceColor',[0, 0, 0],'MarkerSize',1.25)
end
plot([0 0], [0, 1], '-k',[0, 60*1.04], [0 0], '-k','LineWidth',1)
ylim([0 setLim]);
set(gca,'Visible','off')
hold off
set(gcf, 'Units','inches','position',[4 4 2.75 3]);
set(gcf, 'PaperPosition', [4 4 2.75 3]);
print(fig,'T2314_traces_pbs.png', '-r900','-dpng');


T2289_ROIs = data_names{29,4};
T2289_ROIs = T2289_ROIs(~isnan(T2289_ROIs));
T2289_useROI = 1:length(T2289_ROIs);
ROI_plot_general_purpose(T2289_ROIs+1, 'TMEV 2 DPI')

locs2289 = peak_data{29,2}.locs;
T2289_F = data_names{29,3};
time2289 = (1:940)./1.04;
fig = figure();
hold on
for k = 1:length(T2289_ROIs)
    tempLocs = locs2289{T2289_ROIs(k)+1};
    tempF = T2289_F(T2289_ROIs(k)+1,:);
    plot(time2289, (tempF+k*2),'LineWidth',1.0,'Color',[0.9882, 0.8, 0])
    plot(tempLocs./1.04, (tempF(tempLocs)+k*2)*1.05,'o','Color',[0, 0, 0],'MarkerFaceColor',[0, 0, 0],'MarkerSize',1.25)
end
plot([0 0], [0, 1], '-k',[0, 60*1.04], [0 0], '-k','LineWidth',1)
limits = ylim();
setLim = limits(2);
set(gca,'Visible','off')
hold off
set(gcf, 'Units','inches','position',[4 4 2.75 3]);
set(gcf, 'PaperPosition', [4 4 2.75 3]);
print(fig,'T2289_traces_tmev2dpi.png', '-r900','-dpng');


T2235_ROIs = data_names{26,4};
T2235_ROIs = T2235_ROIs(~isnan(T2235_ROIs));
T2235_useROI = 1:length(T2235_ROIs);
T2235plot = figure();
ROI_plot_general_purpose(T2235_ROIs+1, 'TMEV 5 DPI')

locs2235 = peak_data{26,2}.locs;
T2235_F = data_names{26,3};
time2235 = (1:940)./1.04;
fig = figure();
hold on
for k = 1:length(T2235_ROIs)
    tempLocs = locs2235{T2235_ROIs(k)+1};
    tempF = T2235_F(T2235_ROIs(k)+1,:);
    plot(time2235, tempF+k*2,'LineWidth',1.0,'Color',[0.9882, 0.8, 0])
    plot(tempLocs./1.04, (tempF(tempLocs)+k*2)*1.15,'o','Color',[0, 0, 0],'MarkerFaceColor',[0, 0, 0],'MarkerSize',1.25)
end
plot([0 0], [0, 1], '-k',[0, 60*1.04], [0 0], '-k','LineWidth',1)
ylim([0 setLim]);
set(gca,'Visible','off')
hold off
set(gcf, 'Units','inches','position',[4 4 2.75 3]);
set(gcf, 'PaperPosition', [4 4 2.75 3]);
print(fig,'T2235_traces_tmev5dpi.png', '-r900','-dpng');


T2130_ROIs = [1; 13];
T2130_useROI = 1:length(T2130_ROIs);
ROI_plot_general_purpose(T2130_ROIs+1, 'TMEV 15 DPI')

locs2130 = peak_data{13,2}.locs;
T2130_F = data_names{13,3};
time2130 = (1:940)./1.04;
fig = figure();
hold on
for k = 1:length(T2130_ROIs)
    tempLocs = locs2130{T2130_ROIs(k)+1};
    tempF = T2130_F(T2130_ROIs(k)+1,:);
    plot(time2130, tempF+k*2,'LineWidth',1.0,'Color',[0.9882, 0.8, 0])
    plot(tempLocs./1.04, (tempF(tempLocs)+k*2).*1.15,'o','Color',[0, 0, 0],'MarkerFaceColor',[0, 0, 0],'MarkerSize',1.25)
end
plot([0 0], [0, 1], '-k',[0, 60*1.04], [0 0], '-k','LineWidth',1)
ylim([0 setLim]);
set(gca,'Visible','off')
hold off
set(gcf, 'Units','inches','position',[4 4 2.75 3]);
set(gcf, 'PaperPosition', [4 4 2.75 3]);
print(fig,'T2130_traces_tmev15dpi.png', '-r900','-dpng');