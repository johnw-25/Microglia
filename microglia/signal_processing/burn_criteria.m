%% The big question: what is a signal? % USE THIS SCRIPT FOR BURNS
% might be useful to look at F0 over many signals from spont. suite2p
% outputs to determine how far above the baseline a signal can be a signal
%
%% Load spont. suite2p data
[data_names, MyData, iscell_list] = loads2p();
%% discriminate peaks/signals based on just prominence standard deviations
peak_data = cell(21,1);
peak_data(:,1) = data_names(:,1);
SD = 4.5;
Fs = 1.04;
for j = 1:length(MyData)
    tempF = data_names{j,2}.F;
    F0 = mean(tempF(1:25));
    tempF = (tempF - F0)./F0;
    
    % insert row numbers and convert to suite2p roi identifiers in 4th
    % columns
    ROI = [1:size(tempF,1)]-1;
    ROI = ROI';
    temp_iscell = iscell_list{1,j};
    tempF = tempF(logical(temp_iscell(:,1)),:);
    ROI = ROI(logical(temp_iscell(:,1)),:);
    temp_locs = cell(size(tempF,1), 1);
    temp_proms = cell(size(tempF,1), 1);
    temp_peaks = cell(size(tempF,1), 1);
    testskew = skewness(tempF,0,2);
    for k = 1:size(tempF,1)
        % Dr. Pat Parker's script for finding peaks based on SD above
        % baseline
        [pks, locs, w, p] = FindPeaks_Stim_v2(tempF(k,:), Fs, SD, 'Frames');
        if isempty(pks) == 0 % check if any peaks were detected
        temp_locs{k,1} = locs;
        temp_proms{k,1} = p;
        temp_peaks{k,1} = pks;
        temp_locs{k,2} = ROI(k);
        % create structure in the correct cell in peak_data
        else % if no peaks were found, turn data into NaN
            tempF(k,:) = NaN;
            ROI(k) = NaN;
        end
    end
    peak_data{j,2} = struct('locs',{temp_locs});
    peak_data{j,2}.proms = temp_proms;
    peak_data{j,2}.peaks = temp_peaks;
    peak_data{j,2}.ROI = ROI;
    peak_data{j,2}.skew = testskew;
%     data_names{j,3} = tempF(logical(temp_iscell(:,1)),:);
%     data_names{j,4} = ROI(logical(temp_iscell(:,1)),:);
end
allskew = [];
% for 1:size(peak_data,2)
%     allskew = 
%% plotting
all_skew = []; % init. step
all_prom = [];
for jj = 1:size(peak_data,1)
    ROIs = peak_data{jj,2}.ROI;
    ROIs = ROIs(~isnan(ROIs));
    plot_F = data_names{jj,3};
%     temp_skew = peak_data{jj,2}.skew;
%     temp_prom = peak_data{jj,2}.proms;
%     all_skew = [all_skew, temp_skew];
%     all_prom = [all_prom, temp_prom];
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
%% testing individual traces
test_F = data_names{5,3};
F = test_F(34,:);
F0 = mean(F(1:25));
dF = (F - F0)./F0;
figure()
plot(test_F)
[pks, locs, w, p] = findpeaks(dF,'MinPeakDistance', 25, 'MinPeakWidth', 5, 'MaxPeakWidth', 20);
%% open and import .npy file
[filename, pathname] = uigetfile('F.npy', 'Select your .npy F data file');
%cancel if user clicks cancel
if isequal(filename,0) || isequal(pathname,0)
    disp('action canceled')
return
end
    %set current directory to pathname
    cd(pathname);
    %set file to path string
   file = [pathname filename];
   %Open Fluorescence info and neuropil data
   ROIdata_all = readNPY('F.npy')';
 
   %open neuropil data
   Fneu = readNPY('Fneu.npy')';
   
   %open iscell
   iscell=readNPY('iscell.npy');
   %Only analyze neuropil and ROIs that user selected as cells (1=chosen,
   %0=discarded).
   ROIdata_all=ROIdata_all(:,logical(iscell(:,1)));
   Fneu=Fneu(:,logical(iscell(:,1)));