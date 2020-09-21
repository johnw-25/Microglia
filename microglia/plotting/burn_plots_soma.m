function [soma_data, all_F, peak_data, binary] = burn_plots_soma(soma_data)

%% discriminate peaks
peak_data = cell(size(soma_data,1),1);
peak_data(:,1) = soma_data(:,1);
Fs = 1.04;
lag = 30; %window
threshold = 5.0; % threshold in standard deviations
influence = 0.99;
videos = size(soma_data,1);
frames = max(cellfun('size',soma_data(:,3),2));
all_F = zeros(videos, frames);
binary = zeros(videos, frames);
for j = 1:size(soma_data,1)
    tempF = soma_data{j,3};    
    % insert row numbers and convert to suite2p roi identifiers in 4th
    % columns
%     temp_locs = cell(size(all_F,1), 1);
%     temp_proms = cell(size(all_F,1), 1);
%     temp_peaks = cell(size(all_F,1), 1);
%     temp_widths = cell(size(all_F,1), 1);
%     F0 = mean(tempF(1,end-50:end)); % burn F0
    F0 = mean(tempF(1,1:50));
    tempF(1,:) = (tempF(1,:) - F0)./F0;
    all_F(j,1:length(tempF)) = tempF;
%     [~, ~, ~, raw_p] = FindPeaks_Stim_v2(tempF(1,:) ,Fs, SDs, 'Frames');
    [signals,avgFilter,stdFilter] = peak_detect(tempF,lag,threshold,influence);
    binary(j,1:length(signals)) = signals;
    [pks, locs, w, p] = findpeaks(signals);
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

%% Break up into PBS/TMEV ( somaes )
peak_data(:,3) = soma_data(:,4);
peak_data(:,4) = soma_data(:,5);
p_TMEV_logicals = contains(soma_data(:,4),'TMEV'); 
% rows where TMEV occurs
TMEV_rows = find(p_TMEV_logicals);
p_TMEV_dF = all_F(p_TMEV_logicals,:);
% extract main table data
p_TMEV_extracted = soma_data(p_TMEV_logicals,:);
for k = 1:size(p_TMEV_dF,1)
p_TMEV_extracted{k,3} = p_TMEV_dF(k,:);
end
p_TMEV_extracted = sortrows(p_TMEV_extracted,5); % finished

% % % split TMEV into 2 dpi, 5/6 dpi, 14/16 dpi
idx_p_TMEV_2_dpi = find(cell2mat(p_TMEV_extracted(:,5)) == 2);
p_TMEV_2_dpi = p_TMEV_extracted(1:idx_p_TMEV_2_dpi(end),:);

idx_p_TMEV_5_dpi = find(cell2mat(p_TMEV_extracted(:,5)) == 5);
idx_p_TMEV_6_dpi = find(cell2mat(p_TMEV_extracted(:,5)) == 6);
p_TMEV_5_dpi = p_TMEV_extracted(idx_p_TMEV_5_dpi(1):idx_p_TMEV_6_dpi(end),:);

idx_p_TMEV_14_dpi = find(cell2mat(p_TMEV_extracted(:,5)) == 14);
p_TMEV_14_dpi = p_TMEV_extracted(idx_p_TMEV_14_dpi(1):end,:);

% repeat for PBS group
p_PBS_logicals = contains(soma_data(:,4),'PBS'); 
PBS_rows = find(p_PBS_logicals);
p_PBS_dF = all_F(p_PBS_logicals,:);
p_PBS_extracted = soma_data(p_PBS_logicals,:);
for k = 1:size(p_PBS_dF,1)
p_PBS_extracted{k,3} = p_PBS_dF(k,:);
end
p_PBS_extracted = sortrows(p_PBS_extracted,5); % finished

% % % break up into PBS DPI (2-5 and 6-14/16)
idx_p_PBS_2_5_dpi = find(cell2mat(p_PBS_extracted(:,5)) == 5);
p_PBS_2_5_dpi = p_PBS_extracted(1:idx_p_PBS_2_5_dpi(end),:);

idx_p_PBS_6_14_dpi = find(cell2mat(p_PBS_extracted(:,5)) == 6);
p_PBS_6_14_dpi = p_PBS_extracted(idx_p_PBS_6_14_dpi(1):end,:); 

% % % PBS axis labels % % %
p_PBS_videos_2_5_dpi = p_PBS_2_5_dpi(:,1);
p_PBS_dpi = num2str(cell2mat(p_PBS_2_5_dpi(:,5)));
p_PBS_ROI = num2str(cell2mat(p_PBS_2_5_dpi(:,2)));
p_PBS_2_5_labels = cell(size(p_PBS_videos_2_5_dpi,1),1);
for i = 1:length(p_PBS_2_5_labels)
p_PBS_2_5_labels(i) = strcat(p_PBS_videos_2_5_dpi(i),' , ', p_PBS_dpi(i,:),' , ', p_PBS_ROI(i,:));
end

p_PBS_videos_6_14_dpi = p_PBS_6_14_dpi(:,1);
p_PBS_dpi = num2str(cell2mat(p_PBS_6_14_dpi(:,5)));
p_PBS_ROI = num2str(cell2mat(p_PBS_6_14_dpi(:,2)));
p_PBS_6_14_labels = cell(size(p_PBS_videos_6_14_dpi,1),1);
for i = 1:length(p_PBS_6_14_labels)
p_PBS_6_14_labels(i) = strcat(p_PBS_videos_6_14_dpi(i),' , ', p_PBS_dpi(i,:),' , ', p_PBS_ROI(i,:));
end
% % % TMEV axis labels % % %

p_TMEV_videos = p_TMEV_2_dpi(:,1);
p_TMEV_dpi = num2str(cell2mat(p_TMEV_2_dpi(:,5)));
p_TMEV_ROI = num2str(cell2mat(p_TMEV_2_dpi(:,2)));
p_TMEV_2_dpi_labels = cell(size(p_TMEV_videos,1),size(p_TMEV_videos,2));
for i = 1:length(p_TMEV_2_dpi_labels)
p_TMEV_2_dpi_labels(i) = strcat(p_TMEV_videos(i),' , ', p_TMEV_dpi(i,:),' , ', p_TMEV_ROI(i,:));
end

p_TMEV_videos = p_TMEV_5_dpi(:,1);
p_TMEV_dpi = num2str(cell2mat(p_TMEV_5_dpi(:,5)));
p_TMEV_ROI = num2str(cell2mat(p_TMEV_5_dpi(:,2)));
p_TMEV_5_dpi_labels = cell(size(p_TMEV_videos,1),size(p_TMEV_videos,2));
for i = 1:length(p_TMEV_5_dpi_labels)
p_TMEV_5_dpi_labels(i) = strcat(p_TMEV_videos(i),' , ', p_TMEV_dpi(i,:),' , ', p_TMEV_ROI(i,:));
end

p_TMEV_videos = p_TMEV_14_dpi(:,1);
p_TMEV_dpi = num2str(cell2mat(p_TMEV_14_dpi(:,5)));
p_TMEV_ROI = num2str(cell2mat(p_TMEV_14_dpi(:,2)));
p_TMEV_14_dpi_labels = cell(size(p_TMEV_videos,1),size(p_TMEV_videos,2));
for i = 1:length(p_TMEV_14_dpi_labels)
p_TMEV_14_dpi_labels(i) = strcat(p_TMEV_videos(i),' , ', p_TMEV_dpi(i,:),' , ', p_TMEV_ROI(i,:));
end

%% PBS/TMEV/DPI Peak data ( somaes )
% % % TMEV % % %
p_TMEV_peak = peak_data(p_TMEV_logicals,:);
p_TMEV_peak = sortrows(p_TMEV_peak,4);
p_T_2_peak = p_TMEV_peak(1:idx_p_TMEV_2_dpi(length(idx_p_TMEV_2_dpi)),:);

p_T_5_peak = p_TMEV_peak(idx_p_TMEV_5_dpi(1):idx_p_TMEV_6_dpi(length(idx_p_TMEV_6_dpi)),:);

p_T_14_peak = p_TMEV_peak(idx_p_TMEV_14_dpi(1):end,:);

% % % PBS % % %
p_PBS_peak = peak_data(p_PBS_logicals,:);
p_PBS_peak = sortrows(p_PBS_peak,4);
p_P_2_5_peak = p_PBS_peak(1:idx_p_PBS_2_5_dpi(length(idx_p_PBS_2_5_dpi)),:);

p_P_6_14_peak = p_PBS_peak(idx_p_PBS_6_14_dpi(1):end,:);
%% Plotting PBS and TMEV %%
% figure()
% for jj = 1:size(p_PBS_2_5_dpi,1)
%     signal = cell2mat(p_PBS_2_5_dpi(jj,7)) + jj;
%     plot(signal);
%     hold on
%     
%     pks = p_P_2_5_peak{jj,2}.peaks + jj;  
%     locs = p_P_2_5_peak{jj,2}.locs;
%     plot(locs,pks,'o');
%     hold on
% end
% xlabel('Frames'); ylabel('Video,dpi,roi'); title('Soma: 2-5 DPI PBS');
% yticks([1:size(p_PBS_2_5_dpi,1)]);
% yticklabels(p_PBS_2_5_labels)
% % set(gca,'YDir','reverse');
% set(gca,'TickLabelInterpreter','none')
% 
% figure()
% for jj = 1:size(p_PBS_6_14_dpi,1)
%     signal = cell2mat(p_PBS_6_14_dpi(jj,7)) + jj;
%     plot(signal);
%     hold on
%     
%     pks = p_P_6_14_peak{jj,2}.peaks + jj;  
%     locs = p_P_6_14_peak{jj,2}.locs;
%     plot(locs,pks,'o');
%     hold on
% end
% xlabel('Frames'); ylabel('Video,dpi,roi'); title('Soma: 6-16 DPI PBS');
% yticks([1:size(p_PBS_6_14_dpi,1)]);
% yticklabels(p_PBS_6_14_labels)
% % set(gca,'YDir','reverse');
% set(gca,'TickLabelInterpreter','none')
% 
% figure()
% for jj = 1:size(p_TMEV_2_dpi,1)
%     signal = cell2mat(p_TMEV_2_dpi(jj,7)) + jj;
%     plot(signal);
%     hold on
%     
%     pks = p_T_2_peak{jj,2}.peaks + jj;  
%     locs = p_T_2_peak{jj,2}.locs;
%     plot(locs,pks,'o');
%     hold on
% end
% xlabel('Frames'); ylabel('Video,dpi,roi'); title('Soma: 2 DPI TMEV');
% yticks([1:size(p_TMEV_2_dpi,1)]);
% yticklabels(p_TMEV_2_dpi_labels)
% % set(gca,'YDir','reverse');
% set(gca,'TickLabelInterpreter','none')
% 
% figure()
% for jj = 1:size(p_TMEV_5_dpi,1)
%     signal = cell2mat(p_TMEV_5_dpi(jj,7)) + jj;
%     plot(signal);
%     hold on
%     
%     pks = p_T_5_peak{jj,2}.peaks + jj;  
%     locs = p_T_5_peak{jj,2}.locs;
%     plot(locs,pks,'o');
%     hold on
% end
% xlabel('Frames'); ylabel('Video,dpi,roi'); title('Soma: 5-6 DPI TMEV');
% yticks([1:size(p_TMEV_5_dpi,1)]);
% yticklabels(p_TMEV_5_dpi_labels)
% % set(gca,'YDir','reverse');
% set(gca,'TickLabelInterpreter','none')
% 
% 
% figure()
% for jj = 1:size(p_TMEV_14_dpi,1)
%     signal = cell2mat(p_TMEV_14_dpi(jj,7)) + jj;
%     plot(signal);
%     hold on
%     
%     pks = p_T_14_peak{jj,2}.peaks + jj;  
%     locs = p_T_14_peak{jj,2}.locs;
%     plot(locs,pks,'o');
%     hold on
% end
% xlabel('Frames'); ylabel('Video,dpi,roi'); title('Soma: 14-16 DPI TMEV');
% yticks([1:size(p_TMEV_14_dpi,1)]);
% yticklabels(p_TMEV_14_dpi_labels)
% % set(gca,'YDir','reverse');
% set(gca,'TickLabelInterpreter','none')