function [ALL_BINARY, process_data, process_peak_data] = plot_peak_locs_rev2(process_data, soma_data)
% This function plots the shifted peak location data of burn calcium
% responses.

%% generate process and peak data
[process_data, ~, process_peak_data, process_binary] = plot_burns_processes(process_data);
fs = 1.04; % frames per second (sampling rate)
process_timehitR = cell2mat(process_data(:,6)); % extract timehitR to shift peak locs
process_timehitR = round(process_timehitR.*fs);
process_frames = max(cellfun('size',process_data(:,3),2));
N = 2;
% % Initialize matrices for dF/F and binary % %
process_binary_shifted = zeros(size(process_data,1), process_frames*N);
process_dF_shifted = zeros(size(process_data,1), process_frames*N);
for ii = 1:size(process_data,1)
    k = process_timehitR(ii); % temporary shifting index
    if k > process_frames
        process_peak_data{ii,5} = (process_peak_data{ii,2}.locs - k + 1740 - 903)./fs; % shifting
    else
        process_peak_data{ii,5} = (process_peak_data{ii,2}.locs - k - 903)./fs; % shifting
    end
    process_peak_data{ii,6} = process_peak_data{ii,2}.peaks; % move peak heights out of struct
    idx = process_peak_data{ii,2}.locs; % locations where binary = 1
    process_binary(ii,idx) = process_peak_data{ii,5}; % replace 1's with shifted peak locs
    % Call user-defined function for shifting binary trace padded w/ zeros.
    [temp_b_shift, ~] = traceShift(process_binary(ii,:), k, N, process_frames);
    process_binary_shifted(ii,1:length(temp_b_shift)) = temp_b_shift; 
    % Create new dF/F
    tempF = process_data{ii,3};    
    F0 = mean(tempF(1,1:50));
    tempF(1,:) = (tempF(1,:) - F0)./F0;
    [temp_dF_shift, ~] = traceShift(tempF, k, N, process_frames);
    process_dF_shifted(ii,1:length(temp_dF_shift)) = temp_dF_shift;
end

for k = 1:size(process_binary_shifted, 1)
    process_data{k,7} = process_dF_shifted(k,:);
    process_data{k,8} = process_binary_shifted(k,:);
end

% extract logical indices for TMEV/PBS
TMEV_logicals = contains(process_data(:,4), 'TMEV');
PBS_logicals = contains(process_data(:,4), 'PBS');

% sort data by TMEV/PBS
TMEV_process_data = process_peak_data(TMEV_logicals, :);
TMEV_process = process_data(TMEV_logicals, :);
PBS_process_data = process_peak_data(PBS_logicals, :);
PBS_process = process_data(PBS_logicals, :);

% sort data by dpi
PBS_process = sortrows(PBS_process, 5);
PBS_process_data = sortrows(PBS_process_data, 4);
% sort data by dpi
TMEV_process = sortrows(TMEV_process, 5);
TMEV_process_data = sortrows(TMEV_process_data, 4);

% extract dF/F and binarized/thresholded traces
TMEV_process_dF = TMEV_process(:, 7);
PBS_process_dF = PBS_process(:, 7);

% binary
TMEV_process_binary = zeros(size(TMEV_process, 1), process_frames*N);
for jj = 1:size(TMEV_process, 1)
TMEV_process_binary(jj,:) = cell2mat(TMEV_process(jj,8));
end

PBS_process_binary = zeros(size(PBS_process, 1), process_frames*N);
for jj = 1:size(PBS_process, 1)
PBS_process_binary(jj,:) = cell2mat(PBS_process(jj,8));
end

TMEV_process_binary = process_binary_shifted(TMEV_logicals, :);
PBS_process_binary = process_binary_shifted(PBS_logicals, :);

% % % split PBS into time points 2 DPI, 5/6, 14/16 % % %
PBS_process_data = sortrows(PBS_process_data, 4);

idx_p_PBS_2_dpi = find(cell2mat(PBS_process_data(:,4)) == 2);
p_PBS_2_dpi = PBS_process_binary(1:idx_p_PBS_2_dpi(end),:);
p_PBS_2_NAMES = PBS_process_data(1:idx_p_PBS_2_dpi(end),1)';
p_PBS_2_dF = PBS_process_dF(1:idx_p_PBS_2_dpi(end),:);
p_PBS_2_peaks = PBS_process_data(1:idx_p_PBS_2_dpi(end),5:6);

ALL_BINARY.PBS.Process_2_DPI = p_PBS_2_dpi;
ALL_BINARY.PBS.Process_2_DPI_PEAKS = p_PBS_2_peaks;
test_time = (-process_frames*(N/2):process_frames*(N/2))./fs;
test_time = test_time(2:end);
% set(gca,'YDir','reverse');
% plot traces and peak locs
figure()
for jj = 1:size(p_PBS_2_dF, 1)
    L = length(p_PBS_2_dF{jj})/(2*N);
    locs = p_PBS_2_peaks{jj, 1};
    pks = p_PBS_2_peaks{jj, 2};
    plot(test_time, p_PBS_2_dF{jj} + jj)
    hold on
    plot(locs+L, pks + jj, 'o')
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('PBS Process 2 DPI Shifted Traces');
yticks(1:size(p_PBS_2_dpi,1));
yticklabels(p_PBS_2_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on
figure()
for jj = 1:size(p_PBS_2_dpi,1)
    L = length(p_PBS_2_dF{jj})/(2*N);
    % % Get peak data % %
    locs = p_PBS_2_peaks{jj, 1};
    binary_pks = ones(size(locs));
    % % Replace non-zero elements with 1's % %
    idx = p_PBS_2_dpi(jj,:) ~= 0;
    p_PBS_2_dpi(jj,idx) = 1;
    hold on
    plot(test_time, p_PBS_2_dpi(jj,:) + jj);
    hold on
    plot(locs+L, binary_pks + jj, 'd');
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('PBS Process 2 DPI Shifted Binary Traces');
yticks(1:size(p_PBS_2_dpi,1));
yticklabels(p_PBS_2_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on

idx_p_PBS_5_dpi = find(cell2mat(PBS_process_data(:,4)) == 5);
idx_p_PBS_6_dpi = find(cell2mat(PBS_process_data(:,4)) == 6);
p_PBS_5_dpi = PBS_process_binary(idx_p_PBS_5_dpi(1):idx_p_PBS_6_dpi(end),:);
p_PBS_5_NAMES = PBS_process_data(idx_p_PBS_5_dpi(1):idx_p_PBS_6_dpi(end),1);
p_PBS_5_dF = PBS_process_dF(idx_p_PBS_5_dpi(1):idx_p_PBS_6_dpi(end),:);
p_PBS_5_peaks = PBS_process_data(idx_p_PBS_5_dpi(1):idx_p_PBS_6_dpi(end),5:6);
ALL_BINARY.PBS.Process_5_DPI_PEAKS = p_PBS_5_peaks;
ALL_BINARY.PBS.Process_5_DPI = p_PBS_5_dpi;
figure()
for jj = 1:size(p_PBS_5_dF, 1)
    L = length(p_PBS_5_dF{jj})/(2*N);
    locs = p_PBS_5_peaks{jj, 1};
    pks = p_PBS_5_peaks{jj, 2};
    plot(test_time, p_PBS_5_dF{jj} + jj)
    hold on
    plot(locs+L, pks + jj, 'o')
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('PBS Process 5-6 DPI Shifted Traces');
yticks(1:size(p_PBS_5_dpi,1));
yticklabels(p_PBS_5_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on

figure()
for jj = 1:size(p_PBS_5_dpi,1)
    L = length(p_PBS_5_dF{jj})/(2*N);
    locs = p_PBS_5_peaks{jj, 1};
    binary_pks = ones(size(locs));
    idx = p_PBS_5_dpi(jj,:) ~= 0;
    p_PBS_5_dpi(jj,idx) = 1;
    hold on
    plot(test_time, p_PBS_5_dpi(jj,:) + jj);
    hold on
    plot(locs+L, binary_pks + jj, 'd');
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('PBS Process 5-6 DPI Shifted Binary Traces');
yticks(1:size(p_PBS_5_dpi,1));
yticklabels(p_PBS_5_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on

idx_p_PBS_14_dpi = find(cell2mat(PBS_process_data(:,4)) == 14);
p_PBS_14_dpi = PBS_process_binary(idx_p_PBS_14_dpi(1):end,:);
p_PBS_14_NAMES = PBS_process_data(idx_p_PBS_14_dpi(1):end,1)';
p_PBS_14_dF = PBS_process_dF(idx_p_PBS_14_dpi(1):end,:);
p_PBS_14_peaks = PBS_process_data(idx_p_PBS_14_dpi(1):end,5:6);
ALL_BINARY.PBS.Process_14_DPI = p_PBS_14_dpi;
ALL_BINARY.PBS.Process_14_DPI_PEAKS = p_PBS_14_peaks;
figure()
for jj = 1:size(p_PBS_14_dF, 1)
    L = length(p_PBS_14_dF{jj})/(2*N);
    locs = p_PBS_14_peaks{jj, 1};
    pks = p_PBS_14_peaks{jj, 2};
    plot(test_time, p_PBS_14_dF{jj} + jj)
    hold on
    plot(locs+L, pks + jj, 'o')
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('PBS Process 14-16 DPI Shifted Traces');
yticks(1:size(p_PBS_14_dpi,1));
yticklabels(p_PBS_14_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on

figure()
for jj = 1:size(p_PBS_14_dpi,1)
    L = length(p_PBS_14_dF{jj})/(2*N);
    locs = p_PBS_14_peaks{jj, 1};
    binary_pks = ones(size(locs));
    idx = p_PBS_14_dpi(jj,:) ~= 0;
    p_PBS_14_dpi(jj,idx) = 1;
    hold on
    plot(test_time, p_PBS_14_dpi(jj,:) + jj);
    hold on
    plot(locs+L, binary_pks + jj, 'd');
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('PBS Process 14-16 DPI Shifted Binary Traces');
yticks(1:size(p_PBS_14_dpi,1));
yticklabels(p_PBS_14_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on

% % % split TMEV into 2 dpi, 5/6 dpi, 14/16 dpi
idx_p_TMEV_2_dpi = find(cell2mat(TMEV_process_data(:,4)) == 2);
p_TMEV_2_dpi = TMEV_process_binary(1:idx_p_TMEV_2_dpi(end),:);
p_TMEV_2_NAMES = TMEV_process_data(1:idx_p_TMEV_2_dpi(end),1)';
p_TMEV_2_dF = TMEV_process_dF(1:idx_p_TMEV_2_dpi(end),:);
p_TMEV_2_peaks = TMEV_process_data(1:idx_p_TMEV_2_dpi(end),5:6);
ALL_BINARY.TMEV.Process_2_DPI = p_TMEV_2_dpi;
ALL_BINARY.TMEV.Process_2_DPI_PEAKS = p_TMEV_2_peaks;
figure()
for jj = 1:size(p_TMEV_2_dF, 1)
    L = length(p_TMEV_2_dF{jj})/(2*N);
    locs = p_TMEV_2_peaks{jj, 1};
    pks = p_TMEV_2_peaks{jj, 2};
    plot(test_time, p_TMEV_2_dF{jj} + jj)
    hold on
    plot(locs+L, pks + jj, 'o')
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('TMEV Process 2 DPI Shifted Traces');
yticks(1:size(p_TMEV_2_dpi,1));
yticklabels(p_TMEV_2_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on
figure()
for jj = 1:size(p_TMEV_2_dpi,1)
    L = length(p_TMEV_2_dF{jj})/(2*N);
    locs = p_TMEV_2_peaks{jj, 1};
    binary_pks = ones(size(locs));
    idx = p_TMEV_2_dpi(jj,:) ~= 0;
    p_TMEV_2_dpi(jj,idx) = 1;
    hold on
    plot(test_time, p_TMEV_2_dpi(jj,:) + jj);
    hold on
    plot(locs+L, binary_pks + jj, 'd');
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('TMEV Process 2 DPI Shifted Binary Traces');
yticks(1:size(p_TMEV_2_dpi,1));
yticklabels(p_TMEV_2_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on

idx_p_TMEV_5_dpi = find(cell2mat(TMEV_process_data(:,4)) == 5);
idx_p_TMEV_6_dpi = find(cell2mat(TMEV_process_data(:,4)) == 6);
p_TMEV_5_dpi = TMEV_process_binary(idx_p_TMEV_5_dpi(1):idx_p_TMEV_6_dpi(end),:);
p_TMEV_5_NAMES = TMEV_process_data(idx_p_TMEV_5_dpi(1):idx_p_TMEV_6_dpi(end),1)';
p_TMEV_5_dF = TMEV_process_dF(idx_p_TMEV_5_dpi(1):idx_p_TMEV_6_dpi(end),:);
p_TMEV_5_peaks = TMEV_process_data(idx_p_TMEV_5_dpi(1):idx_p_TMEV_6_dpi(end),5:6);
ALL_BINARY.TMEV.Process_5_DPI = p_TMEV_5_dpi;
ALL_BINARY.TMEV.Process_5_DPI_PEAKS = p_TMEV_5_peaks;
figure()
for jj = 1:size(p_TMEV_5_dpi,1)
    L = length(p_TMEV_5_dF{jj})/(2*N);
    locs = p_TMEV_5_peaks{jj, 1};
    binary_pks = ones(size(locs));
    idx = p_TMEV_5_dpi(jj,:) ~= 0;
    p_TMEV_5_dpi(jj,idx) = 1;
    hold on
    plot(test_time, p_TMEV_5_dpi(jj,:) + jj);
    hold on
    plot(locs+L, binary_pks + jj, 'd');
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('TMEV Process 4-6 DPI Shifted Binary Traces');
yticks(1:size(p_TMEV_5_dpi,1));
yticklabels(p_TMEV_5_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on
figure()
for jj = 1:size(p_TMEV_5_dF, 1)
    L = length(p_TMEV_5_dF{jj})/(2*N);
    locs = p_TMEV_5_peaks{jj, 1};
    pks = p_TMEV_5_peaks{jj, 2};
    plot(test_time, p_TMEV_5_dF{jj} + jj)
    hold on
    plot(locs+L, pks + jj, 'o')
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('TMEV Process 4-6 DPI Shifted Traces');
yticks(1:size(p_TMEV_5_dpi,1));
yticklabels(p_TMEV_5_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on

idx_p_TMEV_14_dpi = find(cell2mat(TMEV_process_data(:,4)) == 14);
p_TMEV_14_dpi = TMEV_process_binary(idx_p_TMEV_14_dpi(1):end,:);
p_TMEV_14_NAMES = TMEV_process_data(idx_p_TMEV_14_dpi(1):end,1)';
p_TMEV_14_dF = TMEV_process_dF(idx_p_TMEV_14_dpi(1):end,:);
p_TMEV_14_peaks = TMEV_process_data(idx_p_TMEV_14_dpi(1):end,5:6);
ALL_BINARY.TMEV.Process_14_DPI = p_TMEV_14_dpi;
ALL_BINARY.TMEV.Process_14_DPI_PEAKS = p_TMEV_14_peaks;
figure()
for jj = 1:size(p_TMEV_14_dF, 1)
    L = length(p_TMEV_14_dF{jj})/(2*N);
    locs = p_TMEV_14_peaks{jj, 1};
    pks = p_TMEV_14_peaks{jj, 2};
    plot(test_time, p_TMEV_14_dF{jj} + jj)
    hold on
    plot(locs+L, pks + jj, 'o')
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('TMEV Process 14-16 DPI Shifted Traces');
yticks(1:size(p_TMEV_14_dpi,1));
yticklabels(p_TMEV_14_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on
figure()
for jj = 1:size(p_TMEV_14_dpi,1)
    L = length(p_TMEV_14_dF{jj})/(2*N);
    locs = p_TMEV_14_peaks{jj, 1};
    binary_pks = ones(size(locs));    
    idx = p_TMEV_14_dpi(jj,:) ~= 0;
    p_TMEV_14_dpi(jj,idx) = 1;
    hold on
    plot(test_time, p_TMEV_14_dpi(jj,:) + jj);
    hold on
    plot(locs+L, binary_pks + jj, 'd');
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('TMEV Process 14-16 DPI Shifted Binary Traces');
yticks(1:size(p_TMEV_14_dpi,1));
yticklabels(p_TMEV_14_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on
%% generate soma and peak data
[soma_data, ~, soma_peak_data, soma_binary] = burn_plots_soma(soma_data);
fs = 1.04; % frames per second (sampling rate)
soma_timehitR = cell2mat(soma_data(:,6)); % extract timehitR to shift peak locs
soma_timehitR = round(soma_timehitR.*fs);
soma_frames = max(cellfun('size',soma_data(:,3),2));
N = 2;
test_time = (-soma_frames*(N/2):soma_frames*(N/2))./fs;
test_time = test_time(2:end);

% % Initialize matrices for dF/F and binary % %
soma_binary_shifted = zeros(size(soma_data,1), soma_frames*N);
soma_dF_shifted = zeros(size(soma_data,1), soma_frames*N);
for ii = 1:size(soma_data,1)
    k = soma_timehitR(ii); % temporary shifting index
    if k > soma_frames
        soma_peak_data{ii,5} = (soma_peak_data{ii,2}.locs - k + 1740)./fs; % shifting
    else
        soma_peak_data{ii,5} = (soma_peak_data{ii,2}.locs - k)./fs; % shifting
    end
    soma_peak_data{ii,6} = soma_peak_data{ii,2}.peaks; % move peak heights out of struct
    idx = soma_peak_data{ii,2}.locs; % locations where binary = 1
    soma_binary(ii,idx) = soma_peak_data{ii,5}; % replace 1's with shifted peak locs
    % Call user-defined function for shifting binary trace padded w/ zeros.
    [temp_b_shift, ~] = traceShift(soma_binary(ii,:), k, N, soma_frames);
    soma_binary_shifted(ii,1:length(temp_b_shift)) = temp_b_shift; 
    % Create new dF/F
    tempF = soma_data{ii,3};    
    F0 = mean(tempF(1,1:50));
    tempF(1,:) = (tempF(1,:) - F0)./F0;
    [temp_dF_shift, ~] = traceShift(tempF, k, N, soma_frames);
    soma_dF_shifted(ii,1:length(temp_dF_shift)) = temp_dF_shift;
end

for k = 1:size(soma_binary_shifted, 1)
    soma_data{k,7} = soma_dF_shifted(k,:);
    soma_data{k,8} = soma_binary_shifted(k,:);
end

% extract logical indices for TMEV/PBS
TMEV_logicals = contains(soma_data(:,4), 'TMEV');
PBS_logicals = contains(soma_data(:,4), 'PBS');

% sort data by TMEV/PBS
TMEV_soma_data = soma_peak_data(TMEV_logicals, :);
TMEV_soma = soma_data(TMEV_logicals, :);
PBS_soma_data = soma_peak_data(PBS_logicals, :);
PBS_soma = soma_data(PBS_logicals, :);

% sort data by dpi
PBS_soma = sortrows(PBS_soma, 5);
PBS_soma_data = sortrows(PBS_soma_data, 4);
% sort data by dpi
TMEV_soma = sortrows(TMEV_soma, 5);
TMEV_soma_data = sortrows(TMEV_soma_data, 4);

% extract dF/F and binarized/thresholded traces
TMEV_soma_dF = TMEV_soma(:, 7);
PBS_soma_dF = PBS_soma(:, 7);
% binary
TMEV_soma_binary = zeros(size(TMEV_soma, 1), soma_frames*2);
for jj = 1:size(TMEV_soma, 1)
TMEV_soma_binary(jj,:) = cell2mat(TMEV_soma(jj,8));
end

PBS_soma_binary = zeros(size(PBS_soma, 1), soma_frames*2);
for jj = 1:size(PBS_soma, 1)
PBS_soma_binary(jj,:) = cell2mat(PBS_soma(jj,8));
end

% % % split PBS into time points 2 DPI, 5/6, 14/16 % % %

idx_s_PBS_2_dpi = find(cell2mat(PBS_soma_data(:,4)) == 2);
s_PBS_2_dpi = PBS_soma_binary(1:idx_s_PBS_2_dpi(end),:);
s_PBS_2_NAMES = PBS_soma_data(1:idx_s_PBS_2_dpi(end),1)';
s_PBS_2_dF = PBS_soma_dF(1:idx_s_PBS_2_dpi(end),:);
s_PBS_2_peaks = PBS_soma_data(1:idx_s_PBS_2_dpi(end),5:6);
ALL_BINARY.PBS.Soma_2_DPI = s_PBS_2_dpi;
ALL_BINARY.PBS.Soma_2_DPI_PEAKS = s_PBS_2_peaks;
% plot traces and peak locs
figure()
for jj = 1:size(s_PBS_2_dF, 1)
    locs = s_PBS_2_peaks{jj, 1};
    pks = s_PBS_2_peaks{jj, 2};
    plot(test_time, s_PBS_2_dF{jj} + jj)
    hold on
    plot(locs, pks + jj, 'o')
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('PBS Soma 2 DPI Shifted Traces');
yticks(1:size(s_PBS_2_dpi,1));
yticklabels(s_PBS_2_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on
figure()
for jj = 1:size(s_PBS_2_dpi,1)
    locs = s_PBS_2_peaks{jj, 1};
    binary_pks = ones(size(locs));
    idx = s_PBS_2_dpi(jj,:) ~= 0;
    s_PBS_2_dpi(jj,idx) = 1;
    hold on
    plot(test_time, s_PBS_2_dpi(jj,:) + jj);
    hold on
    plot(locs, binary_pks + jj, 'd');
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('PBS Soma 2 DPI Shifted Binary Traces');
yticks(1:size(s_PBS_2_dpi,1));
yticklabels(s_PBS_2_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on

idx_s_PBS_5_dpi = find(cell2mat(PBS_soma_data(:,4)) == 5);
idx_s_PBS_6_dpi = find(cell2mat(PBS_soma_data(:,4)) == 6);
s_PBS_5_dpi = PBS_soma_binary(idx_s_PBS_5_dpi(1):idx_s_PBS_6_dpi(end),:);
s_PBS_5_NAMES = PBS_soma_data(idx_s_PBS_5_dpi(1):idx_s_PBS_6_dpi(end),1);
s_PBS_5_dF = PBS_soma_dF(idx_s_PBS_5_dpi(1):idx_s_PBS_6_dpi(end),:);
s_PBS_5_peaks = PBS_soma_data(idx_s_PBS_5_dpi(1):idx_s_PBS_6_dpi(end),5:6);
ALL_BINARY.PBS.Soma_5_DPI = s_PBS_5_dpi;
ALL_BINARY.PBS.Soma_5_DPI_PEAKS = s_PBS_5_peaks;
% plot traces and peak locs
figure()
for jj = 1:size(s_PBS_5_dF, 1)
    locs = s_PBS_5_peaks{jj, 1};
    pks = s_PBS_5_peaks{jj, 2};
    plot(test_time, s_PBS_5_dF{jj} + jj)
    hold on
    plot(locs, pks + jj, 'o')
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('PBS Soma 5-6 DPI Shifted Traces');
yticks(1:size(s_PBS_5_dpi,1));
yticklabels(s_PBS_5_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on
figure()
for jj = 1:size(s_PBS_5_dpi,1)
    locs = s_PBS_5_peaks{jj, 1};
    binary_pks = ones(size(locs));
    idx = s_PBS_5_dpi(jj,:) ~= 0;
    s_PBS_5_dpi(jj,idx) = 1;
    hold on
    plot(test_time, s_PBS_5_dpi(jj,:) + jj);
    hold on
    plot(locs, binary_pks + jj, 'd');
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('PBS Soma 5-6 DPI Shifted Binary Traces');
yticks(1:size(s_PBS_5_dpi,1));
yticklabels(s_PBS_5_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on

idx_s_PBS_14_dpi = find(cell2mat(PBS_soma_data(:,4)) == 14);
s_PBS_14_dpi = PBS_soma_binary(idx_s_PBS_14_dpi(1):end,:);
s_PBS_14_NAMES = PBS_soma_data(idx_s_PBS_14_dpi(1):end,1)';
s_PBS_14_dF = PBS_soma_dF(idx_s_PBS_14_dpi(1):end,:);
s_PBS_14_peaks = PBS_soma_data(idx_s_PBS_14_dpi(1):end,5:6);
ALL_BINARY.PBS.Soma_14_DPI = s_PBS_14_dpi;
ALL_BINARY.PBS.Soma_14_DPI_PEAKS = s_PBS_14_peaks;
% plot traces and peak locs
figure()
for jj = 1:size(s_PBS_14_dF, 1)
    locs = s_PBS_14_peaks{jj, 1};
    pks = s_PBS_14_peaks{jj, 2};
    plot(test_time, s_PBS_14_dF{jj} + jj)
    hold on
    plot(locs, pks + jj, 'o')
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('PBS Soma 14-16 DPI Shifted Traces');
yticks(1:size(s_PBS_14_dpi,1));
yticklabels(s_PBS_14_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on
figure()
for jj = 1:size(s_PBS_14_dpi,1)
    locs = s_PBS_14_peaks{jj, 1};
    binary_pks = ones(size(locs));
    idx = s_PBS_14_dpi(jj,:) ~= 0;
    s_PBS_14_dpi(jj,idx) = 1;
    hold on
    plot(test_time, s_PBS_14_dpi(jj,:) + jj);
    hold on
    plot(locs, binary_pks + jj, 'd');
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('PBS Soma 14-16 DPI Shifted Binary Traces');
yticks(1:size(s_PBS_14_dpi,1));
yticklabels(s_PBS_14_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on

% % % split TMEV into 2 dpi, 5/6 dpi, 14/16 dpi
idx_s_TMEV_2_dpi = find(cell2mat(TMEV_soma_data(:,4)) == 2);
s_TMEV_2_dpi = TMEV_soma_binary(1:idx_s_TMEV_2_dpi(end),:);
s_TMEV_2_NAMES = TMEV_soma_data(1:idx_s_TMEV_2_dpi(end),1)';
s_TMEV_2_dF = TMEV_soma_dF(1:idx_s_TMEV_2_dpi(end),:);
s_TMEV_2_peaks = TMEV_soma_data(1:idx_s_TMEV_2_dpi(end),5:6);
ALL_BINARY.TMEV.Soma_2_DPI = s_TMEV_2_dpi;
ALL_BINARY.TMEV.Soma_2_DPI_PEAKS = s_TMEV_2_peaks;
% plot traces and peak locs
figure()
for jj = 1:size(s_TMEV_2_dF, 1)
    locs = s_TMEV_2_peaks{jj, 1};
    pks = s_TMEV_2_peaks{jj, 2};
    plot(test_time, s_TMEV_2_dF{jj} + jj)
    hold on
    plot(locs, pks + jj, 'o')
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('TMEV Soma 2 DPI Shifted Traces');
yticks(1:size(s_TMEV_2_dpi,1));
yticklabels(s_TMEV_2_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on
figure()
for jj = 1:size(s_TMEV_2_dpi,1)
    locs = s_TMEV_2_peaks{jj, 1};
    binary_pks = ones(size(locs));
    idx = s_TMEV_2_dpi(jj,:) ~= 0;
    s_TMEV_2_dpi(jj,idx) = 1;
    hold on
    plot(test_time, s_TMEV_2_dpi(jj,:) + jj);
    hold on
    plot(locs, binary_pks + jj, 'd');
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('TMEV Soma 2 DPI Shifted Binary Traces');
yticks(1:size(s_TMEV_2_dpi,1));
yticklabels(s_TMEV_2_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on

idx_s_TMEV_5_dpi = find(cell2mat(TMEV_soma_data(:,4)) == 5);
idx_s_TMEV_6_dpi = find(cell2mat(TMEV_soma_data(:,4)) == 6);
s_TMEV_5_dpi = TMEV_soma_binary(idx_s_TMEV_5_dpi(1):idx_s_TMEV_6_dpi(end),:);
s_TMEV_5_NAMES = TMEV_soma_data(idx_s_TMEV_5_dpi(1):idx_s_TMEV_6_dpi(end),1);
s_TMEV_5_dF = TMEV_soma_dF(idx_s_TMEV_5_dpi(1):idx_s_TMEV_6_dpi(end),:);
s_TMEV_5_peaks = TMEV_soma_data(idx_s_TMEV_5_dpi(1):idx_s_TMEV_6_dpi(end),5:6);
ALL_BINARY.TMEV.Soma_5_DPI = s_TMEV_5_dpi;
ALL_BINARY.TMEV.Soma_5_DPI_PEAKS = s_TMEV_5_peaks;
% plot traces and peak locs
figure()
for jj = 1:size(s_TMEV_5_dF, 1)
    locs = s_TMEV_5_peaks{jj, 1};
    pks = s_TMEV_5_peaks{jj, 2};
    plot(test_time, s_TMEV_5_dF{jj} + jj)
    hold on
    plot(locs, pks + jj, 'o')
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('TMEV Soma 5-6 DPI Shifted Traces');
yticks(1:size(s_TMEV_5_dpi,1));
yticklabels(s_TMEV_5_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on
figure()
for jj = 1:size(s_TMEV_5_dpi,1)
    locs = s_TMEV_5_peaks{jj, 1};
    binary_pks = ones(size(locs));
    idx = s_TMEV_5_dpi(jj,:) ~= 0;
    s_TMEV_5_dpi(jj,idx) = 1;
    hold on
    plot(test_time, s_TMEV_5_dpi(jj,:) + jj);
    hold on
    plot(locs, binary_pks + jj, 'd');
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('TMEV Soma 4-6 DPI Shifted Binary Traces');
yticks(1:size(s_TMEV_5_dpi,1));
yticklabels(s_TMEV_5_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on

idx_s_TMEV_14_dpi = find(cell2mat(TMEV_soma_data(:,4)) == 14);
s_TMEV_14_dpi = TMEV_soma_binary(idx_s_TMEV_14_dpi(1):end,:);
s_TMEV_14_NAMES = TMEV_soma_data(idx_s_TMEV_14_dpi(1):end,1)';
s_TMEV_14_dF = TMEV_soma_dF(idx_s_TMEV_14_dpi(1):end,:);
s_TMEV_14_peaks = TMEV_soma_data(idx_s_TMEV_14_dpi(1):end,5:6);
ALL_BINARY.TMEV.Soma_14_DPI = s_TMEV_14_dpi;
ALL_BINARY.TMEV.Soma_14_DPI_PEAKS = s_TMEV_14_peaks;
% plot traces and peak locs
figure()
for jj = 1:size(s_TMEV_14_dF, 1)
    locs = s_TMEV_14_peaks{jj, 1};
    pks = s_TMEV_14_peaks{jj, 2};
    plot(test_time, s_TMEV_14_dF{jj} + jj)
    hold on
    plot(locs, pks + jj, 'o')
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('TMEV Soma 14-16 DPI Shifted Traces');
yticks(1:size(s_TMEV_14_dpi,1));
yticklabels(s_TMEV_14_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on
figure()
for jj = 1:size(s_TMEV_14_dpi,1)
    locs = s_TMEV_14_peaks{jj, 1};
    binary_pks = ones(size(locs));
    idx = s_TMEV_14_dpi(jj,:) ~= 0;
    s_TMEV_14_dpi(jj,idx) = 1;
    hold on
    plot(test_time, s_TMEV_14_dpi(jj,:) + jj);
    hold on
    plot(locs, binary_pks + jj, 'd');
end
xlabel('Relative Time (s)'); ylabel('Trace'); title('TMEV Soma 14-16 DPI Shifted Binary Traces');
yticks(1:size(s_TMEV_14_dpi,1));
yticklabels(s_TMEV_14_NAMES)
set(gca,'TickLabelInterpreter','none')
grid on

end
