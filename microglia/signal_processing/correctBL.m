%% This will be a script to test baseline correction in our signals
% First gotta load some data
load process_data.mat
% test_trace = process_data{18,3}; % RAW TRACE 2043_3
test_trace = soma_data{106,3};      % RAW_TRACE 2167_1
% test_trace = process_data{109,3}; % RAW TRACE 2102_10
%%
% % % Start with detrend()
detrended_BL = detrend(test_trace);

figure()
plot(test_trace)
hold on
plot(detrended_BL)
legend('Raw trace', 'Detrended data');
grid on

Raw_stdev = std(test_trace);
DT_stdev = std(detrended_BL);

%%
lag = 30; %window
threshold = 5.0; % threshold in standard deviations
influence = 0.99;
% % % Try polyfit()
x = 1:length(test_trace); % get time vector
[p, S] = polyfit(x, test_trace, 2); % quadratic fit

poly_corrected = polyval(p, test_trace);
[quadratic_signal,avgFilter,stdFilter] = peak_detect(poly_corrected,lag,threshold,influence);

% Filter the polynomial-fitted signal
%%% Amplitude already attenuated so use Savitzky Golay filter
poly_filtered = sgolayfilt(double(poly_corrected), 2, 3);
%%% Try moving avg filter
window = 2; % arbitrary window width
kernel = ones(1,window)/window;
avg_filter = filter(kernel, 1, poly_corrected);

figure()
plot(x, test_trace)
hold on
plot(x, poly_corrected)
hold on
% plot(x,quadratic_signal*100)
plot(x, poly_filtered);
hold on
plot(x,avg_filter);
legend('Raw trace', 'Polyfit n=2','SG Filtered', 'Mean Filter');
grid on



[p, S] = polyfit(x, test_trace, 4); % 4th order fit

poly_corrected = polyval(p, test_trace);
[fourthO_signal,avgFilter,stdFilter] = peak_detect(poly_corrected,lag,threshold,influence);
figure()
plot(x, test_trace)
hold on
plot(x, poly_corrected)
hold on
plot(x, fourthO_signal*100)
legend('Raw trace', 'Polyfit n=4','Binarized');
grid on
%%
% % % Time to put on the big kid pants and revisit fourier transforms
% *gulp*
Fs = 1.04; % FPS / sampling 
Fn = Fs/2; % Nyquist frequency - might try 2.5x sampling interval
time = x./Fs; % time (s)
L = numel(time); % length

FT_trace = fft((test_trace-mean(test_trace)))/L; % fast FT / num elements
Fv = linspace(0, 1, fix(L/2)+1)*Fn; % freq. vector
Iv = 1:numel(Fv);

% Visualize fequency spectrum
figure()
plot(Fv, abs(FT_trace(Iv))*2)
grid on
set(gca, 'XMinorTick','on')
xlabel('Frequency')


Wp = [0.0036 0.01]/Fn; % Passband Frequency
Ws = [0.0030 0.025]/Fn; % Stopband frequency
Rs = 0.002; % Passband Ripple
Rs = 0.5; % stopband ripple
[n, Ws] = cheb2ord(Wp,Ws,Rp,Rs); % Filter order
[z, p, k] = cheby2(abs(n),Rs,Ws);
[sos, g] = zp2sos(z,p,k);

figure()
freqz(sos, 2^16, Fs) % bode plot

trace_filtered = filtfilt(sos, g, double(test_trace));

figure()
plot(time, trace_filtered)
grid on
hold on
plot(time, test_trace)
legend('Filtered', 'Raw');
    
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