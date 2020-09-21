%%
%Function  to find location of peaks with a prominence above 
% desired standard deviations

%%

%input variables
    % data = the trace to be analyzed
    % Fs = sampling rate or frames per second

%output variables
    %pks = peak amplitude of the events
    %locs = index in seconds
    %widths = temporal full-width at half the pks
    %proms = amplitude of events usig prominence instead of absolute peak
    
% Example to call function
    % [Peak,Index_s,FWHM_s,Prominence] = FindPeaks_Stim(Trace_50,500)
    
%%
function [pks,locs,widths,proms] = FindPeaks_Stim_v2(data ,Fs, Input_SDs, timescale)

% % Commenting this out so we can fix SD in the main script - remove
% % intervention step
%     %prompt for # of ROIs
%     prompt = {'Prominence of Peak (No. SDs)'};
%     dlg_title = 'Input Variables';
%     num_lines = 1;
%     defaultans = {'10'};
%     InputVars = str2double(inputdlg(prompt,dlg_title,num_lines,defaultans));
% 
%     %pull out input variables
%     Input_SDs = InputVars(1,1);%Maximum % of peak desired 
    %find peaks
    [~, ~, ~, proms] = findpeaks(data); %find all peaks for individual trace
%     figure,plot(data);
    MeanProms = mean(proms);%mean prominence of all peaks for individual trace
    StdProms = std(proms);%standard deviation of prominence of all peaks for individual trace
    UseProms = MeanProms + StdProms*Input_SDs;%determine N standard deviations
%     findpeaks(data,Fs,'MinPeakProminence',UseProms,'Annotate','extents') %annotate peaks on individual figures
    clear pks locs widths proms 
    if strcmp(timescale, 'Frames') == 1
    [pks,locs,widths,proms] = findpeaks(data,'MinPeakProminence',UseProms,'MinPeakDistance', 15, 'MinPeakWidth', 2.5); %find peaks for individual trace
    elseif strcmp(timescale, 'Seconds') == 1
    [pks,locs,widths,proms] = findpeaks(data,Fs,'MinPeakProminence',UseProms,'MinPeakDistance', 25); %find peaks for individual trace
    else
        error('Please select a valid time scale to find peaks.');
    end  
end