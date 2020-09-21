function [raw_auc, filtered_auc, peakDurations, peakWidth, interEventIntervals, pks, locs, filtFig, removed,removedValleys,addedValleys] = peakMeasure(raw_signal, detrended_signal, filtered_signal, Valley_SDs,Peak_SDs,traceName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% This function extracts characteristics about peaks in a single %
%%%%%%%%%% column or row vector signal. Intended for use in identifying   %
%%%%%%%%%% peaks and valleys in microglial calcium signaling data.        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Use detrended signal to find peak locations and filtered signal to
%%%%%% find valley locations.
[~, ~, ~, proms] = findpeaks(detrended_signal);

[valleys, ~, ~, ~] = findpeaks(-filtered_signal);
MeanProms = mean(proms);%mean prominence of all peaks for individual trace
StdProms = std(proms);%standard deviation of prominence of all peaks for individual trace
UseProms = MeanProms + StdProms*Peak_SDs;%determine N standard deviations

MeanValleys = mean(valleys);%mean prominence of all peaks for individual trace
StdValleys = std(valleys);%standard deviation of prominence of all peaks for individual trace
UseValleys = MeanValleys + StdValleys*Valley_SDs;%determine N standard deviations

[valleys, vly_locs, ~, ~] = findpeaks(-filtered_signal);
% [valleys, vly_locs] = peakfinder(-filtered_signal);
[pks,locs,peakWidth,proms] = findpeaks(filtered_signal, 'MinPeakProminence', UseProms, 'MinPeakDistance', 10);
% % % Take a peak location then attempt to find the two nearest valley locs
% interEventIntervals = diff(locs);
peakDurations = zeros(size(pks));
filtered_auc = zeros(size(pks));
raw_auc = zeros(size(pks));

% % % Throw out valleys that exceed a baseline threshold; we prefer to keep
% % % valleys on the baseline not raised on a peak's plateau
raw_valleys = raw_signal(vly_locs);
valleyMean = mean(raw_valleys);
valleyStd = std(raw_valleys);
flag_valleys = valleyMean + valleyStd*1.5;
vly_locs = vly_locs(raw_valleys < flag_valleys);
valleys = valleys(raw_valleys < flag_valleys);

% visualize results
filtFig = figure();
plot(filtered_signal)
hold on
plot(vly_locs, -valleys, 'r*')
hold on
plot(locs, filtered_signal(locs), 'm*');
linkdata on

h=figure;
axh = axes('Parent',h); % parent axes
true_valleys = raw_signal(vly_locs);
figure(h)
plot(raw_signal)
hold on
valleyLine = plot(axh, vly_locs, true_valleys, 'r*','XDataSource','vly_locs','YDataSource','true_valleys');
hold on
rawLine = plot(axh, locs, raw_signal(locs), 'm*');
legend('signal','valleys','peaks');
title(traceName, 'interpreter', 'none');
linkdata on
% % % Maximize figure w/ robot % % %
robot = java.awt.Robot; 
robot.keyPress(java.awt.event.KeyEvent.VK_ALT);      %// send ALT
robot.keyPress(java.awt.event.KeyEvent.VK_SPACE);    %// send SPACE
robot.keyRelease(java.awt.event.KeyEvent.VK_SPACE);  %// release SPACE
robot.keyRelease(java.awt.event.KeyEvent.VK_ALT);    %// release ALT
robot.keyPress(java.awt.event.KeyEvent.VK_X);        %// send X
robot.keyRelease(java.awt.event.KeyEvent.VK_X);      %// release X

%%% User selected valleys will be removed %%%
valley_removals = str2double(inputdlg('How many valleys would you like to remove?'));
while isnan(valley_removals)
    valley_removals = str2double(inputdlg('Please input a valid integer to remove valleys.'));
end
removedValleys = zeros(valley_removals, 1);
if valley_removals > 0 && isnumeric(valley_removals)
    for i = 1:valley_removals
        confirmation = 0;
        while confirmation == 0
            [x, y] = ginput(1);
            [~, removeIdx] = min(abs(vly_locs-x));
            ANSWER = questdlg('Did you select the correct point?', ...
                'Confirmation Menu', ...
                'Yes', 'No', 'STOP', 'Yes');
            switch ANSWER
                case 'Yes'
                    x = round(x);
                    confirmation = 1;
                    % update valley info
                    vly_locs(removeIdx) = [];
                    valleys(removeIdx) = [];
                    true_valleys(removeIdx) = [];
                    removedValleys(i) = x;
                    set(valleyLine, 'XData', vly_locs,'YData',true_valleys);
%                     refreshdata(valleyLine,'base');
                    
                case 'No'
                    % Repeat removal selection
                    confirmation = 0;
                case 'STOP'
                    lastCall = questdlg('Are you sure you want to exit the program?', 'Quit','Yes','No','Yes');
                    switch lastCall
                        case 'Yes'
                            % exit the function and return to the main
                            % program
                            return
                        case 'No'
                            continue
                    end
                    
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% End of Valley Removal %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% User selected points will be added to valleys %%%
valley_additions = str2double(inputdlg('How many valleys would you like to add?'));
while isnan(valley_additions)
    valley_additions = str2double(inputdlg('Please input a valid integer to add valleys.'));
end
addedValleys = zeros(valley_additions,1);
if valley_additions > 0 && isnumeric(valley_additions)   
    for i = 1:valley_additions
        confirmation = 0;
        while confirmation == 0
            [x, y] = ginput(1);
            ANSWER = questdlg('Did you select the correct point?', ...
                'Confirmation Menu', ...
                'Yes', 'No', 'STOP', 'Yes');
            switch ANSWER
                case 'Yes'
                    x = round(x);
                    confirmation = 1;
                    % update valley info
                    if ~isempty(x) || ~isempty(y)
                        vly_locs = [vly_locs; x];
                        valleys = [valleys;y];
                        true_valleys = [true_valleys,y];
                        addedValleys(i) = x;
                        figure(h)
                        hold on
                        plot(axh, x,y, 'g*')
                    end
                case 'No'
                    % Repeat removal selection
                    confirmation = 0;
                case 'STOP'
                    lastCall = questdlg('Are you sure you want to exit the program?', 'Quit','Yes','No','Yes');
                    switch lastCall
                        case 'Yes'
                            continue
                        case 'No'
                            break
                    end
                    
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% End of Valley Addition %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:length(pks) % iterate thru peaks
    current_peak = locs(i);
    % % Leftward valley location will always be smaller than current peak
    test_valley_locs = current_peak - vly_locs;
    % % Check if there are locations surrounding the peak
    if (any(test_valley_locs < 0) && any(test_valley_locs > 0)) == 1
        left_valley_locs = test_valley_locs;
        left_valley_locs(left_valley_locs < 0) = Inf;
        [~, leftIdx] = min(left_valley_locs); % min value will be closest to peak
        if ~isempty(leftIdx) == 1
            leftLoc = vly_locs(leftIdx);
        end
        % % isolate values <0 because those are to the right of the peak
        right_valley_locs = test_valley_locs;
        right_valley_locs(right_valley_locs > 0) = Inf;
        [~,rightIdx] = min(abs(right_valley_locs)); % min val will be 'smallest' negative
        if ~isempty(rightIdx) == 1
            rightLoc = vly_locs(rightIdx);
        end
        % first check if right and left locs exist to avoid errors
        if (exist('leftLoc', 'var') == 1 && exist('rightLoc', 'var')) == 1
            % % % Integrate between right and left valley locs % % %
            peakRAW = raw_signal(leftLoc:rightLoc);
            raw_auc(i) = measureAUC(peakRAW);
            
            peakFILT = filtered_signal(leftLoc:rightLoc);
            filtered_auc(i) = measureAUC(peakFILT);
            
            peakDurations(i) = rightLoc - leftLoc;
        else
            continue
        end
        % clear valley locations for next loop iteration
        clear rightLoc;
        clear leftLoc;
    else
        continue
    end
end


%%% User selected peaks will be removed %%%
num_removals = str2double(inputdlg('How many peaks would you like to remove?'));
while isnan(num_removals)
    num_removals = str2double(inputdlg('Please input a valid integer to remove peaks.'));
end
removed = zeros(num_removals, 1);
if num_removals > 0 && isnumeric(num_removals)  
    for i = 1:num_removals
        confirmation = 0;
        % instantiate condition for misclick
        while confirmation == 0
            [x, y] = ginput(1);
            [val, removeIdx] = min(abs(locs-x));
            ANSWER = questdlg('Did you select the correct point?', ...
                'Confirmation Menu', ...
                'Yes', 'No', 'STOP', 'Yes');
            switch ANSWER
                case 'Yes'
                    confirmation = 1;
                    % update peak info
                    raw_auc(removeIdx) = NaN;
                    filtered_auc(removeIdx) = NaN;
                    peakDurations(removeIdx) = NaN;
                    locs(removeIdx) = NaN;
                    pks(removeIdx) = NaN;
                    peakWidth(removeIdx) = NaN;
                    removed(i) = x;
                case 'No'
                    % Repeat removal selection
                    confirmation = 0;
                case 'STOP'
                    lastCall = questdlg('Are you sure you want to exit the program?', 'Quit','Yes','No','Yes');
                    switch lastCall
                        case 'Yes'
                            continue
                        case 'No'
                            break
                    end
                    
            end
        end
    end
end
% calculate IEI
interEventIntervals = diff(locs(~isnan(locs)));
end