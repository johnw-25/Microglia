function [signals,avgFilter,stdFilter] = peak_detect(y, lag, threshold, influence)
% Robust peak detection algorithm via filtering and thresholding signal
% Inputs are -
% y: input signal vector
% lag: window size for moving average and moving std filters
% threshold: standard devs. above moving avg
% influence: 0 to 1 scalar for weighting importance of current and past
% values
% # Initialize variables
% set signals to vector 0,...,0 of length of y;   # Initialize signal results
signals = zeros(size(y));
% set filteredY to y(1),...,y(lag)                # Initialize filtered series
filteredY = y(1:lag+1);
%  Initialize average filter
%  Initialize std. filter
avgFilter(lag+1,1) = mean(y(1:lag+1));
stdFilter(lag+1,1) = std(y(1:lag+1));
% Loop over all datapoints y(lag+2),...,y(t)
for i=lag+2:length(y)
    % If new value is a specified number of deviations away
    if abs(y(i)-avgFilter(i-1)) > threshold*stdFilter(i-1)
        if y(i) > avgFilter(i-1)
            % Positive signal
            signals(i) = 1;
        else
            % Negative signal
            signals(i) = 0;
        end
        % Make influence lower
        filteredY(i) = influence*y(i)+(1-influence)*filteredY(i-1);
    else
        % No signal
        signals(i) = 0;
        filteredY(i) = y(i);
    end
    % Adjust the filters
    avgFilter(i) = mean(filteredY(i-lag:i));
    stdFilter(i) = std(filteredY(i-lag:i));
end
% Done, now return results
end