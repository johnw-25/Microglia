function [shifted_trace,loc_shift] = traceShift(RAW_trace, k, N, Z)
% This function takes four inputs:
% % RAW_trace: a signal represented as a row vector
% % k: single integer that shifts RAW_trace leftward
% % N: positive integer determines number of zeros are used to pad vector.
% % Z: initial vector length
% Two output:
% % shifted_trace: trace shifted by -k units padded with excess zeros
% The function is dependent on circshift() as the main operator.
% % loc_shift:
% The linear shift between RAW_trace and shifted_trace.

% Create before_time0 vector of uniform length
% Create uniform length first shifted vector
if length(RAW_trace) < Z
    RAW_trace = [RAW_trace, zeros(1,Z-length(RAW_trace))];
end
before_time0 = zeros(1,Z);
first_shifted_trace = zeros(1,Z);
% Find positions to insert raw trace into padded vector.
% insert_idx = length(padded_trace)-length(RAW_trace)+1:length(padded_trace);
% % Insert raw trace into padded vector.
% padded_trace(1, insert_idx) = RAW_trace;
% Compute leftward shift by factor of k.
first_shifted_trace(1,1:length(RAW_trace)) = circshift(RAW_trace, -k+1);
% Update k if k is a time shift longer than the duration of the trace.
if k > length(first_shifted_trace)
    k = k - length(first_shifted_trace);
end
temp_before = first_shifted_trace(end-k:end);
before_time0(1,end-length(temp_before)+1:end) = temp_before;
first_shifted_trace(end-k:end) = 0;
shifted_trace = [before_time0, first_shifted_trace];
% Use maxima to find lag between signals - maxima won't change.
loc_shift = find(max(shifted_trace)) - find(max(RAW_trace));
end
