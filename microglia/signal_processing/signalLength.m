function signalArcLength = signalLength(signal, Fs)
%%% This function finds the straight line distance of a 1D signal using
%%% diff().
time = (1:length(signal))./Fs;
dt = diff(time);
dFsegments = diff(signal);
ind_lengths = sqrt(dFsegments.^2 + dt.^2);
signalArcLength = sum(ind_lengths);

end