function LPF = noiseRemover(filterOrder, Fc, Fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Simple function to create a FIR low-pass filter to remove high%%%%%%%
%%%%% frequency noise in microglia calcium fluorescence data.       %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LPF = designfilt('lowpassfir', 'FilterOrder', filterOrder, 'CutoffFrequency', Fc, 'SampleRate', Fs);

%response = fvtool(LPF);
end