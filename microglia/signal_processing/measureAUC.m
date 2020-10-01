function AUC = measureAUC(cutSignal)
%%% Basic function that finds the area under the curve of an input signal.
posAUC = sum(cutSignal(cutSignal >=0));
negAUC = sum(cutSignal(cutSignal < 0));

AUC = posAUC + abs(negAUC);

end