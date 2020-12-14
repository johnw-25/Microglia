function PlotZonesSqrtAUCVsAmp(TimeZones, Pathology, limits)
% % Zones - sqrt(AUC) vs. Amplitude % %
figure()
plot(TimeZones.Zone1_PeakMags, sqrt(TimeZones.Zone1_AUC),'kd','MarkerFaceColor','g');
hold on
plot(TimeZones.Zone2_PeakMags, sqrt(TimeZones.Zone2_AUC),'kd','MarkerFaceColor','b');
hold on
plot(TimeZones.Zone3_PeakMags, sqrt(TimeZones.Zone3_AUC),'kd','MarkerFaceColor','r');
xlabel('Peak Amplitude (dF/F)'); ylabel('sqrt(AUC) (sqrt(dF/F))');
title(strcat(Pathology,': Zones 1-3 Event sqrt(AUC) vs Amplitude'));
legend('Zone 1','Zone 2','Zone 3');
ylim(limits);
end