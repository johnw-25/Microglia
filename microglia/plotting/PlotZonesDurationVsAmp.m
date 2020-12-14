function PlotZonesDurationVsAmp(TimeZones, Pathology, limits)
% % Zones - Duration vs. Amplitude % %
figure()
plot(TimeZones.Zone1_PeakMags, TimeZones.Zone1_PeakDurs,'kd','MarkerFaceColor','g');
hold on
plot(TimeZones.Zone2_PeakMags, TimeZones.Zone2_PeakDurs,'kd','MarkerFaceColor','b');
hold on
plot(TimeZones.Zone3_PeakMags, TimeZones.Zone3_PeakDurs,'kd','MarkerFaceColor','r');
xlabel('Peak Amplitude (dF/F)'); ylabel('Peak Duration (s)');
title(strcat(Pathology,': Zones 1-3 Peak Duration vs Amplitude'));
legend('Zone 1','Zone 2','Zone 3');
ylim(limits);
end