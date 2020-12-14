function PlotZonesPeakMags(TimeZones, Pathology, limits)
% % Zones - EVENT PeakMags % %
figure()
plot(TimeZones.Zone1_Locs, TimeZones.Zone1_PeakMags,'kd','MarkerFaceColor','g');
hold on
plot(TimeZones.Zone2_Locs, TimeZones.Zone2_PeakMags,'kd','MarkerFaceColor','b');
hold on
plot(TimeZones.Zone3_Locs, TimeZones.Zone3_PeakMags,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('Peak Amplitude (dF/f)');
title(strcat(Pathology,': Zones 1-3 Peak Amplitude'));
legend('Zone 1','Zone 2','Zone 3');
grid on;
ylim(limits);
end