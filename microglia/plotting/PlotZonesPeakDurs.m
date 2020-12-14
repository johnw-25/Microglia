function PlotZonesPeakDurs(TimeZones, Pathology, limits)

% % Zones - EVENT PeakDurs % %
figure()
plot(TimeZones.Zone1_Locs, TimeZones.Zone1_PeakDurs,'kd','MarkerFaceColor','g');
hold on
plot(TimeZones.Zone2_Locs, TimeZones.Zone2_PeakDurs,'kd','MarkerFaceColor','b');
hold on
plot(TimeZones.Zone3_Locs, TimeZones.Zone3_PeakDurs,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('Peak Durations (s)');
title(strcat(Pathology,': Zones 1-3 Peak Durations'));
legend('Zone 1','Zone 2','Zone 3');
grid on;
ylim(limits);
end