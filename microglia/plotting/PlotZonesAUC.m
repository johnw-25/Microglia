function PlotZonesAUC(TimeZones, Pathology, limits)

% % Zones - EVENT AUC % %
figure()
plot(TimeZones.Zone1_Locs, TimeZones.Zone1_AUC,'kd','MarkerFaceColor','g');
hold on
plot(TimeZones.Zone2_Locs, TimeZones.Zone2_AUC,'kd','MarkerFaceColor','b');
hold on
plot(TimeZones.Zone3_Locs, TimeZones.Zone3_AUC,'kd','MarkerFaceColor','r');
xlabel('Peak locations (frames)'); ylabel('AUC (dF/F)');
title(strcat(Pathology,': Zones 1-3 AUC'));
legend('Zone 1','Zone 2','Zone 3');
grid on;
ylim(limits);
end