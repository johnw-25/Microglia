function ScatterRMS(allRMS, Pathology, faceColor)
A_RMS = allRMS(2:end,4);
A_RMS(cellfun(@(A_RMS) any(isnan(A_RMS)),A_RMS)) = [];
A_RMS(cellfun(@(A_RMS) any(isempty(A_RMS)),A_RMS)) = [];
A_RMS = cell2mat(A_RMS);
B_RMS = allRMS(2:end,5);
B_RMS(cellfun(@(B_RMS) any(isnan(B_RMS)),B_RMS)) = [];
B_RMS(cellfun(@(B_RMS) any(isempty(B_RMS)),B_RMS)) = [];
B_RMS = cell2mat(B_RMS);
C_RMS = allRMS(2:end,6);
C_RMS(cellfun(@(C_RMS) any(isnan(C_RMS)),C_RMS)) = [];
C_RMS(cellfun(@(C_RMS) any(isempty(C_RMS)),C_RMS)) = [];
C_RMS = cell2mat(C_RMS);

MEAN_A_RMS = mean(A_RMS);
MEAN_B_RMS = mean(B_RMS);
MEAN_C_RMS = mean(C_RMS);

[~, p12] = ttest2(A_RMS, B_RMS);
[~, p23] = ttest2(B_RMS, C_RMS);
[~, p13] = ttest2(A_RMS, C_RMS);

y = [MEAN_A_RMS, MEAN_B_RMS, MEAN_C_RMS];

A_SEM = std(A_RMS)/sqrt(length(A_RMS));
B_SEM = std(A_RMS)/sqrt(length(B_RMS));
C_SEM = std(A_RMS)/sqrt(length(C_RMS));

allSEM = [A_SEM, B_SEM, C_SEM];

x = [0 1 2];

figure()
errorbar(x,y,allSEM,'k');
hold on
plot(x,y,'kd','MarkerFaceColor','k');
xlabel('Zones'); ylabel('Average RMS'); title(strcat(Pathology,': Average RMS Between Zones (1-200, 201-timehitR-1, timehitR-end)'));

boxGroups = [ones(size(A_RMS)); ones(size(B_RMS)) * 2; ones(size(C_RMS)) * 3];
allData = [A_RMS; B_RMS; C_RMS];
figure()
boxplot(allData,boxGroups, 'Widths', 0.3, 'Colors', faceColor);
title(strcat(Pathology,': RMS Boxplots of Zones')); xlabel('Zones'); ylabel('RMS');
sigstar({[1,2],[2,3],[1,3]},[p12,p23,p13]);
end