function velocityCell = ExtractVelocityData(pathName,fileName)
test_data = readtable(fullfile(pathName,fileName),'Sheet',2);
id = test_data.TrackID;
points = test_data.Points;
meanV = test_data.Maxvumsec;
pathology = test_data.Groups;
dpi = test_data.DPI;
mice = test_data.Mouse;
meanDtoR = test_data.MeanD2Rum;
maxDtoR = test_data.MaxD2Rum;
lenMicrons = test_data.Len__um_;
meanY = test_data.Mean_y__um_;
maxY = test_data.Max_y__um_;
sdY = test_data.SD_y__um_;
meanS = test_data.Mean_D2S__um_;
maxS = test_data.Max_D2S__um_;
sdS = test_data.SD_D2S__um_;

velocityCell = [id, num2cell(points), num2cell(meanV), pathology, num2cell(dpi), num2cell(mice), num2cell(meanDtoR), ...
    num2cell(maxDtoR), num2cell(lenMicrons), num2cell(meanY), num2cell(maxY), num2cell(sdY), num2cell(meanS), num2cell(maxS), num2cell(sdS)];

pbsCell = velocityCell(strcmp(pathology,'PBS'),:);
tmevCell = velocityCell(strcmp(pathology,'TMEV'),:);

pbsVel = meanV(strcmp(pathology,'PBS'));
tmevVel = meanV(strcmp(pathology,'TMEV'));

pbsDist = meanDtoR(strcmp(pathology,'PBS'));
tmevDist = meanDtoR(strcmp(pathology,'TMEV'));

tmevDpi = dpi(strcmp(pathology,'TMEV'));
groups = [ones(size(pbsVel)); ones(size(tmevVel))*2];
figure()
boxplot([pbsVel;tmevVel], groups)

tmev2Vel = tmevVel(tmevDpi == 2);
tmev15Vel = tmevVel(tmevDpi >= 14);
tmev5Vel = tmevVel(tmevDpi > 2 & tmevDpi < 14); %&& tmevVel(tmevDpi < 14);

tmev2Dist = tmevDist(tmevDpi == 2);
tmev15Dist = tmevDist(tmevDpi >= 14);
tmev5Dist = tmevDist(tmevDpi > 2 & tmevDpi < 14); %&& tmevDist(tmevDpi < 14);

p12 = ranksum(pbsVel, tmev2Vel);
p13 = ranksum(pbsVel, tmev5Vel);
p14 = ranksum(pbsVel, tmev15Vel);
p23 = ranksum(tmev2Vel, tmev5Vel);
p24 = ranksum(tmev2Vel, tmev15Vel);
p34 = ranksum(tmev5Vel, tmev15Vel);

p12D = ranksum(pbsDist, tmev2Dist);
p13D = ranksum(pbsDist, tmev5Dist);
p14D = ranksum(pbsDist, tmev15Dist);
p23D = ranksum(tmev2Dist, tmev5Dist);
p24D = ranksum(tmev2Dist, tmev15Dist);
p34D = ranksum(tmev5Dist, tmev15Dist);

set(0,'defaultfigurecolor',[1 1 1])
boxMarkerColorTMEV = [1, 0, 0];
boxMarkerColorPBS = [0, 0, 0];
boxColors = [1, 0, 0; 1, 0, 0; 1, 0, 0; 0.1, 0.1, 0.1];

grps = [ones(size(pbsVel)); ones(size(tmev2Vel))*2; ones(size(tmev5Vel))*3; ones(size(tmev15Vel))*4];
allVelData = [pbsVel;tmev2Vel;tmev5Vel;tmev15Vel];
figure()
boxplot(allVelData,grps,'Widths',0.35,'OutlierSize',6,'Symbol','o', 'Colors', [0 0 0])
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),boxColors(j,:),'FaceAlpha',.5);
end
hold on
% Pbs = scatter(ones(size(pbsVel)).*(1+(rand(size(pbsVel))-0.5)/2),pbsVel,4,'o','MarkerEdgeColor',boxMarkerColorPBS,'MarkerFaceColor',boxMarkerColorPBS);
% Tmev2 = scatter((ones(size(tmev2Vel)).*2).*(1+(rand(size(tmev2Vel))-0.5)/10),tmev2Vel,4,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
% Tmev5 = scatter((ones(size(tmev5Vel)).*3).*(1+(rand(size(tmev5Vel))-0.5)/10),tmev5Vel,4,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
% Tmev15 = scatter((ones(size(tmev15Vel)).*4).*(1+(rand(size(tmev15Vel))-0.5)/10),tmev15Vel,4,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
% Pbs.MarkerFaceAlpha = 0.5;
% Pbs.MarkerEdgeAlpha = 0.5;
% Tmev2.MarkerFaceAlpha = 0.5;
% Tmev2.MarkerEdgeAlpha = 0.5;
% Tmev5.MarkerFaceAlpha = 0.5;
% Tmev5.MarkerEdgeAlpha = 0.5;
% Tmev15.MarkerFaceAlpha = 0.5;
% Tmev15.MarkerEdgeAlpha = 0.5;
hold off
set(gcf, 'Units','inches','position',[4 4 2.5 2.5]);
set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'Units','inches','position',[4 4 2.5 2.5]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 2.5 2.5]);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-','LineWidth',2)
ylimit = ylim;
yTickL = ChangeTextColor(num2str(round(ylimit(2),2)), [0, 0, 0]);
yticks(ylimit(2));
set(gca,'YTickLabel',yTickL);
set(gca, 'TickLabelInterpreter', 'tex');
% title(strcat(Pathology,': Path Length'));
print(gcf,'AvgVelocityBoxPlots.png', '-r900','-dpng');
% sigstar({[1,2],[1,3],[1,4]}, [p12,p13,p14])

allDistData = [pbsDist;tmev2Dist;tmev5Dist;tmev15Dist];
figure()
boxplot(allDistData,grps,'Widths',0.35,'OutlierSize',6,'Symbol','o', 'Colors', [0 0 0])
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),boxColors(j,:),'FaceAlpha',.5);
end
set(gcf, 'Units','inches','position',[4 4 2.5 2.5]);
set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'Units','inches','position',[4 4 2.5 2.5]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 2.5 2.5]);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-','LineWidth',2)
ylimit = ylim;
yTickL = ChangeTextColor(num2str(round(ylimit(2),2)), [0, 0, 0]);
yticks(ylimit(2));
set(gca,'YTickLabel',yTickL);
set(gca, 'TickLabelInterpreter', 'tex');
% title(strcat(Pathology,': Path Length'));
print(gcf,'AvgDistanceBoxPlots.png', '-r900','-dpng');
end