function IEIBoxPlotScatter(PBS, TMEV2, TMEV5, TMEV15)
set(0,'defaultfigurecolor',[1 1 1])
boxMarkerColorTMEV = [1, 0, 0];
boxMarkerColorPBS = [0, 0, 0];
boxColors = [1, 0, 0; 1, 0, 0; 1, 0, 0; 0.1, 0.1, 0.1];
xDim = 2.5;
yDim = 2.5;
yMax = 0.03;
zoneGroups = [ones(length(PBS),1); ones(length(TMEV2),1)*2; ones(length(TMEV5),1)*3; ones(length(TMEV15),1)*4];
zoneData = [PBS; TMEV2; TMEV5; TMEV15];
figure()
boxplot(zoneData,zoneGroups, 'Widths',0.35, 'Colors', [0 0 0],'OutlierSize',3,'Symbol','');
h = findobj(gca,'Tag','Box');
upperWhiskers = findobj(gca,'Tag','Upper Whisker');
upperAdjVal = findobj(gca,'Tag','Upper Adjacent Value');
lowerWhiskers = findobj(gca,'Tag','Lower Whisker');
lowerAdjVal = findobj(gca,'Tag','Lower Adjacent Value');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),boxColors(j,:),'FaceAlpha',.5);
end
hold on
% set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
ylimit = ylim;
ylim([0 1600]);
xlimit = xlim;
xlim([-0.5, xlimit(2)]);
% ylim([-4, -1.5]);
% yticks(-4:0.5:-1.5);
% yTickL = {ChangeTextColor(num2str(-4), [0, 0, 0]),ChangeTextColor(num2str(-3.5), [0, 0, 0]),...
%     ChangeTextColor(num2str(-3.0), [0, 0, 0]),ChangeTextColor(num2str(-2.5), [0, 0, 0]),...
%     ChangeTextColor(num2str(-2.0), [0, 0, 0]),ChangeTextColor(num2str(-1.5), [0, 0, 0])};
% set(gca,'YTickLabel',yTickL');
set(gca, 'TickLabelInterpreter', 'tex');
set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 xDim yDim]);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-','LineWidth',2)
set(findobj(gca,'tag','Median'),'LineWidth',1.5,'Color',[0,0,0])
hold on
z1Pbs = scatter((ones(size(PBS)).*0.5).*(1+(rand(size(PBS))-0.5)/1),PBS,24,'o','MarkerEdgeColor',boxMarkerColorPBS,'MarkerFaceColor',boxMarkerColorPBS);
z1Tmev2 = scatter((ones(size(TMEV2)).*1.5).*(1+(rand(size(TMEV2))-0.5)/8),TMEV2,24,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
z1Tmev5 = scatter((ones(size(TMEV5)).*2.5).*(1+(rand(size(TMEV5))-0.5)/8),TMEV5,24,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
z1Tmev15 = scatter((ones(size(TMEV15)).*3.5).*(1+(rand(size(TMEV15))-0.5)/8),TMEV15,24,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
z1Pbs.MarkerFaceAlpha = 0.5;
z1Pbs.MarkerEdgeAlpha = 0;
z1Tmev2.MarkerFaceAlpha = 0.5;
z1Tmev2.MarkerEdgeAlpha = 0;
z1Tmev5.MarkerFaceAlpha = 0.5;
z1Tmev5.MarkerEdgeAlpha = 0;
z1Tmev15.MarkerFaceAlpha = 0.5;
z1Tmev15.MarkerEdgeAlpha = 0;
hold off