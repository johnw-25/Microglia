function PVALUES = AllDerivZones(PBS_DATA, TMEV2_DATA, TMEV5_DATA, TMEV15_DATA,Pathology, yMax, xDim, yDim)
set(0,'defaultfigurecolor',[1 1 1])
boxMarkerColorTMEV = [1, 0, 0];
boxMarkerColorPBS = [0, 0, 0];
boxColors = [1, 0, 0; 1, 0, 0; 1, 0, 0; 0.1, 0.1, 0.1];
pbsLengthA = zeros(size(PBS_DATA,1),1);
pbsLengthB = zeros(size(PBS_DATA,1),1);
pbsLengthC = zeros(size(PBS_DATA,1),1);
LPF = noiseRemover(100, 0.012, 0.97);
for ii = 1:size(PBS_DATA,1)
    
    timehitR = double(PBS_DATA{ii,6});
    timehitR_marker = round(timehitR);
    currentTrace = PBS_DATA{ii,3};
    F0 = mean(currentTrace(1:50));
    currentTrace = (currentTrace-F0)./F0;
    currentTrace = filtfilt(LPF, double(currentTrace));
    if timehitR_marker > 200 && timehitR_marker < length(currentTrace)
        pbsLengthA(ii) = sum(abs(diff(currentTrace(1:200))))/length(currentTrace(1:200));
        pbsLengthB(ii) = sum(abs(diff(currentTrace(201:timehitR_marker-1))))/length(currentTrace((201:timehitR_marker-1)));
        pbsLengthC(ii) = sum(abs(diff(currentTrace(timehitR_marker:end))))/length(currentTrace(timehitR_marker:end));
    else
        pbsLengthA(ii) = NaN;
        pbsLengthB(ii) = NaN;
        pbsLengthC(ii) = NaN;
    end
end

tmev2LengthA = zeros(size(TMEV2_DATA,1),1);
tmev2LengthB = zeros(size(TMEV2_DATA,1),1);
tmev2LengthC = zeros(size(TMEV2_DATA,1),1);
LPF = noiseRemover(100, 0.012, 0.97);
for ii = 1:size(TMEV2_DATA,1)
    
    timehitR = double(TMEV2_DATA{ii,6});
    timehitR_marker = round(timehitR);
    currentTrace = TMEV2_DATA{ii,3};
    F0 = mean(currentTrace(1:50));
    currentTrace = (currentTrace-F0)./F0;
    currentTrace = filtfilt(LPF, double(currentTrace));
    if timehitR_marker > 200 && timehitR_marker < length(currentTrace)
        tmev2LengthA(ii) = sum(abs(diff(currentTrace(1:200))))/length(currentTrace(1:200));
        tmev2LengthB(ii) = sum(abs(diff(currentTrace(201:timehitR_marker-1))))/length(currentTrace((201:timehitR_marker-1)));
        tmev2LengthC(ii) = sum(abs(diff(currentTrace(timehitR_marker:end))))/length(currentTrace(timehitR_marker:end));
    else
        tmev2LengthA(ii) = NaN;
        tmev2LengthB(ii) = NaN;
        tmev2LengthC(ii) = NaN;
    end
end

tmev5LengthA = zeros(size(TMEV5_DATA,1),1);
tmev5LengthB = zeros(size(TMEV5_DATA,1),1);
tmev5LengthC = zeros(size(TMEV5_DATA,1),1);
LPF = noiseRemover(100, 0.012, 0.97);
for ii = 1:size(TMEV5_DATA,1)
    
    timehitR = double(TMEV5_DATA{ii,6});
    timehitR_marker = round(timehitR);
    currentTrace = TMEV5_DATA{ii,3};
    F0 = mean(currentTrace(1:50));
    currentTrace = (currentTrace-F0)./F0;
    currentTrace = filtfilt(LPF, double(currentTrace));
    if timehitR_marker > 200 && timehitR_marker < length(currentTrace)
        tmev5LengthA(ii) = sum(abs(diff(currentTrace(1:200))))/length(currentTrace(1:200));
        tmev5LengthB(ii) = sum(abs(diff(currentTrace(201:timehitR_marker-1))))/length(currentTrace((201:timehitR_marker-1)));
        tmev5LengthC(ii) = sum(abs(diff(currentTrace(timehitR_marker:end))))/length(currentTrace(timehitR_marker:end));
    else
        tmev5LengthA(ii) = NaN;
        tmev5LengthB(ii) = NaN;
        tmev5LengthC(ii) = NaN;
    end
end

tmev15LengthA = zeros(size(TMEV15_DATA,1),1);
tmev15LengthB = zeros(size(TMEV15_DATA,1),1);
tmev15LengthC = zeros(size(TMEV15_DATA,1),1);
LPF = noiseRemover(100, 0.012, 0.97);
for ii = 1:size(TMEV15_DATA,1)
    
    timehitR = double(TMEV15_DATA{ii,6});
    timehitR_marker = round(timehitR);
    currentTrace = TMEV15_DATA{ii,3};
    F0 = mean(currentTrace(1:50));
    currentTrace = (currentTrace-F0)./F0;
    currentTrace = filtfilt(LPF, double(currentTrace));
    if timehitR_marker > 200 && timehitR_marker < length(currentTrace)
        tmev15LengthA(ii) = sum(abs(diff(currentTrace(1:200))))/length(currentTrace(1:200));
        tmev15LengthB(ii) = sum(abs(diff(currentTrace(201:timehitR_marker-1))))/length(currentTrace((201:timehitR_marker-1)));
        tmev15LengthC(ii) = sum(abs(diff(currentTrace(timehitR_marker:end))))/length(currentTrace(timehitR_marker:end));
    else
        tmev15LengthA(ii) = NaN;
        tmev15LengthB(ii) = NaN;
        tmev15LengthC(ii) = NaN;
    end
end

PVALUES.Zone1.Pbs_Tmev2 = ranksum(pbsLengthA, tmev2LengthA);
PVALUES.Zone1.Pbs_Tmev5 = ranksum(pbsLengthA, tmev5LengthA);
PVALUES.Zone1.Pbs_Tmev15 = ranksum(pbsLengthA, tmev15LengthA);
PVALUES.Zone1.Tmev2_Tmev5 = ranksum(tmev2LengthA, tmev5LengthA);
PVALUES.Zone1.Tmev2_Tmev15 = ranksum(tmev2LengthA, tmev15LengthA);
PVALUES.Zone1.Tmev5_Tmev15 = ranksum(tmev5LengthA, tmev15LengthA);

PVALUES.Zone2.Pbs_Tmev2 = ranksum(pbsLengthB, tmev2LengthB);
PVALUES.Zone2.Pbs_Tmev5 = ranksum(pbsLengthB, tmev5LengthB);
PVALUES.Zone2.Pbs_Tmev15 = ranksum(pbsLengthB, tmev15LengthB);
PVALUES.Zone2.Tmev2_Tmev5 = ranksum(tmev2LengthB, tmev5LengthB);
PVALUES.Zone2.Tmev2_Tmev15 = ranksum(tmev2LengthB, tmev15LengthB);
PVALUES.Zone2.Tmev5_Tmev15 = ranksum(tmev5LengthB, tmev15LengthB);

PVALUES.Zone3.Pbs_Tmev2 = ranksum(pbsLengthC, tmev2LengthC);
PVALUES.Zone3.Pbs_Tmev5 = ranksum(pbsLengthC, tmev5LengthC);
PVALUES.Zone3.Pbs_Tmev15 = ranksum(pbsLengthC, tmev15LengthC);
PVALUES.Zone3.Tmev2_Tmev5 = ranksum(tmev2LengthC, tmev5LengthC);
PVALUES.Zone3.Tmev2_Tmev15 = ranksum(tmev2LengthC, tmev15LengthC);
PVALUES.Zone3.Tmev5_Tmev15 = ranksum(tmev5LengthC, tmev15LengthC);

zone1Groups = [ones(length(pbsLengthA),1); ones(length(tmev2LengthA),1)*2; ones(length(tmev5LengthA),1)*3; ones(length(tmev15LengthA),1)*4];
zone1Data = log10([pbsLengthA; tmev2LengthA; tmev5LengthA; tmev15LengthA]);
figure()
boxplot(zone1Data,zone1Groups, 'Widths',0.35, 'Colors', [0 0 0],'OutlierSize',3,'Symbol','');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),boxColors(j,:),'FaceAlpha',.5);
end
hold on
% set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
ylimit = ylim;

xlimit = xlim;
xlim([-0.5, xlimit(2)]);
ylim([-4, -1.5]);
yticks(-4:0.5:-1.5);
yTickL = {ChangeTextColor(num2str(-4), [0, 0, 0]),ChangeTextColor(num2str(-3.5), [0, 0, 0]),...
    ChangeTextColor(num2str(-3.0), [0, 0, 0]),ChangeTextColor(num2str(-2.5), [0, 0, 0]),...
    ChangeTextColor(num2str(-2.0), [0, 0, 0]),ChangeTextColor(num2str(-1.5), [0, 0, 0])};
set(gca,'YTickLabel',yTickL');
set(gca, 'TickLabelInterpreter', 'tex');
set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 xDim yDim]);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-','LineWidth',2)
set(findobj(gca,'tag','Median'),'LineWidth',1.5,'Color',[0,0,0])
hold on
z1Pbs = scatter((ones(size(pbsLengthA)).*0.5).*(1+(rand(size(pbsLengthA))-0.5)/1),log10(pbsLengthA),8,'o','MarkerEdgeColor',boxMarkerColorPBS,'MarkerFaceColor',boxMarkerColorPBS);
z1Tmev2 = scatter((ones(size(tmev2LengthA)).*1.5).*(1+(rand(size(tmev2LengthA))-0.5)/8),log10(tmev2LengthA),8,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
z1Tmev5 = scatter((ones(size(tmev5LengthA)).*2.5).*(1+(rand(size(tmev5LengthA))-0.5)/8),log10(tmev5LengthA),8,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
z1Tmev15 = scatter((ones(size(tmev15LengthA)).*3.5).*(1+(rand(size(tmev15LengthA))-0.5)/8),log10(tmev15LengthA),8,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
z1Pbs.MarkerFaceAlpha = 0.5;
z1Pbs.MarkerEdgeAlpha = 0;
z1Tmev2.MarkerFaceAlpha = 0.5;
z1Tmev2.MarkerEdgeAlpha = 0;
z1Tmev5.MarkerFaceAlpha = 0.5;
z1Tmev5.MarkerEdgeAlpha = 0;
z1Tmev15.MarkerFaceAlpha = 0.5;
z1Tmev15.MarkerEdgeAlpha = 0;
hold off
% title(strcat(Pathology,': Zone 1 Path Length'),'FontWeight','normal','FontSize',6);
print(gcf,strcat(Pathology,'Zone1PathLengthBoxPlots.png'), '-r900','-dpng');

zone2Groups = [ones(length(pbsLengthB),1); ones(length(tmev2LengthB),1)*2; ones(length(tmev5LengthB),1)*3; ones(length(tmev15LengthB),1)*4];
zone2Data = log10([pbsLengthB; tmev2LengthB; tmev5LengthB; tmev15LengthB]);
figure()
boxplot(zone2Data,zone2Groups, 'Widths',0.35, 'Colors', [0 0 0],'OutlierSize',3,'Symbol','');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),boxColors(j,:),'FaceAlpha',.5);
end
hold on
z1Pbs = scatter((ones(size(pbsLengthB)).*0.5).*(1+(rand(size(pbsLengthB))-0.5)/1),log10(pbsLengthB),8,'o','MarkerEdgeColor',boxMarkerColorPBS,'MarkerFaceColor',boxMarkerColorPBS);
z1Tmev2 = scatter((ones(size(tmev2LengthB)).*1.5).*(1+(rand(size(tmev2LengthB))-0.5)/8),log10(tmev2LengthB),8,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
z1Tmev5 = scatter((ones(size(tmev5LengthB)).*2.5).*(1+(rand(size(tmev5LengthB))-0.5)/8),log10(tmev5LengthB),8,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
z1Tmev15 = scatter((ones(size(tmev15LengthB)).*3.5).*(1+(rand(size(tmev15LengthB))-0.5)/8),log10(tmev15LengthB),8,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
z1Pbs.MarkerFaceAlpha = 0.5;
z1Pbs.MarkerEdgeAlpha = 0.5;
z1Tmev2.MarkerFaceAlpha = 0.5;
z1Tmev2.MarkerEdgeAlpha = 0.5;
z1Tmev5.MarkerFaceAlpha = 0.5;
z1Tmev5.MarkerEdgeAlpha = 0.5;
z1Tmev15.MarkerFaceAlpha = 0.5;
z1Tmev15.MarkerEdgeAlpha = 0.5;
hold off
ylimit = ylim;
xlimit = xlim;
xlim([-0.5, xlimit(2)]);
ylim([-4, -1.5]);
yticks(-4:0.5:-1.5);
yTickL = {ChangeTextColor(num2str(-4), [0, 0, 0]),ChangeTextColor(num2str(-3.5), [0, 0, 0]),...
    ChangeTextColor(num2str(-3.0), [0, 0, 0]),ChangeTextColor(num2str(-2.5), [0, 0, 0]),...
    ChangeTextColor(num2str(-2.0), [0, 0, 0]),ChangeTextColor(num2str(-1.5), [0, 0, 0])};
set(gca,'YTickLabel',yTickL');
set(gca, 'TickLabelInterpreter', 'tex');
set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 xDim yDim]);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-','LineWidth',2)
set(findobj(gca,'tag','Median'),'LineWidth',1.5,'Color',[0,0,0])
% set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
% title(strcat(Pathology,': Zone 2 Path Length'),'FontWeight','normal','FontSize',6);
print(gcf,strcat(Pathology,'Zone2PathLengthBoxPlots.png'), '-r900','-dpng');

zone3Groups = [ones(length(pbsLengthC),1)*2; ones(length(tmev2LengthC),1)*4; ones(length(tmev5LengthC),1)*6; ones(length(tmev15LengthC),1)*8];
zone3Data = log10([pbsLengthC; tmev2LengthC; tmev5LengthC; tmev15LengthC]);
figure()
boxplot(zone3Data,zone3Groups, 'Widths',0.35, 'Colors', [0 0 0],'OutlierSize',3,'Symbol','');
h = findobj(gca,'Tag','Box');
ax = gca;
ax.YRuler.Exponent = 0;
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),boxColors(j,:),'FaceAlpha',.5);
end
hold on
z1Pbs = scatter((ones(size(pbsLengthC)).*0.5).*(1+(rand(size(pbsLengthC))-0.5)/1),log10(pbsLengthC),8,'o','MarkerEdgeColor',boxMarkerColorPBS,'MarkerFaceColor',boxMarkerColorPBS);
z1Tmev2 = scatter((ones(size(tmev2LengthC)).*1.5).*(1+(rand(size(tmev2LengthC))-0.5)/8),log10(tmev2LengthC),8,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
z1Tmev5 = scatter((ones(size(tmev5LengthC)).*2.5).*(1+(rand(size(tmev5LengthC))-0.5)/8),log10(tmev5LengthC),8,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
z1Tmev15 = scatter((ones(size(tmev15LengthC)).*3.5).*(1+(rand(size(tmev15LengthC))-0.5)/8),log10(tmev15LengthC),8,'o','MarkerEdgeColor',boxMarkerColorTMEV,'MarkerFaceColor',boxMarkerColorTMEV);
z1Pbs.MarkerFaceAlpha = 0.5;
z1Pbs.MarkerEdgeAlpha = 0.5;
z1Tmev2.MarkerFaceAlpha = 0.5;
z1Tmev2.MarkerEdgeAlpha = 0.5;
z1Tmev5.MarkerFaceAlpha = 0.5;
z1Tmev5.MarkerEdgeAlpha = 0.5;
z1Tmev15.MarkerFaceAlpha = 0.5;
z1Tmev15.MarkerEdgeAlpha = 0.5;
hold off
ylimit = ylim;
xlimit = xlim;
xlim([-0.5, xlimit(2)]);
ylim([-4, -1.5]);
yticks(-4:0.5:-1.5);
yTickL = {ChangeTextColor(num2str(-4), [0, 0, 0]),ChangeTextColor(num2str(-3.5), [0, 0, 0]),...
    ChangeTextColor(num2str(-3.0), [0, 0, 0]),ChangeTextColor(num2str(-2.5), [0, 0, 0]),...
    ChangeTextColor(num2str(-2.0), [0, 0, 0]),ChangeTextColor(num2str(-1.5), [0, 0, 0])};
set(gca,'YTickLabel',yTickL');
set(gca, 'TickLabelInterpreter', 'tex');
set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 xDim yDim]);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-','LineWidth',2)
set(findobj(gca,'tag','Median'),'LineWidth',1.5,'Color',[0,0,0])
% set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
% title(strcat(Pathology,': Zone 3 Path Length'),'FontWeight','normal','FontSize',6);
print(gcf,strcat(Pathology,'Zone3PathLengthBoxPlots.png'), '-r900','-dpng');