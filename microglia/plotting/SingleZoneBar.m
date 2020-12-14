function SingleZoneBar(barTitle, xDim, yDim, DATA, faceColor)
if ~isstruct(DATA)
    error('input must be struct')
end
% % Event AUC % %
zone1AUC = DATA.Zone1_AUC;
zone2AUC = DATA.Zone2_AUC;
zone3AUC = DATA.Zone3_AUC;
aucBarData = [max(zone1AUC), max(zone2AUC), max(zone3AUC)];
aucSEM = [sem(zone1AUC), sem(zone2AUC), sem(zone3AUC)];
[~,p12] = ttest2(zone1AUC,zone2AUC);
[~,p23] = ttest2(zone2AUC,zone3AUC);
[~,p13] = ttest2(zone1AUC,zone3AUC);
figure()
aucBar = bar(aucBarData,'FaceColor', 'flat','BarWidth',0.4);
set(aucBar(1), 'FaceColor', faceColor)
set(gca, 'XTickLabel', {})
box off
title(strcat(barTitle,' Signal Area'),'FontWeight','Normal');
nbars = size(aucBarData,2);
x = aucBar.XEndPoints;
hold on
errorbar(x, aucBarData, aucSEM, 'k','linestyle','none','CapSize',3,'LineWidth',1);
sigstar({[1,2],[2,3],[1,3]},[p12,p23,p13]);
ylimit = ylim;
yticks(ylimit(2));
yTickL = ChangeTextColor(num2str(ylimit(2)), [0, 0, 0]);
set(gca,'YTickLabel',yTickL);
set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 xDim yDim]);
end

