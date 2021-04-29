function Phase01HeatMap(X, Y, MATRIX, yMax,lowestValue, highestValue,titleName, myColorMap)
%%%% Custom function to display matrix data as heatmap with our 'unique'
%%%% figure specifications.
%%%% INPUTS:
%%%%% X: Defines horizontal axis limits and imaging dimensions
%%%%% Y: Defines vertical axis limits and image dimensions
%%%%% MATRIX: 2D numerical matrix
%%%%% yMax: Determines size of y-axis (might be == Y)
%%%%% highest/lowestValue: determines color axis scale limits
%%%%% titleName: figure title
%%%%% myColorMap: colormap to apply to MATRIX
set(0,'defaultfigurecolor',[1 1 1])
f = figure();
h1 = imagesc([X(1) X(2)], [1 Y], MATRIX);
axisLines = gca;
set(axisLines,'Parent',f,'YDir','normal', ...
    'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255], ...
    'FontSize',10,'FontName','Arial','Color','none');
xlim([X(1) X(2)]);
xl = xlim(); % Get existing range.
yl = ylim();
yTickMax = ChangeTextColor(num2str(round(yl(2))), [0, 0, 0]);
roiTicks = {yTickMax};
ylim([0 yMax]);
xTickLocations = [round(xl(1)) 60 round(xl(2))];
yTickLocations = [round(Y)];
minTimeTick = ChangeTextColor(num2str(0), [0, 0, 0]);
time0Tick = ChangeTextColor(num2str(1), [0, 0, 0]);
maxTimeTick = ChangeTextColor(num2str(1.5), [0, 0, 0]);
timeTicks = {minTimeTick, time0Tick, maxTimeTick};
set(axisLines,'XTick', xTickLocations, ...
    'YTick',yTickLocations);
set(axisLines, 'XTickLabels',timeTicks,'YTickLabels',roiTicks);
% Make new labels for the new tick marks
title(titleName,'FontWeight','Normal'); %xlabel(labels,'Time (min)'); ylabel('ROI');
set(h1,'AlphaData',~isnan(MATRIX));
set(gcf, 'Units','inches','position',[4 4 2 4.5]);
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 2 4.5]);
caxis(gca,[lowestValue, highestValue]);
colormap(myColorMap);
colorbar
hold on
line(axisLines,[60, 60], [0 90], 'Color', [0.8, 0.8, 0.8],'LineStyle','-','LineWidth',1);
box off
print(f,strcat(titleName,'_DF.png'), '-r900','-dpng');
end