function SimpleBars(data, yMax, graphColor)
% Column oriented data input

numGroups = size(data,2);
groupSize = ones(size(data,1),1);
figure()
hold on
for i = 1:numGroups
    tempData = data(:,i);
    tempMean = mean(tempData);
    tempSem = sem(tempData);
    
    xBar = [i-0.3, i+0.3]; yBar = [tempMean, tempMean];
    xLine = [i, i]; yLine = [tempMean-tempSem, tempMean+tempSem];
    xPosSem = [i-0.15, i+0.15]; yPosSem = [tempMean+tempSem, tempMean+tempSem];
    xNegSem = [i-0.15, i+0.15]; yNegSem = [tempMean-tempSem, tempMean-tempSem];
    
    %tempGroup = groupSize.*i;
    
    yl = ylim;
    
    hold on
    yTickMax = ChangeTextColor(num2str(yMax), [0, 0, 0]);
    yTickMin = ChangeTextColor(num2str(0), [0, 0, 0]);
    yTicks = {yTickMin, yTickMax};
    ylim([0 yMax])
    yticks([0,yMax]);
%     xlim(xlimit);
%     xticks([xlimit(1) 0 xlimit(2)]); 
    
    line(xLine,yLine,'Color',graphColor);
    line(xBar,yBar,'Color',graphColor,'LineWidth',1);
    line(xPosSem,yPosSem,'Color',graphColor,'LineWidth',1);
    line(xNegSem,yNegSem,'Color',graphColor,'LineWidth',1);
    tempHandle = scatter((i+(rand(size(tempData,1),1)-0.5)/6),tempData,8,'o',...
        'MarkerEdgeColor',graphColor,'MarkerFaceColor',graphColor);
    tempHandle.MarkerFaceAlpha = 0.5;
    tempHandle.MarkerEdgeAlpha = 0.25;
    
    
    set(gcf, 'Units','inches','position',[4 4 4 4]);
    set(gca,'FontSize',10,'XColor', [1, 1, 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
    set(gca, 'FontName','Arial');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [4 4 4 4]);
    set(gca,'YTickLabel', yTicks);
end