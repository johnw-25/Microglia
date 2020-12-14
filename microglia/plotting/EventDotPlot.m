function EventDotPlot(DATA, Pathology, xlimit, yMax)
set(0,'defaultfigurecolor',[1 1 1])
if isstruct(DATA)
    fields = fieldnames(DATA);
    fig = figure();
    t = title(Pathology,'FontWeight','Normal');
    titlePos = get(t, 'position');
    titlePos(2) = titlePos(2) + 115;
    set(t, 'position',titlePos,'FontSize',8);
    for k = 1:numel(fields)
        temp_timehitR = DATA.(fields{k}).timehitR;
        R = round(temp_timehitR*0.97, 1);
        tempLocs = DATA.(fields{k}).PeakLocs;
        tempLocs = tempLocs(~isnan(tempLocs));
        frames = length(DATA.(fields{k}).trace);
        if temp_timehitR > frames
            shiftedLocs = (tempLocs - R)./0.97;
        else
            shiftedLocs = (tempLocs - R)./0.97;
        end
        shiftedLocs = shiftedLocs(~isnan(shiftedLocs));
        shiftedLocs = shiftedLocs./60;
        K = ones(size(shiftedLocs))*k;
        hold on
        plot(shiftedLocs, K, '.k', 'MarkerFaceColor','k','MarkerSize',3);
    end
    ylim([0 numel(fields)]);
    yl = ylim;
    yticks(yl(2));
    hold on
    line([0, 0], [0, yl(2)],'Color','blue','LineStyle','-','LineWidth',1);    
%     yTickMin = ChangeTextColor(num2str(round(yl(2)/4)), [0, 0, 0]);
    yTickMax = ChangeTextColor(num2str(round(yl(2))), [0, 0, 0]);
    roiTicks = {yTickMax};
    ylim([0 yMax])
    xlim(xlimit);
    xticks([xlimit(1) 0 xlimit(2)]); 
    set(gcf, 'Units','inches','position',[4 4 1.75 1.75]);
    set(gca,'FontSize',10,'XColor', [137/255 137/255 137/255],'YColor', [137/255 137/255 137/255],'FontName','Arial');
    set(gca, 'FontName','Arial');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [4 4 1.5 1.5]);
    timeLabel = sprintf('\\color[rgb]{%f, %f, %f}%s', [0 0 0], 'Time (min)');
    minTimeTick = sprintf('\\color[rgb]{%f, %f, %f}%s', [0 0 0], num2str(xlimit(1)));
    time0Tick = sprintf('\\color[rgb]{%f, %f, %f}%s', [0 0 0], num2str(0));
    maxTimeTick = sprintf('\\color[rgb]{%f, %f, %f}%s', [0 0 0], num2str(xlimit(2)));
    timeTicks = {minTimeTick, time0Tick, maxTimeTick};
    set(gca,'XTickLabel',timeTicks,'YTickLabel', roiTicks);
    xlabel(timeLabel);
    set(t,'FontSize',8);
%     saveas(fig, strcat(Pathology,'_DotPlot.png'), '-r600');
    print(fig, strcat(Pathology,'_DotPlot.png'), '-r600','-dpng');
else
    error('Function not defined for DATA input being non-structure.')
end
end