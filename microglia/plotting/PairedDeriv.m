function [ZONES, PVALUES,videos] = PairedDeriv(DATA, Pathology, xDim, yDim, faceColor, yMax)
set(0,'defaultfigurecolor',[1 1 1])
sLengthA = zeros(size(DATA,1),1);
sLengthB = zeros(size(DATA,1),1);
sLengthC = zeros(size(DATA,1),1);
lineTransparency = 0.15;
increaseColor = [71, 91, 249]./256; % Orange: [254/256, 127/256, 0/256]
decreaseColor = [71, 91, 249]./256;
% LPF = noiseRemover(100, 0.012, 0.97);
LPF = noiseRemover(100, 0.485, 0.97);
videos = [];
for ii = 1:size(DATA,1)
    currentVideo = DATA{ii,1};
    currentVideo = currentVideo(2:5);
    currentVideo = str2double(currentVideo);
    videos = [videos; currentVideo];
    timehitR = double(DATA{ii,6});
    timehitR_marker = round(timehitR);
    currentTrace = DATA{ii,3};
    F0 = mean(currentTrace(1:50));
    currentTrace = (currentTrace-F0)./F0;
    currentTrace = filtfilt(LPF, double(currentTrace));
    if timehitR_marker > 200 && timehitR_marker < length(currentTrace)
        sLengthA(ii) = sum(abs(diff(currentTrace(1:200))))/length(currentTrace(1:200));
        sLengthB(ii) = sum(abs(diff(currentTrace(201:timehitR_marker-1))))/length(currentTrace((201:timehitR_marker-1)));
        sLengthC(ii) = sum(abs(diff(currentTrace(timehitR_marker:end))))/length(currentTrace(timehitR_marker:end));
    else
        sLengthA(ii) = NaN;
        sLengthB(ii) = NaN;
        sLengthC(ii) = NaN;
    end
end

% remove NaN
sLengthA = log10(sLengthA(~isnan(sLengthA)));
sLengthB = log10(sLengthB(~isnan(sLengthB)));
sLengthC = log10(sLengthC(~isnan(sLengthC)));

q95A = quantile(sLengthA,0.95);
q75A = quantile(sLengthA,0.75);
q25A = quantile(sLengthA,0.25);

q95B = quantile(sLengthB,0.95);
q75B = quantile(sLengthB,0.75);
q25B = quantile(sLengthB,0.25);

q95C = quantile(sLengthC,0.95);
q75C = quantile(sLengthC,0.75);
q25C = quantile(sLengthC,0.25);

% w95 = (q95 - q75) / (q75 - q25)
w95A = (q95A - q75A)/(q75A - q25A);
w95B = (q95B - q75B)/(q75B - q25B);
w95C = (q95C - q75C)/(q75C - q25C);

traces = DATA(:,1);
traces = traces(~isnan(sLengthA));
mice = DATA(:,7);
mice = mice(~isnan(sLengthA));
% sLengthA = sLengthA(~isnan(sLengthA));
% sLengthB = sLengthB(~isnan(sLengthB));
% sLengthC = sLengthC(~isnan(sLengthC));

allZones = [sLengthA,sLengthB,sLengthC];
ZONES = cell(length(sLengthA),5);
ZONES = [num2cell(sLengthA),num2cell(sLengthB),num2cell(sLengthC),traces, mice];
PVALUES.p12 = signrank(sLengthA, sLengthB);
PVALUES.p23 = signrank(sLengthB, sLengthC);
PVALUES.p13 = signrank(sLengthA, sLengthC);

meanLengthA = mean(sLengthA);
meanLengthB = mean(sLengthB);
meanLengthC = mean(sLengthC);

barData = [max(sLengthA), max(sLengthB), max(sLengthC)];
semA = std(sLengthA)/sqrt(length(sLengthA));
semB = std(sLengthB)/sqrt(length(sLengthB));
semC = std(sLengthC)/sqrt(length(sLengthC));

x = [0 1 2];
y = [meanLengthA, meanLengthB, meanLengthC];
semLength = [semA, semB, semC];

% figure()
% errorbar(x,y,semLength,'k');
% hold on
% plot(x,y,'kd','MarkerFaceColor','k');
% xlabel('Zones'); ylabel('Average sum of deriv'); title(strcat(Pathology,': Average sum of deriv Between Zones (1-200, 201-timehitR-1, timehitR-end)'));
groups = [ones(length(sLengthA),1); ones(length(sLengthB),1)*2; ones(length(sLengthC),1)*3];
boxData = [sLengthA, sLengthB, sLengthC];
allData = [sLengthA; sLengthB; sLengthC];
figure()
boxplot(allData,groups, 'Widths',0.35, 'Colors', [0 0 0],'OutlierSize',3,'Symbol','o','whisker',w95A);
h = findobj(gca,'Tag','Box');
upperWhiskers = findobj(gca,'Tag','Upper Whisker');
upperAdjVal = findobj(gca,'Tag','Upper Adjacent Value');
lowerWhiskers = findobj(gca,'Tag','Lower Whisker');
lowerAdjVal = findobj(gca,'Tag','Lower Adjacent Value');
outliers = findobj(gca,'Tag','Outliers');
% hold on
% plot(ones(size(sLengthA)).*(1+(rand(size(sLengthA))-0.5)/10),sLengthA,'o','MarkerEdgeColor',[0.2, 0.2, 0.2],'MarkerFaceColor',[0.2, 0.2, 0.2],'MarkerSize',2)
% plot((ones(size(sLengthB)).*2).*(1+(rand(size(sLengthB))-0.5)/10),sLengthB,'o','MarkerEdgeColor',[0.2, 0.2, 0.2],'MarkerFaceColor',[0.2, 0.2, 0.2],'MarkerSize',2)
% plot((ones(size(sLengthC)).*3).*(1+(rand(size(sLengthC))-0.5)/10),sLengthC,'o','MarkerEdgeColor',[0.2, 0.2, 0.2],'MarkerFaceColor',[0.2, 0.2, 0.2],'MarkerSize',2)
for j=1:length(h)
    groupData = boxData(:,j);
    [data_highpercent(j)] = quantile(groupData,0.95);
    [data_lowpercent(j)] = quantile(groupData,1.00-0.95);
    currentOutliers = get(outliers(length(h)-j+1),'ydata');
    newOutliers = currentOutliers(outliers > data_highpercent(j));
    patch(get(h(j),'XData'),get(h(j),'YData'),faceColor,'FaceAlpha',.5);
    upperLims = get(upperWhiskers(length(h)-j+1),'ydata');
    lowerLims = get(lowerWhiskers(length(h)-j+1),'ydata');
    set(upperWhiskers(length(h)-j+1),'ydata',[upperLims(1), data_highpercent(j)]);
    set(upperAdjVal(length(h)-j+1),'ydata',[data_highpercent(j), data_highpercent(j)]);
    set(lowerWhiskers(length(h)-j+1),'ydata',[data_lowpercent(j), lowerLims(2)]);
    set(lowerAdjVal(length(h)-j+1),'ydata',[data_lowpercent(j), data_lowpercent(j)]);
    set(outliers(length(h)-j+1),'ydata',newOutliers);
end
set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 xDim yDim]);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-','LineWidth',2)
set(findobj(gca,'tag','Median'),'LineWidth',1.5,'Color',[0,0,0])
ylimit = ylim;
% ylim([-4, -1.5]);
% yticks(-4:0.5:-1.5);
% yTickL = {ChangeTextColor(num2str(-4), [0, 0, 0]),ChangeTextColor(num2str(-3.5), [0, 0, 0]),...
%     ChangeTextColor(num2str(-3.0), [0, 0, 0]),ChangeTextColor(num2str(-2.5), [0, 0, 0]),...
%     ChangeTextColor(num2str(-2.0), [0, 0, 0]),ChangeTextColor(num2str(-1.5), [0, 0, 0])};
% set(gca,'YTickLabel',yTickL');
set(gca, 'TickLabelInterpreter', 'tex');
% title(strcat(Pathology,': Path Length'));
print(gcf,strcat(Pathology,'_PathLengthBoxPlots.png'), '-r900','-dpng');
figure()
hold on
for k = 1:length(sLengthA)
    if sLengthA(k) > sLengthB(k)
        plot1 = plot([1,2],[sLengthA(k),sLengthB(k)],'o','MarkerEdgeColor',decreaseColor,'MarkerFaceColor',decreaseColor,'MarkerSize',3);
        plot1.Color(4) = 0.9;
        plot2 = plot([1,2],[sLengthA(k),sLengthB(k)],'Color',decreaseColor,'LineWidth',2);
        plot2.Color(4) = lineTransparency;
    else
        plot3 = plot([1,2],[sLengthA(k),sLengthB(k)],'Color', increaseColor,'LineWidth',2);
        plot3.Color(4) = lineTransparency;
        plot4 = plot([1,2],[sLengthA(k),sLengthB(k)],'o','MarkerEdgeColor',increaseColor,'MarkerFaceColor',increaseColor,'MarkerSize',3);
        plot4.Color(4) = 0.9;
    end
    
    if sLengthB(k) > sLengthC(k)
        plot5 = plot([2,3],[sLengthB(k),sLengthC(k)],'o','MarkerEdgeColor',decreaseColor,'MarkerFaceColor',decreaseColor,'MarkerSize',3);
        plot5.Color(4) = 0.9;
        plot6 = plot([2,3],[sLengthB(k),sLengthC(k)],'Color',decreaseColor,'LineWidth',2);
        plot6.Color(4) = lineTransparency;
    else
        plot7 = plot([2,3],[sLengthB(k),sLengthC(k)],'o','MarkerEdgeColor',increaseColor,'MarkerFaceColor',increaseColor,'MarkerSize',3);
        plot7.Color(4) = 0.9;
        plot8 = plot([2,3],[sLengthB(k),sLengthC(k)],'Color',increaseColor,'LineWidth',2);
        plot8.Color(4) = lineTransparency;
    end
end
% plot9 = plot([1,2], [mean(sLengthA), mean(sLengthB)],'s','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',9);
% plot10 = plot([2,3], [mean(sLengthB), mean(sLengthC)],'s','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',9);
plot11 = plot([1,2], [mean(sLengthA), mean(sLengthB)],'Color',[0 0 0],'LineWidth',3);
plot12 = plot([2,3], [mean(sLengthB), mean(sLengthC)],'Color',[0 0 0],'LineWidth',3);
ylimit = ylim;
% ylim([-4, -1.5]);
% yticks(-4:0.5:-1.5);
% yTickL = {ChangeTextColor(num2str(-4), [0, 0, 0]),ChangeTextColor(num2str(-3.5), [0, 0, 0]),...
%     ChangeTextColor(num2str(-3.0), [0, 0, 0]),ChangeTextColor(num2str(-2.5), [0, 0, 0]),...
%     ChangeTextColor(num2str(-2.0), [0, 0, 0]),ChangeTextColor(num2str(-1.5), [0, 0, 0])};
% set(gca,'YTickLabel',yTickL');
set(gca, 'TickLabelInterpreter', 'tex');
set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 xDim yDim]);
box off
print(gcf,strcat(Pathology,'_PathLengthScatter.png'), '-r900','-dpng');
% title(strcat(Pathology,': Path Length'),'FontWeight','Normal','FontSize',7);
% sigstar({[1,2], [2,3], [1,3]}, [p12, p23, p13]);

figure()
aucBar = bar(barData,'FaceColor', 'flat','BarWidth',0.4);
ylimit = ylim;
ylim([0 yMax]);
yticks(ylimit(2));
yTickL = ChangeTextColor(num2str(ylimit(2)), [0, 0, 0]);
aucBar.EdgeColor = 'none';
set(aucBar(1), 'FaceColor', faceColor)
set(gca, 'XTickLabel', {})
box off
title(Pathology,'FontWeight','Normal','FontSize',7);
nbars = size(barData,2);
x = aucBar.XEndPoints;
hold on
errorbar(x, barData, semLength, 'k','linestyle','none','CapSize',3,'LineWidth',1,'Color',[0.4, 0.4, 0.4]);
% sigstar({[1,2],[2,3],[1,3]},[p12,p23,p13]);
set(gca,'YTickLabel',yTickL);
set(gca, 'XTick',[])
set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 xDim yDim]);
% print(gcf,strcat(Pathology,'_PathLengthBars.png'), '-r900','-dpng');
end