function PairedZonesDeriv(barTitle, xDim, yDim, PBS, TMEV2, tmev5, tmev15)
set(0,'defaultfigurecolor',[1 1 1])
pbsLengthA = zeros(size(PBS,1),1);
pbsLengthB = zeros(size(PBS,1),1);
pbsLengthC = zeros(size(PBS,1),1);
LPF = noiseRemover(100, 0.012, 0.97);
for ii = 1:size(PBS,1)
    
    timehitR = double(PBS{ii,6});
    timehitR_marker = round(timehitR);
    currentTrace = PBS{ii,3};
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

% remove NaN
pbsLengthA = pbsLengthA(~isnan(pbsLengthA));
pbsLengthB = pbsLengthB(~isnan(pbsLengthB));
pbsLengthC = pbsLengthC(~isnan(pbsLengthC));

set(0,'defaultfigurecolor',[1 1 1])
tmev2LengthA = zeros(size(TMEV2,1),1);
tmev2LengthB = zeros(size(TMEV2,1),1);
tmev2LengthC = zeros(size(TMEV2,1),1);
LPF = noiseRemover(100, 0.012, 0.97);
for ii = 1:size(TMEV2,1)
    
    timehitR = double(TMEV2{ii,6});
    timehitR_marker = round(timehitR);
    currentTrace = TMEV2{ii,3};
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

% remove NaN
tmev2LengthA = tmev2LengthA(~isnan(tmev2LengthA));
tmev2LengthB = tmev2LengthB(~isnan(tmev2LengthB));
tmev2LengthC = tmev2LengthC(~isnan(tmev2LengthC));

set(0,'defaultfigurecolor',[1 1 1])
tmev5LengthA = zeros(size(tmev5,1),1);
tmev5LengthB = zeros(size(tmev5,1),1);
tmev5LengthC = zeros(size(tmev5,1),1);
LPF = noiseRemover(100, 0.012, 0.97);
for ii = 1:size(tmev5,1)
    
    timehitR = double(tmev5{ii,6});
    timehitR_marker = round(timehitR);
    currentTrace = tmev5{ii,3};
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

% remove NaN
tmev5LengthA = tmev5LengthA(~isnan(tmev5LengthA));
tmev5LengthB = tmev5LengthB(~isnan(tmev5LengthB));
tmev5LengthC = tmev5LengthC(~isnan(tmev5LengthC));

set(0,'defaultfigurecolor',[1 1 1])
tmev15LengthA = zeros(size(tmev15,1),1);
tmev15LengthB = zeros(size(tmev15,1),1);
tmev15LengthC = zeros(size(tmev15,1),1);
LPF = noiseRemover(100, 0.012, 0.97);
for ii = 1:size(tmev15,1)
    
    timehitR = double(tmev15{ii,6});
    timehitR_marker = round(timehitR);
    currentTrace = tmev15{ii,3};
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

% remove NaN
tmev15LengthA = tmev15LengthA(~isnan(tmev15LengthA));
tmev15LengthB = tmev15LengthB(~isnan(tmev15LengthB));
tmev15LengthC = tmev15LengthC(~isnan(tmev15LengthC));

zoneGroups = [ones(3,1), ones(3,1)*2, ones(3,1)*3];

zoneAuc = [max(pbsLengthA),max(pbsLengthB),max(pbsLengthC);
    max(tmev2LengthA),max(tmev2LengthB),max(tmev2LengthC);
    max(tmev5LengthA),max(tmev5LengthB),max(tmev5LengthC);
    max(tmev15LengthA),max(tmev15LengthB),max(tmev15LengthC)];

allSEM = [sem(pbsLengthA),sem(pbsLengthB),sem(pbsLengthC);
    sem(tmev2LengthA),sem(tmev2LengthB),sem(tmev2LengthC);
    sem(tmev5LengthA),sem(tmev5LengthB),sem(tmev5LengthC);
    sem(tmev15LengthA),sem(tmev15LengthB),sem(tmev15LengthC)];

figure()
aucBar = bar(zoneAuc,'FaceColor', 'flat','BarWidth',0.6,'EdgeColor','none');
% sigstar({[1,2],[1,3]},[0.05,0.05]);
% aucBar(1).CData(1,:) = [0 0 0];
% aucBar(1).CData(2,:) = [0 0 0];
% aucBar(1).CData(3,:) = [0 0 0];
% 
% aucBar(2).CData(1,:) = [1 0 0];
% aucBar(2).CData(2,:) = [1 0 0];
% aucBar(2).CData(3,:) = [1 0 0];
% 
% aucBar(3).CData(1,:) = [1 0 0];
% aucBar(3).CData(2,:) = [1 0 0];
% aucBar(3).CData(3,:) = [1 0 0];
% 
% aucBar(4).CData(1,:) = [1 0 0];
% aucBar(4).CData(2,:) = [1 0 0];
% aucBar(4).CData(3,:) = [1 0 0];

set(gca, 'XTickLabel', {})
box off
nbars = size(zoneAuc,2);
x = [];
for i = 1:nbars
    x = [x; aucBar(i).XEndPoints];
end
hold on
errorbar(x', zoneAuc, allSEM, 'k','linestyle','none','CapSize',3,'LineWidth',1,'Color',[0.4, 0.4, 0.4]);
ylimit = ylim;
yticks(ylimit(2));
yTickL = ChangeTextColor(num2str(ylimit(2)), [0, 0, 0]);
set(gca,'YTickLabel',yTickL);
set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
set(gca, 'FontName','Arial');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [4 4 xDim yDim]);