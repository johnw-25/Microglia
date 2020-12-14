function PairedBarGraph(barTitle, xDim, yDim, varargin)
set(0,'defaultfigurecolor',[1 1 1])
if length(varargin) == 4 % full comparison across PBS and all TMEV
    for i = 1:length(varargin)
        if isempty(varargin{i}.Zone1_AUC)
            varargin{i}.Zone1_AUC = 0;
            varargin{i}.Zone1_Locs = 0;
            varargin{i}.Zone1_PeakDurs = 0;
            varargin{i}.Zone1_PeakMags = 0;
        elseif isempty(varargin{i}.Zone2_AUC)
            varargin{i}.Zone2_AUC = 0;
            varargin{i}.Zone2_Locs = 0;
            varargin{i}.Zone2_PeakDurs = 0;
            varargin{i}.Zone2_PeakMags = 0;
        elseif isempty(varargin{i}.Zone3_AUC)
            varargin{i}.Zone3_AUC = 0;
            varargin{i}.Zone3_Locs = 0;
            varargin{i}.Zone3_PeakDurs = 0;
            varargin{i}.Zone3_PeakMags = 0;
        else
            %nothing
        end
    end
            

    zoneGroups = [ones(3,1), ones(3,1)*2, ones(3,1)*3];
    % % % Do AUC First % % %
    zoneAuc = [max(varargin{1}.Zone1_AUC),max(varargin{2}.Zone1_AUC),max(varargin{3}.Zone1_AUC),max(varargin{4}.Zone1_AUC);
        max(varargin{1}.Zone2_AUC),max(varargin{2}.Zone2_AUC),max(varargin{3}.Zone2_AUC),max(varargin{4}.Zone2_AUC);
        max(varargin{1}.Zone3_AUC),max(varargin{2}.Zone3_AUC),max(varargin{3}.Zone3_AUC),max(varargin{4}.Zone3_AUC)];
    
    allSEM = [sem(varargin{1}.Zone1_AUC),sem(varargin{2}.Zone1_AUC),sem(varargin{3}.Zone1_AUC),sem(varargin{4}.Zone1_AUC),
        sem(varargin{1}.Zone2_AUC),sem(varargin{2}.Zone2_AUC),sem(varargin{3}.Zone2_AUC),sem(varargin{4}.Zone2_AUC),
        sem(varargin{1}.Zone3_AUC),sem(varargin{2}.Zone3_AUC),sem(varargin{3}.Zone3_AUC),sem(varargin{4}.Zone3_AUC)];
    figure()
    aucBar = bar(zoneAuc,'FaceColor', 'flat','BarWidth',0.6,'EdgeColor','none');
    sigstar({[1,2],[1,3]},[0.05,0.05]);
    set(aucBar(1), 'FaceColor', [0.2 0.2 0.2])
    set(aucBar(2), 'FaceColor', [1 0 0])
    set(aucBar(3), 'FaceColor', [1 0 0])
    set(aucBar(4), 'FaceColor', [1 0 0])
    set(gca, 'XTickLabel', {})
    box off
    title(strcat(barTitle,' Signal Area'),'FontWeight','Normal');
    nbars = size(zoneAuc,2);
    x = [];
    for i = 1:nbars
        x = [x; aucBar(i).XEndPoints];
    end
    hold on
    errorbar(x', zoneAuc, allSEM, 'k','linestyle','none','CapSize',3,'LineWidth',1);
    ylimit = ylim;
    yticks(ylimit(2));
    yTickL = ChangeTextColor(num2str(round(ylimit(2))), [0, 0, 0]);
    set(gca,'YTickLabel',yTickL);
    set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
    set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
    set(gca, 'FontName','Arial');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [4 4 xDim yDim]);

    saveas(gcf, strcat(barTitle,' Signal Area Bar Graph.jpg'));
    % Peak Duration %
    zonePeakDurs = [max(varargin{1}.Zone1_PeakDurs),max(varargin{2}.Zone1_PeakDurs),max(varargin{3}.Zone1_PeakDurs),max(varargin{4}.Zone1_PeakDurs);
        max(varargin{1}.Zone2_PeakDurs),max(varargin{2}.Zone2_PeakDurs),max(varargin{3}.Zone2_PeakDurs),max(varargin{4}.Zone2_PeakDurs);
        max(varargin{1}.Zone3_PeakDurs),max(varargin{2}.Zone3_PeakDurs),max(varargin{3}.Zone3_PeakDurs),max(varargin{4}.Zone3_PeakDurs)];
    
    allSEM = [sem(varargin{1}.Zone1_PeakDurs),sem(varargin{2}.Zone1_PeakDurs),sem(varargin{3}.Zone1_PeakDurs),sem(varargin{4}.Zone1_PeakDurs),
        sem(varargin{1}.Zone2_PeakDurs),sem(varargin{2}.Zone2_PeakDurs),sem(varargin{3}.Zone2_PeakDurs),sem(varargin{4}.Zone2_PeakDurs),
        sem(varargin{1}.Zone3_PeakDurs),sem(varargin{2}.Zone3_PeakDurs),sem(varargin{3}.Zone3_PeakDurs),sem(varargin{4}.Zone3_PeakDurs)];
    figure()
    PeakDursBar = bar(zonePeakDurs,'FaceColor', 'flat','BarWidth',0.6,'EdgeColor','none');
    set(PeakDursBar(1), 'FaceColor', [0.2 0.2 0.2])
    set(PeakDursBar(2), 'FaceColor', [1 0 0])
    set(PeakDursBar(3), 'FaceColor', [1 0 0])
    set(PeakDursBar(4), 'FaceColor', [1 0 0])
    set(gca, 'XTickLabel', {})
    box off
    title(strcat(barTitle,' Duration'),'FontWeight','Normal');
    
    nbars = size(zonePeakDurs,2);
    x = [];
    for i = 1:nbars
        x = [x; PeakDursBar(i).XEndPoints];
    end
    hold on
    errorbar(x', zonePeakDurs, allSEM, 'k','linestyle','none','CapSize',3);
    ylimit = ylim;
    yticks(round(ylimit(2)));
    yTickL = ChangeTextColor(num2str(round(ylimit(2))), [0, 0, 0]);
    set(gca,'YTickLabel',yTickL);
    set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
    set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
    set(gca, 'FontName','Arial');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [4 4 xDim yDim]);
    saveas(gcf, strcat(barTitle,' Peak duration Bar Graph.jpg'));
    % Peak Amp %
    zonePeakMags = [max(varargin{1}.Zone1_PeakMags),max(varargin{2}.Zone1_PeakMags),max(varargin{3}.Zone1_PeakMags),max(varargin{4}.Zone1_PeakMags);
        max(varargin{1}.Zone2_PeakMags),max(varargin{2}.Zone2_PeakMags),max(varargin{3}.Zone2_PeakMags),max(varargin{4}.Zone2_PeakMags);
        max(varargin{1}.Zone3_PeakMags),max(varargin{2}.Zone3_PeakMags),max(varargin{3}.Zone3_PeakMags),max(varargin{4}.Zone3_PeakMags)];
    
    allSEM = [sem(varargin{1}.Zone1_PeakMags),sem(varargin{2}.Zone1_PeakMags),sem(varargin{3}.Zone1_PeakMags),sem(varargin{4}.Zone1_PeakMags),
        sem(varargin{1}.Zone2_PeakMags),sem(varargin{2}.Zone2_PeakMags),sem(varargin{3}.Zone2_PeakMags),sem(varargin{4}.Zone2_PeakMags),
        sem(varargin{1}.Zone3_PeakMags),sem(varargin{2}.Zone3_PeakMags),sem(varargin{3}.Zone3_PeakMags),sem(varargin{4}.Zone3_PeakMags)];
    figure()
    PeakMagsBar = bar(zonePeakMags,'FaceColor', 'flat','BarWidth',0.6,'EdgeColor','none');
    set(PeakMagsBar(1), 'FaceColor', [0.2 0.2 0.2])
    set(PeakMagsBar(2), 'FaceColor', [1 0 0])
    set(PeakMagsBar(3), 'FaceColor', [1 0 0])
    set(PeakMagsBar(4), 'FaceColor', [1 0 0])
    set(gca, 'XTickLabel', {})
    box off
    title(strcat(barTitle,' Amplitude'),'FontWeight','Normal');
    
    nbars = size(zonePeakMags,2);
    x = [];
    for i = 1:nbars
        x = [x; PeakMagsBar(i).XEndPoints];
    end
    hold on
    errorbar(x', zonePeakMags, allSEM, 'k','linestyle','none','CapSize',3);
    ylimit = ylim;
    yticks(ylimit(2));
    yTickL = ChangeTextColor(num2str(ylimit(2)), [0, 0, 0]);
    set(gca,'YTickLabel',yTickL);
    set(gcf, 'Units','inches','position',[4 4 xDim yDim]);
    set(gca,'FontSize',10,'XColor', [1 1 1],'YColor', [137/255 137/255 137/255],'FontName','Arial');
    set(gca, 'FontName','Arial');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [4 4 xDim yDim]);
    saveas(gcf, strcat(barTitle,' Peak amplitude Bar Graph.jpg'));

else
    error('input all groups please')
end
end