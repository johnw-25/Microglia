function [R_square, Ypred, mdl, CI] = LinearAnalysis(x, y, aTitle, varNames, yMax,xMax)
%%% Robust function that attempts to generate a model for input int,
%%% single, double vector, x.
set(0,'defaultfigurecolor',[1 1 1])
if sum(size(x) > 1 & size(y) > 1) == 2
    error('input must be vectors')
end

if size(x,2) > 1
    x = x';
end

if size(y,2) > 1
    y = y';
end

mdl = fitlm(x, y);                                  % Fit Data
B = mdl.Coefficients.Estimate;                      % Coefficients
CI = coefCI(mdl);                                   % Coefficient Confidence Intervals                     % Fitted Regression Line & Confidence Intervals
[~,YCI] = predict(mdl, x); 
SEM = sem(y);
ts = tinv([0.05 0.95],length(y)-1);
CI = mean(y)+ts*SEM;
% yOuts = y(y<CI(1) | y>CI(2));
% xOuts = x(y<CI(1) | y>CI(2));
% yNew = y(y>=CI(1) & y<=CI(2));
% xNew = x(y>=CI(1) & y<=CI(2));
% newMdl = fitlm(xNew, yNew);
R_square = mdl.Rsquared.Ordinary;
[Ypred,~] = predict(mdl, x); 
figure()
plot(x, y, 'ob')
hold on
plot(x, Ypred,'-r')
hold off
t = title(strcat(aTitle,strcat(' R^2: ',num2str(R_square))));
grid
yl = ylim;
yticks(yMax);
yTickMax = ChangeTextColor(num2str(round(yMax)), [0, 0, 0]);
roiTicks = {yTickMax};
ylim([-2 yMax])
xlim([-2 xMax])
% xlimit = xlim();
% xlim(xlimit);
print(gcf,strcat(aTitle,'.png'), '-r900','-dpng');
% xticks([xlimit(1) 0 xlimit(2)]);
% timeLabel = sprintf('\\color[rgb]{%f, %f, %f}%s', [0 0 0], 'Time (min)');
% minTimeTick = sprintf('\\color[rgb]{%f, %f, %f}%s', [0 0 0], num2str(xlimit(1)));
% time0Tick = sprintf('\\color[rgb]{%f, %f, %f}%s', [0 0 0], num2str(0));
% maxTimeTick = sprintf('\\color[rgb]{%f, %f, %f}%s', [0 0 0], num2str(xlimit(2)));
% timeTicks = {minTimeTick, time0Tick, maxTimeTick};
% set(gca,'XTickLabel',timeTicks,'YTickLabel', roiTicks);
% xlabel(timeLabel);
set(t,'FontSize',8);
dF = mdl.DFE;
rmse = mdl.RMSE;
slope = mdl.Coefficients.Estimate(2);
slopePvalue = mdl.Coefficients.pValue(2);
statsTable = table(R_square,dF, rmse, slope, slopePvalue);
resultsTable = table(x, y, Ypred,'VariableNames',varNames);
fileName = strcat(aTitle,'.xlsx');
writetable(resultsTable,fileName,'Sheet','data')
writetable(statsTable,fileName,'Sheet','model stats')

end
