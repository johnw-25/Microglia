function [R_square, Ypred, mdl, CI] = LinearAnalysis(x, y,aTitle)
%%% Robust function that attempts to generate a model for input int,
%%% single, double vector, x.

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
yOuts = y(y<CI(1) | y>CI(2));
xOuts = x(y<CI(1) | y>CI(2));
yNew = y(y>=CI(1) & y<=CI(2));
xNew = x(y>=CI(1) & y<=CI(2));
% yOuts = y(y < YCI(:,1) | y > YCI(:,2));
% xOuts = x(y < YCI(:,1) | y > YCI(:,2));
% yNew = y(y >= YCI(:,1) & y <= YCI(:,2));
% xNew = x(y >= YCI(:,1) & y <= YCI(:,2));
newMdl = fitlm(xNew, yNew);
R_square = mdl.Rsquared.Ordinary;
[Ypred,~] = predict(mdl, x); 
figure()
plot(xNew, yNew, 'ob')
hold on
plot(x, Ypred,'-r')
% plot(xNew, Ypred,'-r',x, YCI(:,1), '.k',x, YCI(:,2), '.k') %  
% patch(x, YCI(:,1),[1 0 0])
% hold on
plot(xOuts, yOuts,'sr');
hold off
title(strcat(aTitle,strcat(' R^2: ',num2str(R_square))));
grid
end
