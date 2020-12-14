function [Rsquare, fitLine, beta] = StandardizedRegression(x, y)

if sum(size(x) > 1 & size(y) > 1) == 2
    error('input must be vectors')
end

if size(x,2) > 1
    x = x';
end

if size(y,2) > 1
    y = y';
end

% first check for normality using lillietest()
h = lillietest(x);


% if h
%     % turn data into normal distribution
%     z = zscore(x);
%     hz = lillietest(x);
%     if hz
%         xAvg = mean(z);
%         xSD = std(z);
%         xNew = x(z < xAvg+(2*xSD));
%         xOuts = x(x >= xAvg+(2*xSD));
%     end
% else % data is already normal
%     xAvg = mean(x);
%     xSD = std(x);
%     xNew = x(x < xAvg+(2*xSD));
%     xOuts = x(x >= xAvg+(2*xSD));
% end

hy = lillietest(y);

if hy
    % turn data into normal distribution
    z = zscore(y);
    hzy = lillietest(y);
    if hzy
        yAvg = mean(z);
        ySD = std(z);
        yNew = y(z < yAvg+(2*ySD));
        yOuts = y(z >= yAvg+(2*ySD));
        xNew = x(z < yAvg+(2*ySD));
        xOuts = x(z >= yAvg+(2*ySD));
    end
else
    yAvg = mean(y);
    ySD = std(y);
    yNew = y(y < yAvg+(2*ySD));
    yOuts = y(y >= yAvg+(2*ySD));
    xNew = x(y < yAvg+(2*ySD));
    xOuts = x(y >= yAvg+(2*ySD));
end

X = [ones(length(xNew),1) xNew];
b = X\yNew; % b(1) is intercept, b(2) is slope

fitLine = X*b;
R_square = 1 - sum((yNew - fitLine).^2)/sum((yNew - mean(yNew)).^2);
% try polyfit
p = polyfit(xNew, yNew, 1); % P(1) is slope, p(2) is intercept
fitLinePoly = polyval(p,xNew); % y = p(1)*x + p(2) is the same thing
polyResiduals = sum((yNew - fitLinePoly).^2);
totalResiduals = (length(yNew)-1 * var(yNew));
polyRsq = 1 - polyResiduals/totalResiduals;
if polyRsq > R_square
    R_square = polyRsq;
    fitLine = fitLinePoly;
end

h = figure();
hold on
plot(xNew,yNew, 'o');
plot(xNew, fitLine,'--');
if ~isempty(xOuts) || ~isempty(yOuts)
    plot(xOuts, yOuts, 'or', 'MarkerFaceColor','r');
end
legend('Data','Linear model','Outliers');

end