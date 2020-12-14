function PairedSignalLength(DATA, Pathology)
set(0,'defaultfigurecolor',[1 1 1])
sLengthA = zeros(size(DATA,1),1);
sLengthB = zeros(size(DATA,1),1);
sLengthC = zeros(size(DATA,1),1);
LPF = noiseRemover(100, 0.02, 0.97);

for ii = 1:size(DATA,1)
    
    timehitR = double(DATA{ii,6});
    timehitR_marker = round(timehitR);
    currentTrace = DATA{ii,3};
    F0 = mean(currentTrace(1:50));
    currentTrace = (currentTrace-F0)./F0;
    currentTrace = filtfilt(LPF, double(currentTrace));
    if timehitR_marker > 200 && timehitR_marker < length(currentTrace)
        sLengthA(ii) = sum(currentTrace(1:200))/length(currentTrace(1:200));
        sLengthB(ii) = sum(currentTrace(201:timehitR_marker-1))/length(currentTrace((201:timehitR_marker-1)));
        sLengthC(ii) = sum(currentTrace(timehitR_marker:end))/length(currentTrace(timehitR_marker:end));
    else
        sLengthA(ii) = NaN;
        sLengthB(ii) = NaN;
        sLengthC(ii) = NaN;
    end
end

% remove NaN
sLengthA = sLengthA(~isnan(sLengthA));
sLengthB = sLengthB(~isnan(sLengthB));
sLengthC = sLengthC(~isnan(sLengthC));

meanLengthA = mean(sLengthA);
meanLengthB = mean(sLengthB);
meanLengthC = mean(sLengthC);

semA = std(sLengthA)/sqrt(length(sLengthA));
semB = std(sLengthB)/sqrt(length(sLengthB));
semC = std(sLengthC)/sqrt(length(sLengthC));

x = [0 1 2];
y = [meanLengthA, meanLengthB, meanLengthC];
semLength = [semA, semB, semC];

figure()
errorbar(x,y,semLength,'k');
hold on
plot(x,y,'kd','MarkerFaceColor','k');
xlabel('Zones'); ylabel('Average Signal Length'); title(strcat(Pathology,': Average Signal Length Between Zones (1-200, 201-timehitR-1, timehitR-end)'));

groups = [ones(length(sLengthA),1); ones(length(sLengthB),1)*2; ones(length(sLengthC),1)*3];
allData = [sLengthA; sLengthB; sLengthC];
figure()
boxplot(allData,groups);
title(strcat(Pathology,': Signal Length Boxplots of Zones')); xlabel('Zones'); ylabel('RMS');
end