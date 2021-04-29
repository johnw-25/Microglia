function [CENTROIDS,peak_data, timehitr] = SomaCentroids(excelFile, fallList, videos, excelList)
pbsRadius = 71*1.5; tmev2Radius = 72*1.5; tmev5Radius = 60*1.5; tmev15Radius = 60*1.5;

% Loop on each folder
MyData = cell(size(fallList,1),1);
for i = 1:length(fallList) % 1 to number of Fall.mat files
    filetoread = fullfile(fallList(i).folder,'Fall.mat'); % extract suite2p data from video subfolders
    % store each Fall.mat file into MyData
    fileID = fopen(filetoread);
    MyData{i} = load(filetoread); % the format depends of your files
    fclose(fileID);
end
Peak_SDs = 1.8; filterOrder = 100; Fc = 0.06; bp = 400;
% perform automatic event detection
peak_data = AutoEventDetect(MyData,videos,Peak_SDs,filterOrder,Fc,bp);
timehitr = cell(length(fallList),2);
th = 0:pi/50:2*pi;
excelData = readtable(excelFile, 'Sheet', 2);
for i = 1:length(fallList)
    % get data from burn tracks
    excelToRead = fullfile(excelList(i).folder,excelList(i).name);
    fileID = fopen(excelToRead);
    currentExcel = readtable(excelToRead);
    currentVideo = str2double(videos(i).name);
    subData = excelData(excelData.Image == currentVideo,:);
    groups = subData.Groups;
    dpi = subData.DPI;
    timeHitR = subData.TimeHitRadius1*1.04;
    timeHitR = mean(timeHitR);
    
    timehitr{i,1} = currentVideo;
    timehitr{i,2} = timeHitR;
    fclose(fileID);
    % retrieve registration shifts from suite2p output
    xoff = MyData{i,1}.ops.xoff;
    yoff = MyData{i,1}.ops.yoff;
    [~,xIdx] = max(abs(xoff));
    [~,yIdx] = max(abs(yoff));
    xShift = 2.62*xoff(xIdx);
    yShift = 2.62*yoff(yIdx);
    xM = round(currentExcel.XM.*2.6242);
    yM = round(currentExcel.YM.*2.6242);
    % exclude burn region from imagej COM coordinates
    logicals = [xM > 240 & xM < 265, yM > 230 & yM < 300];
    removeBurn = ~(logicals(:,1) == 1 & logicals(:,2) == 1);
    xM = xM(removeBurn);
    yM = yM(removeBurn);
    xM = xM-xShift;
    yM = yM-yShift;
    
    % find out which radius to use
    if strcmp(groups(1),'PBS')
        r = pbsRadius;
    else %TMEV
        if dpi(1) == 2
            r = tmev2Radius;
        elseif dpi(1) >=14
            r = tmev15Radius;
        else
            r = tmev5Radius;
        end
    end
    
    % get data from current video
    stat = MyData{i,1}.stat;
    events = peak_data{i,2}.locs;
    tempRoi = (1:length(stat))'-1;
    locs = events(:,1);
    %Open mean image image (1 plane)
    [filename, pathname] = uigetfile('*.tif', 'select mean image');
    cd(pathname);
    %set file to path string
    file = [pathname filename];
    %Open Fluorescence info
    image = imread(file);
    %measure size of image
    [ypixels,xpixels]=size(image);
    cROI=zeros(ypixels,xpixels);
    figure();
    hold on
    imshow(image)%,'InitialMagnification','fit')
    truesize([512,512])
    unassignedROIs = ones(1,size(stat,2));
    unassignedROIs = logical(unassignedROIs);
    allLocs = peak_data{i,2}.locs;
    for aa = 1:length(xM) % assume every center of mass is a cell
        tempCell = [xM(aa), yM(aa)];
        cellStr = strcat('Cell',num2str(aa));
        CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).circle = [xM(aa),yM(aa),r];
        for jj=1:length(stat)
            center = stat{1,jj}.med;
            roiStr = strcat('ROI',num2str(jj-1));
            tempLocs = allLocs{jj};
            for ii=1:length(stat{1,jj}.lam)
                cROI(stat{1,jj}.ypix(ii),stat{1,jj}.xpix(ii))=stat{1,jj}.lam(ii);
            end
            x = center(2); y = center(1);
            if (x-tempCell(1))^2 + (y-tempCell(2))^2 <= r^2
                unassignedROIs(jj) = 0;
                % if the roi matches with a cell, count the events in that
                % roi with fields we matched earlier
                % initialize number of events
                phase1tempLocs = tempLocs(tempLocs <=90);
                phase2tempLocs = tempLocs(tempLocs > 90 & tempLocs < timeHitR);
                phase3tempLocs = tempLocs(tempLocs >= timeHitR);
                if tempRoi(jj)+1 == jj
                    if ~isfield(CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr),roiStr)
                        CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).(roiStr).Phase1.NumEvents = numel(phase1tempLocs);
                        CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).(roiStr).Phase2.NumEvents = numel(phase2tempLocs);
                        CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).(roiStr).Phase3.NumEvents = numel(phase3tempLocs);
                        %                                             CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).(roiStr).NumEvents = numel(tempLocs);
                        CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).(roiStr).circle = [x,y,r];
                        CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).(roiStr).Mask = cROI;
                    else
                        old1Count = CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).(roiStr).Phase1.NumEvents;
                        old2Count = CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).(roiStr).Phase2.NumEvents;
                        old3Count = CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).(roiStr).Phase3.NumEvents;
%                         oldCount = CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).(roiStr).NumEvents;
                        CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).(roiStr).Phase1.NumEvents = old1Count + numel(phase1tempLocs);
                        CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).(roiStr).Phase2.NumEvents = old2Count + numel(phase2tempLocs);
                        CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).(roiStr).Phase3.NumEvents = old3Count + numel(phase3tempLocs);
                        %                                             CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).(roiStr).NumEvents = oldCount + numel(tempLocs);
                    end
                end
                
                CENTROIDS.(strcat('T',num2str(currentVideo))).(cellStr).circle = [tempCell(1),tempCell(2),r];
            else
                
            end
        end
    end
    unassignedStat = stat(unassignedROIs);
    [~,idx] = find(unassignedROIs == 1);
    unassignedRoiNums = idx-1;
    CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.temp = 0;
    CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned = rmfield(CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned,'temp');
    for w = 1:length(unassignedRoiNums)
        for ii=1:length(unassignedStat{1,w}.lam)
            cROI(unassignedStat{1,w}.ypix(ii),unassignedStat{1,w}.xpix(ii))=unassignedStat{1,w}.lam(ii);
        end
        for xx = 1:length(stat)
            if tempRoi(xx) == unassignedRoiNums(w)
                roiStr = strcat('ROI',num2str(unassignedRoiNums(w)));
                tempLocs = locs{xx};
                                phase1tempLocs = tempLocs(tempLocs <=90);
                                phase2tempLocs = tempLocs(tempLocs > 90 & tempLocs < timeHitR);
                                phase3tempLocs = tempLocs(tempLocs >= timeHitR);
                if ~isfield(CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned,roiStr)
                                        CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiStr).Phase1.NumEvents = numel(phase1tempLocs);
                                        CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiStr).Phase2.NumEvents = numel(phase2tempLocs);
                                        CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiStr).Phase3.NumEvents = numel(phase3tempLocs);
%                     CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiStr).NumEvents = numel(tempLocs);
                    CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiStr).circle = [x,y,r];
                    CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiStr).Mask = cROI;
                else
                                        old1Count = CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiStr).Phase1.NumEvents;
                                        old2Count = CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiStr).Phase2.NumEvents;
                                        old3Count = CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiStr).Phase3.NumEvents;
                                        CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiStr).Phase1.NumEvents = old1Count + numel(phase1tempLocs);
                                        CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiStr).Phase2.NumEvents = old2Count + numel(phase2tempLocs);
                                        CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiStr).Phase3.NumEvents = old3Count + numel(phase3tempLocs);
%                     oldCount = CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiStr).NumEvents;
%                     CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiStr).NumEvents = oldCount + numel(tempLocs);
                end
            else
            end
        end
    end
    
    hold on
    roiColorMask=zeros(ypixels,xpixels,3);
    cellFields = fieldnames(CENTROIDS.(strcat('T',num2str(currentVideo))));
    % check for overlapped ROIs and assign to only one cell
    for a = 1:numel(cellFields)
        firstRoiFields = fieldnames(CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{a}));
        if strcmp(cellFields{a},'Unassigned')
            continue
        end
        firstCenter = CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{a}).circle;
        for b = 1:numel(cellFields)%length(cellFields)-a+1;
            if strcmp(cellFields{b},'Unassigned')
                continue
            end
            secondRoiFields = fieldnames(CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{b}));
            secondCenter = CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{b}).circle;
            for c = 1:length(firstRoiFields)
                if strcmp(firstRoiFields{c},'circle')
                    continue
                end
                for d = 1:length(secondRoiFields)
                    if strcmp(firstRoiFields{c},secondRoiFields{d})
                        if strcmp(secondRoiFields{d},'circle')
                            continue
                        end
                        %do stuff
                        if ~isfield(CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{a}),firstRoiFields{c})
                            continue
                        end
                        if ~isfield(CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{b}),secondRoiFields{d})
                            continue
                        end
                        firstRoiC = CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{a}).(firstRoiFields{c}).circle;
                        secondRoiC = CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{b}).(secondRoiFields{d}).circle;
                        d1 = sqrt((firstRoiC(1)-firstCenter(1))^2+(firstRoiC(2)-firstCenter(2))^2);
                        d2 = sqrt((secondRoiC(1)-secondCenter(1))^2+(secondRoiC(2)-secondCenter(2))^2);
                        if d1 < d2
                            CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{b}) = rmfield(CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{b}), secondRoiFields{d});
                        elseif d1 > d2
                            CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{a}) = rmfield(CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{a}), firstRoiFields{c});
                        else
                            %nothing
                        end
                    end
                end
            end
        end
    end
    
    
    for x = 1:numel(cellFields)
        if strcmp(cellFields{x},'Unassigned')
            continue
        end
        cellLoc = CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{x}).circle;
        cellID = cellFields{x};
        cellID = cellID(5:end);
        sX = cellLoc(3) * cos(th) + cellLoc(1);
        sY = cellLoc(3) * sin(th) + cellLoc(2);
        h = fill(sX,sY,'r');
        h.FaceAlpha = 0.2;
        cellText = text(cellLoc(1),cellLoc(2),cellID,'Color',[1,1,1]);
        uistack(cellText,'top');
        roiFields = fieldnames(CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{x}));
        for y = 1:numel(roiFields)
            if isstruct(CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{x}).(roiFields{y}))
                roiColorMask(:,:,1)=rand(1,1);
                roiColorMask(:,:,2)=rand(1,1);
                roiColorMask(:,:,3)=rand(1,1);
                cellEx=imshow(roiColorMask);
                roiMask = CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{x}).(roiFields{y}).Mask;
                roiLoc = CENTROIDS.(strcat('T',num2str(currentVideo))).(cellFields{x}).(roiFields{y}).circle;
                set(cellEx,'AlphaData',(cROI)./(0.9*max(max(roiMask))))
                roiText = text(roiLoc(1),roiLoc(2),cellID,'Color',[1,1,1],'BackgroundColor',[0,0,0]);
                uistack(roiText,'top');
            end
        end
    end
    if isfield(CENTROIDS.(strcat('T',num2str(currentVideo))),'Unassigned')
        roiFields = fieldnames(CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned);
        for y = 1:numel(roiFields)
            if isstruct(CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiFields{y}))
                roiColorMask(:,:,1)=rand(1,1);
                roiColorMask(:,:,2)=rand(1,1);
                roiColorMask(:,:,3)=rand(1,1);
                cellEx=imshow(roiColorMask);
                roiMask = CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiFields{y}).Mask;
                roiLoc = CENTROIDS.(strcat('T',num2str(currentVideo))).Unassigned.(roiFields{y}).circle;
                set(cellEx,'AlphaData',(roiMask)./(0.9*max(max(roiMask))))
%                 text(roiLoc(1),roiLoc(2),cellID,'Color',[1,1,1]);
            end
        end
    end
    set(gcf,'MenuBar','none')
    set(gca,'DataAspectRatio',[1,1,1])
    scaleBar = 2.6241 * 20;
    line([400, 400+scaleBar],[480,480],'Color',[1,1,1],'LineWidth',1.25)
    hold off
    print(gcf, strcat(videos(i).name,'_ROIs.tif'), '-r900','-dtiff');
end



end