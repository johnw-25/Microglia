function CENTROIDS = SomaCentroids(excelFile, fallList, videos)
pbsRadius = 71; tmev2radius = 72; tmev5Radius = 60; tmev15Radius = 60;

% Loop on each folder
MyData = [];
for i = 1:length(fallList) % 1 to number of Fall.mat files
  filetoread = fullfile(fallList(i).folder,'Fall.mat'); % extract suite2p data from video subfolders
  % store each Fall.mat file into MyData
  fileID = fopen(filetoread);
  MyData{end+1} = load(filetoread); % the format depends of your files
  fclose(fileID);
end

test_data = readtable(excelFile,'Sheet',2);
somaRoi = test_data.Final_calcium_soma;

for i = 1:length(fallList)
    currentVideo = str2double(videos(i).name);    
    subTable = test_data(test_data.Image == currentVideo,:);
    somaRoi = subTable.Final_calcium_soma;
    subTable = subTable(~isnan(somaRoi),:);
    somaRoi = somaRoi(~isnan(somaRoi))+1;
    groups = subTable.Groups;
    dpi = subTable.DPI;
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
        
    stat = MyData{1, i}.stat;  
    %Open mean image image (1 plane)
    [filename, pathname] = uigetfile('*.tif', 'select mean image');
    cd(pathname);
    %set file to path string
    file = [pathname filename];
    %Open Fluorescence info
    image = imread(file);   
    %measure size of image
    [ypixels,xpixels]=size(image);
    somaColorMask=zeros(ypixels,xpixels,3);
    somaColorMask(:,:,1)=0;
    somaColorMask(:,:,2)=0;
    somaColorMask(:,:,3)=1;
    cROI=zeros(ypixels,xpixels);
    fig = figure();
    hold on
    imshow(image)%,'InitialMagnification','fit')
    truesize([512,512])
    for jj=1:length(somaRoi)
        center = stat{1,somaRoi(jj)}.med;
        x = center(2); y = center(1);
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;
        for ii=1:length(stat{1,somaRoi(jj)}.lam)
            cROI(stat{1,somaRoi(jj)}.ypix(ii),stat{1,somaRoi(jj)}.xpix(ii))=stat{1,somaRoi(jj)}.lam(ii);
        end
        CENTROIDS.(strcat('T',num2str(currentVideo))).(strcat('ROI',num2str(somaRoi(jj)))).Mask = cROI;
        CENTROIDS.(strcat('T',num2str(currentVideo))).(strcat('ROI',num2str(somaRoi(jj)))).circle = [x,y,r];
        fieldNames = fields(CENTROIDS.(strcat('T',num2str(currentVideo))));
%         if jj==1
%             CENTROIDS.(strcat('T',num2str(currentVideo))).(strcat('ROI',num2str(somaRoi(jj)))).cellID = 1;
%         else
%             for x = 1:numel(fields(CENTROIDS.(strcat('T',num2str(currentVideo)))))
%                 
                
        end
        hold on
        somaEx=imshow(somaColorMask);
        set(somaEx,'AlphaData',(cROI)./(0.75*max(max(cROI))))
        h = fill(xunit,yunit,'r');
        h.FaceAlpha = 0.2;
        %     hold on
        %     processEx=imshow(processColorMask);
        %     set(processEx,'AlphaData',(cpROI)./(0.75*max(max(cpROI))))
%         set(gcf, 'Units','inches','position',[4 4 2.5 2.5]);
%         set(gcf, 'PaperPosition', [4 4 2.5 2.5]);
        set(gcf,'MenuBar','none')
        set(gca,'DataAspectRatioMode','auto')
        hold off
    end

    
    
end