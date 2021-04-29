function PlotSomaProcess(somaROI,processROI, figTitle,somaColor,processColor,imageFile,file,pathname)

%%  Script is to plot ROI(s) of interest
%ROI should be scalar or vector of ROIs
%This relies on the fall.mat output from Suite2P. 

%select Fall.mat file
% [filename, pathname] = uigetfile('*.mat', 'select Fall.mat file');
%  %set current directory to pathname
    cd(pathname);
%     %set file to path string
%    file = [pathname filename]; 
%    %%%Open Fluorescence info 
  load(file);

    
%update iscell with the user curated output from iscell
% iscell=readNPY('iscell.npy');
         
%Open mean image image (1 plane)
% [filename, pathname] = uigetfile('*.tif', 'select mean image');
 
    %set current directory to pathname
%     cd(pathname);
    %set file to path string
%    file = [pathname filename]; 
   %Open Fluorescence info 
   image = imread(imageFile);

%measure size of image
[ypixels,xpixels]=size(image);
blackImage = zeros(512,512,3,'uint8');
%make blank ROI field
csROI=zeros(ypixels,xpixels);
cpROI=zeros(ypixels,xpixels);

%Remove discarded ROIs from from stat based on your iscell selections
% stat=stat(1,logical(iscell(:,1)));
stat=stat;
%collect ROI pixel positions for all ROIs you wanted
for jj=1:length(somaROI)
    tempROI = zeros(ypixels,xpixels);
    tempColor = zeros(ypixels,xpixels,3);
    tempColor(:,:,1) = somaColor(jj,1);
    tempColor(:,:,2) = somaColor(jj,2);
    tempColor(:,:,3) = somaColor(jj,3);
for ii=1:length(stat{1,somaROI(jj)}.lam)
csROI(stat{1,somaROI(jj)}.ypix(ii),stat{1,somaROI(jj)}.xpix(ii))=stat{1,somaROI(jj)}.lam(ii);
tempROI(stat{1,somaROI(jj)}.ypix(ii),stat{1,somaROI(jj)}.xpix(ii))=stat{1,somaROI(jj)}.lam(ii);
end
SMASKS.(strcat('Soma',num2str(somaROI(jj)))).mask = tempROI;
SMASKS.(strcat('Soma',num2str(somaROI(jj)))).colorMap = tempColor;
end

for jj=1:length(processROI)
    tempROI = zeros(ypixels,xpixels);
    tempColor = zeros(ypixels,xpixels,3);
    tempColor(:,:,1) = processColor(jj,1);
    tempColor(:,:,2) = processColor(jj,2);
    tempColor(:,:,3) = processColor(jj,3);
    for ii=1:length(stat{1,processROI(jj)}.lam)
        cpROI(stat{1,processROI(jj)}.ypix(ii),stat{1,processROI(jj)}.xpix(ii))=stat{1,processROI(jj)}.lam(ii);
        tempROI(stat{1,processROI(jj)}.ypix(ii),stat{1,processROI(jj)}.xpix(ii))=stat{1,processROI(jj)}.lam(ii);
    end
    PMASKS.(strcat('Process',num2str(processROI(jj)))).mask = tempROI;
    PMASKS.(strcat('Process',num2str(processROI(jj)))).colorMap = tempColor;
end

%plot
fig = figure();
image = imclearborder(image);
% imshow(image)%,'InitialMagnification','fit')
imshow(blackImage)%,'InitialMagnification','fit')
truesize([512,512])
% hold on
% somaEx=imshow(somaColorMask);
% set(somaEx,'AlphaData',(csROI)./(0.75*max(max(csROI))))
hold on
pFields = fieldnames(PMASKS);
for k = 1:numel(pFields)
    mask = PMASKS.(pFields{k}).mask;
    color = PMASKS.(pFields{k}).colorMap;
    center = stat{1,processROI(k)}.med;
    x = center(2); y = center(1);
    processEx=imshow(color);
    set(processEx,'AlphaData',(mask)./(0.75*max(max(mask))))
    text(x,y,num2str(processROI(k)-1),'Color',[1,1,1],'BackgroundColor',[0,0,0],'FontSize',6);
end
sFields = fieldnames(SMASKS);
for k = 1:numel(sFields)
    mask = SMASKS.(sFields{k}).mask;
    color = SMASKS.(sFields{k}).colorMap;
    center = stat{1,somaROI(k)}.med;
    x = center(2); y = center(1);
    somaEx=imshow(color);
    set(somaEx,'AlphaData',(mask)./(0.75*max(max(mask))))
    text(x,y,num2str(somaROI(k)-1),'Color',[1,1,1],'BackgroundColor',[0,0,0],'FontSize',6);
end
% processEx=imshow(processColorMask);
% set(processEx,'AlphaData',(cpROI)./(0.75*max(max(cpROI))))
set(gcf, 'Units','inches','position',[4 4 2.5 2.5]);
set(gcf, 'PaperPosition', [4 4 5 5]);
set(gcf,'MenuBar','none')
set(gca,'DataAspectRatio',[1,1,1])
scaleBar = 2.6241 * 20;
line([400, 400+scaleBar],[480,480],'Color',[1,1,1],'LineWidth',1.25)
hold off
print(fig, strcat(figTitle,'_ROIs.tif'), '-r900','-dtiff');



