function PlotSomaProcess(somaROI,processROI, figTitle,somaColor,processColor)

%%  Script is to plot ROI(s) of interest
%ROI should be scalar or vector of ROIs
%This relies on the fall.mat output from Suite2P. 

%select Fall.mat file
[filename, pathname] = uigetfile('*.mat', 'select Fall.mat file');
 %set current directory to pathname
    cd(pathname);
    %set file to path string
   file = [pathname filename]; 
   %Open Fluorescence info 
  load(file);

    
%update iscell with the user curated output from iscell
iscell=readNPY('iscell.npy');
         
%Open mean image image (1 plane)
[filename, pathname] = uigetfile('*.tif', 'select mean image');
 
    %set current directory to pathname
    cd(pathname);
    %set file to path string
   file = [pathname filename]; 
   %Open Fluorescence info 
   image = imread(file);

%measure size of image
[ypixels,xpixels]=size(image);

%make blank ROI field
csROI=zeros(ypixels,xpixels);
cpROI=zeros(ypixels,xpixels);

%Make background colors
somaColorMask=zeros(ypixels,xpixels,3);
somaColorMask(:,:,1)=somaColor(1);
somaColorMask(:,:,2)=somaColor(2);
somaColorMask(:,:,3)=somaColor(3);

processColorMask=zeros(ypixels,xpixels,3);
processColorMask(:,:,1)=processColor(1);
processColorMask(:,:,2)=processColor(2);
processColorMask(:,:,3)=processColor(3);

%Remove discarded ROIs from from stat based on your iscell selections
% stat=stat(1,logical(iscell(:,1)));
stat=stat;
%collect ROI pixel positions for all ROIs you wanted
for jj=1:length(somaROI)   
for ii=1:length(stat{1,somaROI(jj)}.lam)
csROI(stat{1,somaROI(jj)}.ypix(ii),stat{1,somaROI(jj)}.xpix(ii))=stat{1,somaROI(jj)}.lam(ii);
end
end

for jj=1:length(processROI)   
for ii=1:length(stat{1,processROI(jj)}.lam)
cpROI(stat{1,processROI(jj)}.ypix(ii),stat{1,processROI(jj)}.xpix(ii))=stat{1,processROI(jj)}.lam(ii);
end
end

%plot
fig = figure();
imshow(image)%,'InitialMagnification','fit')
truesize([512,512])
hold on
somaEx=imshow(somaColorMask);
set(somaEx,'AlphaData',(csROI)./(0.75*max(max(csROI))))
hold on
processEx=imshow(processColorMask);
set(processEx,'AlphaData',(cpROI)./(0.75*max(max(cpROI))))
set(gcf, 'Units','inches','position',[4 4 2.5 2.5]);
set(gcf, 'PaperPosition', [4 4 2.5 2.5]);
set(gcf,'MenuBar','none')
set(gca,'DataAspectRatioMode','auto')
hold off
print(fig, strcat(figTitle,'_ROIs.png'), '-r900','-dpng');



