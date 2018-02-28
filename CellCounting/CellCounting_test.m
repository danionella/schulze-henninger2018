%% Counting neurons in brain slices: test example using Slice14.tif file
clear all
close all

%% load
% load image
img=TIFFStack('./Slice14.tif'); % TIFFStack function can be downloaded from https://github.com/DylanMuir/TIFFStack

% create 2D gaussian function for smoothing
sigma=3; % if smaller, you get several maxima within one cell, if larger you miss some cells
[X Y]=meshgrid(1:100);
g=exp(-((X-50).^2+(Y-50).^2)/(2*sigma^2));

% threshold values (previously manually adjusted in Fiji), works better than automatic thresholding
thrMan=[166 229 223 208 183 206 171 199 214 206 192 167 192 173 194 201 195 168 166 208 155 169 219 216 192 148 199 173];

% manual roi selection of brain region (also done in Fiji)
brainArea=ReadImageJROI('./RoiSetStackRGB.zip'); % function from https://github.com/DylanMuir/ReadImageJROI

%% image processing (threshold, smooth, detect maxima)

% threshold the images using the manually determined values
disp('thresholding')
img.bInvert = false; % the threshold values have been determined with light background
maskThr=double((img(:,:,1)<thrMan(14))); %(1-->cell,0-->no cell)

% select brain region
disp('selecting')
maskSelec=0.*img(:,:,1);
x=brainArea{14}.mnCoordinates(:,1);
y=brainArea{14}.mnCoordinates(:,2);
maskSelec=roipoly(maskSelec,x,y); %(1-->brain,0-->no brain)

% crop image around ROI
disp('cropping')
y1=brainArea{14}.vnRectBounds(1)-100;
y2=brainArea{14}.vnRectBounds(3)+100;
x1=brainArea{14}.vnRectBounds(2)-100;
x2=brainArea{14}.vnRectBounds(4)+100;

img.bInvert = true;  % change to black background
clear im
im=img(y1:y2,x1:x2,1);

% smooth the cropped around brain but not thresholded image by convolving with 2D gaussian
disp('convolving')
clear gImg
gImg=conv2(double(im),g,'same');

% count local maxima in smoothed image
disp('counting')
A=imregionalmax(gImg);

% count only the maxima which are within cells (maskThr=1) and within brain (maskSelec=1)
bw=maskThr(y1:y2,x1:x2).*maskSelec(y1:y2,x1:x2).*A;
count=length(find(bw));

% plot results
[i j]=find(bw);
figure(666), imagesc(im), axis image, colormap gray
hold on, scatter(j,i,30,'r','+'),hold off
title(['There are ' num2str(length(i)) ' maxima here'])
drawnow

Ntotal=sum(count);

