clear all
close all
clc

%function WaveletSeeding_ActiveContour_v1()
%SHOWALLPOSITIVESAMPLES Summary of this function goes here
%   Detailed explanation goes here

%addpath('../tight_subplot');
addpath('../thresh_tool');
addpath('../spotDetector_101410');
addpath('../cutFirstHistMode');
feature('UseOldFileDialogs', 1);

GUI = 1;
merge = 75;

%% Ratio of the pixel to real value
eachPixelLengthInMicro = 1.29;
safetyFactor = 0;

%% Size Constraints
MaxPSMAaxisInMicro = 30;
MaxDAPIaxisInMicro = 25;
minPSMAaxisInMicro = 5;
minDAPIaxisInMicro = 3;

PSMA.maxAxisLength = MaxPSMAaxisInMicro*(1+safetyFactor) / eachPixelLengthInMicro;
DAPI.maxAxisLength = MaxDAPIaxisInMicro*(1+safetyFactor) / eachPixelLengthInMicro;
PSMA.minAxisLength = minPSMAaxisInMicro*(1-safetyFactor) / eachPixelLengthInMicro;
DAPI.minAxisLength = minDAPIaxisInMicro*(1-safetyFactor) / eachPixelLengthInMicro;

PSMA.maxArea = pi * PSMA.maxAxisLength^2;
PSMA.minArea = pi * PSMA.minAxisLength^2;

DAPI.maxArea = pi * DAPI.maxAxisLength^2;
DAPI.minArea = pi * DAPI.minAxisLength^2;


histJump = 1;
% inputPath = 'E:\Dropbox\Matlab\Prostate_Cancer\New_folder\Cells';
%
% PSMA.im = imread(fullfile(inputPath,'9740-101Ficoll02whole_c1_ORG.tif'));
% DAPI.im = imread(fullfile(inputPath,'9740-101Ficoll02whole_c2_ORG.tif'));
% CD45.im = imread(fullfile(inputPath,'9740-101Ficoll02whole_c3_ORG.tif'));

inputPath = '';
inputPath = 'E:\MATLAB\CTC\Alex Dropbox\10xNewFilter\Flash Camera\10x';

%% Reads the files and crops them
if GUI
    [FileName,PathName,FilterIndex] = uigetfile({'*.tif';'*.jpg';'*.tiff';'*.*'},'Select PSMA image',inputPath);
    PSMA.origIm = double(imread(fullfile(PathName, FileName)));
    %PSMA.origIm = PSMA.origIm(2000:2300,2000:2300);
    PSMA.im = PSMA.origIm;
    
    [FileName,PathName,FilterIndex] = uigetfile({'*.tif';'*.jpg';'*.tiff';'*.*'},'Select DAPI image',PathName);
    DAPI.origIm = double(imread(fullfile(PathName, FileName)));
    %DAPI.origIm = DAPI.origIm(2000:2300,2000:2300);
    DAPI.im = DAPI.origIm;
    
    [FileName,PathName,FilterIndex] = uigetfile({'*.tif';'*.jpg';'*.tiff';'*.*'},'Select CD45 image',PathName);
    CD45.origIm = double(imread(fullfile(PathName, FileName)));
    %CD45.origIm = CD45.origIm(2000:2300,2000:2300);
    CD45.im = CD45.origIm;
end

if GUI
    [FileName,PathName_xls,FilterIndex] = uigetfile({'*.*';'*.xls'; '*.xlsx'},'Cancel if you don''t have annotation!',PathName);
    if FileName
        [x,y] = readAnnotation(fullfile(PathName_xls, FileName));
    else
        x = [];
        y = [];
    end
else
    [x,y] = readAnnotation(fullfile(inputPath,'Manual Coords9740-103Ficoll01whole_with_new-filter.xlsx'));
end
% inputPath = 'E:\MATLAB\CTC\DATA\9740-103-02';
%
% PSMA.im = imread(fullfile(inputPath,'9740-103-02_c1_ORG.tif'));
% DAPI.im = imread(fullfile(inputPath,'9740-103-02_c2_ORG.tif'));
% CD45.im = imread(fullfile(inputPath,'9740-103-02_c3_ORG.tif'));
%
% [x,y] = readAnnotation(fullfile(inputPath,'Manual Coords9740-103Ficoll02whole_with_new-filter.xlsx'));

%% Fixing the Brightness
Gfilt = fspecial('gaussian',250,35);

PSMA.brightness = imfilter(PSMA.im,Gfilt,'replicate','same');

PSMA.im = 100*PSMA.im ./ max(PSMA.brightness,1);

%% Calculate the wavelet parallel
waveletSteps = 1003;
waveletBorder = 50;
partsI = ceil(size(PSMA.im,1)/waveletSteps);
partsJ = ceil(size(PSMA.im,2)/waveletSteps);
PSMAwavelet = cell(partsI,partsJ);
PSMAim = PSMA.im;
parfor i = 1:partsI
    for j = 1:partsJ
        [~,PSMAwavelet{i,j}] = spotDetector(double(PSMAim(max(1,waveletSteps*(i-1)-waveletBorder+1):min(waveletSteps*i+waveletBorder,size(PSMAim,1)),max(1,waveletSteps*(j-1)-waveletBorder+1):min(waveletSteps*j+waveletBorder,size(PSMAim,2)))));
        PSMAwavelet{i,j} = (PSMAwavelet{i,j}>0);
        if i>1 && i<partsI
            PSMAwavelet{i,j} = PSMAwavelet{i,j}(waveletBorder+1:end-waveletBorder,:);
        elseif i>1
            PSMAwavelet{i,j} = PSMAwavelet{i,j}(waveletBorder+1:end,:);
        elseif i<partsI
            PSMAwavelet{i,j} = PSMAwavelet{i,j}(1:end-waveletBorder,:);
        end
        if j>1 && j<partsJ
            PSMAwavelet{i,j} = PSMAwavelet{i,j}(:,waveletBorder+1:end-waveletBorder);
        elseif j>1
            PSMAwavelet{i,j} = PSMAwavelet{i,j}(:,waveletBorder+1:end);
        elseif j<partsJ
            PSMAwavelet{i,j} = PSMAwavelet{i,j}(:,1:end-waveletBorder);
        end
        [~,DAPIwavelet{i,j}] = spotDetector(double(DAPIim(max(1,waveletSteps*(i-1)-waveletBorder+1):min(waveletSteps*i+waveletBorder,size(DAPIim,1)),max(1,waveletSteps*(j-1)-waveletBorder+1):min(waveletSteps*j+waveletBorder,size(DAPIim,2)))));
        DAPIwavelet{i,j} = (DAPIwavelet{i,j}>0);
        if i>1 && i<partsI
            DAPIwavelet{i,j} = DAPIwavelet{i,j}(waveletBorder+1:end-waveletBorder,:);
        elseif i>1
            DAPIwavelet{i,j} = DAPIwavelet{i,j}(waveletBorder+1:end,:);
        elseif i<partsI
            DAPIwavelet{i,j} = DAPIwavelet{i,j}(1:end-waveletBorder,:);
        end
        if j>1 && j<partsJ
            DAPIwavelet{i,j} = DAPIwavelet{i,j}(:,waveletBorder+1:end-waveletBorder);
        elseif j>1
            DAPIwavelet{i,j} = DAPIwavelet{i,j}(:,waveletBorder+1:end);
        elseif j<partsJ
            DAPIwavelet{i,j} = DAPIwavelet{i,j}(:,1:end-waveletBorder);
        end
% % %         [~,CD45wavelet{i,j}] = spotDetector(double(CD45im(max(1,waveletSteps*(i-1)-waveletBorder+1):min(waveletSteps*i+waveletBorder,size(CD45im,1)),max(1,waveletSteps*(j-1)-waveletBorder+1):min(waveletSteps*j+waveletBorder,size(CD45im,2)))));
% % %         CD45wavelet{i,j} = (CD45wavelet{i,j}>0);
% % %         if i>1 && i<partsI
% % %             CD45wavelet{i,j} = CD45wavelet{i,j}(waveletBorder+1:end-waveletBorder,:);
% % %         elseif i>1
% % %             CD45wavelet{i,j} = CD45wavelet{i,j}(waveletBorder+1:end,:);
% % %         elseif i<partsI
% % %             CD45wavelet{i,j} = CD45wavelet{i,j}(1:end-waveletBorder,:);
% % %         end
% % %         if j>1 && j<partsJ
% % %             CD45wavelet{i,j} = CD45wavelet{i,j}(:,waveletBorder+1:end-waveletBorder);
% % %         elseif j>1
% % %             CD45wavelet{i,j} = CD45wavelet{i,j}(:,waveletBorder+1:end);
% % %         elseif j<partsJ
% % %             CD45wavelet{i,j} = CD45wavelet{i,j}(:,1:end-waveletBorder);
% % %         end
    end
end

PSMA.wavelet = cell2mat(PSMAwavelet);
DAPI.wavelet = cell2mat(DAPIwavelet);
% % % CD45.wavelet = cell2mat(CD45wavelet);

clear PSMAwavelet;
clear DAPIwavelet;
% % % clear CD45wavelet;


%% Active Contour and Watershed

% % % DAPI.activeCont = activecontour(DAPI.im,DAPI.wavelet);
% % % DAPI.waveletWatershed = watershed(1-(DAPI.wavelet>0));
% % % DAPI.borderIm = DAPI.im - min(DAPI.im(:));
% % % DAPI.borderIm = DAPI.borderIm / max(DAPI.im(:));
% % % DAPI.borderIm = 1 - DAPI.borderIm;
% % % DAPI.borderIm(DAPI.waveletWatershed == 0) = 1;
% % % DAPI.activeContSeparate = activecontour(DAPI.borderIm, DAPI.wavelet);
% % % DAPI.segments = (DAPI.activeContSeparate>0) & (DAPI.activeCont>0);

PSMA.processedIm = PSMA.im;
sorted = sort(PSMA.processedIm(:));
p = round(997*length(PSMA.processedIm(:))/1000);
PSMA.processedIm(PSMA.processedIm>sorted(p)) = sorted(p);
p = round(3*length(PSMA.processedIm(:))/1000);
PSMA.processedIm(PSMA.processedIm<sorted(p)) = sorted(p);
PSMA.processedIm(PSMA.processedIm<0) = 0;

PSMA.processedIm = PSMA.processedIm - min(PSMA.processedIm(:));
PSMA.processedIm = PSMA.processedIm / max(PSMA.processedIm(:))*256;

waveletSteps = 1005;
waveletBorder = 50;
partsI = ceil(size(PSMA.im,1)/waveletSteps);
partsJ = ceil(size(PSMA.im,2)/waveletSteps);
activeCont = cell(partsI,partsJ);
waveletWatershed = activeCont;
borderIm = activeCont;
activeContSeparate = activeCont;
segments = activeCont;
warning('Careful! Here I am using DAPI instead of PSMA in the next line!');
PSMAwavelet = PSMA.wavelet;
DAPIwaveletWatershed = DAPI.waveletWatershed;
PSMAwaveletWatershed = PSMA.waveletWatershed;
PSMAim = PSMA.im;
PSMAprocessedIm = PSMA.processedIm;
parfor i = 1:partsI
    for j = 1:partsJ
        wavelet = PSMAwavelet(max(1,waveletSteps*(i-1)-waveletBorder+1):min(waveletSteps*i+waveletBorder,size(PSMAim,1)),max(1,waveletSteps*(j-1)-waveletBorder+1):min(waveletSteps*j+waveletBorder,size(PSMAim,2)));
        imToProcess = PSMAprocessedIm(max(1,waveletSteps*(i-1)-waveletBorder+1):min(waveletSteps*i+waveletBorder,size(PSMAim,1)),max(1,waveletSteps*(j-1)-waveletBorder+1):min(waveletSteps*j+waveletBorder,size(PSMAim,2)));
        waveletWatershed{i,j} = DAPIwaveletWatershed(max(1,waveletSteps*(i-1)-waveletBorder+1):min(waveletSteps*i+waveletBorder,size(PSMAim,1)),max(1,waveletSteps*(j-1)-waveletBorder+1):min(waveletSteps*j+waveletBorder,size(PSMAim,2)));
        activeCont{i,j} = activecontour(imToProcess,wavelet,1000);
        waveletWatershed{i,j} = uint16(watershed(1-(wavelet>0)));
        borderIm{i,j} = imToProcess - min(imToProcess(:));
        borderIm{i,j} = borderIm{i,j} / max(imToProcess(:));
        borderIm{i,j} = 1 - borderIm{i,j};
        borderIm{i,j}(waveletWatershed{i,j} == 0) = 1;
        activeContSeparate{i,j} = activecontour(borderIm{i,j}, wavelet,1000);
        segments{i,j} = (activeContSeparate{i,j}>0) & (activeCont{i,j}>0);
    end
end

PSMA.activeCont = cutBorders(activeCont,waveletBorder);
PSMA.waveletWatershed = cutBorders(waveletWatershed,waveletBorder);
PSMA.borderIm = cutBorders(borderIm,waveletBorder);
PSMA.activeContSeparate = cutBorders(activeContSeparate,waveletBorder);
PSMA.segments = cutBorders(segments,waveletBorder);
se = [0 1 0; 1 1 1; 0 1 0];
tmp = imfilter(double(PSMA.segments),se);
tmp = (PSMA.segments & tmp>2);
PSMA.segments = tmp;


keyboard;

% % % CD45.activeCont = activecontour(CD45.im,CD45.wavelet);
% % % CD45.waveletWatershed = watershed(1-(CD45.wavelet>0));
% % % CD45.borderIm = CD45.im - min(CD45.im(:));
% % % CD45.borderIm = CD45.borderIm / max(CD45.im(:));
% % % CD45.borderIm = 1 - CD45.borderIm;
% % % CD45.borderIm(CD45.waveletWatershed == 0) = 1;
% % % CD45.activeContSeparate = activecontour(CD45.borderIm, CD45.wavelet);
% % % CD45.segments = (CD45.activeContSeparate>0) & (CD45.activeCont>0);

save('167_all_data','-v7.3');

% load('107_all_data');

CD45.filter = generateFilter([40,40],5,3);
CD45.filtered = conv2(double(CD45.im),CD45.filter,'same');
CD45.filteredHist = hist(CD45.filtered(:), min(min(CD45.filtered)):histJump:max(max(CD45.filtered)));
[~, preThre, ~, ~] = cutFirstHistMode(CD45.filtered(3000:3400,3000:3400),0);
[CD45.thre,~] = thresh_tool_v2(CD45.filtered(3000:3400,3000:3400), CD45.origIm(3000:3400,3000:3400),[],preThre);
CD45.binary = (CD45.filtered > CD45.thre);
% % % CD45.segments = CD45.segments | CD45.binary;
CD45.segments = CD45.binary;



PSMA.activeCont = [];
PSMA.waveletWatershed = [];
PSMA.borderIm = [];
PSMA.activeContSeparate = [];
PSMA.brightness = [];

DAPI.activeCont = [];
DAPI.waveletWatershed = [];
DAPI.borderIm = [];
DAPI.activeContSeparate = [];
DAPI.brightness = [];

CD45.activeCont = [];
CD45.waveletWatershed = [];
CD45.borderIm = [];
CD45.activeContSeparate = [];
CD45.brightness = [];

save('167_all_data_short','-v7.3');
%load('107_all_data_short');

h1 = figure;
AX(1) = subplot(3,3,1);
myimshow(log(double(PSMA.origIm)),[],[0.005,0.005]);
%imshow(log(double(PSMA.origIm)),[]);
hold on
title('PSMA original Input');
for i = 1:length(x)
    rectangle('Position',[x(i)-10, y(i) - 10, 21,21],'EdgeColor','g');
end


AX(2) = subplot(3,3,2);
%imshow(log(double(DAPI.im)),[]);
myimshow(DAPI.origIm,[],[0.005,0.001]);
for i = 1:length(x)
    rectangle('Position',[x(i)-10, y(i) - 10, 21,21],'EdgeColor','g');
end
title('DAPI original Input');

AX(3) = subplot(3,3,3);
myimshow(log(double(CD45.origIm)),[],[0.01 0.001]);
for i = 1:length(x)
    rectangle('Position',[x(i)-10, y(i) - 10, 21,21],'EdgeColor','g');
end
title('CD45 original Input');



PSMA.reginoprops = regionprops(logical(PSMA.segments),'Area','BoundingBox','PixelIdxList','MajorAxisLength','MinorAxisLength','PixelList','PixelIdxList','Eccentricity');
for i = 1:length(PSMA.reginoprops)
    if PSMA.reginoprops(i).Area > PSMA.maxArea
        PSMA.segments(PSMA.reginoprops(i).PixelIdxList) = 0;
    elseif PSMA.reginoprops(i).Area < PSMA.minArea
        PSMA.segments(PSMA.reginoprops(i).PixelIdxList) = 0;
    elseif PSMA.reginoprops(i).MajorAxisLength > PSMA.maxAxisLength
        PSMA.segments(PSMA.reginoprops(i).PixelIdxList) = 0;
    elseif PSMA.reginoprops(i).MinorAxisLength < PSMA.minAxisLength
        PSMA.segments(PSMA.reginoprops(i).PixelIdxList) = 0;
    elseif PSMA.reginoprops(i).MinorAxisLength*2<PSMA.reginoprops(i).MajorAxisLength
        PSMA.segments(PSMA.reginoprops(i).PixelIdxList) = 0;
    elseif PSMA.reginoprops(i).Eccentricity > 0.4
        PSMA.binary(PSMA.reginoprops(i).PixelIdxList) = 0;
    end
end
AX(4) = subplot(3,3,4);
imshow(PSMA.segments,[]);
hold on
for i = 1:length(x)
    rectangle('Position',[x(i)-10, y(i) - 10, 21,21],'EdgeColor','g');
end
title('PSMA Processed');
se = strel('disk',2);
[yInd, xInd] = find(imerode(PSMA.segments == 1,se,'same')==1);
KeepInd = (((xInd > merge) & (xInd < size(PSMA.segments,1)-merge)) & ((yInd > merge) & (yInd < size(PSMA.segments,2)-merge)));
xInd = xInd(KeepInd);
yInd = yInd(KeepInd);
K_PSMA = convhull(xInd, yInd);
hold on;
plot(xInd(K_PSMA), yInd(K_PSMA), '-b', 'LineWidth',2);
hold off;
BW_tmp = poly2mask(xInd(K_PSMA), yInd(K_PSMA), size(PSMA.segments,1),size(PSMA.segments,2));
pointsMask = BW_tmp;


DAPI.reginoprops = regionprops(DAPI.segments,'Area','BoundingBox','PixelIdxList','MajorAxisLength','MinorAxisLength','PixelList','PixelIdxList','Eccentricity');
for i = 1:length(DAPI.reginoprops)
    if DAPI.reginoprops(i).Area > DAPI.maxArea
        DAPI.segments(DAPI.reginoprops(i).PixelIdxList) = 0;
    elseif DAPI.reginoprops(i).Area < DAPI.minArea
        DAPI.segments(DAPI.reginoprops(i).PixelIdxList) = 0;
    elseif DAPI.reginoprops(i).MajorAxisLength > DAPI.maxAxisLength
        DAPI.segments(DAPI.reginoprops(i).PixelIdxList) = 0;
    elseif DAPI.reginoprops(i).MinorAxisLength < DAPI.minAxisLength
        DAPI.segments(DAPI.reginoprops(i).PixelIdxList) = 0;
    elseif DAPI.reginoprops(i).Eccentricity > 0.4
        DAPI.binary(DAPI.reginoprops(i).PixelIdxList) = 0;
    end
end

AX(5) = subplot(3,3,5);
imshow(DAPI.segments,[]);
for i = 1:length(x)
    rectangle('Position',[x(i)-10, y(i) - 10, 21,21],'EdgeColor','g');
end
title('DAPI Processed');
se = strel('disk',2);
[yInd, xInd] = find(imerode(DAPI.segments,se,'same') == 1);
KeepInd = (((xInd > merge) & (xInd < size(DAPI.segments,1)-merge)) & ((yInd > merge) & (yInd < size(DAPI.segments,2)-merge)));
xInd = xInd(KeepInd);
yInd = yInd(KeepInd);
K_DAPI = convhull(xInd, yInd);
hold on;
plot(xInd(K_DAPI), yInd(K_DAPI), '-b', 'LineWidth',2);
hold off;
BW_tmp = poly2mask(xInd(K_DAPI), yInd(K_DAPI), size(DAPI.segments,1),size(DAPI.segments,2));
pointsMask = pointsMask .* BW_tmp;


AX(6) = subplot(3,3,6);
imshow(CD45.segments,[]);
for i = 1:length(x)
    rectangle('Position',[x(i)-10, y(i) - 10, 21,21],'EdgeColor','g');
end
title('CD45 Processed');

se = strel('disk',2);
tmp_CD45.segments = imerode(CD45.segments,se,'same');

for test = 1:10
    [yInd, xInd] = find(tmp_CD45.segments == 1);
    KeepInd = (((xInd > merge) & (xInd < size(CD45.segments,1)-merge)) & ((yInd > merge) & (yInd < size(CD45.segments,2)-merge)));
    xInd = xInd(KeepInd);
    yInd = yInd(KeepInd);
    K_CD45 = convhull(xInd, yInd);
    L = bwlabel(tmp_CD45.segments);
    for i = 1:length(K_CD45)
        tmp_CD45.segments(L == L(yInd(K_CD45(i)),xInd(K_CD45(i)))) = 0;
    end
end
[yInd, xInd] = find(tmp_CD45.segments == 1);
KeepInd = (((xInd > merge) & (xInd < size(CD45.segments,1)-merge)) & ((yInd > merge) & (yInd < size(CD45.segments,2)-merge)));
xInd = xInd(KeepInd);
yInd = yInd(KeepInd);
K_CD45 = convhull(xInd, yInd);

hold on;
plot(xInd(K_CD45), yInd(K_CD45), '-b', 'LineWidth',2);
hold off;
BW_tmp = poly2mask(xInd(K_CD45), yInd(K_CD45), size(CD45.segments,1),size(CD45.segments,2));
pointsMask = pointsMask .* BW_tmp;


PSMA.segments = PSMA.segments .* pointsMask;
DAPI.segments = DAPI.segments .* pointsMask;
CD45.segments = CD45.segments .* pointsMask;


AX(7) = subplot(3,3,7);
bDA_PS = DAPI.segments & PSMA.segments;
imshow(bDA_PS,[]);
for i = 1:length(x)
    rectangle('Position',[x(i)-10, y(i) - 10, 21,21],'EdgeColor','g');
end
title('Positive PSMA & DAPI');

out = (DAPI.segments) & (PSMA.segments) & (1-CD45.segments);
AX(8) = subplot(3,3,8);
imshow(out,[]);
for i = 1:length(x)
    rectangle('Position',[x(i)-10, y(i) - 10, 21,21],'EdgeColor','g');
end
title('Positive PSMA & DAPI ,Negative CD45');

CC = regionprops(bDA_PS,'Area','BoundingBox','Centroid','PixelIdxList','MajorAxisLength','MinorAxisLength');
L = bwlabel(bDA_PS);
%hL = hist(L(:), 0:max(L(:)));
LandCD45 = L .* (1-CD45.segments);
hLandCD45 = hist(LandCD45(:),0:max(L(:)));
%hL = hL(2:end);
hLandCD45 = hLandCD45(2:end);
detections = (abs([CC.Area] - hLandCD45)<0.3*[CC.Area]);
detections([CC.Area]<20) = 0;
%detections([CC.MajorAxisLength]>1.5*[CC.MinorAxisLength]) = 0;
tmp = reshape([CC(detections==1).Centroid],2,sum(detections));
[PSMA_Reg, PSMA_indReg] = bwselect(PSMA.segments,tmp(1,:),tmp(2,:));
PSMA_Reg_Prop = regionprops(PSMA_Reg,'Area','BoundingBox','Centroid','PixelIdxList','MajorAxisLength','MinorAxisLength');
%PSMA_Reg_Prop([PSMA_Reg_Prop(:).MajorAxisLength]>1.5*[PSMA_Reg_Prop(:).MinorAxisLength]) = [];
tmp = reshape([PSMA_Reg_Prop(:).Centroid],2,length(PSMA_Reg_Prop));


AX(9) = subplot(3,3,9);
imshow(out,[]);
hold on;
for i = 1:length(x)
    rectangle('Position',[x(i)-10, y(i) - 10, 21,21],'EdgeColor','g');
end
hold on;
%keyboard;
detectionArea = [PSMA_Reg_Prop(:).Area];
detectionPixelIdxList = {PSMA_Reg_Prop(:).PixelIdxList};
detectionMajorAxis = [PSMA_Reg_Prop(:).MajorAxisLength];
detectionMinorAxis = [PSMA_Reg_Prop(:).MinorAxisLength];

for i = 1:length(detectionPixelIdxList)
    detectionMax(i) = max(PSMA.origIm(detectionPixelIdxList{i}));
    detectionMedian(i) = median(PSMA.origIm(detectionPixelIdxList{i}));
    detectionmin(i) = min(PSMA.origIm(detectionPixelIdxList{i}));
end
plot(tmp(1,:),tmp(2,:),'Marker','*','Color','r','lineStyle','none');

title('Detected CTCs (Red Marks)');

linkaxes(AX,'xy');
keyboard;

AX(10) =figure;
imshow(out,[]);
hold on;
for i = 1:length(x)
    rectangle('Position',[x(i)-10, y(i) - 10, 21,21],'EdgeColor','g');
end
hold on;
plot(tmp(1,:),tmp(2,:),'Marker','*','Color','r','lineStyle','none');

outputExcel = cell(size(tmp,2)+1, 5);
outputExcel{1,1} = 'x in Pixels';
outputExcel{1,2} = 'y in Pixels';
outputExcel{1,4} = 'x in Micro Meter';
outputExcel{1,5} = 'y in Micro Meter';

outputExcel(2:1+size(tmp,2),1:2) = mat2cell(tmp',ones(1,size(tmp,2)),ones(1,size(tmp,1)));
outputExcel(2:1+size(tmp,2),4:5) = mat2cell(tmp'*eachPixelLengthInMicro,ones(1,size(tmp,2)),ones(1,size(tmp,1)));

uiwait(msgbox('Press OK if you are ready to save the data!'));

prompt = {'Enter excel file name (Leave empty if you do not need it)','Enter experiment name to save results (Leave empty if you do not need it)'};
dlg_title = 'Save Results (Cancel if you do not want anything to be saved)';
num_lines = 1;
def = {sprintf('DetectionCoordinates_P%0.3f_D%0.3f_C%0.3f.xlsx',PSMA.thre,DAPI.thre,CD45.thre),sprintf('Experiment_Name_P%0.3f_D%0.3f_C%0.3f',PSMA.thre,DAPI.thre,CD45.thre)};
answer = inputdlg(prompt,dlg_title, num_lines, def);

if isempty(answer)
    return;
end
if ~isempty(answer{1})
    [pathstr,name,ext] = fileparts(answer{1});
    if isempty(ext)
        answer{1} = [answer{1} '.xlsx'];
    end
    if isempty(pathstr)
        answer{1} = fullfile(PathName,'Results', answer{1});
    end
    [pathstr,name,ext] = fileparts(answer{1});
    if ~exist(pathstr,'dir')
        mkdir(pathstr);
    end
    xlswrite(answer{1}, outputExcel);
end

if ~isempty(answer{2})
    
    [pathstr,name,ext] = fileparts(answer{2});
    if isempty(pathstr)
        answer{2} = fullfile(PathName,'Results', answer{2});
    end
    
    if ~exist(fullfile(answer{2},'Annotations__1'),'dir')
        mkdir(fullfile(answer{2},'Annotations__1'));
    end
    title('Detections2');
    for i = 1:length(x)
        figure(h1);
        subplot(3,3,9);
        xlim([x(i)-70,x(i)+70]);
        ylim([y(i)-70,y(i)+70]);
        title(sprintf('x = %d, y = %d', round(x(i)),round(y(i))));
        saveas(h1,fullfile(answer{2},'Annotations__1',sprintf('%d',i)),'tif');
    end
    
    if ~exist(fullfile(answer{2},'Detections__1'),'dir')
        mkdir(fullfile(answer{2},'Detections__1'));
    end
    for i = 1:size(tmp,2)
        figure(h1);
        subplot(3,3,9);
        xlim([tmp(1,i)-70,tmp(1,i)+70]);
        ylim([tmp(2,i)-70,tmp(2,i)+70]);
        title(sprintf('x = %d, y = %d', round(tmp(1,i)),round(tmp(2,i))));
        subplot(3,3,1);
        xlabel(sprintf('Area=%i,\nIntensity: Max=%i,median=%i,min=%i',detectionArea(i),detectionMax(i),detectionMedian(i),detectionmin(i)));
        subplot(3,3,4);
        xlabel(sprintf('MajorAxis = %2.2f ,minorAxis = %2.2f',detectionMajorAxis(i),detectionMinorAxis(i)));
        saveas(h1,fullfile(answer{2},'Detections__1',sprintf('%d',i)),'tif');
    end
    
    x = round(ceil(x));
    y = round(ceil(y));
    
    if ~exist(fullfile(answer{2},'Annotations__2'),'dir')
        mkdir(fullfile(answer{2},'Annotations__2'));
    end
    for i = 1:length(x)
        figure(h1);
        subplot(3,3,1);
        hold off;
        m = min(min(PSMA.origIm((y(i)-70):(y(i)+70),(x(i)-70):(x(i)+70))));
        M = max(max(PSMA.origIm((y(i)-70):(y(i)+70),(x(i)-70):(x(i)+70))));
        imshow(PSMA.origIm,[m,M]);
        title('PSMA High Contrast (not compareable with other patches)');
        subplot(3,3,9);
        xlim([x(i)-70,x(i)+70]);
        ylim([y(i)-70,y(i)+70]);
        title(sprintf('x = %d, y = %d', round(x(i)),round(y(i))));
        saveas(h1,fullfile(answer{2},'Annotations__2',sprintf('%d',i)),'tif');
    end
    
    
    if ~exist(fullfile(answer{2},'Detections__2'),'dir')
        mkdir(fullfile(answer{2},'Detections__2'));
    end
    tmp = round(tmp);
    for i = 1:size(tmp,2)
        figure(h1);
        subplot(3,3,1);
        hold off;
        m = min(min(PSMA.origIm((tmp(2,i)-70):(tmp(2,i)+70),(tmp(1,i)-70):(tmp(1,i)+70))));
        M = max(max(PSMA.origIm((tmp(2,i)-70):(tmp(2,i)+70),(tmp(1,i)-70):(tmp(1,i)+70))));
        imshow(PSMA.origIm,[m,M]);
        title('PSMA High Contrast (not compareable with other patches)');
        subplot(3,3,1);
        xlabel(sprintf('Area=%i,\nIntensity: Max=%i,median=%i,min=%i',detectionArea(i),detectionMax(i),detectionMedian(i),detectionmin(i)));
        
        subplot(3,3,9);
        xlim([tmp(1,i)-70,tmp(1,i)+70]);
        ylim([tmp(2,i)-70,tmp(2,i)+70]);
        title(sprintf('x = %d, y = %d', round(tmp(1,i)),round(tmp(2,i))));
        subplot(3,3,4);
        xlabel(sprintf('MajorAxis = %2.2f ,minorAxis = %2.2f',detectionMajorAxis(i),detectionMinorAxis(i)));
        saveas(h1,fullfile(answer{2},'Detections__2',sprintf('%d',i)),'tif');
    end
end

if ~isempty(answer{1}) || ~isempty(answer{2})
    uiwait(msgbox('Save Procedure Complete! Press OK to see the results'));
    winopen(fullfile(PathName,'Results'));
end

return;

bound = 100;
for i = 1:length(x)
    Yboundry = (x(i)-bound):(x(i)+bound);
    Xboundry = (y(i)-bound):(y(i)+bound);
    figure;
    subplot(3,3,1);
    imshow(PSMA.im(Xboundry,Yboundry),[]);
    hold on
    rectangle('Position',[bound-10, bound - 10, 21,21],'EdgeColor','g');
    title('PSMA original Input');
    
    subplot(3,3,2);
    imshow(DAPI.im(Xboundry,Yboundry),[]);
    rectangle('Position',[bound-10, bound - 10, 21,21],'EdgeColor','g');
    title('DAPI original Input');
    
    subplot(3,3,3);
    imshow(CD45.im(Xboundry,Yboundry),[]);
    rectangle('Position',[bound-10, bound - 10, 21,21],'EdgeColor','g');
    title('CD45 original Input');
    
    subplot(3,3,4);
    imshow(PSMA.segments(Xboundry,Yboundry),[]);
    hold on
    rectangle('Position',[bound-10, bound - 10, 21,21],'EdgeColor','g');
    title('PSMA Processed');
    
    subplot(3,3,5);
    imshow(DAPI.segments(Xboundry,Yboundry),[]);
    rectangle('Position',[bound-10, bound - 10, 21,21],'EdgeColor','g');
    title('DAPI Processed');
    
    subplot(3,3,6);
    imshow(CD45.segments(Xboundry,Yboundry),[]);
    rectangle('Position',[bound-10, bound - 10, 21,21],'EdgeColor','g');
    title('CD45 Processed');
    
    subplot(3,3,8);
    imshow(out(Xboundry,Yboundry),[]);
    rectangle('Position',[bound-10, bound - 10, 21,21],'EdgeColor','g');
    title('Detections');
    
end

%end