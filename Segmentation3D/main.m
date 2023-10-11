% Center for Research in Computer Vision (Mubarak Shah, PhD)
%% Cleanup
clear all
close all
clc

%% Input images folders
c3.imDir = 'E:\MATLAB\CTC\Alex Dropbox\63xNewFilter\3D\9740-103-02-1822-5846';
c3.imageNames = dir(fullfile(c3.imDir,'*c3*.tif'));
c4.imDir = c3.imDir;
c4.imageNames = dir(fullfile(c4.imDir,'*c4*.tif'));

%% Debug == 1 if you are debugging the code, Debugg == 0 if you need to just use the code
debug = 0;

%% Cut First Mode Implimentation
addpath('cutFirstHistMode\');

%% Histogram bin steps
octave = 0.005;

%% Setting up the sample name for being used in results figure
UnderLineSpot = strfind(c3.imDir,'_');
sampleName = c3.imDir;
sampleName(UnderLineSpot) = '-';
if ~isempty(UnderLineSpot)
    sampleName = sampleName(UnderLineSpot(1)+1:end);
end
backSlashSpot = strfind(c3.imDir, '\');
if ~isempty(backSlashSpot)
    sampleName = sampleName(backSlashSpot(end)+1:end);
end
%% Reading input images
for i = 1:length(c3.imageNames)
    c3.im(:,:,i) = imread(fullfile(c3.imDir, c3.imageNames(i).name));
    c4.im(:,:,i) = imread(fullfile(c4.imDir, c4.imageNames(i).name));
end

%% Normalizing the input images
c3.im = double(c3.im);
c4.im = double(c4.im);

c3.im = c3.im - min(min(min(c3.im)));
c3.im = c3.im / max(max(max(c3.im)));

c4.im = c4.im - min(min(min(c4.im)));
c4.im = c4.im / max(max(max(c4.im)));

%% Mean of each pixel in 3D
c3.mean = mean(c3.im,3);
c4.mean = mean(c4.im,3);

%% User selects the cell area manually
h3 = figure;
imshow(c3.mean,[]);
r3 = getrect(h3);
rectangle('Position',r3,'LineWidth',2,'EdgeColor','r');
pause(0.2);

%% Cropping input images
c3.cropped = c3.im(round(r3(2):r3(2)+r3(4)),round(r3(1):r3(1)+r3(3)),:);
c4.cropped = c4.im(round(r3(2):r3(2)+r3(4)),round(r3(1):r3(1)+r3(3)),:);

ann = AnnotationTool(c3.cropped);
ann = AnnotationTool(c4.cropped);

%% Calculating the histograms
c3.imHist = hist(c3.im(:),[0:octave:1]);
c3.imHistNorm = c3.imHist / sum(c3.imHist);
c4.imHist = hist(c4.im(:),[0:octave:1]);
c4.imHistNorm = c4.imHist / sum(c4.imHist);

c3.croppedHist = hist(c3.cropped(:),[0:octave:1]);
c3.croppedHistNorm = c3.croppedHist / sum(c3.croppedHist);
c4.croppedHist = hist(c4.cropped(:),[0:octave:1]);
c4.croppedHistNorm = c4.croppedHist / sum(c4.croppedHist);

% [cutoffIndex, c3.cutoffValue, sp, axesH] = cutFirstHistMode(c3.im,0);
% [cutoffIndex, c4.cutoffValue, sp, axesH] = cutFirstHistMode(c4.im,0);


%% Setting the probability map
c3.probMap = c3.croppedHistNorm ./ max(c3.imHistNorm,eps);
c4.probMap = c4.croppedHistNorm ./ max(c4.imHistNorm,eps);

%% Showing the histograms and the probability map
if debug
    figure
    plot(c3.imHistNorm); hold on; plot(c3.croppedHistNorm,'r');
    figure
    plot(c3.probMap);
    
    figure
    plot(c4.imHistNorm); hold on; plot(c4.croppedHistNorm,'r');
    figure
    plot(c4.probMap);
end


%% Segmentation using the Histogram method
c3.croppedSegment_ours = zeros(size(c3.cropped));
c3.croppedSegment_ours(c3.probMap(floor(c3.cropped / octave)+1)>1) = 1;
se = strel('disk',2);
for i = 1:size(c3.croppedSegment_ours,3)
    c3.croppedSegment_ours(:,:,i) = imopen(c3.croppedSegment_ours(:,:,i),se);
end
c3.segment_ours = zeros(size(c3.im));
c3.connComp = regionprops(bwlabeln(c3.croppedSegment_ours,6),'Area','PixelIdxList');
c3.croppedSegment_ours(:,:,:) = 0;
areas = cat(1,c3.connComp.Area);
[a b] = max(areas);
c3.croppedSegment_ours(c3.connComp(b).PixelIdxList) = 1;
c3.segment_ours(round(r3(2):r3(2)+r3(4)),round(r3(1):r3(1)+r3(3)),:) = c3.croppedSegment_ours;

%% Segmantation using Otsu's method
c3.level = graythresh(c3.cropped);
c3.croppedSegment_t1 = zeros(size(c3.cropped));
c3.croppedSegment_t1(c3.cropped>c3.level) = 1;
c3.segment_t1 = zeros(size(c3.im));
c3.segment_t1(round(r3(2):r3(2)+r3(4)),round(r3(1):r3(1)+r3(3)),:) = c3.croppedSegment_t1;

%% Segmentation using the first cut mode
[cutoffIndex, c3.level2, sp, axesH] = cutFirstHistMode(c3.im,0);
c3.croppedSegment_t2 = zeros(size(c3.cropped));
c3.croppedSegment_t2(c3.cropped>c3.level2) = 1;
c3.segment_t2 = zeros(size(c3.im));
c3.segment_t2(round(r3(2):r3(2)+r3(4)),round(r3(1):r3(1)+r3(3)),:) = c3.croppedSegment_t2;


%% Segmentation using the Histogram method
c4.croppedSegment_ours = zeros(size(c4.cropped));
c4.croppedSegment_ours(c4.probMap(floor(c4.cropped / octave)+1)>1) = 1;
c4.segment_ours = zeros(size(c4.im));
c4.connComp = regionprops(bwlabeln(c4.croppedSegment_ours,6),'Area','PixelIdxList');
areas = cat(1,c4.connComp.Area);
[a b] = max(areas);
c4.croppedSegment_ours(:,:,:) = 0;
c4.croppedSegment_ours(c4.connComp(b).PixelIdxList) = 1;
c4.segment_ours(round(r3(2):r3(2)+r3(4)),round(r3(1):r3(1)+r3(3)),:) = c4.croppedSegment_ours;

%% Segmantation using Otsu's method
c4.level = graythresh(c4.cropped);
c4.croppedSegment_t1 = zeros(size(c4.cropped));
c4.croppedSegment_t1(c4.cropped>c4.level) = 1;
c4.segment_t1 = zeros(size(c4.im));
c4.segment_t1(round(r3(2):r3(2)+r3(4)),round(r3(1):r3(1)+r3(3)),:) = c4.croppedSegment_t1;

%% Segmentation using the first cut mode
[cutoffIndex, c4.level2, sp, axesH] = cutFirstHistMode(c4.cropped,0);
c4.croppedSegment_t2 = zeros(size(c4.cropped));
c4.croppedSegment_t2(c4.cropped>c4.level2) = 1;
c4.croppedSegment_t2 = imfill(c4.croppedSegment_t2,'holes');
c4.segment_t2 = zeros(size(c4.im));
c4.segment_t2(round(r3(2):r3(2)+r3(4)),round(r3(1):r3(1)+r3(3)),:) = c4.croppedSegment_t2;
c4.connComp = regionprops(bwlabeln(c4.segment_ours,6),'Area','PixelIdxList');

%% Calculating the sum of intensity values
c4.sum_ours = sum(sum(sum(c3.im(c4.segment_ours>0))))/100;
c4.sum_t1 = sum(sum(sum(c3.im(c4.segment_t1>0))))/100;
c4.sum_t2 = sum(sum(sum(c3.im(c4.segment_t2>0))))/100;

c3.sum_ours = sum(sum(sum(c3.im(c3.segment_ours>0))))/100;
c3.sum_t1 = sum(sum(sum(c3.im(c3.segment_t1>0))))/100;
c3.sum_t2 = sum(sum(sum(c3.im(c3.segment_t2>0))))/100;

pause(2);
close all
h = figure('Position',[300,300,800,800]);
for i = 1:size(c3.im,3)
    subplot(2,2,1);
    imagesc(c3.im(:,:,i),[0,1]);
    ylabel('AR');
    title(sprintf('input %s', sampleName));
    subplot(2,2,2);
    imagesc(c3.segment_ours(:,:,i),[0,1]);
    title('Mask of Cellular Area');
    xlabel(sprintf('%3.0f',c3.sum_ours));
    subplot(2,2,3);
    imagesc(c4.im(:,:,i),[0,1]);
    ylabel('DAPI');
    subplot(2,2,4);
    imagesc(c4.segment_t2(:,:,i),[0,1]);
    xlabel(sprintf('%3.0f',c4.sum_t2));
    title('Mask of Nucleus Area');
    pause(0.3);
    %print('-djpeg100',sprintf('%.3i.jpg',i));
end
disp(c4.sum_t2/c3.sum_ours);

keyboard;

if debug
h = figure;
for i = 1:size(c3.im,3)
    subplot(2,4,1);
    imagesc(c3.im(:,:,i),[0,1]);
    ylabel('AR');
    title('input');
    subplot(2,4,2);
    imagesc(c3.segment_ours(:,:,i),[0,1]);
    title('ours');
    xlabel(sprintf('%3.0f',c3.sum_ours));
    subplot(2,4,3);
    imagesc(c3.segment_t1(:,:,i),[0,1]);
    title('Otsu');
    xlabel(sprintf('%3.0f',c3.sum_t1));
    subplot(2,4,4);
    imagesc(c3.segment_t2(:,:,i),[0,1]);
    title('Cut mode');
    xlabel(sprintf('%3.0f',c3.sum_t2));
    subplot(2,4,5);
    imagesc(c4.im(:,:,i),[0,1]);
    ylabel('DAPI');
    subplot(2,4,6);
    imagesc(c4.segment_ours(:,:,i),[0,1]);
    xlabel(sprintf('%3.0f',c4.sum_ours));
    subplot(2,4,7);
    imagesc(c4.segment_t1(:,:,i),[0,1]);
    xlabel(sprintf('%3.0f',c4.sum_t1));
    subplot(2,4,8);
    imagesc(c4.segment_t2(:,:,i),[0,1]);
    xlabel(sprintf('%3.0f',c4.sum_t2));
    print('-djpeg100',sprintf('%.3i.jpg',i));
end
end
