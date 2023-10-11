frame = double(imread('testImage.tif'));
[detectionResults, detectionMask] = spotDetector(frame);

figure; 
h(1) = subplot(1,2,1); imagesc(frame); colormap(gray(256)); myimshow(frame,[],0.01); axis image; title('Input');
h(2) = subplot(1,2,2); imagesc(detectionMask); axis image; title('Detection');
linkaxes(h,'xy');
