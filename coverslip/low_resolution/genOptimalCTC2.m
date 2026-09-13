function genOptimalCTC2

% Generate figures en masse for the CTC detection algorithm
% Analyzes the performance of the CTC detection algorithm
% Nik Mihaylov (July, 2014) - Coverslip analysis

[keys,sections,subsections] = inifile('datasets.ini','readall');

function ret = getParams(section, subsection)
  ret = {};
  pos = 1;
  for ii = 1:size(keys,1)
    if strcmp(keys{ii,1}, section) && strcmp(keys{ii,2}, subsection)
      ret{pos,1} = keys{ii,3};
      ret{pos,2} = keys{ii,4};
      pos = pos + 1;
    end
  end
  ks = ret(:,1);
  values = ret(:,2);
  ret = cell2struct(values, ks, 1);
end

globalParams = getParams('','');

databasedir = globalParams.databasedir;
xlsbasedir = globalParams.xlsbasedir;
resbasedir = globalParams.resbasedir;
mkdir(resbasedir);

grandTotals = zeros(1, 4);

tic


for dataset = 1:1 %size(sections,1)

params = getParams(sections{dataset}, '');


% initialize the parameters
labelCKimg = imread(params.labelckimg);
labelDAPIimg = imread(params.labeldapiimg);
labelCD45img = imread(params.labelcd45img);
name = params.name;
datadir = fullfile(databasedir, params.datadir);
coordinates = xlsread(fullfile(xlsbasedir, params.coordinates));
resdir = fullfile(resbasedir, params.resdir);
nameCK = params.nameck;
nameDAPI = params.namedapi;
nameCD45 = params.namecd45;
positions = [1 : str2num(params.positions)];
coeffCD45 = str2num(params.coeffcd45);
coeffDAPIarea = str2num(params.coeffdapiarea);
coeffCKarea = str2num(params.coeffckarea);
maxDist = str2num(params.maxdist);
rectSize = str2num(params.rectsize);
minCK = str2num(params.minck);
maxCK = str2num(params.maxck);
minDAPI = str2num(params.mindapi);
maxDAPI = str2num(params.maxdapi);
minCD45 = str2num(params.mincd45);
maxCD45 = str2num(params.maxcd45);
maxCD45pixels = str2num(params.maxcd45pixels);


mkdir(resdir);

imgdir = fullfile(resdir, 'images');
mkdir(imgdir);

out = fopen(fullfile(resdir, 'parameters.txt'), 'w');
fprintf(out, 'Data set: %s\n\n', name);

%[status, version] = system('svnversion');
%if status == 0
%    fprintf(out, 'Code Revision: %s\n', version);
%end

fprintf(out, 'Data dir = %s\n\n', datadir);

% Detection algorithm parameters

fprintf(out, 'CD45 cut off coefficient = %d\n', coeffCD45);
fprintf(out, 'DAPI area limit = %d\n', coeffDAPIarea);
fprintf(out, 'CK area limit = %d\n', coeffCKarea);
fprintf(out, 'Max Dist = %f\n', maxDist);
fprintf(out, 'Max CD45 Pixels = %d\n\n', maxCD45pixels);

fprintf(out, 'CK intensity stretch = [%d, %d]\n', minCK, maxCK);
fprintf(out, 'DAPI intensity stretch = [%d, %d]\n', minDAPI, maxDAPI);
fprintf(out, 'CD45 intensity stretch = [%d, %d]\n', minCD45, maxCD45);

% Run the detection algorithm for all (position, z) combinations

res = zeros(0, 3);

totalTP = 0;
totalFP = 0;
totalTN = 0;
totalFN = 0;

maxDAPIarea = 0;
maxCKarea = 0;

allstats = [];

for pos = positions
    manCoords = coordinates(coordinates(:, 1) == pos, 2:3);
    numManualCTCs = size(manCoords, 1);
    %if numManualCTCs == 0
    %    continue
    %end

    try
        imgCK   = imread(fullfile(datadir, sprintf(nameCK, pos)), 'tiff');
        imgDAPI = imread(fullfile(datadir, sprintf(nameDAPI, pos)), 'tiff');
        imgCD45 = imread(fullfile(datadir, sprintf(nameCD45, pos)), 'tiff');
    catch err
        fprintf('Error: at pos = %d: %s\n', pos, err.message);
        % We get here when some file doesn't exist; then just go on with the next one
        continue;
    end

    % Segment the channels
    [~, spotsDAPI] = spotDetector(double(imgDAPI));
    [~, spotsCK] = spotDetector(double(imgCK));
    % Segment the CD45 channel
    dimgCD45 = double(imgCD45);
    [~, cutCD45] = cutFirstHistMode(dimgCD45, 0);
    maskCD45 = dimgCD45 > cutCD45 * coeffCD45;


    [numCTCs, ctcs, stats] = optimalDetectCTC2(manCoords, spotsDAPI, spotsCK, maskCD45, coeffCD45, maxCD45pixels, coeffDAPIarea, coeffCKarea, maxDist);
    ctcs = round(ctcs);

    for j = 1 : numCTCs
        res = [res; [pos, int32(ctcs(j, 1)), int32(ctcs(j, 2))]];
    end
    [tp, fp, tn, fn] = analyzeCTC2(manCoords, ctcs, maxDist);
    totalTP = totalTP + tp;
    totalFP = totalFP + fp;
    totalTN = totalTN + tn;
    totalFN = totalFN + fn;

    % make the figures
    figCK = mkFigure(imgCK, minCK, maxCK, manCoords, ctcs, maxDist, rectSize);
    figDAPI = mkFigure(imgDAPI, minDAPI, maxDAPI, manCoords, ctcs, maxDist, rectSize);
    figCD45 = mkFigure(dimgCD45, minCD45, maxCD45, manCoords, ctcs, maxDist, rectSize);

    gap = uint8(255 * ones(size(imgCK, 1) + size(labelCKimg,1), 1, 3));
    fig = cat(2, cat(1, labelCKimg, figCK), gap, cat(1, labelDAPIimg, figDAPI), gap, cat(1, labelCD45img, figCD45));

    figCKspots = mkFigure(spotsCK, minCK, maxCK, manCoords, ctcs, maxDist, rectSize);
    figDAPIspots = mkFigure(spotsDAPI, minDAPI, maxDAPI, manCoords, ctcs, maxDist, rectSize);
    figCD45spots = mkFigure(maskCD45 .* dimgCD45, minCD45, maxCD45, manCoords, ctcs, maxDist, rectSize);

    gap = uint8(255 * ones(size(imgCK, 1), 1, 3));
    figSpots = cat(2, figCKspots, gap, figDAPIspots, gap, figCD45spots);

    gap = uint8(255 * ones(1, size(figSpots,2), 3));
    fig = cat(1, fig, gap, figSpots);

    imwrite(fig, fullfile(imgdir, sprintf('p%03d_1.png', pos)), 'png');

    maxDAPIarea = max(maxDAPIarea, stats.DAPIarea(3));
    maxCKarea = max(maxCKarea, stats.CKarea(3));
    stats.pos = pos;
    stats.numCTCs = numCTCs;
    stats.numManualCTCs = numManualCTCs;
    allstats = [allstats; stats];
end % pos

fprintf(out, '\nTotal:\nTrue Positives:  %d\nFalse Positives: %d\nTrue Negatives:  %d\nFalse Negatives: %d\n', totalTP, totalFP, totalTN, totalFN);

fprintf(out, 'Max DAPI area found = %d\nMax CK area found = %d\n', maxDAPIarea, maxCKarea);
fclose(out);

% save the results
csvwrite(fullfile(resdir, 'ctc.csv'), res);

out = fopen(fullfile(resdir, 'thresholds.csv'), 'w');
% print the header line
fprintf(out, 'Position,Num_Manual_CTCs,Num_CTCs');
fprintf(out, ',DAPI_Area_Low,DAPI_Area_High,DAPI_Area_Max,DAPI_Area_Uni_Thr,DAPI_Area_Coef');
fprintf(out, ',DAPI_Mean_Low,DAPI_Mean_High,DAPI_Mean_Max,DAPI_Mean_Uni_Thr,DAPI_Mean_Coef');
fprintf(out, ',DAPI_Total_Low,DAPI_Total_High,DAPI_Total_Max,DAPI_Total_Uni_Thr,DAPI_Total_Coef');
fprintf(out, ',DAPI_Min_Low,DAPI_Min_High,DAPI_Min_Max,DAPI_Min_Uni_Thr,DAPI_Min_Coef');
fprintf(out, ',DAPI_Max_Low,DAPI_Max_High,DAPI_Max_Max,DAPI_Max_Uni_Thr,DAPI_Max_Coef');
fprintf(out, ',DAPI_Ecc_Low,DAPI_Ecc_High,DAPI_Ecc_Max,DAPI_Ecc_Uni_Thr,DAPI_Ecc_Coef');
fprintf(out, ',CK_Area_Low,CK_Area_High,CK_Area_Max,CK_Area_Uni_Thr,CK_Area_Coef');
fprintf(out, ',CK_Mean_Low,CK_Mean_High,CK_Mean_Max,CK_Mean_Uni_Thr,CK_Mean_Coef');
fprintf(out, ',CK_Total_Low,CK_Total_High,CK_Total_Max,CK_Total_Uni_Thr,CK_Total_Coef');
fprintf(out, ',CK_Min_Low,CK_Min_High,CK_Min_Max,CK_Min_Uni_Thr,CK_Min_Coef');
fprintf(out, ',CK_Max_Low,CK_Max_High,CK_Max_Max,CK_Max_Uni_Thr,CK_Max_Coef');
fprintf(out, '\n');

for i = 1:size(allstats,1)
    stats = allstats(i);
    fprintf(out, '%d,%d,%d', stats.pos, stats.numManualCTCs, stats.numCTCs);
    fprintf(out, ',%d,%d,%d,%f,%f', stats.DAPIarea(1), stats.DAPIarea(2), stats.DAPIarea(3), stats.DAPIarea(4), stats.DAPIarea(2)/stats.DAPIarea(4));
    fprintf(out, ',%f,%f,%f,%f,%f', stats.DAPImean(1), stats.DAPImean(2), stats.DAPImean(3), stats.DAPImean(4), stats.DAPImean(2)/stats.DAPImean(4));
    fprintf(out, ',%f,%f,%f,%f,%f', stats.DAPItotal(1), stats.DAPItotal(2), stats.DAPItotal(3), stats.DAPItotal(4), stats.DAPItotal(2)/stats.DAPItotal(4));
    fprintf(out, ',%d,%d,%d,%f,%f', stats.DAPImin(1), stats.DAPImin(2), stats.DAPImin(3), stats.DAPImin(4), stats.DAPImin(2)/stats.DAPImin(4));
    fprintf(out, ',%d,%d,%d,%f,%f', stats.DAPImax(1), stats.DAPImax(2), stats.DAPImax(3), stats.DAPImax(4), stats.DAPImax(2)/stats.DAPImax(4));
    fprintf(out, ',%f,%f,%f,%f,%f', stats.DAPIecc(1), stats.DAPIecc(2), stats.DAPIecc(3), stats.DAPIecc(4), stats.DAPIecc(2)/stats.DAPIecc(4));
    fprintf(out, ',%d,%d,%d,%f,%f', stats.CKarea(1), stats.CKarea(2), stats.CKarea(3), stats.CKarea(4), stats.CKarea(2)/stats.CKarea(4));
    fprintf(out, ',%f,%f,%f,%f,%f', stats.CKmean(1), stats.CKmean(2), stats.CKmean(3), stats.CKmean(4), stats.CKmean(2)/stats.CKmean(4));
    fprintf(out, ',%f,%f,%f,%f,%f', stats.CKtotal(1), stats.CKtotal(2), stats.CKtotal(3), stats.CKtotal(4), stats.CKtotal(2)/stats.CKtotal(4));
    fprintf(out, ',%d,%d,%d,%f,%f', stats.CKmin(1), stats.CKmin(2), stats.CKmin(3), stats.CKmin(4), stats.CKmin(2)/stats.CKmin(4));
    fprintf(out, ',%d,%d,%d,%f,%f', stats.CKmax(1), stats.CKmax(2), stats.CKmax(3), stats.CKmax(4), stats.CKmax(2)/stats.CKmax(4));
    fprintf(out, '\n');
end
fclose(out);

grandTotals = grandTotals + [totalTP, totalFP, totalTN, totalFN];

end % imageType

toc

fprintf('Total TP: %d\n', grandTotals(1));
fprintf('Total FP: %d\n', grandTotals(2));
fprintf('Total TN: %d\n', grandTotals(3));
fprintf('Total FN: %d\n', grandTotals(4));

end
