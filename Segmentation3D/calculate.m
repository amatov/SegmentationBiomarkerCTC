sumAR = 0;
sumDAPI = 0;

c3.cropped = ARAnn.cuboid;
c4.cropped = DAPI.cuboid;

for i = 3:size(ARAnn.cuboid,3)
    ARmask = poly2mask(ARAnn.x{i},ARAnn.y{i},size(c3.cropped,1),size(c3.cropped,2));
    DAPImask = poly2mask(DAPIAnn.x{i},DAPIAnn.y{i},size(c3.cropped,1),size(c3.cropped,2));
%     figure;
    myimshow(c3.cropped(:,:,i),[],0.001);
    hold on;
    plot(ARAnn.x{i},ARAnn.y{i},'.','markersize',12);
    plot(ARAnn.x{i}(:),ARAnn.y{i}(:),'g');
%     figure;
    myimshow(c4.cropped(:,:,i),[],0.001);
    hold on;
    plot(DAPIAnn.x{i},DAPIAnn.y{i},'.','markersize',12);
    plot(DAPIAnn.x{i}(:),DAPIAnn.y{i}(:),'g');
    
    sumAR = sumAR + sum(sum(ARmask.*c3.cropped(:,:,i)))/1000;
    sumDAPI = sumDAPI + sum(sum(DAPImask.*c3.cropped(:,:,i)))/1000;
    
    
    
end