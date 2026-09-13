function annotations = AnnotationTool(cuboid)

h = figure;
for i = 1:size(cuboid,3)
    figure(h);
    myimshow(cuboid(:,:,i),[],0.001);
    [x{i},y{i}] = getline(h,'closed');    
end


keyboard;

return;









