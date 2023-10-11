function [ outputMat ] = cutBorders(inputCell, borderSize )
%CUTBORDERS Summary of this function goes here
%   Detailed explanation goes here

for i = 1:size(inputCell,1)
    for j = 1:size(inputCell,2)
        if i>1 && i<size(inputCell,1)
            inputCell{i,j} = inputCell{i,j}(borderSize+1:end-borderSize,:);
        elseif i>1
            inputCell{i,j} = inputCell{i,j}(borderSize+1:end,:);
        elseif i<size(inputCell,1)
            inputCell{i,j} = inputCell{i,j}(1:end-borderSize,:);
        end
        
        if j>1 && j<size(inputCell,2)
            inputCell{i,j} = inputCell{i,j}(:,borderSize+1:end-borderSize);
        elseif j>1
            inputCell{i,j} = inputCell{i,j}(:,borderSize+1:end);
        elseif j<size(inputCell,2)
            inputCell{i,j} = inputCell{i,j}(:,1:end-borderSize);
        end
    end
end

outputMat = cell2mat(inputCell);

end

