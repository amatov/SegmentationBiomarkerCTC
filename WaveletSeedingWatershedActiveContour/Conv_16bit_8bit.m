function [ outImg ] = Conv_16bit_8bit( inpImg )
%CONV_16BIT_8BIT Summary of this function goes here
%   Detailed explanation goes here
I = inpImg;
se = strel('disk', 3);
Io = imopen(I, se);
figure, imshow(Io), title('Opening (Io)')

Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);

Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);

fgm = imregionalmax(Iobrcbr);


[s,ind] = sort(inpImg(:));
L = length(s);
outImg = zeros(size(inpImg),'uint8');
grid = (linspace(1,L,257));

for i = 0:255
    outImg(inpImg>=s(grid(i+1)) & inpImg<s(grid(i+2))) = i;
end

end

