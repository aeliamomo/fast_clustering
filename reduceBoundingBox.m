function reducedImage = reduceBoundingBox(ImagePath)
%REDUCEBOUNDINGBOX 
%Image = imread('/Users/chingyukao/Documents/MATLAB/Multi-Modal-Similarity/Multi-Modal-Similarity/Dataset/T2_05.TIFF');
%imshow(Image);

Image = imread(ImagePath);
%simshow(Image);

row = size(Image,1)
col = size(Image,2) 

getRows = find(sum(Image,2)>0);
newRowStart = min(getRows);
newRowEnd = max(getRows);

getCols = find(sum(Image,1)>0);
newColStart = min(getCols);
newColEnd = max(getCols);

reducedImage = Image(newRowStart:newRowEnd,newColStart:newColEnd);
size(reducedImage)
%figure, imshow(reducedImage, []);


end