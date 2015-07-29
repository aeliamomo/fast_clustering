%function patches = ExtractBlocksPerPixel(file_1,file_2, patchSize, i, j)
patchSize = 15

halfPatchSize = floor(patchSize/2);
totalPointsInPatch = patchSize*patchSize;

a = imread('/Users/chingyukao/Documents/MATLAB/Multi-Modal-Similarity-till-08062015/Multi-Modal-Similarity/Dataset/T1_11.TIFF');
b = padarray(a,[halfPatchSize,halfPatchSize]);

result = [];
for i = 1+halfPatchSize:size(a,1)
    for j = 1+halfPatchSize:size(a,2)
        temp = b(i-halfPatchSize:i+halfPatchSize,j-halfPatchSize:j+halfPatchSize);
        temp = reshape(temp',[1,patchSize*patchSize]);
        pixels = numel(temp);
        result = [result temp];
    end
end
patches = ones(size(result,2)/pixels,totalPointsInPatch);
offset=1;
for i=1:size(result,2)/pixels
  patches(i,:)=result(offset:(offset+totalPointsInPatch-1));
  offset =offset+totalPointsInPatch;
end
patches;
%end