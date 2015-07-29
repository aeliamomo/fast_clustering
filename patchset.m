function [ patch patchtemp patchreshaped ] = patchset( I,row,cols )



%PATCHSET Summary of this function goes here

%   Detailed explanation goes here
I = imread('/Users/chingyukao/Documents/MATLAB/Multi-Modal-Similarity-till-08062015/Multi-Modal-Similarity/Dataset/T1_11.TIFF');
patchSize = row;
halfPatchSize = floor(patchSize/2);
I = padarray(I,[halfPatchSize,halfPatchSize]) ;
row1 = 1;
row2=row;
column1 = 1;
column2 = cols;
patchnum=1;

while(~(column2==size(I,2) && (row2==size(I,1))))
 if column2==size(I,2)
 patch{patchnum}=I(row1:row2,column1:column2);
 patchnum=patchnum+1;
 column1=1;
 column2=cols;
 row2=row2+1;
 row1=row1+1;
 else
 patch{patchnum}=I(row1:row2,column1:column2);
 patchnum=patchnum+1;
 column2 = column2+1;
 column1 = column1+1;
 end
 patch{patchnum}=I(row1:row2,column1:column2);
end
 patchtemp=cellfun(@(x) reshape(x',[1,patchSize*patchSize]) ,patch,'UniformOutput', false);
 patchreshaped = vertcat(patchtemp{:});
end