originalImage  

transformedImage = ImageTransformations(imagePath, transformation, boolTranslation)
patchSize = 15 
noOfsample = 20

similarity_patches =  extractSimilarPosPatches(originalImage, transformedImage, patchSize, noOfsample)

% A = 1:1:size(X,2)
% c = linspace(1,10,length(A));
% a = 25;
% scatter(A,X,a,c,'filled')

% data = dlmread('CLUSTER_ASSIGNATION_sorted');
% noOfCentroids = data(end,1)
% %data = [1 1;2 0;2 1; 2 1;3 0 ;3 0]
% %noOfCentroids = data(end,1)
% 
% X = []
% 
% for i = 1:noOfCentroids
%     
%     %%fractions / numerator
%     fraction = sum(data(data(:,1)==i,2))
%     numerator = size(data(data(:,1)==i,2),1)
%     temp_result = fraction / numerator;
%     X = [X;temp_result]
%     
% end