
%% MLMI Project on Multi-Modal-Similarity.

% Steps to extract patches out of each modality image, represent them into Distance matrix, label them and compare with other modality.

%% Step 1: Read the grayscale image.
rng(16);
CreateTIFF();
fig_index = 3;

%% Step 2: Divide the image into number of patches using blockproc function in Matlab.
%path = strcat(pwd,'/Multi-Modal-Similarity/Dataset/');

% Call initial_config function to initialize file name, path & patch size.
[file, pathFile, patchSize, applyTranslationRotation] = initial_config();

patchPairMatrix = extractAppendedPatches(applyTranslationRotation);
%[optional] using reduced image
%reducedImg = reduceBoundingBox(pathFile);

% Changes to pass Patch Size dynamically as an input parameter.
%patches = ExtractBlocks(reducedImg,patchSize);

patches = ExtractBlocks(pathFile, patchSize);
appendedPatches = patchPairMatrix(:,1:patchSize*patchSize*2);
%[a,b,patches] = patchset(pathFile, patchSize,patchSize);

%write patches to file 
%dlmwrite('patches.dat',patches);

%% Step 3a: Compare patch by patch using SSD or NCC as distance measures.
% SSD - Make sure to use float, as raising integer to power of 2 will go out of range.
% NCC - Use Normalized 2-D cross-correlation (normxcorr2) function in Matlab.

%%using appened patches instead of patches 
patches = appendedPatches;
noOfPatches = size(patches, 1);

%save('noOfPatches.mat','noOfPatches');




distanceMatrix = zeros(noOfPatches,noOfPatches);
% t = cputime
% for i = 1:noOfPatches
%     for j = 1:noOfPatches
%         %distanceMatrix(i,j) = mean(normxcorr2(patches(i,:),patches(j,:))); % for NCC
%         %distanceMatrix(i,j) = sum((patches(i,:)-patches(j,:)).^2); % for SSD
%         distanceMatrix(i,j) = sum((bsxfun(@minus, patches(i,:), patches(j,:)).^2));
%         %distanceMatrix(i,j) = sqrt(sum((patches(i,:) - patches(j,:)).^2)./patchSize);
%         %distanceMatrix(i,j) = norm(patches(i,:) - patches(j,:));
%         %dotProduct = dot(patches(i,:),patches(j,:));
%         %distanceMatrix(i,j) = dotProduct/(norm(patches(i,:))*norm(patches(j,:)));
%         distanceMatrix(j,i) = distanceMatrix(i,j);
%     end
% end
% withTwoForLoops = cputime-t

patches = single(patches);
t2 = cputime
IP = dot(patches',patches');
distanceMatrix = bsxfun(@plus,IP',IP)-2*(patches*patches');
e2 = cputime-t2
% Step 3b: Plot the patch by patch co-relations on a Distance Matrix.

[rows,cols] = size(distanceMatrix);
totalRows = rows*(cols-1)/2; % Excluding the elements (j,i) for every (i,j) as it is a Symmetric Matrix.
distMatFile = zeros(totalRows,3);

% Computes the Euclidean distance between pairs of objects in m-by-n data matrix using pdist(distanceMatrix) Matlab function.
% D is a row vector of length m(m?1)/2, corresponding to pairs of observations in distanceMatrix.
% D = pdist(distanceMatrix);

maxValue = max(max(distanceMatrix));
k=1;
for i = 1:rows
   for j = i+1:cols
       distMatFile(k,:) = [i,j,(distanceMatrix(i,j)/maxValue)];
       %distMatFile(k,:) = [i,j,D(k)];
       k = k+1;
   end
end

extDat = '.dat';
datfile = strcat(file,extDat);

%save(datfile, 'distMatFile', '-ascii');
dlmwrite(datfile,distMatFile)

%% Step 4: Compute the clustering according to the idea proposed in the paper using the Distance matrix computed above.
% Step 4a: Compute the Density matrix using
%          dc as some cut-off or threshold value.
%          For all points di in one row across all coulmns j of the distance matrix.
%          X(dij - dc)=1 if (dij - dc)<0 and X(dij - dc)=0 if(dij - dc)>0
%          Take the sum which is the density for point i.

% Step 4b: Compute the Delta matrix using
%          Minimum distance between point i and any other point with higher density value.
%          Take all points having density greater than the current point and take the least one out of them.

% Step 4c: Plot Density as X-Axis and Delta as Y-Axis matrices.
%          Select points with higher density and delta values as cluster centroids.


% Step 4d: Plot Density as X-Axis and Delta as Y-Axis matrices.
%          Select points with higher density and delta values as cluster centroids.

% Above Step 4 is already implemented in paper in the file cluster_dp.m
% Use it by wrapping a function around it.

%cluster_dp;
noOfCutOffPerVal = 1;
stepCutOffPerVal = 0.2;
cutOffPercentage = ones(noOfCutOffPerVal,1);
cutOffPercentage(:,1)=[stepCutOffPerVal:stepCutOffPerVal:noOfCutOffPerVal*stepCutOffPerVal];
cutOffPercentage = [10];
%noOfCutOffPerVal = 2;
%cutOffPercentage = [1,1.3,1.5,1.8,2,5,10,50,98,99];
% lstClusterPatches = {};
% for c=1:noOfCutOffPerVal
%     clusterPatch = cluster_dp(cutOffPercentage(c));
%     lstClusterPatches = [lstClusterPatches; clusterPatch];
% end

%% Step 5: Label the patches in the same modality using the cluster centroids.

% Map clusters to patches and display it.
%load('noOfPatches.mat');
%sortedClusterPatches = displayClusteredPatches(noOfPatches);

% Visualize the patches and clusters.
%visualizePatches(sortedClusterPatches);

for c=1:noOfCutOffPerVal
    clusterPatch = cluster_dp(distMatFile, cutOffPercentage(c));
    %sortedClusterPatches = displayClusteredPatches(file,clusterPatch, patchSize, cutOffPercentage(c));
    %visualizePatches(patches, sortedClusterPatches, cutOffPercentage(c));
end

assignedCluster = dlmread('CLUSTER_ASSIGNATION');

assignedCluster = [assignedCluster, patchPairMatrix(:,end)];%adding labels 
dlmwrite('CLUSTER_ASSIGNATION_withIndex',assignedCluster);
assignedCluster = [assignedCluster(:,2),assignedCluster(:,end)]%taking only two columns: namely no. of cluster and the labels
assignedCluster = sortrows(assignedCluster,1)%sorting according no. of cluster
dlmwrite('CLUSTER_ASSIGNATION_sorted',assignedCluster);

data = dlmread('CLUSTER_ASSIGNATION_sorted');
noOfCentroids = data(end,1);
%data = [1 1;2 0;2 1; 2 1;3 0 ;3 0]
%noOfCentroids = data(end,1)

X = [];

%convert tranlation and rotation labels to zero
if(applyTranslationRotation==1)
    data(data(:,2)==2,2) = 0
    data(data(:,2)==3,2) = 0
end
%computer probability of each cluster
for i = 1:noOfCentroids
    %%fractions / numerator
    fraction = sum(data(data(:,1)==i,2));
    numerator = size(data(data(:,1)==i,2),1);
    temp_result = fraction / numerator;
    %X is similarity vector
    X = [X, temp_result];
    
end

%%according the centroid to get coresponding patch pairs
centroids = dlmread('CLUSTER_ASSIGNATION');
assignedcentroids = sortrows(centroids,2);
uniqueCentroids = unique(centroids(:,4));
%%append similarity and centroids

final = [X', uniqueCentroids];
dlmwrite('similarityAndCentroids',final);

%taking translated similarity of all 12 images(ohne no.8) 
k_value = 10;
translation_similarity_nn = [];
rotation_similarity_nn = [];
translation_similarity_temp = [];
rotation_similarity_temp = [];
translation_roc = [];
rotation_roc = [];
target = zeros(13,1);%13 is because 13 different translation scales.
target(ceil(size(target,1)/2),:)=1;
for p = 1:12

isTransformed = 1;
if p<10
    if p == 8
    continue;
    end
    image_number = strcat(num2str(0),num2str(p))
end

misaligned_patches = extractTranslatedPatches(isTransformed,image_number);
rotated_patches = extractRotatedPatches(isTransformed,image_number);

dis = []
dis_r = []
for i = 1:size(misaligned_patches,1)
    d = [];
    d_r = [];
    misaligned_patches(i,:)';
    rotated_patches(i,:)';
    for k = 1:noOfCentroids
        centroid_intensity = appendedPatches(final(k,2),:);
        %d is the distance between centroid intensity and misaligned_patch
        %intensity
        d = [d, sqrt(sum((misaligned_patches(i,:)'-centroid_intensity').^2))];
        d_r = [d_r, sqrt(sum((rotated_patches(i,:)'-centroid_intensity').^2))];
    end
    dis = [dis;d];
    dis_r = [dis_r; d_r];
end
%%make translation labels
%if dis = 143*14 means that there are 11images *14 translated scales = 143,
%and 14 means there 14 centroids

dis(:,end+1) = 99;
dis_r(:,end+1) = 99;
index = 1;
for i = 1:size(dis,1)
    dis(i,end) = index;
    dis_r(i,end) = index;
    if mod(i,11) == 0
        index = index +1 ;
    end
end



translation_similarity = AssignToClusterKNN(dis,k_value,final);
rotation_similarity = AssignToClusterKNN(dis_r,k_value,final);

translation_similarity_temp = [translation_similarity_temp,translation_similarity];
rotation_similarity_temp = [rotation_similarity_temp,rotation_similarity];

translation_label = [translation_similarity,target];
translation_roc = [translation_roc;translation_label];
rotation_label = [rotation_similarity,target];
rotation_roc = [rotation_roc;rotation_label];

end

translation_similarity_nn = sum(translation_similarity_temp')/size(translation_similarity_temp',1);
rotation_similarity_nn = sum(rotation_similarity_temp')/size(rotation_similarity_temp',1);
% 
% translation_similarity_nn = similarity;
% rotation_similarity_nn = similarity;
% translation_similarity_nn= [];
% rotation_similarity_nn = [];
% translation_similarity_knn = [];
% rotation_similarity_knn = [];
% temp_similarity = [];
% temp_similarity_r = [];
% min_matrix(:,end+1)= dis(:,end);
% min_matrix_r(:,end+1) = dis_r(:,end);
% index = 1
% for i = 1:size(min_matrix(:,1))
%     temp = min_matrix(min_matrix(:,3)==index,2);
%     temp_r = min_matrix_r(min_matrix_r(:,3)==index,2);
%     %m = mode(temp);
%     %misaligned patch should be classified to the clusters which contain low
%     %similarity.. not high similarity.
%     m = temp(mod(i,11)+1,1);
%     mr = temp_r(mod(i,11)+1,1);
%     temp_similarity = [temp_similarity ,final(m,1)];
%     temp_similarity_r = [temp_similarity_r ,final(m,1)];
%     if mod(i,11) == 0
%         average_similarity =  sum(temp_similarity)/11;
%         average_similarity_r =  sum(temp_similarity_r)/11;
%         translation_similarity_nn = [translation_similarity_nn, average_similarity];
%         rotation_similarity_nn = [rotation_similarity_nn, average_similarity_r];
%         temp_similarity = [];
%         temp_similarity_r = [];
%         index = index +1;
%     end
% end

%% Capture-Range Plot

translations = [[-60,-60];[-50,-50];[-40,-40];[-30,-30];[-20,-20];[-10,-10];[0,0];[10,10];[20,20];[30,30];[40,40];[50,50];[60,60]];
totalTranlations = size(translations,1);

rotations = [-90;-75;-60;-45;-30;-15;0;15;30;45;60;75;90];
figure(3);
%plot3(translations(:,1),translations(:,2),translationSimilarity);
% Plot a 3D Surface Plot for Translation Similarity
%X = reshape(translations(:,1),[sqrt(totalTranlations),sqrt(totalTranlations)]);
%Y = reshape(translations(:,2),[sqrt(totalTranlations),sqrt(totalTranlations)]);
%Z = reshape(translationSimilarity,[sqrt(totalTranlations),sqrt(totalTranlations)]);
plot3(translations(:,1),translations(:,2),translation_similarity_nn,'-o');
%surf(X,Y,Z);
xlabel('Translation x')
ylabel('Translation y')
zlabel('Similarity')
title('Capture range plot for fast clustering')

figure(4);
plot(rotations(:,1),rotation_similarity_nn,'-o');
xlabel('Rotation')
ylabel('Similarity')
title('Capture range plot for fast clustering')


%%computing ROC and PR curve : 

figure(5);
true_label = translation_roc(:,2)';
data = translation_roc(:,1)';
[tpr, fpr, thresholds] = roc(true_label,data);
%[prec, tpr, fpr, thresh] = prec_rec(translation_similarity_nn, target,'plotPR',1)

plot(fpr,tpr,'-o')
xlabel('False positive rate');
ylabel('True positive rate');
title('ROC for Classification by Clustering');


%ROC curve 
% target = zeros(1,size(translation_similarity_nn,2));
% target(1,ceil(size(target,2)/2))=1;
% figure(5);
% [tpr, fpr, thresholds] = roc(target, translation_similarity_nn);
% %[prec, tpr, fpr, thresh] = prec_rec(translation_similarity_nn, target,'plotPR',1)
% 
% plot(fpr,tpr,'-o')
% xlabel('False positive rate')
% ylabel('True positive rate')
% title('ROC for Classification by Clustering')

%PR courve
figure(6)
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(true_label, data, 1, 'xCrit', 'reca', 'yCrit', 'prec');
plot(Xpr,Ypr,'-o')
xlabel('Recall'); ylabel('Precision')
title(['Precision-recall curve (AUC: ' num2str(AUCpr) ')'])


%%plotting roc and AUC curves using different k values.
%AUC
x =[1 ,3,5,6,7,8,9,10];
y = [0.024772,0.0627, 0.0727,0.076,0.078,0.071,0.071,0.059];%receive from different k values.
plot(x,y,'-*')
xlabel('k value'); ylabel('AUC');
title('AUC according to different k-values for kNN');

y1 = [0    0.6364];
x1 = [0    0.8409];

y3 = [0    0.2727    0.6364    0.8182];
x3 = [ 0    0.2727    0.4545    0.8636];

y5 = [0    0.1818    0.1818    0.2727    0.8182    0.8182];
x5 = [0    0.0606    0.1212    0.1970    0.7955    0.9242];

y6 = [0    0    0.0909    0.1818    0.2727    0.5455    0.8182];
x6 = [0    0    0.0530    0.1212    0.1818    0.6591    0.9015];

y7 = [0         0    0.1818    0.3636    0.7273    0.8182    1.0000    1.0000];
x7 = [0         0    0.1212    0.2273    0.6742    0.8712    0.9470    1.0000];


y9 = [0         0    0.2727    0.7273    0.8182    1.0000    1.0000    1.0000];
x9 = [0         0    0.2652    0.6364    0.8333    0.9621    0.9697    1.0000];
figure(7)
plot(x1,y1,x3,y3,x5,y5,x7,y7,x9,y9)
xlabel('False positive rate');ylabel('True positive rate');
title('ROC by using different k-values with kNN');
legend('k = 1','k = 3', 'k = 5', 'k = 7','k = 9');
