function sortedClusterPatches = displayClusteredPatches(file,clusterPatch, patchSize,percent)
%% Map the patches to its clusters and display it.

%load('clusteredPatches.mat');
clusterPatch = dlmread('CLUSTER_ASSIGNATION');
totalPatches = size(clusterPatch,1);
%patchSize = sqrt(noOfPatches);

Image  = [];
sortedClusterPatches = sortrows(clusterPatch,1);

% uniqueClusters = unique(sortedClusterPatches(:,2));
% totalClusters = size(uniqueClusters,2);


% for k = 1:totalPatches
%     clusterNo = sortedClusterPatches(k,2);
%     subPatch = ones(patchSize,patchSize);
%     clusteredImage = clusterNo * subPatch;
%     Image = [Image ,clusteredImage];
% end
% Image = reshape(Image,[patchSize*sqrt(totalPatches),patchSize*sqrt(totalPatches)]);


patchNoPerRow = sqrt(totalPatches);

finalImage = [];
for k = 0:totalPatches-1
    
    clusterNo = sortedClusterPatches(k+1,2);
    subPatch = ones(patchSize,patchSize);
    clusteredImage = clusterNo * subPatch;
    if mod(k,patchNoPerRow) == 0
       finalImage = [finalImage;Image];
       Image = [];
    end
    Image = [Image ,clusteredImage];
end


%% Plot the points belonging to each cluster with different colours.
figure;
dlmwrite(sprintf('%s_finalImage',file),finalImage);
imagesc(finalImage);
savefig(sprintf('tuningParameters/displayClusteredPatches_%f.fig',percent));
end

