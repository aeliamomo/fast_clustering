function visualizePatches(patches, sortedClusterPatches,percent)

%patches = dlmread('patches.dat');
sortedClusterPatches = dlmread('CLUSTER_ASSIGNATION');
noOfClusters = size(unique(sortedClusterPatches(:,2)),1);

img_size = sqrt(size(patches,2));

new_cluster = [];
for k = 1:noOfClusters
    new_cluster = sortedClusterPatches(sortedClusterPatches(:,2) == k);
    no_of_figure = size(new_cluster,1);
    imgs = cell(no_of_figure,1);
    for i=1:no_of_figure
        new_cluster(i);
        imgs{i} = mat2gray(reshape(patches(new_cluster(i),:),[img_size,img_size]));
    end
    figure(2+k);
    for i=1:no_of_figure
        if(mod(no_of_figure,10)>0 )  
            subplot(round(no_of_figure/10)+1,10,i);   
        else
            subplot(round(no_of_figure/10),10,i);
        end
        h = imshow(imgs{i}, 'InitialMag',100, 'Border','tight');
        title(num2str(i));
    end
    savefig(sprintf('tuningParameters/patchesVisualize_%f_cluser_%d.fig',percent,k));
end

end
