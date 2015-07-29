function [ similarity ] = AssignToClusterKNN( distance_diff , k_value ,final_centroid_similarity)
final = final_centroid_similarity;
dis = distance_diff;
k = k_value;

knn_results = [];%saving result for translation
%knn_results_r = [];%saving result for rotation
min_matrix = [];
dis_idx_length = size(dis,2)-1;

idx_matrix = zeros(1,dis_idx_length);
for i = 1:dis_idx_length
    idx_matrix(1,i) = i;
end
for i =  1:size(dis,1)
    
    [minval, minidx] = min(dis(i,1:end-1));
%    [minvalr, minidxr] = min(dis_r(i,1:end-1));
    %index_check = [index_check, ;minidx];
    min_matrix = [min_matrix;minval, minidx];
%    min_matrix_r = [min_matrix_r;minvalr,minidxr];
    
    %using KNN
    knn_sorted = sortrows([dis(i,1:end-1);idx_matrix]');
%    knn_sorted_r = sortrows([dis_r(i,1:end-1);idx_matrix]');
    
    %take k = 3 for example
    knn_results = [knn_results;knn_sorted(1:k,2)];
%    knn_results_r = [knn_results_r;knn_sorted_r(1:k,2)];
    
end


similarity = [];

begin = 1;
tail = k;
iterationEnd = size(knn_results,1)/(k*11);
for i= 1:iterationEnd
    
    nearest_nb = knn_results(begin:tail);

    final_inverse = final(:,2)';
    final_sim = sum(final(nearest_nb))/k;
    similarity = [similarity; final_sim];
    
    begin = begin +(k) ;
    tail = tail +(k);

end

% knn_results = [];%saving result for translation
% %knn_results_r = [];%saving result for rotation
% min_matrix = []
% dis_idx_length = size(dis,2)-1
% 
% idx_matrix = zeros(1,dis_idx_length)
% for i = 1:dis_idx_length
%     idx_matrix(1,i) = i;
% end
% for i =  1:size(dis,1)
%     
%     [minval, minidx] = min(dis(i,1:end-1));
% %    [minvalr, minidxr] = min(dis_r(i,1:end-1));
%     %index_check = [index_check, ;minidx];
%     min_matrix = [min_matrix;minval, minidx];
% %    min_matrix_r = [min_matrix_r;minvalr,minidxr];
%     
%     %using KNN
%     knn_sorted = sortrows([dis(i,1:end-1);idx_matrix]');
% %    knn_sorted_r = sortrows([dis_r(i,1:end-1);idx_matrix]');
%     
%     %take k = 3 for example
%     knn_results = [knn_results;knn_sorted(1:k,2)];
% %    knn_results_r = [knn_results_r;knn_sorted_r(1:k,2)];
%     
% end
% 
% 
% similarity = [];
% 
% 
% for i= 1:size(knn_results,1)
%     
%     nearest_nb = knn_results(k-(k-1):k)
% 
%     final = final(:,2)';
%     final_sim = sum(final(nearest_nb))/k;
%     similarity = [similarity, final_sim];
%     
%     k = k+k;
%     
% end

end

