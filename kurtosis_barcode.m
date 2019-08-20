function [kurt] = kurtosis_barcode(MST)
% [kurt] = kurtosis_barcode(MST)
%
% The function estimates the kurtosis of barcode-shape of the connected components. 
% 
% MST  : the minimum spanning tree
%             Each row vector consists of [node1 node2 edge_weight].
%             The number of row vector is same to the number of edges in MST. 
% kurt : kurtosis which represents the peakedness and heavy-tailness of the shape. 
%
%
% EXAMPLE: 
% [kurt] = kurtosis_barcode(MST);
%
%
% The code is downloaed from https://sites.google.com/site/hkleebrain/home/persistent-homology-2
%
% (C) Hyekyoung Lee, hklee.brain@gmail.com
%  Department of Nuclear Medicine
%  Seoul National University, South Korea
%
% 
% If you use this code, please reference one of the following papers. The details on
% the mathematical basis of of the algorithm can be found in these papers.
%
% [1] Lee, H., Kang, H.K., Chung, M.K., Kim, B.-N., Lee, D.S. 2012. 
%     Persistent brain network homology from the perspective of dendrogram, 
%     IEEE Transactions on Medical Imaging. in press. 
%     https://docs.google.com/file/d/0BzqCeYxcOj3bOFRXT3lZY0YwQUU/preview
%
% [2] Lee, H., Kang, H.K., Chung, M.K., Kim, B.-N., Lee, D.S. 2011. 
%     Computing the shape of brain network using graph filtration and Gromov-Haudorff metric,
%     MICCAI2011.
%     https://docs.google.com/file/d/0BzqCeYxcOj3bYWU5YzhhYzYtM2Q4NC00YTFmL
%     WE4YWItZjkwZjI5ZThkNjlk/edit?pli=1
%
% Update history: Sep 12, 2012




x = sort(MST,'ascend');
sample = [];
for i = 1:length(MST)
    sample = [sample; x(i)*ones(i,1); -x(i)*ones(i,1);];
end
kurt = kurtosis(sample);
