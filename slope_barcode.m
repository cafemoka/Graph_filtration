function [s] = slope_barcode(MST)
% [s] = slope_barcode(MST)
%
% The function estimates the slope of barcode-shape of the connected components [1]. 
% 
% MST  : the minimum spanning tree
%             Each row vector consists of [node1 node2 edge_weight].
%             The number of row vector is same to the number of edges in MST. 
% s    : slope of the barcode of the connected components
%
%
% EXAMPLE: 
% [s] = slope_barcode(MST);
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


p = length(MST);
y = [p:-1:1]';
x = sort(MST,'ascend');
s = (p*sum(x.*y) - sum(x)*sum(y))/(p*sum(x.^2) - sum(x)*sum(x));
