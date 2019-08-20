function [Z] = distance_between_mtx(X,Y)
% [Z] = distance_between_mtx(X,Y)
%
% The function estimates the distance between the point in X and one in Y. 
% 
% X : the point set 
%     Each row vector consists of [node1 node2]. 
%     The number of row vector is same to the number of nodes in set X. 
% Y : the point set 
%     Each row vector consists of [node1 node2]. 
%     The number of row vector is same to the number of nodes in set X. 
% Z : the distance matrix between the point in X and one in Y 
%     If the number of points in X is p1 and one in Y is p2, the dimension
%     of Z is p1*p2. 
%
%
% EXAMPLE: 
% [Z] = distance_between_mtx(X,Y);
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


m = size(X,1);
n = size(Y,1);
p = size(X,2); % == size(Y,2)

Xdiag = sum(X.^2,2);
Ydiag = sum(Y.^2,2);

Z = repmat(Xdiag,[1 n]) + repmat(Ydiag',[m 1]) - 2*X*Y';