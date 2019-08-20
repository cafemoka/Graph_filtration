function [dx,MST] = single_linkage_matrix(cx,threshold)
% [dx,MST] = single_linkage_matrix(cx,threshold)
%
% The function calculates the single linkage matrix and minimum spanning tree. 
% 
% cx        : distance matrix (= number of nodes-by-number of nodes)
% threshold : the maximum weight that you guess that all nodes are connected. 
% dx        : single linkage matrix [1] (= number of nodes-by-number of nodes)
% MST      : the minimum spanning tree
%            Each row vector consists of [node1 node2 edge_weight].
%            The number of row vector is same to the number of edges in MST. 
%
%
% EXAMPLE: 
% [dx,MST] = single_linkage_matrix(cx);
% [dx,MST] = single_linkage_matrix(cx,1);
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


p = size(cx,1);
ind_tril = find(tril(ones(p,p),-1));
ind_diag = find(eye(p,p));

x = repmat([1:p]',[1 p]);
y = repmat([1:p],[p 1]);
E = [x(ind_tril) y(ind_tril) cx(ind_tril)];

[tval tind] = sort(E(:,3),'ascend');
E = E(tind,:);

if nargin < 2
    imax = size(E,1);
else
    imax = max(find(E(:,3) <= threshold));
end

CC = mat2cell([1:p],1,ones(1,p));

dx = zeros(p,p);
ck_mtx = zeros(p,p);
ind_mtx = reshape([1:(p*p)],p,p);

ck_mtx(ind_diag) = 1;

MST = [];

for i = 1:imax
    t1 = find(cellfun(@(x)sum(x == E(i,1)),CC)>0);
    t2 = find(cellfun(@(x)sum(x == E(i,2)),CC)>0);
    
    if t1 ~= t2
        if t1 < t2
            CC{t1} = [CC{t1} CC{t2}]; 
            ind_sel = CC{t1}; 
            tind = setdiff([1:length(CC)],t2);
            CC = CC(tind);
        else 
            CC{t2} = [CC{t2} CC{t1}]; 
            ind_sel = CC{t2}; 
            tind = setdiff([1:length(CC)],t1);
            CC = CC(tind);
        end
    
        ck = ck_mtx(ind_sel,ind_sel);
        tind = find(ck == 0);
        temp = ind_mtx(ind_sel,ind_sel);
        tind = temp(tind);
        
        dx(tind) = E(i,3);
        MST = [MST; E(i,:)];
        ck_mtx(tind) = 1;
    end
%    display([num2str(i) '/' num2str(size(E,1))]);
end
   