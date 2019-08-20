function [pair] = correspondence(X,Y)
% [pair] = correspondence(X,Y)
%
% The function estimates the correspondence between the point in X and one in Y. 
% 
% X    : the point set 
%        Each row vector consists of [node1 node2]. 
%        The number of row vector is same to the number of nodes in set X. 
% Y    : the point set 
%        Each row vector consists of [node1 node2]. 
%        The number of row vector is same to the number of nodes in set X. 
% pair : Each row vector consists of [node1 node2] 
%        where node1 is the index of the point in X and node2 is one in Y. 
%        node1 corresponds to node2. 
%
%
% EXAMPLE: 
% [pair] = correspondence(X,Y);
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


Z = distance_between_mtx(X,Y);
xind = repmat([1:size(X,1)]',[1 size(Y,1)]);
yind = repmat([1:size(Y,1)],[size(X,1) 1]);

tval = reshape(Z,size(X,1)*size(Y,1),1);
xind = reshape(xind,size(X,1)*size(Y,1),1);
yind = reshape(yind,size(X,1)*size(Y,1),1);

[tval tind] = sort(tval,'ascend');

tlist = [xind(tind) yind(tind) tval];

xtag = zeros(size(X,1),1);
ytag = zeros(size(Y,1),1);

pt = 1;
pair = [];
while prod(xtag) == 0 & prod(ytag) == 0
    if xtag(tlist(pt,1)) == 0 & ytag(tlist(pt,2)) == 0
        pair = [pair; tlist(pt,1:2)];
        xtag(tlist(pt,1)) = 1; 
        ytag(tlist(pt,2)) = 1;
    end
    pt = pt + 1;
end

% color_list = colormap(jet(103));
% color_list = color_list(randperm(103),:);
% figure;
% for i = 1:103
%     hold on;
%     plot(X(i,1),X(i,2),'o','MarkerEdgeColor',color_list(i,:),'LineWidth',3);
%     plot(Y(i,1),Y(i,2),'x','MarkerEdgeColor',color_list(i,:),'LineWidth',3,'MarkerSize',10);
% end
    

