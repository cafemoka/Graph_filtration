function [d] = bottleneck_distance(X,Y)
% [d] = bottleneck_distance(X,Y);
%
% The function estimates the bottleneck distance between two networks. 
% 
% X : Each row vector consists of [birth_time death_time] of connected components.
%     The birth time of connected components in the Rips complex is always 0. 
%     Whenever the edge is added to MST by increasing the threshold, the
%     connected component disappear. 
%     The number of row vector is same to the number of connected components in the network. 
% Y : Each row vector consists of [birth_time death_time] of connected components.
%     The birth time of connected components in the Rips complex is always 0. 
%     Whenever the edge is added to MST by increasing the threshold, the
%     connected component disappear. 
%     The number of row vector is same to the number of connected components in another network. 
% d : the bottleneck distance [1]
%
%
% EXAMPLE: 
% [d] = bottleneck_distance(X,Y);
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


D = distance_between_mtx(X,Y);
Dxy = zeros(size(X,1)+size(Y,1),size(X,1)+size(Y,1));
Dxy(1:size(X,1),(size(X,1)+1):(size(X,1)+size(Y,1))) = D;
Dxy = Dxy + Dxy';
[Tree, pred] = graphminspantree(sparse(Dxy));
d = full(max(max(Tree)));


% % old version
% D = distance_between_mtx(X,Y);
% 
% [Dval Dind] = sort(D,'ascend');
% 
% mappingXY = zeros(size(X,1),1);
% mappingYX = zeros(size(Y,1),1);
% %count = 1;
% while ~(sum(mappingXY==0)==0 & sum(mappingYX==0)==0)
%     [tval tind] = sort(Dval(1,:),'ascend');
%     x = Dind(1,tind(1));
%     y = tind(1);
%     if mappingXY(x) == 0 & mappingYX(y) == 0
%         mappingXY(x) = y;
%         mappingYX(y) = x;
%     end 
%     Dval(:,tind(1)) = [Dval(2:end,tind(1)); inf];
%     Dind(:,tind(1)) = [Dind(2:end,tind(1)); inf];
%     
%     %display(num2str(count));
%     %count = count + 1;
% end
% 
% D = distance_between_mtx(X,Y(mappingXY,:));
% d = max(diag(D));