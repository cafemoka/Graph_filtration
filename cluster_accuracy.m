function [acc] = cluster_accuracy(D,ytrue,nc)
% [acc] = cluster_accuracy(D,ytrue,nc)
%
% The function estimates the clustering accuracy by comparing the estimated
% cluster label based on the distance matrix 'D' to the true label 'ytrue'. 
% 
% D       : (Gromov-Hasdorff) distance matrix between samples/networks
%            = number of networks-by-number of networks
% ytrue   : true label of samples/networks, for example, 1,2,...
%           vector with the dimension, number of networks 
% nc      : the number of predetermined number of clusters 
%           This is used for estimating the cluster label based on 'D'
% acc     : clustering accuracy, [0 1]
%
%
% EXAMPLE: 
% D = [0 100 100; 100 0 100; 100 100 0];
% ytrue = [1; 2; 3];
% nc = 3;
% [acc] = cluster_accuracy(D,ytrue,nc);
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


n = size(D,1);
ind_tril = find(tril(ones(n,n),-1));
Y = D(ind_tril)';
Z = linkage(Y,'ward');
yhat = cluster(Z,'maxclust',nc);

nn = [];
for i = 1:nc
    [nn(i,:) x] = hist(yhat(find(ytrue == i)),[1:nc]);
end

[tval tind] = sort(nn');
[ttval ttind] = sort(tval(end,:),'descend');
tag = zeros(nc,1);
estimated_label = zeros(nc,1);

i = 1;
while prod(tag) == 0
    if tag(tind(end,ttind(i))) == 0 & estimated_label(ttind(i)) == 0
        estimated_label(ttind(i)) = tind(end,ttind(i));
        tag(tind(end,ttind(i))) = 1;
        tval(:,ttind(i)) = [tval(end,ttind(i)); tval(1:(end-1),ttind(i))];
        tind(:,ttind(i)) = [tind(end,ttind(i)); tind(1:(end-1),ttind(i))];
        i = i + 1;
    else
        tval(:,ttind(i)) = [tval(end,ttind(i)); tval(1:(end-1),ttind(i))];
        tind(:,ttind(i)) = [tind(end,ttind(i)); tind(1:(end-1),ttind(i))];
        [ttval ttind] = sort(tval(end,:),'descend');
        if unique(ttval) == 0
            estimated_label(find(estimated_label==0)) = find(tag == 0);
            tag(find(tag==0)) = 1;
        else
            i = 1;
        end
    end
end

temp = [];
for i = 1:nc
    temp(yhat == estimated_label(i)) = i;
end
yhat = temp; 

acc = sum(reshape(yhat,prod(size(yhat)),1) == reshape(ytrue,prod(size(ytrue)),1))/length(ytrue);
