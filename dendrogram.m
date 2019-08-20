function [max_dist] = dendrogram(MST)
% [max_dist] = dendrogram(MST)
%
% The function displays the single linkage dendrogram. 
% The dendrogram is colorized according to the distance from the giant component. 
% 
% MST      : the minimum spanning tree
%             Each row vector consists of [node1 node2 edge_weight].
%             The number of row vector is same to the number of edges in MST. 
% max_dist : the maximum distance from the giant component in [1]. 
%
%
% EXAMPLE: 
% MST = [1 2 0.5; 2 3 0.2];
% dendrogram(MST)
% [max_dist] = dendrogram(MST);
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

p = size(MST,1) + 1;


% connected components. 
cc_count = 1;
bar = [];
for i = 1:p
    bar(cc_count).x1 = 0;
    bar(cc_count).y1 = p-i+1;
    bar(cc_count).cc = i;
    bar(cc_count).from = 0;
    cc_count = cc_count + 1;
end

for i = 1:size(MST,1)
    m = max(find(arrayfun(@(x)(sum(x.cc == MST(i,1))>0),bar)));
    n = max(find(arrayfun(@(x)(sum(x.cc == MST(i,2))>0),bar)));
    
    if m ~= n % always m < n
        bar(m).x2 = MST(i,3);
        bar(m).to = cc_count;
        bar(n).x2 = MST(i,3);
        bar(n).to = cc_count;
        
        bar(cc_count).x1 = MST(i,3);
        bar(cc_count).cc = [bar(m).cc bar(n).cc];
        bar(cc_count).from = [m n];
        bar(cc_count).y1 = (bar(m).y1 + bar(n).y1)/2;
        
        bar(m).y2 = bar(cc_count).y1;
        bar(n).y2 = bar(cc_count).y1;

        cc_count = cc_count + 1;
    end 
end
bar(cc_count-1).x2 = bar(bar(cc_count-1).from(1)).x2*1.1;
bar(cc_count-1).y2 = bar(cc_count-1).y1;
bar(cc_count-1).to = cc_count;
tx = bar(bar(cc_count-1).from(1)).x2*1.1;
ty = p + 0.5; 

bar(end).dist = 0;
for i = (length(bar)-1):-1:1
    bar(i).dist = bar(bar(i).to).dist + 1;
    % display(num2str(bar(i).dist));
end

max_dist = max(arrayfun(@(x)(x.dist),bar));

if nargin == 1
    color_list = colormap(jet(max_dist+1));
    % color_list = color_list(1:end,:);
    % color_list = color_list(round(1:size(color_list,1)/(max_dist+1):size(color_list,1)),:);
    color_list = color_list(end:-1:1,:);
%     S = 255;                        % saturation
%     V = 220;                        % brightness
%     H = [255/(max_dist+1):255/(max_dist+1):255]; % hue
%     color_list = hsv2rgb([H' S*ones(max_dist+1,1) V*ones(max_dist+1,1)]/255);
%     color_list = color_list(end:-1:1,:);

    for i = length(bar):-1:1
        line([bar(i).x1 bar(i).x2],[bar(i).y1 bar(i).y1],'color',color_list(bar(i).dist+1,:),'LineWidth',2); hold on;
        line([bar(i).x2 bar(i).x2],[bar(i).y1 bar(i).y2],'color',color_list(bar(i).dist+1,:),'LineWidth',2);
    end
    
%     figure;
%     for i = length(bar):-1:1
%         i
%         if bar(i).dist <max_dist
%         line([bar(i).x1 bar(i).x2],[bar(i).y1 bar(i).y1],'color',color_list(bar(i).dist+1,:),'LineWidth',2); hold on;
%         line([bar(i).x2 bar(i).x2],[bar(i).y1 bar(i).y2],'color',color_list(bar(i).dist+1,:),'LineWidth',2);
%         end;
%     end
     
    xlim([0 tx]);
    ylim([0 ty]);
    colormap(color_list);
    % xlim([0 3.5]);
    % ylim([0 8.5]);
end