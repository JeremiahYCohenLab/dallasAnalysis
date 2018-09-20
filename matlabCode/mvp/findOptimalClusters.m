function findOptimalClusters(range)

if nargin < 1
    range = 10;
end

for i = 1:range
    [tbl, idx, D] = mvpClassification(0,0,i);
    tmp = [];
    for j = 1:length(idx)
        tmp(j) = (D(j, idx(j)))^2; 
    end
    sumSq(i) = sum(tmp);
end
        
figure; plot(sumSq, '-b', 'LineWidth', 2)
ylabel('total within-clusters sum of squares')
xlabel('number of clusters')
title('elbow method')

end