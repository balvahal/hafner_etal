function [] = scatterCorrelation(x,y)
scatter(x,y,10,[0.4,0.4,0.4]);
title(sprintf('Spearman r = %0.2f', corr(x,y,'type', 'spearman')));
end