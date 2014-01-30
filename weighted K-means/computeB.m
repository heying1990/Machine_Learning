function result=computeB(xi, belongCenteri, j, center)
result=(pdist2(xi, center(belongCenteri, :),'Euclidean'))^2/(pdist2(xi, center(j, :),'Euclidean'));