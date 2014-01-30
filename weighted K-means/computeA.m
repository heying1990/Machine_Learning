function result=computeA(xi, j, center, numCenter)
result=pdist2(xi, center(j,:),'Euclidean');
buffer = zeros(numCenter,1);
for i=1:numCenter
    buffer(i)=2*pdist2(xi, center(i, :),'Euclidean');
end
result = result+sum(buffer);



