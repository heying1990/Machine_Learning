function result=computeDistW(data, center, belongCenter)
centerNum=size(center, 1);
result=0;
for k=1:centerNum
   result=result+sqrt(dist2(data, center(k, :)));
end
result=result*dist2(data, center(belongCenter,:));