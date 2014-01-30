function [center, belongCenter]=WKA(data, numCenter)
[dataNum,d]=size(data);
tmp=randperm(dataNum);
index=tmp(1:2*numCenter);
center =zeros(numCenter,d);
for i = 1:numCenter
    center(i,:)=(data(index(i),:)+data(index(i+1),:))/2;
end
%for i = 1:numCenter
 %   center(i,:)=(test0(i,:)+test0(i+1,:))/2;
%end

belongCenter=zeros(dataNum,1);
a=zeros(dataNum,numCenter);
b=zeros(dataNum,numCenter);
newC = zeros(dataNum,1);
oldC = ones(dataNum,1);
while(newC~= oldC)
   oldC=newC;
   %distMatrix=dist2(data, center);
   distMatrix = pdist2(data,center,'Euclidean');
   for i=1:dataNum
       [minVal,minIndex]=min(distMatrix(i,:));
       belongCenter(i,1)=minIndex;
   end
   for i=1:dataNum
       for j=1:numCenter
           a(i,j)=computeA(data(i, :), j, center, numCenter);
           b(i,j)=computeB(data(i, :), belongCenter(i,1), j, center);
       end
   end
   for j=1:numCenter
       indexRow=find(belongCenter==j);
       center(j,:)=computeM(data, indexRow, a, b, j);
   end
   newC= belongCenter;
end