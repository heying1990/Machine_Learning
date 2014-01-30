function [center,belongCenter,centralizedData,C]=kmeansCenter(data,numCenter,epsilon)
% compute the numCenter centroid of the data, termination criteria is
% epsilon, belongCenter is the final result of the index number of centroids
% that each vector is divided into
% centralized data is the compressed version, replacing each block with the
% centroid it belongs to
% innitialize the centroid
[dataNum,d]=size(data);
tmp=randperm(dataNum);
index=tmp(1:numCenter);
center=data(index,:);
belongCenter=zeros(dataNum,1);
distanceVec=zeros(dataNum,1);%record the distance of each sample wrt its center
difference=1000;
newC=0;
while(difference>=epsilon)
    oldC=newC;
    distMatrix=dist2(data,center);
    for i=1:dataNum
        [minVal,minIndex]=min(distMatrix(i,:));
        belongCenter(i,1)=minIndex;
    end
    for i=1:numCenter
        indexRow=find(belongCenter==i);
        center(i,:)=mean(data(indexRow,:));
    end
    for i=1:dataNum
        distanceVec(i,1)=dist2(data(i,:),center(belongCenter(i,1),:));
    end
    newC=sum(distanceVec);
    difference=newC-oldC;
end
C=newC;
centralizedData=zeros(dataNum,d);
for i=1:dataNum
    centralizedData(i,:)=center(belongCenter(i,1),:);
end

