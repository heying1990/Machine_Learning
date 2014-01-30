function [ Cluster,center_a ] = kernel_kmean( K,x,center )
% UNTITLED Summary of this function goes here
%   Detailed explanation goes here    
    [n,d] = size(x);
 %   center_a = zeros(K,d);
    buffer = zeros(n,1);
    for i = 1:n
       buffer(i,:) = (pdist2(x(i,:),mean(x),'Euclidean'))^2;
    end
    sum_temp = sum(buffer);
    beta = sum_temp/n;
 %   r = randperm(n);
  %  for i = 1:K
   %    center_a(i,:) = x(r(i),:); 
   % end 
    center_a = center;
    newC = zeros(n,1);
    oldC = ones(n,1);
    while(newC ~= oldC)
        oldC = newC;
        Cluster = zeros(n,1);
        distanceMat = pdist2(x,center_a,'Euclidean');
        fenmu = zeros(K,1);
        fenzi = zeros(K,d);
        for j = 1:n
            [~,Cluster(j)] = min(distanceMat(j,:));
            Kg = kernelize(x(j,:),center_a(Cluster(j),:),beta);
            fenmu(Cluster(j)) = fenmu(Cluster(j)) + Kg;
            fenzi(Cluster(j),:) = fenzi(Cluster(j),:) + Kg*x(j,:);
            center_a(Cluster(j),:) = fenzi (Cluster(j),:)/fenmu(Cluster(j));
        end
        newC = Cluster;
    end
end