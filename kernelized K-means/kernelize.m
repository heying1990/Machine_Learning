function [ Kg ] = kernelize( x,y,beta )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    dis2xy = (pdist2(x,y,'Euclidean'))^2;
    Kg = exp(-dis2xy/beta);

end

