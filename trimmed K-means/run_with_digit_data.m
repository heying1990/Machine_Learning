clear;
clc;

load mnist_all;

% data = [test0(1:500,:)' test1(1:500,:)' ...
%         test2(1:500,:)' test3(1:500,:)' ...
%         test4(1:500,:)' test5(1:500,:)' ...
%         test6(1:500,:)' test7(1:500,:)' ...
%         test8(1:500,:)' test9(1:500,:)'];
% data = [test0(1:100,:)' test1(1:100,:)' ...
%         test3(1:100,:)' test4(1:100,:)' ...
%         test5(1:100,:)' test6(1:100,:)'];
% data = [test1(1:100,:)' test3(1:100,:)' ...
%         test9(1:100,:)'];
data = [test0' test1' ...
        test3' test4' ...
        test5' test6'];



[d n] = size(data);
k = 6;
n


% [IDX,C] = kmeans(double(data)',k);
[ Cluster,x_bar ] = soltn_kmeans( 6, double(data'));
size(x_bar)

alpha = 0.1;
lass = 5;


% [Res, Num_C, Weights, Clusters, Centers] = RSK_Means(double(data), k, alpha, lass);

% imagesc(reshape(Centers(:,i),28,28)')