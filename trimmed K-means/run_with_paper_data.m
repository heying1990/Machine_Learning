clear;
clc;

n=300;
data1 = zeros(5, 100);
data2 = zeros(5, 100);
data3 = zeros(5, 100);

mu1 = [-4; -4; 0; 0; 0];
mu2 = [0; 0; 0; 0; 0];
mu3 = [4; 4; 0; 0; 0];
sigma = eye(5);

for i=1:100
    data1(:,i) = mvnrnd(mu1, sigma)';
    data2(:,i) = mvnrnd(mu2, sigma)';
    data3(:,i) = mvnrnd(mu3, sigma)';
end

data3(4:5,1) = [1000; 1000];
data3(4:5,2) = [1000; 1000];
data3(4:5,3) = [1000; 1000];


xs = [data1(1,:) data2(1,:) data3(1,:)];
ys = [data1(2,:) data2(2,:) data3(2,:)];


d = 5;
k = 3;
alpha = 0.1;
lass = 1.25;
data = [data1 data2 data3];


[Res, Num_C, Weights] = RSK_Means(data, k, alpha, lass);

hold on;
for i=1:k
    scatter(Res(1,1:Num_C(i),i), Res(2,1:Num_C(i),i));
end
hold off;