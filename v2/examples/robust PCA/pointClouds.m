clear; clc; close all;

d = 2;

X = pointCloudClass.GaussianMixture(d);

scatter(X(:,1), X(:,2))
axis([-10 10 -10 10])
