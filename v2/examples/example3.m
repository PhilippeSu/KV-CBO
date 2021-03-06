clear; clc; close all;

%% 
% Numerical test for robust PCA, where the point cloud is generated by the
% Haystack model (see setUpClass.m)

d = 100;  % euclidean dimension (>= 2)

%% 
sum_error = 0;

for i = 1:10

    [problem, param] = setUpClass.Haystack(d);
    param.verbose = 0;
    
    [Va, info] = KVCBO(problem, param);
    
    error = info(end).error;
    sum_error = sum_error + error;
    fprintf('%i. Error: %.2e\n', i, error);
    
end

fprintf('\nAverage Error: %.2e\n', sum_error/i);
