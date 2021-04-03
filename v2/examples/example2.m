clear; clc; close all;

%%
% Numerical tests for the success rate of various functions (e.g. Ackley,
% Rastrigin, ...)

d = 20;     % euclidean dimension (>=2)

[problem, param] = setUpClass.Ackley(d);
param.verbose = 0;

test_param.runs = 10;          % number of runs
test_param.const = 0.05;        % success if error <= const

[rate, final_info] = test_kvcbo(problem, param, test_param);

%% average error

sum = 0;
for i = 1:size(final_info,2)
   sum = sum + final_info(i).error(end); 
end
avg_err = sum/test_param.runs;

fprintf('Rate: %.2f, Average error: %.2e\n', rate, avg_err);

