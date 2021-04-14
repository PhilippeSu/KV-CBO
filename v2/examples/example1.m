clear; clc; close all;

%%
% Numerical tests for various test function (e.g. Ackley, Rastrigin, ...)

d = 20;  % euclidean dimension (>= 2)

% define the problem
[problem, param] = setUpClass.Ackley(d);

% call KV-CBO
[Va, info] = KVCBO(problem, param);

%% 
figure;
subplot(1,3,1)
semilogy([info.iter], [info.cost], '.-')
xlabel('Iterations')
ylabel('Cost')

subplot(1,3,2)
semilogy([info.iter], [info.error], '.-');
xlabel('Iterations')
ylabel('Error')

subplot(1,3,3)
semilogy([info.iter], [info.variance], '.-');
xlabel('Iterations')
ylabel('Variance')
