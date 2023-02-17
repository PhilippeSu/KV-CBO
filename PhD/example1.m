clear; clc; close all;

%%

% Numerical tests for various test function:
% 
% 1. Rastrigin, 2. Ackley, 3. Alpine, 4. Schaffer, 5. Solomon, 6. Lévi,
% 7. XSY random, 8. Griewank

d = 3;  % euclidean dimension (>= 2)

[costfunction, KVCBOparam] = setUpClass.Ackley(d);

[Va, info] = KVCBO(costfunction, KVCBOparam);

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
