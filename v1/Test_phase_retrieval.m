clear; clc; close all;

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultLegendFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth',2)

% Implementation of the KV-CBO method for the phase retrieval problem

d = 32;                                         

%% Parameters for phase retrieval
energy_functional = @Energy_Class.PhaseRetrieval; 

% Set up of the problem
M = 6*d;                                        % number of frame vectors
X = Generate_Points.UniformOnSphere(d, M);      % columns of X form a frame
A = LowerFrameBound(X);                         % optimal lower frame bound
z = randn(d,1);                                 % vector to recover
Y = ((dot(z*ones(1,M),X)).^2)';                 % vector in R^M
C = sqrt(norm(Y,1)/A);
X = [X; zeros(1,M)];
parameters = struct;
parameters.trueSol = z;
parameters.Recon = @Reconstruct;

%% Parameters for KV-CBO method
parameters.alpha = 1e15;
parameters.sigma =  0.2;                                   % initial sigma, sigma = 0 is the deterministic Kuramoto-Vicsek model
parameters.lambda = 1;
parameters.dt = 0.10;          
parameters.N = 5000;                                       % number of particles
parameters.N2 = parameters.N;                              % number of particles to be evolved (randomly chosen)
parameters.N3 = parameters.N;                              % number of particles to compute Valpha (randomly chosen)

% Implementation parameters
parameters.T = 1000;                                        % maximal number of iterations
parameters.time_limit = 15;                                % time limit in seconds
parameters.epsilon = 1e-5;                                 % for stopping condition
parameters.mu = 1;                                         % parameter mu to drop particles
parameters.Nmin = 1000;                                    % minimal number of particles
parameters.tau = 1.2; parameters.eta = 2;                  % parameters to adaptively change the sigma
parameters.l = 25;                                         % particles, sigma and alpha are updated every l-th iteration

%% KV-CBO method 
[Valpha, Iter, time, Variance, Error] = KuramotoVicsekMethod(energy_functional, d+1, X, parameters, Y, C);

% Reconstruct z from Valpha
z_KV = C*Valpha(1:end-1,1);
error = min([norm(z_KV-z,2), norm(z_KV+z,2)]);
fprintf('\n\n Final Error: \t \t %.2e\n', error)
fprintf(' Runtime: \t \t %.2f\n', Iter*time);

% Plot empirical variance and error
figure
semilogy(1:Iter-1, Error(1:Iter-1))
set(gca, 'YScale', 'log')
ax = gca;
ax.YColor = [0, 0.4470, 0.7410];
ylabel('Error')
title({'KV-CBO method'; 'Error and Variance'})
hold on
yyaxis right
ax = gca;
ax.YColor = [0.25, 0.25, 0.25];
plot(1:Iter-1,Variance(1:Iter-1))
ylabel('Variance')
legend('Error', 'Variance', 'location', 'southeast')
xlabel('Time Steps')
axis([0 Iter-1 min(Variance(1:Iter-1)) max(Variance(1:Iter-1))])
box off

function v = Reconstruct(Valpha, C)
  v = C*Valpha(1:end-1);
end