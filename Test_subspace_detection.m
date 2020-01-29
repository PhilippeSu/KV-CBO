clear; clc; close all;

% Implementation of the KV-CBO method for the robust subspace detection problem

% Dimension
d = 100;                                        

%% Parameters for subspace detection energy
p = 2;
delta = 0;
energy_functional = @Energy_Class.RobustSubspaceDetection; 

%% Parameters for KV-CBO method
parameters = struct;
parameters.alpha = 2000;
parameters.alphamax = 1e15;                                % comment this ligne out to have a constant alpha
parameters.sigma =  1;                                     % initial sigma, sigma = 0 is the deterministic Kuramoto-Vicsek model
parameters.lambda = 1;
parameters.dt = 0.25;          
parameters.N = 100;                                        % number of particles
parameters.N2 = parameters.N;                              % number of particles to be evolved (randomly chosen)
parameters.N3 = parameters.N;                              % number of particles to compute Valpha (randomly chosen)  

% Implementation parameters
parameters.T = 2000;                                       % maximal number of iterations
parameters.time_limit = 10;                                % time limit in seconds
parameters.epsilon = 1e-5;                                 % for stopping condition
parameters.mu = 1;                                         % parameter mu to drop particles
parameters.Nmin = 100;                                     % minimal number of particles
parameters.tau = 1.2; parameters.eta = 2;                  % parameters to adaptively change the sigma
parameters.l = 25;                                         % particles, sigma and alpha are updated every l-th iteration

% Parameters for point cloud
Nsp = 5;                                                   % number of subspaces
M = 100;                                                   % number of points in each subspace
ds = 1;                                                    % dimension of each subspace
noise = 0.01;                                              % Gaussian noise
MOutliers = 0;                                             % number of outliers

%% Generate Points
X = Generate_Points.NearlyParallel(M, d, ds, Nsp, noise, MOutliers);
% X = Generate_Points.RandomlyChosen(M, d, ds, Nsp, noise, MOutliers);
% X = Generate_Points.LoadDataFromFile(s, ',');
X = bsxfun(@minus, X, mean(X,2));                          % subtract the mean 

% Compute minimum by SVD
if p == 2
    [Vg,D] = eig(X*X');
    [~,i] = max(diag(D));
    Vmin = Vg(:,i);
end
parameters.Recon = @Reconstruct;
parameters.trueSol = Vmin;

%% KV-CBO method 
[Valpha, Iter, time, Variance, Error] = KuramotoVicsekMethod(energy_functional, d, X, parameters, p, delta);

error = min([norm(Valpha-Vmin,2), norm(Valpha+Vmin,2)]);
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

function v = Reconstruct(Valpha, varargin)
  v = Valpha;
end


