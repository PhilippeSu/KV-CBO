clear; clc; close all;

% Implementation of the KV-CBO method for the computation of eigenfaces

% Parameters for robsut subspace detection
p = 2;
delta = 0;
energy_functional = @Energy_Class.RobustSubspaceDetection; 
pixelx = 45;                                                % number of pixels in x direcction
pixely = 64;                                                % number of pixels in y direction

s = 'Faces.txt';
X = Generate_Points.LoadDataFromFile(s, ',');
X = bsxfun(@minus, X, mean(X,2));                           % subtract the mean
d = pixelx*pixely;

parameters = struct;
% Compute minimum by SVD
fprintf(' Compute SVD... \n\n')
if p == 2
    [Vg,D] = eig(X*X');
    [~,i] = max(diag(D));
    Vmin = Vg(:,i);
end
parameters.Recon = @Reconstruct;
parameters.trueSol = Vmin;

%% Parameters for KV-CBO method
parameters.alpha = 2e15;
parameters.sigma =  0.019;                                     % initial sigma, sigma = 0 is the deterministic Kuramoto-Vicsek model
parameters.lambda = 1;
parameters.dt = 0.25;          
parameters.N = 5000;                                        % number of particles
parameters.N2 = parameters.N;                              % number of particles to be evolved (randomly chosen)
parameters.N3 = parameters.N;                              % number of particles to compute Valpha (randomly chosen)

% Implementation parameters
parameters.T = 10000;                                        % maximal number of iterations
parameters.time_limit = 1e5;                               % time limit in seconds
parameters.epsilon = 1e-5;                                 % for stopping condition
parameters.mu = 1;                                         % parameter mu to drop particles
parameters.Nmin = 200;                                     % minimal number of particles
parameters.tau = 1.2; parameters.eta = 2;                  % parameters to adaptively change the sigma
parameters.l = 5;                                         % particles, sigma and alpha are updated every l-th iteration

%% KV-CBO method 
fprintf(' Starting KV-CBO method ...\n\n')
msg = sprintf([' Iterations: \t \t %i\n', ...
    ' max. Iterations: \t %i\n\n', ...
    ' Number of Particles: \t %i'] ...
    , 0, parameters.T, parameters.N);
fprintf(msg);

[Valpha, Iter, time, Variance, Error] = KuramotoVicsekMethod(energy_functional, d, X, parameters, p, delta);

P1 = reshape(Valpha(:,1), pixely, pixelx);
minP1 = min(P1(:));
P1 = double(P1 - minP1) ./ double( max(P1(:)) - minP1 );
imwrite(im2double(P1), 'Eigenface1.png');

P2 = reshape(-Valpha(:,1), pixely, pixelx);
minP2 = min(P2(:));
P2 = double(P2 - minP2) ./ double( max(P2(:)) - minP2 );
imwrite(im2double(P2), 'Eigenface2.png');

error = min([norm(Valpha-Vmin,2), norm(Valpha+Vmin,2)]);
fprintf('\n Final Error: \t \t %.2e\n', error)
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


