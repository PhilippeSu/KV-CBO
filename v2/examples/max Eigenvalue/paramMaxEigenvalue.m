function param = paramMaxEigenvalue()

% TODO

param.N = 100;
param.Nbatch = 60;
param.Nmin = 10;
param.logicVa = 'average'; 
param.alpha = 10000;
param.mode = 'anisotropic';
param.sigma = 5;
param.dt = 0.0025;
param.lambda = 1;
param.mu = 0;
param.l = 10;
param.T = 10000;
param.epsilon = 1e-3;
param.nstall = 0;
param.epsilonstall = 0.05;
param.half = 0;

param.distri = 'uniform';

% param.distri = 'fisher';
% param.meanDir = [zeros(4,1), 1];
% param.k = 200;


% param.seed = randi([1,10000]); 
param.verbose = 1; 
param.visual = 0;
param.useGrad = 0;

end

