function param = paramHaystack()

param.N = 100;
param.Nbatch = 50;
param.Nmin = 10;
param.alpha = 50000;
param.logicVa = 'average'; 
param.mode = 'anisotropic';
param.sigma = 1;
param.dt = 0.05;
param.lambda = 1;
param.mu = 0;
param.l = 10;
param.T = 1000;
param.epsilon = 1e-5;
param.nstall = 0;
param.half = 0;

param.distri = 'uniform';

% param.distri = 'fisher';
% param.meanDir = [zeros(19,1); 1] + 0.001*randn(20,1);
% param.k = 200;

% param.seed = randi([1,10000]); 
param.verbose = 1;
param.visual = 0;
param.useGrad = 1;

end

