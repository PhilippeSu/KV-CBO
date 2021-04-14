function param = paramEF()

param.N = 500;
param.Nbatch = 50;
param.Nmin = 100;
param.alpha = 1e5;
param.logicVa = 'average'; % 'min';  
param.mode = 2;
param.sigma = 1;
param.adaptSigma = 'fix';
param.dt = 0.25;
param.lambda = 1;
param.mu = 0;
param.l = 10;
param.T = 100;
param.epsilon = 1e-5;
param.nstall = 250;
param.epsilonstall = 1e-4;
param.half = 1;

param.distri = 'uniform';

% param.distri = 'fisher';
% param.meanDir = [zeros(19,1); 1] + 0.001*randn(20,1);
% param.k = 200;

% param.seed = randi([1,10000]); 
param.verbose = 1;
param.visual = 0;
param.useGrad = 0;

end

