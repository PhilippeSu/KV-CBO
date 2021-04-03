function param = paramPhaseRetrieval()

param.N = 200;
param.Nbatch = 100;
param.Nmin = 100;
param.alpha = 1e4;
param.logicVa = 'min'; % 'average';  
param.mode = 2;
param.sigma = 1;
param.dt = 0.5;
param.lambda = 1;
param.mu = 0;
param.l = 10;
param.T = 2000;
param.epsilon = 1e-8;
param.nstall = 0;
param.epsilonstall = 0.05;
param.distri = 'uniform'; % 'fisher'; 
param.half = 1; 

% param.seed = randi([1,10000]); 
param.verbose = 1;
param.visual = 0;
param.useGrad = 0;

end

