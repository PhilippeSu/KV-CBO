function param = params_rotatedMin(d)
 
param.N = 50;
param.M = 50; 
param.Nmin = 100; 
param.mu = 0;

param.alpha = 1e6; 
param.sigma = 5;
param.dt = 0.05;
param.mode = 'anisotropic';
param.lambda = 1;

param.logicVa = setlogicVa(d); 
param.initData = initialData(d);
param.stopCond = stoppingCondition();
param.out = output_param();
param.adapt = adapt_param();

end

%%

function logicVa = setlogicVa(d)

logicVa.method = 1;

switch logicVa.method
    case 1
        logicVa.name = 'weighted average';
    case 2
        logicVa.name = 'greedy'; 
    case 3
        logicVa.name = 'Cucker-Smale (multiple consensus points)';
        logicVa.gamma = 2;
        logicVa.C = sqrt(d);
        logicVa.beta = 50;
        logicVa.p = 2;
    case 4
        logicVa.name = 'kNN (multiple consensus points)';
        logicVa.k = 10;
        logicVa.p = 2;
end

end

function initData = initialData(d)

initData.name = 'uniform';

% switch method
%     case 1
%         initData.name = 'uniform';
%     case 2
%         initData.name = 'von Mises-Fisher';
% %         initData.meanDir = [zeros(d-1,1); 1];
% %         initData.k = 100; 
%     case 3
%         initData.name = 'uniform over half sphere';
%     otherwise
% end

end

function stopCond = stoppingCondition()

stopCond.nT = 5000;    % maximum number of iterations
stopCond.method = 1;

switch stopCond.method
    case 1
        stopCond.name = 'consensus reached';
        stopCond.epsilon = 1e-4;
    case 2
        stopCond.name = 'change in Va sufficiently small';
        stopCond.nstall = 100;
        stopCond.epsilonstall = 0.05;
    case 3
        stopCond.name = 'change in E(Va) sufficiently small';
        stopCond.nstall = 100;
        stopCond.epsilonstall = 0.05;
    case 4
        stopCond.name = 'variance sufficiently small';
        stopCond.epsilon = 1e-3;
    otherwise
        stopCond.name = 'Maximum number of iterations reached';
end

end

function adapt = adapt_param()

adapt.l = 10;    
adapt.methodForSigma = 0;
adapt.methodForAlpha = 0;
adapt.useGrad = 0;

switch adapt.methodForSigma
    case 1
        adapt.tau = 1.05;
        adapt.sigmamin = 0.01;
    case 2
        adapt.sigma0 = 10;
        adapt.sigmamin = 0.01;
end

switch adapt.methodForAlpha
    case 1
        adapt.alphamax = 1e8;
end

end

function out = output_param()

out.fixseed = 0;
if out.fixseed == 1
    out.seed = 1234;
end

out.verbose = 0;
out.visual = 0;
out.p = 2;  % power for error norm, e.g. p = 'inf'

end
