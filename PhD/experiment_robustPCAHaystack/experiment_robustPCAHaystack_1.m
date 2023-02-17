%%

clc; clear; close all;

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultLegendFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth',2)

d = 100;
runs = 20;

p = 2;

error_KVCBO = zeros(1,runs);
error_SVD = zeros(1,runs);
error_FMS = zeros(1,runs);
error_GD = zeros(1,runs);

for i = 1:runs

    fprintf('\nRun %i\n', i)
    
    %% Generate point cloud (no contamination)

    P = 200;                     % data points 
    outlier_per = 0;
    cloud_param.sigma_in = 1;
    cloud_param.sigma_out = 0;
    cloud_param.ds = 1;           % dimension of subspace
    cloud_param.N_in = round(P*(1-outlier_per),0);
    cloud_param.N_out = round(P*outlier_per,0);
    cloud_param.noise = 0;
    [X, W] = pointCloudClass.Haystack(d, cloud_param);

    %% KVCBO

    delta = 0;
    costfunction = CostFunction_robustPCAHaystack(X, p, delta, 0, W, d);
    param = params_robustPCAHaystack(d);

    [Va, info] = KVCBO(costfunction, param);

    error_KVCBO(1,i) = min([norm(W-Va), norm(W+Va)]);
    fprintf('KV-CBO error: %.2e, number of iterations: %i\n', error_KVCBO(1,i), info(end).iter)

    %% FMS

    if and(exist('fms.m', 'file'), exist('calc_sdist.m', 'file'))

        options=struct;
        options.scaleopt='normal';
        options.p=p;
        options.initopt='random';
        options.svdopt='normal';
        options.maxiter=100;
        options.epsilon=10^-10;

        L=fms(X',d,options);

        error_FMS(1,i) = min([norm(W-L(:,1)), norm(W+L(:,1))]);
        fprintf('FMS error: %.2e\n', error_FMS(1,i))

    end
    
    %% SVD

    [U,~,~] = svd(X);
    error_SVD(1,i) = min([norm(W-U(:,1)), norm(W+U(:,1))]);
    fprintf('SVD error: %.2e\n', error_SVD(1,i))

    %% GD

    x = randn(d,1);
    x = x/norm(x);
    eta = 0.001;

    for j=1:param.stopCond.nT
        grad_sph = gradient(costfunction, x);
        x = x - eta*grad_sph;
        x = x/norm(x);

        if norm(grad_sph) < 1e-3
            break
        end
    end

    error_GD(1,i) = min([norm(W-x), norm(W+x)]);
    fprintf('GD error: %.2e, number of iterations: %i\n \n', error_GD(1,i), j)

end

fprintf('\n KV-CBO\n')
rate = sum(error_KVCBO < 0.05)/runs;
avg = mean(error_KVCBO);
std_dev = std(error_KVCBO);
fprintf('success rate: %.2f, avg. error: %.2e, std. dev.: %.2e\n', rate, avg, std_dev)

fprintf('\n FMS\n')
rate = sum(error_FMS < 0.05)/runs;
avg = mean(error_FMS);
std_dev = std(error_FMS);
fprintf('success rate: %.2f, avg. error: %.2e, std. dev.: %.2e\n', rate, avg, std_dev)

fprintf('\n SVD\n')
rate = sum(error_SVD < 0.05)/runs;
avg = mean(error_SVD);
std_dev = std(error_SVD);
fprintf('success rate: %.2f, avg. error: %.2e, std. dev.: %.2e\n', rate, avg, std_dev)

fprintf('\n GD\n')
rate = sum(error_GD < 0.05)/runs;
avg = mean(error_GD);
std_dev = std(error_GD);
fprintf('success rate: %.2f, avg. error: %.2e, std. dev.: %.2e\n', rate, avg, std_dev)

