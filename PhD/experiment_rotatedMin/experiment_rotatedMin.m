%%

clc; clear; close all;

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultLegendFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth',2)

d = 5;
runs = 10;          

angles = [0, pi/4];
kparam = [0, 200];

params = params_rotatedMin(d);   

errors = zeros(runs, size(kparam,2), size(angles,2));

for phi = angles  
    for k = kparam
        
        fprintf('phi = %.2f, k = %.2f', phi, k)
        
        costfunction = RastriginCostFunction_rotatedMin(d, phi);
        
        if k == 0
            params.initData.method = 1;
            params.initData.name = 'uniform';
        else
            params.initData.method = 2;
            params.initData.name = 'von Mises-Fisher';
            params.initData.k = k;
            params.initData.meanDir = costfunction.minimizer;
        end

        for i =1:runs      
            [Va, info] = KVCBO(costfunction, params);
            info(end).error
            errors(i, k==kparam, phi==angles) = info(end).error;
        end
    end
end

rate = errors < 0.05;
rate = mean(rate)
std_dev = std(errors)
min_err = min(errors)
avg_err = mean(errors)
max_err = max(errors)
