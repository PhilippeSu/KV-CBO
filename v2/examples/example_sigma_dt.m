clear; clc; close all;

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultLegendFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth',2)

%%
d = 5;  % euclidean dimension (>= 2)
runs = 10;

% define the problem
[problem, param] = setUpClass.Ackley(d);
param.verbose = 0;

s_range = 0:0.5:5;
t_range = .05:.1:.25;

mat = zeros(size(t_range,2), size(s_range,2));

for s = s_range
    
    for t = t_range
        
        fprintf('t = %.2f, s = %.2f\n', t, s)
        
        param.sigma = s;
        param.dt = t;
        
        success = 0;
        
        for i = 1:runs
            [Va, info] = KVCBO(problem, param);
            error = info(end).error;
            success = success + (error < 0.05); 
        end
        
        mat(t == t_range, s == s_range) = success;
    end

end

mat = mat./runs;

%%
imagesc([min(s_range) max(s_range)], [min(t_range) max(t_range)], mat)
xlabel('$\sigma$')
ylabel('$\Delta t$')
colorbar
set(gca,'XTick', s_range);
set(gca,'YTick', t_range);

