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
runs = 20;

% define the problem
costfunction = RastriginCostFunction_isoVsAniso(d);
params = params_isoVsAniso(d);

s_range = 0:0.25:5;
t_range = 0.05:.05:.25;
mode = {'isotropic','anisotropic'};

for j = 1:2
    params.mode = mode{j};
    sprintf(' %s \n ************************** \n\n', params.mode);
    
    mat = zeros(size(t_range,2), size(s_range,2));
    
    for s = s_range    
        for t = t_range
            fprintf('*******************\n')
            fprintf('t = %.2f, s = %.2f\n', t, s)

            params.sigma = s;
            params.dt = t;

            success = 0;
            for i = 1:runs
                [Va, info] = KVCBO(costfunction, params);
                error = info(end).error;
                success = success + (error < 0.05); 
            end
                
            fprintf('success rate = %.2f\n\n', success/ runs)
            mat(t == t_range, s == s_range) = success;
        end
    end

    mat = mat./runs;

    %%
    imagesc([min(s_range) max(s_range)], [min(t_range) max(t_range)], mat)
    caxis([0 1]) 
    xlabel('$\sigma$')
    ylabel('$\Delta t$')
    title(params.mode)
    set(gca,'XTick', s_range);
    set(gca,'YTick', t_range);
    colorbar();
    xtickangle(90)
    print(params.mode, '-depsc')

end

