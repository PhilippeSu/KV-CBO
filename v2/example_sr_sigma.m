clc; clear; close all;

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultLegendFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth',2);

d = 20;

test_param.runs = 10;           % number of runs
test_param.const = 0.05;        % success if error <= const

[problem, param] = setUpClass.Ackley(d);
param.verbose = 0;
param.T = 2000;

sigma_range = [1 2 3 4 5];

avg_error = zeros(size(sigma_range,2),param.T);
success_rate = zeros(1, size(sigma_range, 2));

for sigma = sigma_range
   
    fprintf('\nsigma = %.2f\n\n', sigma)
    
    param.sigma = sigma;
    [rate, final_info] = test_kvcbo(problem, param, test_param);
    
    %% Statistics
    
    sum = zeros(1, param.T);
    for i = 1:test_param.runs
        sum = sum + final_info(i).error;
    end
    avg_error(sigma == sigma_range, :) = sum./test_param.runs;
    
    success_rate(1, sigma == sigma_range) = rate;
    
end

%% Plot

figure('Position', [0 100 800 500])
for j = 1:size(avg_error,1)
    hold on
    % semilogy(1:T, average_error(j,:));
    semilogy(1:param.T, smoothdata(avg_error(j,:)));
    set(gca, 'YScale', 'log')
    ylabel('Average Error');
    xlabel('Time Steps');  
end
legendCell = cellstr(num2str(sigma_range', '$\\sigma$=%.2f'));
legend(legendCell);
% title({'Kuramoto-Vicsek Minimizer'; 
%           sprintf('Dimension = %d, $\\alpha$ = %.2e, N = %i, $\\lambda $ = %i', d, param.alpha, param.N, param.lambda);
%           sprintf('Runs = %i, T=%.2e, dt = %.4f', test_param.runs, param.T, param.dt); 
%           },'interpreter','latex')

%%

sigmamin = min(sigma_range); 
sigmamax = max(sigma_range);

figure
plot(sigma_range, 100*success_rate, '-o')
axis([sigmamin sigmamax 0 110])
xlabel('$\sigma$');
yticks([0 20 40 60 80 100])
yticklabels({'0\%', '20\%', '40\%', '60\%', '80\%', '100\%'});
ylabel('Success Rate');

% title({'Success rate in terms of sigma';
%    sprintf('Dimension = %i, $\\alpha$ = %.2e, N = %i, $\\lambda $ = %i' , d, param.alpha, param.N, param.lambda);
%    sprintf('Runs = %i, T = %.2e, dt = %.4f', test_param.runs, param.T, param.dt);
%    },'interpreter','latex')


