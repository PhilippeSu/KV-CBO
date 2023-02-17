%%
clc; clear; close all;

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultLegendFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth',2);

d = 20;
runs = 100;          

costfunction = RastriginCostFunction_Error_Time(d);
params = params_Error_Time(d);

nT = params.stopCond.nT;

best_final_error = 2;
worst_final_error = 0;

best_run_errors = zeros(1, nT);
worst_run_errors = zeros(1, nT);
errors = zeros(runs, nT);
   
for i = 1:runs
    
    [Va, info] = KVCBO(costfunction, params);
    iter = info(end).iter;
    error = info(iter).error;

    for j = 1:nT
       if j <= iter
           errors(i,j) = errors(i,j) + info(j).error; 
       else
           errors(i,j) = errors(i,j) + error; 
       end
    end  

    if error < best_final_error
       best_run_errors = errors(i,:);
       best_final_error = error;
    end

    if error > worst_final_error
       worst_run_errors = errors(i,:);
       worst_final_error = error;
    end

    if error < 0.05
        fprintf('%i. Success, Error = %.2e\n', i, error)
    else
        fprintf('%i. Fail, Error = %.2e\n', i, error)
    end

end
       
%% Plot
 
figure
semilogy(1:nT, smoothdata(best_run_errors))                            
hold on
semilogy(1:nT, smoothdata(worst_run_errors))
hold on
semilogy(1:nT, smoothdata(mean(errors))) 
grid on
hold on
set(gca, 'YScale', 'log')
ylabel('Error');
xlabel('Time Steps');  
