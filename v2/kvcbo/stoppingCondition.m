function [stop,reason] = stoppingCondition(problem, param, info, last)

    persistent counter_distVa counter_cost;

    if isempty(counter_distVa)
       counter_distVa = 0; 
    end

    if isempty(counter_cost)
       counter_cost = 0; 
    end

    stop = 0;
    reason = ' ';
    stats = info(last);

    %%
    if stats.iter >= param.T
        reason = sprintf('Maximum number of iterations reached');
        stop = 1;
    end

    %%
    if stats.iter > 0 && stats.consensus < param.epsilon
        reason = sprintf('Consensus reached');
        stop = 1;
    end

    %%
    if param.nstall > 0
        if stats.distVa < param.epsilonstall
            counter_distVa = counter_distVa + 1;
        else
            counter_distVa = 0;
        end
        if counter_distVa >= param.nstall
            reason = sprintf('Change in Va sufficently small');
            stop = 1;
        end
    end

    %%
    % if isfield(problem, 'min')
    %     min = problem.min;
    %     
    %     if abs(stats.cost - min) < 1e-5 
    %         counter_cost = counter_cost + 1;
    %     else
    %         counter_cost = 0;
    %     end
    %     if counter_cost >= param.nstall
    %         reason = sprintf('Change in cost sufficently small');
    %         stop = 1;
    %         % return
    %     end
    % end

end