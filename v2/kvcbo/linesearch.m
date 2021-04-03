%%
function [stepsize, newx_pos] = linesearch(problem, x, desc_dir)

    % Linesearch method for the step size
    
    f0 = problem.cost(x);
    df0 = -norm(desc_dir)^2;
    
    norm_d = norm(desc_dir, 2); % descent direction is not normalized
    gamma_pos = 1/norm_d;
    
    cost_evaluations = 1;
    suff_decr = 1e-4; 
        
    Armijo_cond = 0;
    
    % Backtracking
    while ~Armijo_cond 
        
        % GD with new step size
        newx_pos = x + gamma_pos*desc_dir;
        newx_pos = newx_pos/norm(newx_pos, 2);
        newf = problem.cost(newx_pos);
        
        Armijo_cond = newf <= f0 + suff_decr*gamma_pos*df0;
        
        cost_evaluations = cost_evaluations + 1;
        if cost_evaluations >= 25
            break;
        end
        
        % Reduce the step size,
        gamma_pos = 0.5 * gamma_pos;
        
    end
    
    %%
    % If we got here without obtaining a decrease, we reject the step.
    if newf > f0
        gamma_pos = 0;
        newx_pos = x;
    end

     stepsize = gamma_pos * norm_d;
    
end
