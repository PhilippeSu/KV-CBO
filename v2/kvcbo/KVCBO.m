function [Va, info] = KVCBO(problem, param)    

    % clear functions with persistent variables
    clear stoppingCondition visualize

     if isfield(param, 'seed')
        rng(param.seed, 'twister') 
     end

    N = param.N;
    l = param.l;

    iter = 0;
    
    stats = savestats();
    info(1) = stats;
    info(param.T).iter = []; 
    
    V = distributionClass.initialParticles(problem, param);    
    Va = computeVa(problem, param, V);
    
    %% Main Loop
    
    while true       
        
        Va_prev = Va;
        [Va, V] = KuramotoIteration(problem, param, V);
        
        if and(param.useGrad == 1, mod(iter, l) == 0)
            [Va, V] = GradientDescent(problem, param, V);
        end
        
        % Update Statistics
        iter = iter + 1;
        stats = savestats();
        info(iter) = stats; 
        
        [stop, reason] = stoppingCondition(problem, param, info, iter); 
        if stop      
            info(iter).stopCond = reason;
            verbose(info, param, problem, iter)
            break
        end
        
        % Adapt Parameters  
        if mod(iter,l) == 0
            [param, V] = adaptParameters(param, V, info, iter);
            verbose(info, param, problem, iter)
        end
        
        % Visualization
        if param.visual == 1
           visualize(V, Va, problem)
        end

    end
    
    info(iter+1:param.T) = [];
    
    %% Statistics
    
    function stats = savestats()
        stats.iter = iter;
        if iter == 0
            stats.N = N;
            stats.variance = 0;
            stats.consensus = 0;
            stats.distVa = 0;
            stats.gradnorm = 0;
            stats.cost = 0;
            stats.error = 0;
            stats.stopCond = ' ';
            stats.sigma = 0;
        else
            stats.N = param.N;
            stats.variance = sum(vecnorm(V - mean(V,2)*ones(1,stats.N),2).^2)/stats.N;
            stats.consensus = sum(vecnorm(V - Va*ones(1,stats.N),2))/stats.N;
            stats.distVa = min(norm(Va - Va_prev), norm(Va + Va_prev));
            grad_euc = problem.grad(Va);
            grad_sph = (eye(problem.d) - Va*Va')*grad_euc;
            stats.gradnorm = norm(grad_sph);
            stats.cost = problem.cost(Va);
            stats.error = problem.error(Va);
            stats.stopCond = ' ';
            stats.sigma = param.sigma;
        end
    end

end

