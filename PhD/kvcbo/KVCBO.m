function [Va, info] = KVCBO(costfunction, KVCBOparam)    

     if KVCBOparam.out.fixseed
        rng(KVCBOparam.out.seed) 
     end

    N = KVCBOparam.N;
    l = KVCBOparam.adapt.l;
    nT = KVCBOparam.stopCond.nT;

    iter = 0;
    count_stall = 0;
    
    stats = savestats();
    info(1) = stats;
    info(nT).iter = []; 
    
    if KVCBOparam.out.visual == 1
         figure('Position', [10 500 650 400]);
        % figure('Position', [10 500 1400 400]);
        colormap parula;
        [VE, ve] = outputClass.create_mesh(costfunction.d);
    end
    
    V = initialParticles(costfunction, KVCBOparam);    
    Va = computeVa(costfunction, KVCBOparam, V);
    
    % track n randomly selected particles
    n = 0;
    Ntr = randperm(N, n);
    Vtr = V(:, Ntr)';
    
    %% Main Loop  
    
    while true       
        
        Va_prev = Va;
        [Va, V] = KuramotoIteration(costfunction, KVCBOparam, V);        
        Vtr = [Vtr V(:, Ntr)'];
        
        % Update Statistics
        iter = iter + 1;
        stats = savestats();
        info(iter) = stats; 
        
        % check stopping Condition
        [stop, reason, count_stall] = ...
            stoppingCondition(costfunction, KVCBOparam, info, iter, count_stall); 
        
        if stop      
            info(iter).stopCond = reason;
            break
        end
        
        % Adapt Parameters, discard particles, use gradient 
        if mod(iter,l) == 0
            [KVCBOparam, V] = postProcessing(costfunction, KVCBOparam, V, info, iter);
            
            % output
            if KVCBOparam.out.verbose == 1
                outputClass.verbose(info, KVCBOparam, costfunction, iter)
            end

            % Visualization
            if KVCBOparam.out.visual == 1
               outputClass.visualizeCost(costfunction, VE, ve, V, Va, Vtr)
            end
        end

    end
    
    info(iter+1:nT) = [];
    
    %% Statistics
    
    function stats = savestats()
        stats.iter = iter;
        if iter == 0
            stats.N = N;
            stats.variance = 0;
            stats.consensus = 0;
            stats.changeVa = 0;
            stats.gradnorm = 0;
            stats.cost = 0;
            stats.error = 0;
            stats.stopCond = ' ';
            stats.sigma = 0;
        elseif size(Va,2) == 1
            % one consensus point (standard)
            
            stats.N = KVCBOparam.N;
            stats.variance = sum(vecnorm(V - mean(V,2)*ones(1,stats.N),2).^2)/stats.N;
            stats.consensus = sum(vecnorm(V - Va*ones(1,stats.N),2))/stats.N;
            stats.changeVa = min(norm(Va - Va_prev), norm(Va + Va_prev));
            % grad = gradient(costfunction, Va);   % gradient projected onto the manifold       
            % stats.gradnorm = norm(grad); 
            stats.gradnorm = 0;
            stats.cost = costfunction.cost(Va);
            stats.error = costfunction.error(Va, KVCBOparam);
            stats.stopCond = ' ';
            stats.sigma = KVCBOparam.sigma;
        else
            % multiple consensus points
            
            stats.N = KVCBOparam.N;
            stats.variance = sum(vecnorm(V - mean(V,2)*ones(1,stats.N),2).^2)/stats.N;
            stats.consensus = 0;
            stats.changeVa = 0;
            stats.gradnorm = 0;
            stats.cost = costfunction.cost(Va(:,1));   
            stats.error = costfunction.error(Va(:,1), KVCBOparam);
            stats.stopCond = ' ';
            stats.sigma = KVCBOparam.sigma;
        end
    end

end

