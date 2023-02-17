function [stop,reason,counter_stall] = stoppingCondition(costfunction, KVCBOparam, info, iter, counter_stall)

    stop = 0;
    reason = ' ';
    stats = info(iter);

    %% 
    if stats.iter >= KVCBOparam.stopCond.nT
        reason = sprintf('Maximum number of iterations reached');
        stop = 1;
    end

    %%
    switch KVCBOparam.stopCond.method
        case 1 % consensus reached
            if stats.iter > 0 && stats.consensus < KVCBOparam.stopCond.epsilon
                reason = sprintf(KVCBOparam.stopCond.name);
                stop = 1;
            end
            
        case 2 % change in Va sufficiently small
            if stats.changeVa < KVCBOparam.stopCond.epsilonstall
                counter_stall = counter_stall + 1;
            else
                counter_stall = 0;
            end
            if counter_stall >= KVCBOparam.stopCond.nstall
               reason = sprintf(KVCBOparam.stopCond.name);
               stop = 1;
            end
            
         case 3 % change in E(Va) sufficiently small
            if abs(stats.cost - costfunction.minimizer) < KVCBOparam.stopCond.epsilonstall 
                counter_stall = counter_stall + 1;
            else
                counter_stall = 0;
            end
            if counter_stall >= KVCBOparam.stopCond.nstall
               reason = sprintf(KVCBOparam.stopCond.name);
               stop = 1;
            end
            
        case 4 % variance sufficiently small
            if stats.iter > 0 && stats.variance < KVCBOparam.stopCond.epsilon
                reason = sprintf(KVCBOparam.stopCond.name);
                stop = 1;
            end
    end

end
