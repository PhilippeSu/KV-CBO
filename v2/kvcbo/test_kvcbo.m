function [rate, final_info] = test_kvcbo(problem, param, test_param)

    runs = test_param.runs;
    const = test_param.const;

    i = 0;
    count_success = 0;

    stats = savestats();
    final_info(1) = stats;
    final_info(runs).iter = [];

    %% Main Loop

    while true

        [~, info] = KVCBO(problem, param);

        i = i + 1;
        stats = savestats();
        final_info(i) = stats;

        error = final_info(i).error(end);
        final_info(i).error(end:param.T) = error;
        if error < const
            count_success = count_success + 1;
        end
        
        verbose_test_kvcbo(final_info, info, test_param, i)

        if i >= runs
            break
        end

    end

    rate = count_success/runs;

    %% Statistics

     function stats = savestats()

         stats.iter = i;
         if i == 0
             stats.N = [];
             stats.variance = [];
             stats.consensus = [];
             stats.distVa = [];
             stats.gradnorm = [];
             stats.cost = [];
             stats.error = [];
             % stats.stopCond = [];
         else
             siz = size(info,2);

             for j=1:siz
                 stats.N(j) = info(j).N;
                 stats.variance(j) = info(j).variance;
                 stats.consensus(j) = info(j).consensus;
                 stats.distVa(j) = info(j).distVa;
                 stats.gradnorm(j) = info(j).gradnorm;
                 stats.cost(j) = info(j).cost;
                 stats.error(j) = info(j).error;
                 % stats.stopCond(j) = info(j).stopCond;
             end
         end

     end

end
