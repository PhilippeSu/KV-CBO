%%
function verbose_test_kvcbo(final_info, info, test_param, i)

    const = test_param.const;

    error = final_info(i).error(end);
    consensus = final_info(i).consensus(end);
    Navg = mean([info.N]);

    if error < const
        fprintf('%i. Success, Error = %.2e, Consensus = %.2e, N_avg = %.2f, Iter = %i\n' ...
            , i, error, consensus, Navg, size(info,2))
    else
        fprintf('%i. Fail, Error = %.2e, Consensus = %.2e, N_avg = %.2f, Iter = %i\n' ...
            , i, error, consensus, Navg, size(info,2))
    end
       
%     figure;
%     subplot(1,3,1)
%     semilogy([info.iter], [info.cost], '.-')
%     xlabel('Iterations')
%     ylabel('Cost')
%     
%     subplot(1,3,2)
%     semilogy([info.iter], [info.error], '.-');
%     xlabel('Iterations')
%     ylabel('Error')
%     
%     subplot(1,3,3)
%     semilogy([info.iter], [info.gradnorm], '.-');
%     xlabel('Iterations')
%     ylabel('Gradnorm')

end