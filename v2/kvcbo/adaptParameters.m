%%
function [param, V] = adaptParameters(param, V, info, last)

    l = param.l;
    N = param.N;
    Nmin = param.Nmin;
    mu = param.mu;

    % drop particles
    var_curr = info(last).variance;
    var_prev = info(last - l + 1).variance;

    Nnew = fix(N*(1+mu*((var_curr-var_prev)/var_prev)));
    Nnew = max(Nnew, Nmin);
    Nnew = min(Nnew, N);

    if Nnew ~= N
      I = randperm(N,Nnew);
      V = V(:,I);
      param.N = Nnew;
    end
    
    if strcmp(param.adaptSigma, 'method1') 
        param.sigma = param.sigma / 1.05;
    elseif strcmp(param.adaptSigma, 'method2')
        param.sigma = 1/log(last);
    end

end
