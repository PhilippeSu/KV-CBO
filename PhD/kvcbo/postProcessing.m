%%
function [KVCBOparam, V] = postProcessing(costfunction, KVCBOparam, V, info, iter)

    l = KVCBOparam.adapt.l;
    N = KVCBOparam.N;
    Nmin = KVCBOparam.Nmin;
    mu = KVCBOparam.mu;

    %% discard particles
    var_curr = info(iter).variance;
    var_prev = info(iter - l + 1).variance;

    Nnew = fix(N*(1+mu*((var_curr-var_prev)/var_prev)));
    Nnew = max(Nnew, Nmin);
    Nnew = min(Nnew, N);

    if Nnew ~= N
      I = randperm(N,Nnew);
      V = V(:,I);
      KVCBOparam.N = Nnew;
    end
    
    %% addapt parameters
    switch KVCBOparam.adapt.methodForSigma
        case 1
            KVCBOparam.sigma = max(KVCBOparam.sigma/KVCBOparam.adapt.tau, KVCBOparam.adapt.sigmamin);
        case 2
            KVCBOparam.sigma = ...
                max(KVCBOparam.sigma/(KVCBOparam.adapt.sigma0*log(iter+1)), KVCBOparam.adapt.sigmamin);
    end
    
    switch KVCBOparam.adapt.methodForAlpha
        case 1  
            KVCBOparam.alpha = (iter+1)*KVCBOparam.adapt.alphamax/KVCBOparam.nT;
    end
    
    %% use gradient descent step
     if KVCBOparam.adapt.useGrad
        V = GradientDescent(costfunction, KVCBOparam, V);
     end

end
