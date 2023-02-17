%%
function [Va, V] = KuramotoIteration(costfunction, KVCBOparam, V)
        
     M = KVCBOparam.M;
     N = KVCBOparam.N;
     sigma = KVCBOparam.sigma;
     lambda = KVCBOparam.lambda;
     mode = KVCBOparam.mode;
     dt = KVCBOparam.dt;
     d = costfunction.d;
     
     if M < N
         sp = randperm(N,M);
         Vbatch = V(:,sp);
     else
         M = N;
         Vbatch = V;
     end
     
     Va = computeVa(costfunction, KVCBOparam, Vbatch); 
      
     %% 
     if strcmp(mode, 'isotropic') 
         
         % Vectorialization
         Nalphi = ones(d,1)*vecnorm(Vbatch-Va,2);
         eta = randn(d,M);         
         
         % Coefficients
         diff = Nalphi.*(eta-ones(d,1)*sum(eta.*Vbatch).*Vbatch);
         coll = (Va-ones(d,1)*sum(Va.*Vbatch).*Vbatch);
         coll2 = vecnorm(Vbatch-Va,2).^2.*Vbatch;
         
         % Updating the particles
         Vbatch = Vbatch + dt*lambda*coll - dt*(d-1)*coll2*(sigma.^2./2) + sqrt(dt).*diff*sigma;
         Vbatch = Vbatch./vecnorm(Vbatch,2);
         
     elseif strcmp(mode, 'anisotropic')
         
         eta = randn(d,M);
         Nalphi = (Vbatch-Va);
         
         % Coefficients
         Salphi = ones(d,1)*sum(Va.*Vbatch);
         coll = (Va-Salphi.*Vbatch); 
         
         coll2 = (Nalphi.^2 + vecnorm(Nalphi,2).^2 - 2*vecnorm(Nalphi.*Vbatch,2).^2).*Vbatch;
         
         Setai = ones(d,1)*sum(Nalphi.*eta.*Vbatch);
         diff = Nalphi.*eta-Setai.*Vbatch;
         
         % Updating the particles
         Vbatch = Vbatch + dt*lambda*coll + sqrt(dt).*sigma.*diff - dt*(sigma.^2./2).*coll2;
         Vbatch = Vbatch./vecnorm(Vbatch,2);
         
     end

     %% Update
     if N > M
        V(:,sp) = Vbatch;
     else
        V = Vbatch; 
     end

     Va = computeVa(costfunction, KVCBOparam, V);
     
     if KVCBOparam.logicVa.method == 1 || KVCBOparam.logicVa.method == 2
        Va = Va(:,1);
     end
     
end
    