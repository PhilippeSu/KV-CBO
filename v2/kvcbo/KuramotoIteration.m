%%
function [Va, V] = KuramotoIteration(problem, param, V)
        
     Nbatch = param.Nbatch;
     N = param.N;
     sigma = param.sigma;
     lambda = param.lambda;
     mode = param.mode;
     dt = param.dt;
     d = problem.d;
     
     if Nbatch < N
         sp = randperm(N,Nbatch);
         Vbatch = V(:,sp);
     else
         Nbatch = N;
         Vbatch = V;
     end
     
     Va = computeVa(problem, param, Vbatch);
      
     %% 
     if strcmp(mode, 'isotropic')
         
         % Vectorialization
         Nalphi = ones(d,1)*vecnorm(Vbatch-Va*ones(1,Nbatch),2);
         eta = randn(d,Nbatch);
         Va = Va*ones(1,Nbatch);
         
         % Coefficients
         diff = Nalphi.*(eta-ones(d,1)*sum(eta.*Vbatch).*Vbatch);
         coll = (Va-ones(d,1)*sum(Va.*Vbatch).*Vbatch);
         coll2 = vecnorm(Vbatch-Va,2).^2.*Vbatch;
         
         % Updating the particles
         Vbatch = Vbatch + dt*lambda*coll - dt*(d-1)*coll2*(sigma.^2./2) + sqrt(dt).*diff*sigma;
         Vbatch = Vbatch./vecnorm(Vbatch,2);
         
     elseif strcmp(mode, 'anisotropic')
         
         eta = randn(d,Nbatch);
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
     if N > Nbatch
        V(:,sp) = Vbatch;
     else
        V = Vbatch; 
     end

     Va = computeVa(problem, param, V);
     
end
    