function [Valpha, V] = KuramotoIteration(energy_functional, d, Valpha, V, X, ...
                     alpha, sigma, lambda, dt, N, N2, N3, ...
                     varargin)

     % evolve only few randomly choosen particles
     if N > N2
         sp = randperm(N,N2);
         V2 = V(:,sp);
     else
         N2 = N;
         V2 = V;
     end
                 
     % Vectorialization
     Nalphi = ones(d,1)*vecnorm(V2-Valpha*ones(1,N2),2);
     sgn = 1;
     eta = randn(d,N2);
     Valpha = Valpha*ones(1,N2);

     % Coefficients
     diff = Nalphi.*(eta-ones(d,1)*sum(eta.*V2).*V2);
     coll = sgn.*(Valpha-ones(d,1)*sum(Valpha.*V2).*V2);
     coll2 = vecnorm(V2-Valpha,2).^2.*V2;
     
     % Updating the particles
     V2 = V2 + dt*lambda*coll - dt*(d-1)*coll2*(sigma.^2./2) + sqrt(dt).*diff*sigma; 
     V2 = V2./vecnorm(V2,2);

     if N > N2
        V(:,sp) = V2;
     else
        V = V2; 
     end
     
     Valpha = ComputeValpha(energy_functional, X, V, N3, alpha, varargin{:});
      
end
    