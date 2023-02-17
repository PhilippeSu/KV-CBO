classdef CostFunction_robustPCAHaystack
   properties
      name, A, B, C, D, phi, d, 
      minimizer
   end

   methods
       
      function obj = CostFunction_robustPCAHaystack(A,B,C,D,U,d)
          obj.name = 'robust PCA (Haystack experiment)';
          obj.A = A;
          obj.B = B;
          obj.C = C;
          obj.D = D;
          obj.d = d;
          obj.minimizer = U; % 1d subspace we wish to detect
          
      end
       
      function E = cost(obj, V)
            X = obj.A;
            p = obj.B;
            delta = obj.C;
            
            N = size(V,2);
            XN = ones(N,1)*vecnorm(X,2);
            XV = V'*X;
            Z = XN.*XN-XV.*XV;                  % We use the identity||x - <x,v> v ||^2 = ||x||^2 - <x,v>^2
            Z = Z + (delta^2)*ones(size(Z));    % delta = perturbation parameter (for p < 2)
            E = sum((Z)'.^(p/2));
      end
      
     function err = error(obj, V, KVCBOparam)
         p = KVCBOparam.out.p;
         err = min([norm(V - obj.minimizer, p), norm(V + obj.minimizer, p)]); 
      end
      
   end
end
