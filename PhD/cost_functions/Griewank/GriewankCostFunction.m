classdef GriewankCostFunction
   properties
      name, A, B, C, D, phi, d, 
      minimizer
   end

   methods
       
      function obj = GriewankCostFunction(d)
          obj.name = 'Griewank';
          obj.A = 1;
          obj.B = 1/50; % 1/4000;
          obj.C = 0;
          obj.D = 50; % 600;
          obj.phi = 0;
          obj.d = d;
          
          minimizer = [zeros(d-1,1); 1];
          obj.minimizer = rotationMatrix(obj.phi, d) * minimizer;
          
      end
       
      function E = cost(obj, V)
         x = obj.D*(V - obj.minimizer); 
         N = size(V,2);
         a = vecnorm(x,2)/sqrt(obj.d);
         b = prod(cos(x)./sqrt((1:obj.d)'*ones(1,N)));
         
         E = obj.A + obj.B*a - b;
      end
      
       function err = error(obj, V, KVCBOparam)
         p = KVCBOparam.out.p;
         err = min([norm(V - obj.minimizer, p), norm(V + obj.minimizer, p)]); 
      end
      
   end
end
