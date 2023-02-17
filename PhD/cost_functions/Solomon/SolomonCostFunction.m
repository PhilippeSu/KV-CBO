classdef SolomonCostFunction
   properties
      name, A, B, C, D, phi, d, 
      minimizer
   end

   methods

      function obj = SolomonCostFunction(d)
          obj.name = 'Solomon';
          obj.A = -1;
          obj.B = 0.1;
          obj.C = 1;
          obj.D = 5;
          obj.phi = 0;
          obj.d = d;
          
          minimizer = [zeros(d-1,1); 1];
          obj.minimizer = rotationMatrix(obj.phi, d) * minimizer;
          
      end
       
      function E = cost(obj, V)
         x = obj.D*(V - obj.minimizer); 
         b = vecnorm(x,2);
         
         E = obj.A*cos(2*pi*b) + obj.B*b + obj.C;
      end
      
     function err = error(obj, V, KVCBOparam)
         p = KVCBOparam.out.p;
         err = min([norm(V - obj.minimizer, p), norm(V + obj.minimizer, p)]); 
      end
      
   end
end
