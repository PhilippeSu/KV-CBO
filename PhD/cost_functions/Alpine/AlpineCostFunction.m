classdef AlpineCostFunction
   properties
      name, A, B, C, D, phi, d, 
      minimizer
   end

   methods
       
      function obj = AlpineCostFunction(d)
          obj.name = 'Alpine';
          obj.A = 0.1;
          obj.B = 0;
          obj.C = 0;
          obj.D = 6;
          obj.phi = 0;
          obj.d = d;
          
          minimizer = [zeros(d-1,1); 1];
          obj.minimizer = rotationMatrix(obj.phi, d) * minimizer;
          
      end
       
      function E = cost(obj, V)
        x = obj.D*(V - obj.minimizer); 
        
        E = sum(abs(x.*sin(x) - obj.A*x));
      end
      
      function err = error(obj, V, KVCBOparam)
         p = KVCBOparam.out.p;
         err = min([norm(V - obj.minimizer, p), norm(V + obj.minimizer, p)]); 
      end
      
   end
end
