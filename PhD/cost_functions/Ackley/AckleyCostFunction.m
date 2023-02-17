classdef AckleyCostFunction
   properties
      name, A, B, C, D, phi, d, 
      minimizer
   end

   methods
       
      function obj = AckleyCostFunction(d)
          obj.name = 'Ackley';
          obj.A = -20;
          obj.B = -0.2;
          obj.C = 20;
          obj.D = 5;
          obj.phi = 0;
          obj.d = d;
          
          minimizer = [zeros(d-1,1); 1];
          obj.minimizer = rotationMatrix(obj.phi, d) * minimizer;
          
      end
       
      function E = cost(obj, V)
         x = obj.D*(V - obj.minimizer); 
         a = vecnorm(x,2)/sqrt(obj.d);
         b = mean(cos(2*pi*(x)));
             
         E = obj.A*exp(obj.B*a) - exp(b) + exp(1) + obj.C;

      end
      
       function err = error(obj, V, KVCBOparam)
         p = KVCBOparam.out.p;
         err = min([norm(V - obj.minimizer, p), norm(V + obj.minimizer, p)]); 
      end
      
   end
end
