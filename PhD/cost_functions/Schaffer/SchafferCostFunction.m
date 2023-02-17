classdef SchafferCostFunction
   properties
      name, A, B, C, D, phi, d, 
      minimizer
   end

   methods
       
      function obj = SchafferCostFunction(d)
          obj.name = 'Schaffer';
          obj.A = 0.5;
          obj.B = 1;
          obj.C = 0;
          obj.D = 2;
          obj.phi = 0;
          obj.d = d;
          
          minimizer = [zeros(d-1,1); 1];
          obj.minimizer = rotationMatrix(obj.phi, d) * minimizer;
          
      end
       
      function E = cost(obj, V)
          x = obj.D*(V - obj.minimizer); 
          
          a = vecnorm(x,2);
          b = sin(a).^2;
          c = (1+obj.B*a.^2).^2;
          
          E = obj.A + (b-obj.A)./c;
      end
      
       function err = error(obj, V, KVCBOparam)
         p = KVCBOparam.out.p;
         err = min([norm(V - obj.minimizer, p), norm(V + obj.minimizer, p)]); 
      end
      
   end
end
