classdef RastriginCostFunction_isoVsAniso
   properties
      name, A, B, C, D, phi, d, 
      minimizer
   end

   methods

      function obj = RastriginCostFunction_isoVsAniso(d)
          obj.name = 'Rastrigin';
          obj.A = 10;
          obj.B = 0;
          obj.C = 0;
          obj.D = 10;
          obj.phi = 0;
          obj.d = d;
          
          minimizer = [zeros(d-1,1); 1]; % randn(d,1); 
          obj.minimizer = rotationMatrix(obj.phi, d) * minimizer;          
      end
       
      function E = cost(obj, V) 
          x = obj.D*(V - obj.minimizer);
          
          a = sum(x.^2);
          b = sum(cos(2*pi*(x)));
          
          E = obj.A*obj.d + a - obj.A*b;
      end
      
      function err = error(obj, V, KVCBOparam)
         p = KVCBOparam.out.p;
         err = min([norm(V - obj.minimizer, p), norm(V + obj.minimizer, p)]); 
      end
     
   end
end
