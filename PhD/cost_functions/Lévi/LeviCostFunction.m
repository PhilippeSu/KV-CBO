classdef LeviCostFunction
   properties
      name, A, B, C, D, phi, d, 
      minimizer
   end

   methods
       
      function obj = LeviCostFunction(d)
          obj.name = 'Lévi';
          obj.A = 1;
          obj.B = 10;
          obj.C = 0;
          obj.D = 4; 
          obj.phi = 0;
          obj.d = d;
          
          minimizer = [zeros(d-1,1); 1];
          obj.minimizer = rotationMatrix(obj.phi, d) * minimizer;
          
      end
       
      function E = cost(obj, V)
          x = obj.D*(V - obj.minimizer); 
          
          w = ones(obj.d,1) + x./obj.A;
          
          a = sin(pi*w(1,:)).^2;
          b = (w-1).^2.*(1+obj.B*sin(pi*w+1).^2);
          c = (w(obj.d,:)-1).^2.*(1+sin(2*pi*w(obj.d,:)).^2);
          
          E = a + sum(b,1) +c;
      end
      
       function err = error(obj, V, KVCBOparam)
         p = KVCBOparam.out.p;
         err = min([norm(V - obj.minimizer, p), norm(V + obj.minimizer, p)]); 
      end
      
   end
end
