classdef XSYrandomCostFunction
   properties
      name, A, B, C, D, phi, d, 
      minimizer
   end

   methods
  
      function obj = XSYrandomCostFunction(d)
          obj.name = 'XSYrandom';
          obj.A = 0;
          obj.B = 0;
          obj.C = 0;
          obj.D = 0.5;
          obj.phi = 0;
          obj.d = d;
          
          minimizer = [zeros(d-1,1); 1];
          obj.minimizer = rotationMatrix(obj.phi, d) * minimizer;
          
      end
       
      function E = cost(obj, V)
          x = obj.D*(V - obj.minimizer);
          N = size(V,2);
        
          a = rand(size(x));
          b = abs(x).^((1:obj.d)'*ones(1,N));
        
          E = sum(a.*b);
          
      end
      
     function err = error(obj, V, KVCBOparam)
         p = KVCBOparam.out.p;
         err = min([norm(V - obj.minimizer, p), norm(V + obj.minimizer, p)]); 
      end
      
   end
end
