classdef NNreconCostFunction
   properties
      name, Lstar, d, minimizer
   end

   methods

      function obj = NNreconCostFunction(d,Lstar,A)
          obj.name = 'Reconstruction of neural nets';  
          obj.Lstar = Lstar;
          obj.d = d;
          obj.minimizer = obj.cost(A(:,1));
      end
       
      function E = cost(obj, V)
         U = obj.Lstar;
         H = kr(V,V); 
         Z = U'*H;
         Z = U*Z;
         E = -(vecnorm(Z).^2) +1;
      end
      
      function err = error(obj, V, KVCBOparam)
         err = 0;
      end
     
   end
end
