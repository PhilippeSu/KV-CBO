classdef robustPCAClass
   properties
      name
      A
      B
      C
      D
      d
      Vmin
      min
   end

   methods
       
      function obj = robustPCAClass(A,B,C,D,U,d)
          obj.name = 'robust PCA';
          obj.A = A;
          obj.B = B;
          obj.C = C;
          obj.D = D;
          obj.d = d;
          obj.Vmin = U; % das ist ja eigentlich falsch: das Minimum liegt ja nicht mehr hier wenn wir Ausreißer hinzufügen
          
          obj.min = obj.cost(obj.Vmin);
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
      
      function err = error(obj, V)
         err = min([norm(V - obj.Vmin), norm(V + obj.Vmin)]); 
      end
      
      function g = grad(obj, V)
          h = 1e-6;
  
          g = (obj.cost(V*ones(1,obj.d) + h*eye(obj.d)) - obj.cost(V*ones(1,obj.d) - h*eye(obj.d)))./(2*h);
          g = g'; 
      end
      
      function g = grad_sph(obj, V)
         g = (eye(obj.d) - V*V')*grad(obj, V); 
      end
      
   end
end
