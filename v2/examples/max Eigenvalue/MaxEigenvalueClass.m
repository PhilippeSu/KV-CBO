classdef MaxEigenvalueClass
   properties
      name
      A
      B
      C
      D
      phi
      d
      Vmin
      min
   end

   methods
  
      function obj = MaxEigenvalueClass(d)
          obj.name = 'maximum Eigenvalue';
          
          A = randn(d);
          A = 0.5*(A'*A);
          
          obj.A = A;
          obj.B = 0;
          obj.C = 0;
          obj.D = 0;
          obj.phi = 0;
          obj.d = d;
          
          [U,~] = eig(A);
          obj.Vmin = U(:,end);
          
          obj.min = obj.cost(obj.Vmin);
      end
       
      function E = cost(obj, V)
          E = diag(V'*(obj.A*V))';
      end
      
      function err = error(obj, V)
         p = 'inf';
         % p = 2;
         err = min([norm(V - obj.Vmin, p), norm(V + obj.Vmin, p)]); 
      end
      
      function g = grad(obj, V)
          h = 1e-3;
  
          g = (obj.cost(V*ones(1,obj.d) + h*eye(obj.d)) - obj.cost(V*ones(1,obj.d) - h*eye(obj.d)))./(2*h);
          g = g'; 
      end
      
      function g = grad_sph(obj, V)
          g = (eye(obj.d) - V*V')*grad(obj, V);
      end
      
   end
end
