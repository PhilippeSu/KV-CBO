classdef RastriginClass
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

      function obj = RastriginClass(d)
          obj.name = 'Rastrigin';
          obj.A = 0;
          obj.B = 10;
          obj.C = 5;
          obj.D = 5.12;
          obj.phi = 0;
          obj.d = d;
          
          argmin = [zeros(d-1,1); 1]; % randn(d,1); 
          argmin = rotationMatrix(obj.phi, d) * argmin;
          
          obj.Vmin = argmin;
          obj.min = obj.cost(argmin);
      end
       
      function E = cost(obj, V)
          
          x = obj.D*(V - obj.Vmin);
          
          a = sum(x.^2);
          omega = 2*pi;
          b = sum(cos(omega*(x)));
          
          E = obj.B*obj.d + a - obj.B*b;
      end
      
      function err = error(obj, V)
         % p = 'inf';
         p = 2;
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
