classdef GriewankClass
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
       
      function obj = GriewankClass(d)
          obj.name = 'Griewank';
          obj.A = 1;
          obj.B = 1/4000;
          obj.C = 0;
          obj.D = 600;
          obj.phi = 0;
          obj.d = d;
          
          argmin = [zeros(d-1,1); 1];
          argmin = rotationMatrix(obj.phi, d) * argmin;
          
          obj.Vmin = argmin;
          obj.min = obj.cost(argmin);
      end
       
      function E = cost(obj, V)
         x = obj.D*(V - obj.Vmin); 
         N = size(V,2);
         a = vecnorm(x,2)/sqrt(obj.d);
         b = prod(cos(x)./sqrt((1:obj.d)'*ones(1,N)));
         
         E = obj.A + obj.B*a - b;
      end
      
      function err = error(obj, V)
         err = min([norm(V - obj.Vmin), norm(V + obj.Vmin)]); 
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
