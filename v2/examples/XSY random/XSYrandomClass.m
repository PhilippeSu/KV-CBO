classdef XSYrandomClass
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
  
      function obj = XSYrandomClass(d)
          obj.name = 'XSYrandom';
          obj.A = 0;
          obj.B = 0;
          obj.C = 0;
          obj.D = 5;
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
        

          a = rand(size(x));
          b = abs(x).^((1:obj.d)'*ones(1,N));
        
          E = sum(a.*b);
          
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
