classdef PhaseRetrievalClass
   properties
      name
      X
      Y
      C
      z
      d
      Vmin
      min
   end

   methods
       
      function obj = PhaseRetrievalClass(X,Y,C,z,d)
          obj.name = 'Phase Retrieval';
          obj.X = X;
          obj.Y = Y;
          obj.C = C;
          obj.d = d;
          obj.Vmin = z;          
          % obj.min = obj.cost(obj.Vmin);
      end
       
      function E = cost(obj, V)
          nn = size(V,2);
            
          XV = V'*obj.X;
          Z = (XV.*XV - ones(nn,1)*obj.Y').^2;
            
          E = sum(Z');
      end
      
      function err = error(obj, V)
         recon = obj.C*V(1:end-1);
         err = min([norm(recon - obj.Vmin), norm(recon + obj.Vmin)]);
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
