classdef Energy_Class
    
    methods (Static)
        
        function E = Ackley(V,B) 
            [d,~] = size(V);
            x = 3*(V - B);
            C=0;
            
            if d > 2
                a = vecnorm(x,2)/sqrt(d);
                b = mean(cos(2*pi*(x)));
            else
                a = abs(x)/sqrt(d);
                b = cos(2*pi*(x))/(d);
            end
            
            E=-20*exp(-0.2*a)-exp(b)+exp(1)+20+C;
        end
        
        function E = Rastrigin(V,B)
            [d,~]=size(V);
            x = 3*(V - B);
            C = 0 ;  
            
            if d > 2
                a = mean(x.^2);
                b = mean(cos(2*pi*(x)));
            else
                a = abs(x).^2/(d);
                b = cos(2*pi*(x))/d;
            end
            
            E= a - 10*b + 10 + C;
        end
        
        function E = RobustSubspaceDetection(X,V,p,delta)
            nn = size(V,2);
            XN = ones(nn,1)*vecnorm(X,2);
            XV = V'*X;
            Z = XN.*XN-XV.*XV;                  % We use the identity||x - <x,v> v ||^2 = ||x||^2 - <x,v>^2
            Z = Z + (delta^2)*ones(size(Z));    % delta = perturbation parameter (for p < 2)
            E = sum((Z)'.^(p/2));
        end
        
        function E = PhaseRetrieval(X,V,Y,C)  
            nn = size(V,2);
            Y = Y./(C*C);
            
            XV = V'*X;
            Z = (XV.*XV - ones(nn,1)*Y').^2;
            
            E = sum(Z');
        end
        
    end
end

