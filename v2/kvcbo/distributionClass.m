%%
classdef distributionClass
   
    methods (Static)
        
        function  V = initialParticles(problem, param)
            
            if strcmp(param.distri, 'uniform')
                V = distributionClass.uniform(problem, param);
            elseif strcmp(param.distri, 'fisher')
                V = distributionClass.fisher(problem, param);
            end
            
        end
        
        function V = uniform(problem, param)
            N = param.N;
            d = problem.d;
            half = param.half;
     
            V = randn(d,N);
            V = V./vecnorm(V,2);
            
            if half
                phi = InvHyperSphere(V);
                phi = mod(phi,pi);
                V = HyperSphere(phi);
            end
        end
        
        function V = fisher(problem, param)
             % von-Mises-Fisher distribution
            N = param.N; 
            meanDir = param.meanDir;
            k = param.k; 
            V = randVMF(N, meanDir, k)';
        end
        
    end
end

