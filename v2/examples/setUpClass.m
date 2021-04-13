classdef setUpClass
   
    methods (Static)
        
        function [problem, param] = Ackley(d) 
            problem = AckleyClass(d);
            param = paramAckley();
        end
        
        function [problem, param] = Rastrigin(d)
            problem = RastriginClass(d);
            param = paramRastrigin();
        end
        
        function [problem, param] = XSYrandom(d)
            problem = XSYrandomClass(d); 
            param = paramXSYrandom();
        end
        
        function [problem, param] = Alpine(d)
           problem = AlpineClass(d); 
           param = paramAlpine();
        end
        
        function [problem, param] = Griewank(d)
           problem = GriewankClass(d);
           param = paramGriewank();
        end
        
        function [problem, param] = Solomon(d)
           problem = SolomonClass(d);
           param = paramSolomon();
        end
        
        function [problem, param] = MaxEigenvalue(d)
           problem = MaxEigenvalueClass(d);
           param = paramMaxEigenvalue();
        end
        
        function [problem, param] = Haystack(d)
            
            P = 200; % data points 
            outlier_per = 0.75;
            cloud_param.sigma_in = 1;
            cloud_param.sigma_out = 1;
            cloud_param.ds = 1;           % dimension of subspace
            cloud_param.N_in = round(P*(1-outlier_per),0);
            cloud_param.N_out = round(P*outlier_per,0);
            cloud_param.noise = 1e-4;
            
            [X, U] = pointCloudClass.Haystack(d, cloud_param);
            
            p = 1;
            delta = 1e-4;
            problem = robustPCAClass(X, p, delta, 0, U, d);
            
            param = paramHaystack();
        end
        
        function [problem, param] = Eigenfaces(~)
            
            %% Load data
            headerLines = 0;
            field_separator = ',';
            
            s = 'faces.txt';
            % s_benchmark = 'faces_bench.txt';
            if ~exist(s)
                fprintf('\nPlease provide the database faces.txt\n');
            end
            X = importdata(s, field_separator, headerLines);
            X = bsxfun(@minus, X, mean(X,2)); 
                        
            % Benchmark = importdata(s_benchmark, field_separator, headerLines);
            % Benchmark = bsxfun(@minus, Benchmark, mean(Benchmark,2)); 
              
            %% Benchmark method
            [Vg,D] = eig(X*X');
            [~,i] = max(diag(D));
            Vbench = Vg(:,i);
            
            %% set up problem
            d = size(X,1);            
            p = 0.1;
            problem = robustPCAClass(X, p, 0, 0, Vbench, d);
            
            %% set parameters
            param = paramEF(); 
        end
        
        function [problem, param] = PhaseRetrieval(d)           
            M = 2*d; 
            
            % frame
            X = randn(d,M);
            X = X./vecnorm(X,2); 
            
            % lower frame bound
            S = zeros(d,d);
            for i=1:d
                Sd = sum(X*diag(X(i,:)),2);
                S(i,:) = Sd;
            end          
            A = min(eig(S));
            
            % original problem
            z = randn(d,1);                                 % vector to recover
            Y = ((dot(z*ones(1,M),X)).^2)';                 % vector in R^M (measurements)
            
            % problem on the sphere
            C = sqrt(norm(Y,1)/A);
            Y = Y./(C*C);            
            X = [X; zeros(1,M)];
            
            problem = PhaseRetrievalClass(X, Y, C, z, d+1);
            
            param = paramPhaseRetrieval();
        end
        
    end
end

