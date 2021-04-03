
classdef pointCloudClass
     
methods (Static)
    
    function X = Gaussian(d)
      Mu = zeros(d,1); 
      Sigma = eye(d);
      N = 1000;
        
      X = mvnrnd(Mu, Sigma, N);
    end
    
    function X = Gaussian2(d)
      Mu = zeros(d,1);
      v =  randn(d,1);
      v = v/norm(v);
      Sigma = v*v';
      N = 1000;
        
      X = mvnrnd(Mu, Sigma, N);
    end
    
    function X = Tdistribution(d)
      Sigma = eye(d);
      dof = 2;
      N = 1000;
        
      X = mvtrnd(Sigma, dof, N);
    end
    
    function X = Tdistribution2(d)
      v =  randn(d,1);
      v = v/norm(v);
      Sigma = v*v';
      dof = 3;
      N = 1000;
        
      X = mvtrnd(Sigma, dof, N);
    end
    
    function X = LogNormal(d)
      Mu = zeros(d,1); 
      Sigma = eye(d);
      N = 1000;
        
      Y = mvnrnd(Mu, Sigma, N);
      X = exp(Y);
    end
    
    function X = GaussianMixture(d)
      Mu_1 = zeros(d,1);
      Mu_2 = zeros(d,1);
      Sigma_1 = eye(d);
      Sigma_2 = 100.*eye(d);
      N = 1000;
      epsilon = 0.05;
      
      Mu = (1-epsilon).*Mu_1 + epsilon.*Mu_2;
      Sigma = (1-epsilon).*Sigma_1 + epsilon.*Sigma_2;
      
      X = mvnrnd(Mu, Sigma, N);
    end
    
    function [W, w] = RandomWithMinimumSpacing(d,m)
        % TODO: vereinfachen
        W = zeros(d^2, m);
        w = zeros(d,m);

        for i=1:m

            while (true)
                % unknown vector
                wi = randn(d,1);
                wi = wi/norm(wi);

                % wi on upper half sphere
                wi = InvHyperSphere(wi);
                wi = mod(wi,pi);
                wi = HyperSphere(wi);

                w(:,i) = wi;                  
                h = wi*wi';
                h = reshape(h, d^2, 1);
                W(:,i) = h;

                % dist = vecnorm(w(:,1:i-1)-wi*ones(1,i-1),2);
                dist = acos(bsxfun(@dot, w(:,1:i-1),wi*ones(1,i-1)));
                if and(i>1, all(dist > 1.5))
                    break
                elseif i == 1
                    break
                end

            end
        end
    end
    
    function [X, U] = Haystack(d, varargin)
        if nargin == 1
            sigma_in = 1;
            sigma_out = 2;
            ds = 1;           % dimension of subspace
            N_in = 100;
            N_out = 50;
            noise = 0;
        elseif nargin == 2
            cloud_param = varargin{:};
            sigma_in = cloud_param.sigma_in;
            sigma_out = cloud_param.sigma_out;
            ds = cloud_param.ds;           % dimension of subspace
            N_in = cloud_param.N_in;
            N_out = cloud_param.N_out;
            noise = cloud_param.noise;
        end
      
      U = orth(randn(d, ds));   
      
      cov_in = (sigma_in^2/ds).*U*U';
      cov_out = (sigma_out^2/d).*eye(d);
      mu = zeros(d,1);
      
      X_in = mvnrnd(mu, cov_in, N_in) + sqrt(noise/d)*randn(N_in, d);
      X_out = mvnrnd(mu, cov_out, N_out);
      
      X = [X_in; X_out]';  
    end
        
end

end

