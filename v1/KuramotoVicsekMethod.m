function [Valpha, t, varargout] = KuramotoVicsekMethod(energy_functional, d, X, parameters, varargin)
    
    %% Set settings    
    if isfield(parameters,'alpha')
        alpha = parameters.alpha;
    end 
    if isfield(parameters,'alphamax')
        alphamax = parameters.alphamax;
        alpha_adaptive = 1;
    else
        alpha_adaptive = 0;
    end 
    if isfield(parameters,'sigma')
        sigma = parameters.sigma;
    end
    if isfield(parameters,'N')
        N = parameters.N;
    end
    if isfield(parameters,'N2')
        N2 = parameters.N2;
    end
    if isfield(parameters,'N3')
        N3 = parameters.N3;
    end
    if isfield(parameters,'Nmin')
        Nmin = parameters.Nmin;
    end
    if isfield(parameters,'T')
        T = parameters.T;
    end
    if isfield(parameters,'epsilon')
        epsilon = parameters.epsilon;
    end
    if isfield(parameters,'time_limit')
        time_limit = parameters.time_limit;
    end
    if isfield(parameters,'sigma')
        sigma = parameters.sigma;
    end
    if isfield(parameters,'lambda')
        lambda = parameters.lambda;
    end
    if isfield(parameters,'dt')
        dt = parameters.dt;
    end
    if isfield(parameters,'mu')
        mu = parameters.mu;
    end
    if isfield(parameters,'eta')
        eta = parameters.eta;
    end
    if isfield(parameters,'tau')
        tau = parameters.tau;
    end
    if isfield(parameters,'l')
        l = parameters.l;
    end
    if isfield(parameters,'trueSol')
        trueSol = parameters.trueSol;
        computeError = 1;
        Recon = parameters.Recon;
    else
        computeError = 0;
    end

    %% sKV method
    % Initial Data
    V = RandHyperSphere(d,N,1);                         % initial particles in euclidean coordinates
    v = InvHyperSphere(V);                              % initial particles in spherical coordinates 
  
    Ep = energy_functional(X,V,varargin{1},varargin{2});
        
    % Energy rescaling
    Em = min(Ep);
    EMM = max(Ep);
    Ne = EMM-Em;
    Ep = (Ep-Em)/Ne;

    t = 1;
    Valpha = ComputeValpha(energy_functional, X, V, N, alpha, varargin{1}, varargin{2});
    variance = sum(vecnorm(V-mean(V,2)*ones(1,size(V,2)),2).^2);  

    error_hist = zeros(1,T);
    variance_hist = zeros(1,T);
    avg_variance_hist = zeros(1,T);
    var_avg_variance_hist = zeros(1,T);
    stopping_criterion = sum(vecnorm(V-Valpha*ones(1,N),2))/N < epsilon;

    f = @() KuramotoIteration(energy_functional, d, Valpha, V, X, ...
                        alpha, sigma, lambda, dt, N, N2, N3, ...
                        varargin{1}, varargin{2});
                    
    while and(t < T, ~stopping_criterion)

        [Valpha, V] = KuramotoIteration(energy_functional, d, Valpha, V, X, ...
                        alpha, sigma, lambda, dt, N, N2, N3, ...
                        varargin{1}, varargin{2});

         % discard particles every l-th iteration depending on change in variance
        [variance, I] = Variance(V, mu, Nmin, variance);
        variance_hist(t) = variance;
        avg_variance_hist(t) = sum(variance_hist)/t;
        var_avg_variance_hist(t) = var(avg_variance_hist(1:t));
        r = Recon(Valpha, varargin{2});
        error_hist(t) = min([norm(r-trueSol,2), norm(r+trueSol,2)]);
        
        if mod(t,l) == 0
            if var_avg_variance_hist(t) < eta
                sigma = sigma./tau;
            end
            if alpha_adaptive == 1
                alpha = min(2*alpha,alphamax);
            end
            if computeError == 1
                r = Recon(Valpha, varargin{2});
                error = min([norm(r-trueSol,2), norm(r+trueSol,2)]);
                 
                msg = sprintf([' Iterations: \t \t %i\n', ...
                               ' max. Iterations: \t %i\n\n', ...
                               ' Variance: \t \t %.2e\n', ...
                               ' Error: \t \t %.2e\n\n', ...
                               ' Number of Particles: \t %i'] ...
                               , t, T, variance, error, N);
                % pause(0.2)
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
                fprintf([reverseStr, msg]);
            end
            V(:,I) = [];
            Ep(:,I) = [];
            N = size(V,2);
        end

        stopping_criterion = sum(vecnorm(V-Valpha*ones(1,N),2))/N < epsilon;
        t = t + 1;
    end

    varargout{1} = timeit(f); % median runtime per call of KuramotoIteration
    varargout{2} = variance_hist;
    varargout{3} = error_hist;

end
