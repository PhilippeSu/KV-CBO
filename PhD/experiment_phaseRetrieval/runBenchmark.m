function [successRates, results] = runBenchmark(xitem, xvalues, algorithms, params)

runs = params.runs;
successConstant = params.successConstant;

results = zeros(length(xvalues), length(algorithms), runs);

% Loop over the xvalues
for p = 1:length(xvalues)
    
    fprintf('Running trails: %s=%g\n', xitem, xvalues(p))
    [params] = setupTrialParameters(xitem, xvalues(p), params);
        
    for q = 1:runs
        
        if runs>1 && params.verbose==0
            fprintf('%i: ',q);
        end
               
        % Create a random test problem with the right dimensions resp. SNR
        [A, At, z, y] = generateProblem(params.genpar); 
        size(A) % output the size of A, just to be sure
        d = numel(z);
        
        % Loop over the algorithms
        for k = 1:length(algorithms)
            fprintf('  %s:', algorithms{k}.algorithm);
            
            % Load options of kth algorithm
            opts = algorithms{k};
            
            startTime = tic;
            [zstar, ~, opts] = solvePhaseRetrieval(A, At, y, d, opts);
            elapsedTime = toc(startTime);
            results(p, k, q) = min([norm(zstar-z), norm(zstar+z)]);

            if params.verbose
                fprintf('\nTrial: %d, initial method: %s, algorithm: %s, error: %.2e, time: %f\n',...
                    q, opts.initMethod, opts.algorithm, results(p, k, q), elapsedTime);
            end
            
            if p == length(xvalues) && q == runs && k == length(algorithms)
                plot(1:d,z,'Color',[0 0.4470 0.7410])
                hold on
                
                % alignment of the signs of zstar and z
                zstar = zstar*sign(zstar'*z);
                plot(1:d,zstar,'-.','Color', [0 0.4470 0.7410])
                axis([1, d, -3, 3])
            end
            
        end
    end

end

successRates = mean(results < successConstant, 3);
plotComparison(xitem, xvalues, algorithms, successRates);

end

function [params] = setupTrialParameters(xitem, xval, params)
switch xitem
    case 'M'
        params.genpar.M = xval;
        params.genpar.snr = 0;
        params.genpar.addNoise = false; 
        
    case 'SNR'
        params.genpar.M = 6*params.genpar.d;
        params.genpar.snr = xval;
        params.genpar.addNoise = true; 
end
end 

function plotComparison(xitem, xvalues, algorithms, finalResults)
    
    % get the labels to appear in the legend
    algNames = {};
    for k=1:length(algorithms)
        algNames{k} =  algorithms{k}.label;
    end
    
    autoplot(xvalues, finalResults, algNames);
    
    xlabel(xitem) 
    ylabel('Success Rate')
    yticks([0 0.2 0.4 0.6 0.8 1])
    yticklabels({'0\%', '20\%', '40\%', '60\%', '80\%', '100\%'})
    
    legend('show', 'Location', 'southeast')
    grid on
end


function [A, At, z, y] = generateProblem(pars)
    
    d = pars.d;
    M = pars.M;
    snr_val = pars.snr;

    % generate true signal
    z = generateSignal(pars);
    
    %%
       
    % measurement vectors uniform on the sphere
    A = randn(M,d);
    A = (A'./vecnorm(A'))'; % normalize the rows
    At = A';
    
    % generate (linear) measurements
    y = abs(A*z);
    
%     % check: if the signs were available
%     s = sign(A*z);
%     ys = y.*s;    
%     zs = A\ys;
%     norm(zs-z)

    if pars.addNoise
        
        ynl = y;
        noise_power = sqrt(norm(y)^2 / (10^(snr_val/10)));
        
        % Add noise to measurements
        noise = randn(size(y));
        noise = noise_power*noise/norm(noise);
        
        y = ynl + noise;        
        snr(ynl, noise)         % should be equal to snr_val
        y = max(y,0);
    end
    
%     plot(1:6*d,ynl)
%     hold on
%     plot(1:6*d,y)
    
    % y = generateMeasurements(z,A,pars); % hier werden quadratische Messungen erzeugt
    % y = sqrt(y);  % jetzt haben wir lineare Messungen    
end

function z = generateSignal(pars)

d = pars.d;

switch pars.signalType
    case 'gaussian'
        z = randn(d,1);
        
    case 'sinusoid'
        n = 10;               % number of frequencies
        ampl_max = 1;
        z = zeros(d,1);    
        t = (0:d-1)';

        for j =1:n
            ampl = 2*ampl_max*(rand-0.5);       % amplitude randomly in [-ampl_max, ampl_max]
            freq = randi(floor(d/2));
            z = z + ampl.*sin(freq*t/d);   
        end
        % plot(t,z)
end

end
