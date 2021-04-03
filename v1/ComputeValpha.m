function [Valpha] = ComputeValpha(energy_functional, X, V, N3, alpha, varargin)
   N = size(V,2);

    if N3 < N
        % choose a random subset of the particles
        sp = randperm(N,N3);
        Vs = V(:,sp);
    else
        Vs = V;
    end
    
    Es = energy_functional(X, Vs, varargin{1}, varargin{2});   
    
    Esmin = min(Es);    % Vmin is the particles with the lowest energy in this iteration
    Valpha = sum((Vs.*exp(-alpha*(Es-Esmin))),2);
    Valpha = Valpha/sum(exp(-alpha*(Es-Esmin)));
end
