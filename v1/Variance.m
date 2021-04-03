function [variance_new, I] = Variance(V, mu, Nmin, variance)
  variance_new = sum(vecnorm(V-mean(V,2)*ones(1,size(V,2)),2).^2);  
  N = size(V,2);
  
  % Compute new number of particles
  N_new = fix(N*(1+mu*((variance_new-variance)/variance)));
  
  if N_new > Nmin
      N_drop = max(N - N_new, 0);
      I = randperm(N,N_drop);
  else
      % do not drop any particles
      I = [];
  end
end