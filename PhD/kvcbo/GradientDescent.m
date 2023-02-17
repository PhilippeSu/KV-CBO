function V = GradientDescent(costfunction, KVCBOparam, V)

  N = KVCBOparam.N;

  % randomly pick one particle
  sub = randi([1 N], 1, 1);
  Vgrad = V(:,sub);

  % compute descent direction
  grad_sph = gradient(costfunction, Vgrad);
  desc_dir = -grad_sph;

  [~, Vgrad_new] = linesearch(costfunction, Vgrad, desc_dir);

  % update
  V(:,sub) = Vgrad_new;

end