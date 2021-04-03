function [Va, V] = GradientDescent(problem, param, V)

  N = param.N;

  % randomly pick one particle
  sub = randi([1 N], 1, 1);
  Vgrad = V(:,sub);

  % compute descent direction
  grad_sph = problem.grad_sph(Vgrad);
  desc_dir = -grad_sph;

  [~, Vgrad_new] = linesearch(problem, Vgrad, desc_dir);

  % update
  V(:,sub) = Vgrad_new;
  Va = computeVa(problem, param, V);

end