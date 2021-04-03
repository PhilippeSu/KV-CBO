%%
function Va = computeVa(problem, param, V)
    
  alpha = param.alpha;
  logicVa = param.logicVa; 
  
  Es = problem.cost(V);
  
  [Esmin, i] = min(Es);
  
  if strcmp(logicVa, 'min')
    Va = V(:,i);
  else
    Va = sum((V.*exp(-alpha*(Es-Esmin))),2);
    Va = Va./sum(exp(-alpha*(Es-Esmin)));    
  end
  
end
