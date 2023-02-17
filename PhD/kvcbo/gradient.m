function grad = gradient(costfunction, V)

    h = 1e-3;

    d = costfunction.d;
    grad_euc = (costfunction.cost(V*ones(1,d) + h*eye(d)) - costfunction.cost(V*ones(1,d) - h*eye(d)))./(2*h);
    grad_euc = grad_euc'; 

    grad = (eye(d) - V*V')*grad_euc; % gradient projected onto the manifold

end