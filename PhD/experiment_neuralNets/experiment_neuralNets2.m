clear; clc; close all;

d = 3;
runs = 20;
m_range = [1,2,3];

KVCBOparam = KVCBOparamForNNrecon(d);

success_kNN = zeros(1,size(m_range,2));
success_CS = zeros(1,size(m_range,2));

for m = m_range

for i = 1:runs

    fprintf('run=%i\n', i)
    
    % b and theta (weights in network)
    b = 0.3*randn(1,m);
    t = 0.3*randn(1,m);

    % activation function
    sigma = @(x) tanh(x); 

    %% Reconstruction

    h = 0.01;       % step size for disc. of gradient
    mx = 1000;      % samples from the Hessian
    k = m;          

    A_recon = zeros(d,m);

    % A orthogonal
    A = randn(d,m);
    A = A*(A'*A)^(-1/2);

    % make A quasi orthogonal
    noise = 0.01;
    xi = noise*randn(min(d,m),1);
    [U,S,V] = svd(A);
    Sn = padarray(diag(xi),d-m,0,'post') + S;
    A = U*Sn*V';
    A = A./vecnorm(A);

    f = @(x) NN(x,A,b,t,sigma);
    Lstar = approx_subspaceL(f,d,h,mx,k);

    %% KV-CBO

    costfunction = NNreconCostFunction(d,Lstar,A);
    
    % kNN weight
    k = 20;
    p = 1;
    
    KVCBOparam.logicVa = logicVakNN(k,p);
    [Va,info] = KVCBO(costfunction,KVCBOparam);
    
    c = count_foundMinimizers(Va,A)
    if c == m
        success_kNN(m == m_range) = success_kNN(m == m_range) + 1
    end
    
    % CS weight
    b = 100;
    C = 1;
    p = 1;
    
    KVCBOparam.logicVa = logicVaCS(b,C,p);
    [Va,info] = KVCBO(costfunction,KVCBOparam);
    
    c = count_foundMinimizers(Va,A)
     if c == m
        success_CS(m == m_range) = success_CS(m == m_range) + 1
    end

end

end




%% Auxiliary functions

function c = count_foundMinimizers(Va,A)

c = 0;
n = size(Va,2);
m = size(A,2);

for j = 1:m
   a = vecnorm(A(:,j)*ones(1,n) - Va);
   I1 = a < 0.05;
   I2 = a > 1.95;

   if sum(I1) > sum(I2)
       I = I1;
   else
       I = I2;
   end

   if sum(I) > 1
       % A_recon(:,j) = mean(Va(:,I),2);
       c = c + 1;
   end
end
    
end

function logicVa = logicVakNN(k,p)

    logicVa.method = 4;
    logicVa.name = 'kNN (multiple consensus points)';
    logicVa.k = k;
    logicVa.p = p;
    
end

function logicVa = logicVaCS(b,C,p)

    logicVa.method = 3;
    logicVa.name = 'CS (multiple consensus points)';
    logicVa.b = b;
    logicVa.C = C;
    logicVa.p = p;
    
end

function Lstar = approx_subspaceL(f,d,h,mx,k)

    hess_f = zeros(d^2,mx);

    % sample points for Hessian evaluation
    xi = randn(d,mx);
    xi = xi./vecnorm(xi);

    % matrix of vectorized Hessians
    for i = 1:mx
        hess_f(:,i) = reshape(hessNN_disc(f,xi(:,i),h), [d^2 1]); % war [1 d^2]
    end

    % construct Lstar
    [U,~,~] = svd(hess_f);
    Lstar = U(:,1:k);
    
end

function f = NN(X,A,b,t, sigma)

[~,n] = size(X);
f = sigma(A'*X + t'*ones(1,n));
f = dot((ones(n,1)*b)',f);

end

function hessf = hessNN_disc(f,x,h)
    % x is a vector

    d = size(x,1);
    hessf = zeros(d,d);
    I = eye(d,d);

    for i = 1:d
        for j = 1:d
            a1 = f(x + h*I(:,i) + h*I(:,j));
            a2 = f(x + h*I(:,i));
            a3 = f(x + h*I(:,j));
            a4 = f(x);

            hessf(i,j) = a1 - a2 - a3 + a4;
        end
    end
    hessf = hessf./(h^2);

end
