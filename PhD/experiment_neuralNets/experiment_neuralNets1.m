clear; clc; close all;

%% Setup

d = 5;          % input dimension
m = d;          % number of hidden neurons

% generate orthogonal weights
A = randn(d,m);
A = A*(A'*A)^(-1/2); 

% make weights quasi-orthogonal
noise = 0.01;
xi = noise*randn(min(d,m),1); 
[U,S,V] = svd(A);
Sn = padarray(diag(xi),d-m,0,'post') + S;
A = U*Sn*V';
A = A./vecnorm(A);
epsilon = norm(diag(Sn - eye(d,m))); 

% b, theta, sigma (bias, threshold and activation of neural net)
b = 0.3*randn(1,m);
t = 0.3*randn(1,m);
sigma = @(x) tanh(x); 

%% Reconstruct weights with algorithm 6

mx = 10*d; 
h = 0.01; 
k = 2*m;             % in {m,...,d^2}
gamma = 10; 
N = 10000; 

f = @(x) NN(x,A,b,t,sigma);

[A_recon, counter] = algorithm6(f,d,m,mx,h,k,gamma,N,A);

[~,S,~] = svd(A);
[~,Sr,~] = svd(A_recon);

error = norm(diag(S-Sr),1)
counter

%% Reconstruct weights with algorithm 7

[A_recon, counter] = algorithm7(f,d,m,mx,h,k,gamma,N,A);

[~,S,~] = svd(A);
[~,Sr,~] = svd(A_recon);

error = norm(diag(S-Sr),1)
counter


%% auxiliary functions

function [A_recon, counter] = algorithm7(f,d,m,mx,h,k,gamma,N,A)

    Lstar = approx_subspaceL(f,d,h,mx,k);    
    Pstar = Lstar*Lstar';
    
    counter = 0;
    found = 0;
    A_recon = zeros(d,m);

    while found < m

        % starting point
        u = randn(d,1);
        u = u/norm(u);
        
        for i = 1:N
            uu = reshape(u*u', [d^2, 1]);
            Pu = Pstar*uu;
            Pu = reshape(Pu, [d, d]);
            u = u + 2*gamma*Pu*u;
            u = u/vecnorm(u);
        end
        
        [dist, ind] = min([vecnorm(u*ones(1,m) - A), vecnorm(u*ones(1,m) + A)]);
        if ind > m
            ind = ind - m;
        end

        if dist < 0.05 && norm(A_recon(:,ind)) < eps
            found = found + 1;
            A_recon(:,ind) = u;
        end

        counter = counter + 1;
        
        if counter > 500
           % not found minimizer nf
           [~,nf] = min(vecnorm(A_recon));
           break 
        end
    end

end

function [A_recon, counter] = algorithm6(f,d,m,mx,h,k,gamma,N,A)

    Lstar = approx_subspaceL(f,d,h,mx,k);
       
    counter = 0;
    found = 0;
    A_recon = zeros(d,m);

    while found < m

        u1 = projected_gradient_ascent(Lstar,N,gamma);
        
        [dist, ind] = min([vecnorm(u1*ones(1,m) - A), vecnorm(u1*ones(1,m) + A)]);
        if ind > m
            ind = ind - m;
        end

        if dist < 0.05 && norm(A_recon(:,ind)) < eps
            found = found + 1;
            A_recon(:,ind) = u1;
        end

        counter = counter + 1;
        
        if counter > 500
           % not found minimizer nf
           [~,nf] = min(vecnorm(A_recon));
           break 
        end
    end

end

function Lstar = approx_subspaceL(f,d,h,mx,k)

    hess_f = zeros(d^2,mx);

    % sample points for Hessian evaluation
    xi = randn(d,mx);
    xi = xi./vecnorm(xi);

    % matrix of vectorized Hessians
    for i = 1:mx
        hess_f(:,i) = reshape(hessNN_disc(f,xi(:,i),h), [1 d^2]);
    end

    % construct Lstar
    [U,~,~] = svd(hess_f);
    Lstar = U(:,1:k);
    
end

function u1 = projected_gradient_ascent(Lstar,N,gamma)

    Pstar = Lstar*Lstar';

    d = sqrt(size(Lstar,1));
    k = size(Lstar,2);

    W = orth(Lstar); 
    xi = randn(k,1);
    xi = ones(d^2,1)*xi';

    X = xi.*W;      
    X = sum(X,2);
    X = X./norm(X,'fro');

    X = reshape(X, [d d]);

    for i = 1:N
        Pg = P_gamma(X,gamma);
        Pg = reshape(Pg, [d^2, 1]);
        X = Pstar*Pg;
        X = reshape(X, [d d]);
        X = X./norm(X,'fro');
    end

    [U,~,~] = svd(X);
    u1 = U(:,1);

end

function Pg = P_gamma(X,gamma)
    [U,S,V] = svd(X);
    S(1,1) = gamma*S(1,1);

    a = 1/norm(diag(S));
    Pg = a*U*S*V';
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

function f = NN(X,A,b,t,sigma)
    [~,n] = size(X);
    f = sigma(A'*X + t'*ones(1,n));
    f = dot((ones(n,1)*b)',f);
end
