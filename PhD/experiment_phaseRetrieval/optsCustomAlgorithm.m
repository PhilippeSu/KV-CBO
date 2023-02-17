function [x, outs] = optsCustomAlgorithm(A, At, y, z0, pars)

% convert handles to matrices
d = size(z0,1);
M = size(A(z0),1);

Amat = zeros(M,d);

for i=1:d  
  ei = [zeros(i-1,1);1;zeros(d-i,1)];
  Amat(:,i) = A(ei);
end

A = Amat;
At = A';

% recast problem
[d, M] =  size(At);

% KV-CBO method needs quadratic measurements
y = y.^2; 

% lower frame bound
S = zeros(d,d);
for i=1:d
    S(i,:) = sum(At*diag(At(i,:)),2);
end          
C1 = min(eig(S));

checkframebound(A,C1);

% define problem on the sphere
R = sqrt(norm(y,1)/C1);
y = y./(R*R);            
At = [At; zeros(1,M)];
A = At';

%% 

% z0 = initSpectral(A,At,y,d,true,true,false);
% z0 = z0/norm(z0);
% k = 200;

costfunction = PhaseRetrievalCostFunction(At,y,R,zeros(d,1),d+1);
% param = paramForPhaseRetrieval(d,x0,k);
param = paramForPhaseRetrieval(d);
 
[Va, ~] = KVCBO(costfunction,param); 
x = R*Va(1:end-1);
outs = [];

end

function checkframebound(A,C1)

[~,d] = size(A); 
N = 2;
z = randn(d,N);

a = sum(abs(A*z).^2,1);
b = C1*vecnorm(z).^2;

if b <= a
    % everything fine
elseif b > a
    error('frame bound error\n')
end

end

