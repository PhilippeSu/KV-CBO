function [phi,r] = InvHyperSphere(X) 
%
% syntax: [phi,r] = InvHyperSphere(X)
% 
% INPUT ARGUMENTS
%     
%         X       n x N every row of this matrix is a tuple of n
%                 catesian coordinates of the n-dimensional hypersphere
%         
% OUTPUT ARGUMENTS
% 
%         -phi    (n - 1) x N every row of this matrix is a tuple of n-1
%                 parametric coordinates and N, is the number of points
%                 
%         -r      the radius 
% 

n = size(X,1);

% assuming all points are on the sphere

r=sqrt(sum(X(:,1).^2));

% partial sums of radius like terms

U = triu(ones(n-1,n));
R = sqrt(U*(X.*X));

% computation of the angles

phi = acos((X(1:n-1,:)+eps)./(R+eps));

% correcting last component

phi(n-1,:)=phi(n-1,:).*(X(n,:)>=0)+(2*pi-phi(n-1,:)).*(X(n,:)<0); %TODO:

% 
% n = size(X,1);
% tol=1e-6;
% 
% r=sqrt(sum(X(:,1).^2));
% 
% r2 = sqrt(sum(X(n-1:n,:).^2));
% 
% phi(1:n-2,:)=acos(X(1:n-2,:)./r);
% phi(n-1,:)=acos(X(n-1,:)./(r2+eps)).*(r2>tol);
% phi(n-1,:)=phi(n-1,:).*(X(n,:)>=0)+(2*pi-phi(n-1,:)).*(X(n,:)<0);
% 
