function X = HyperSphere(phi,r) 
%
% syntax: X = HyperSphere(phi,r)
%         X = HyperSphere(phi)
% 
% INPUT ARGUMENTS
%     
%         -phi    (n - 1) x N every row of this matrix is a tuple of
%                 n-1 parametric coordinates and N, the number of points
%                 each column of phi stores the angles of one vector v
%
%                 phi(i,:), in [0,pi], i=1,...,n-2 
%                 phi(n-1,:), in [0,2pi)
%
%                 
%         -r      the radius, default value: r = 1
%         
% OUTPUT ARGUMENTS
% 
%         X       n x N every row of this matrix is a tuple of n
%                 catesian coordinates of the n-dimensional hypersphere
%
%         X(1,:)=r*cos(phi(1,:))
%         X(2,:)=r*cos(phi(2,:))*sin(phi(1,:))
%         
%         X(i,:)=r*cos(phi(i,:))*sin(phi(1,:))*...*sin(phi(i-1,:)), i=2,...,n-1
%         X(n,:)=r*sin(phi(1,:))*...*sin(phi(n-1,:))
%
% 
% this function calculates the cartesian coordinates of an n-dimensional hypersphere
% given the n-1 dimensions vector 'phi' of parametric coordinates and the
% radius 'r'
% 


if nargin==1
    r=1;
end

% Number of points

N = size(phi,2); 

% cosine terms 

C = [cos(phi); ones(1,N)];

% sine terms 

S = [ones(1,N); cumprod(sin(phi),1)];

% calculate output 

X = r * C .* S;




