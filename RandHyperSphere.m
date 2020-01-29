function X = RandHyperSphere(m,n,half)
% 
% Projection method of Gaussian samples Muller (1959)
%
% REFERENCES:
% 
% Hicks, J. S. ad Wheeling, R. F. "An Efficient Method for Generating 
% Uniformly Distributed Points on the Surface of an n-Dimensional Sphere." 
% Comm. Assoc. Comput. Mach. 2, 13-15, 1959.
% 
% Marsaglia, G. "Choosing a Point from the Surface of a Sphere." 
% Ann. Math. Stat. 43, 645-646, 1972.
%
% Muller, M. E. "A Note on a Method for Generating Points Uniformly on 
% N-Dimensional Spheres." Comm. Assoc. Comput. Mach. 2, 19-20, Apr. 1959.

if nargin==2
    half=0;
end

X = randn(m,n);
X = X./vecnorm(X,2);

if half
    phi = InvHyperSphere(X);
    phi = mod(phi,pi);
    X = HyperSphere(phi);
end
