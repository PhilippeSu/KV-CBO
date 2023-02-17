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
    r = sqrt(sum(X(:,1).^2));

    % partial sums of radius like terms
    U = triu(ones(n-1,n));
    R = sqrt(U*(X.*X));

    % computation of the angles
    phi = acos((X(1:n-1,:)+eps)./(R+eps));

    % correcting last component
    phi(n-1,:) = phi(n-1,:).*(X(n,:)>=0)+(2*pi-phi(n-1,:)).*(X(n,:)<0); 

end
